"""P3: predict the other measured bioprocess variables in GSE30321 from expression.

Modality: transcriptomics only. Pipeline is cho_poc's, unchanged: nested CV, all
fitting inside the training fold, leak assertion on every split, cluster bootstrap.
Headline split is leave-one-experiment-out; leave-one-expression-family-out is the
strict check; culture-grouped is reported only alongside them.
"""
from __future__ import annotations

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import cho_poc as cp  # noqa: E402

RESULTS = ROOT / "results"
ENDPOINTS = {  # GEO characteristic -> (label, unit)
    "cell viability": ("Viability", "%"),
    "growth rate": ("Growth rate", "1/h"),
    "cell density": ("Viable cell density", "1e6 cells/mL"),
    "qlactate": ("q-lactate", "as deposited"),
    "qammonia": ("q-ammonium", "as deposited"),
    "qproductivity": ("Qp", "pg/cell/day"),
}
SCHEMES = ["leave_one_experiment_out", "leave_one_expression_cluster_out", "culture_grouped_5fold"]
SCHEME_UNIT = {"leave_one_experiment_out": "experiment", "leave_one_expression_cluster_out": "expression family",
               "culture_grouped_5fold": "culture"}


def temperature_baseline(meta, mask, y, cv_name, seed=cp.SEED):
    """Linear model on culture temperature alone, same split."""
    from sklearn.linear_model import LinearRegression

    idx = np.where(mask)[0]
    col, n_splits, scheme, _ = cp.CV_SCHEMES[cv_name]
    cultures = meta["culture_id"].to_numpy()[idx]
    exps = meta["experiment_id"].to_numpy()[idx]
    groups = {"culture": cultures, "experiment": exps, "cluster": meta["expression_cluster"].to_numpy()[idx]}[col]
    T = meta[["temperature"]].to_numpy(dtype=float)[idx]
    yy = y[idx]
    pred, null, fold = np.full(len(idx), np.nan), np.full(len(idx), np.nan), np.full(len(idx), -1)
    for f, (tr, te) in enumerate(cp.make_splits(groups, n_splits, seed)):
        cp.assert_no_group_leak(groups, tr, te)
        pred[te] = LinearRegression().fit(T[tr], yy[tr]).predict(T[te])
        null[te], fold[te] = yy[tr].mean(), f
    return pd.DataFrame({"y": yy, "pred": pred, "null": null, "fold": fold, "group": groups, "culture": cultures,
                         "experiment": exps, "array": meta.index.to_numpy()[idx]})


def run_endpoint(args):
    key, quick, n_boot = args
    t0 = time.time()
    D = cp.prepare()
    meta = D["meta"]
    label, unit = ENDPOINTS[key]
    y = meta[key].to_numpy(dtype=float)
    mask = ~np.isnan(y) & ~meta["qc_label_excluded"].to_numpy()
    rows, preds = [], []

    def add(df, cv, encoder, nfeat, model, stats, note=""):
        unit_col = cp.CV_SCHEMES[cv][3]
        ci = cp.cluster_bootstrap(df, stats, cluster_col=unit_col, n_boot=n_boot)
        for metric, (pt, lo, hi) in ci.items():
            rows.append({"phase": "P3", "modality": "transcriptomics (GSE30321 microarray)", "endpoint": label,
                         "endpoint_unit": unit, "holdout_scheme": cv, "held_out_unit": SCHEME_UNIT[cv],
                         "encoder": encoder, "n_features": nfeat, "model": model, "metric": metric,
                         "value": round(pt, 4), "ci_lo": round(lo, 4), "ci_hi": round(hi, 4), "ci_level": 0.95,
                         "n": len(df), "n_unit": "arrays", "n_cultures": df["culture"].nunique(),
                         "n_experiments": df["experiment"].nunique(),
                         "n_holdout_groups": int(df["fold"].nunique()), "bootstrap_unit": unit_col, "note": note})
        preds.append(df.assign(endpoint=label, holdout_scheme=cv, encoder=encoder, model=model))

    reg_models = ["pls"] if quick else ["pls", "enet", "gbm"]
    for cv in SCHEMES:
        for enc in ["raw", "pathway"]:
            for model in reg_models:
                if cv == "culture_grouped_5fold" and not (enc == "raw" and model == "pls"):
                    continue  # culture-level shown only as the leakage contrast for the primary model
                df = cp.run_cv(D, mask, y, model, enc, cv)
                add(df, cv, enc, cp.n_features_of(D, enc), model, cp.REG_STATS)
        if cv != "culture_grouped_5fold":
            add(temperature_baseline(meta, mask, y, cv), cv, "temperature only", 1, "linear", cp.REG_STATS,
                note="reference: culture temperature alone, no expression")
    # high vs low tercile
    lo, hi = np.nanquantile(y[mask], [1 / 3, 2 / 3])
    cmask = mask & ((y <= lo) | (y >= hi))
    ycls = (y >= hi).astype(float)
    for cv in SCHEMES:
        for enc in ["raw", "pathway"]:
            if cv == "culture_grouped_5fold" and enc != "raw":
                continue
            df = cp.run_cv(D, cmask, ycls, "logit", enc, cv)
            add(df, cv, enc, cp.n_features_of(D, enc), "logit", cp.CLF_STATS,
                note=f"top vs bottom tercile; cut-offs {lo:.4g} / {hi:.4g}")
    print(f"  {label}: n={int(mask.sum())} arrays, {meta.loc[mask, 'experiment_id'].nunique()} experiments, "
          f"{meta.loc[mask, 'expression_cluster'].nunique()} families ({time.time() - t0:.0f}s)", flush=True)
    return rows, pd.concat(preds)


def main(quick=False, n_boot=2000):
    keys = [k for k in ENDPOINTS if k != "qproductivity"]
    with Pool(min(4, len(keys))) as pool:
        out = pool.map(run_endpoint, [(k, quick, 200 if quick else n_boot) for k in keys])
    rows = [r for o in out for r in o[0]]
    m = pd.DataFrame(rows)
    m.to_csv(RESULTS / "metrics_endpoints.csv", index=False)
    pd.concat([o[1] for o in out]).to_csv(RESULTS / "oof_predictions_endpoints.csv", index=False)
    show = m[(m.encoder == "raw") & (m.model.isin(["pls", "logit"])) & (m.metric.isin(["rmse", "r2_cv", "auc", "rmse_null"]))]
    pd.set_option("display.width", 220)
    print(show.pivot_table(index=["endpoint", "holdout_scheme"], columns="metric", values="value").round(3).to_string())
    return m


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
