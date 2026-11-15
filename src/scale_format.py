"""P5 (opportunistic): culture-format contrast on GSE104787 (transcriptomics).

Shake flask (uncontrolled) vs bioreactor (controlled), one cell line, 18 arrays,
6 culture units. A format contrast, NOT a scale contrast; descriptive only.
Pre-registered in DECISIONS.md.
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
import cho_poc as cp  # noqa: E402

RAW, RESULTS = ROOT / "raw", ROOT / "results"
SEED = 20261116
MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104787/matrix/GSE104787_series_matrix.txt.gz"
GPL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL24076&targ=self&form=text&view=data"


def load():
    import yao

    mp = yao._get(MATRIX, RAW / "geo2/GSE104787_series_matrix.txt.gz")
    gp = yao._get(GPL, RAW / "geo2/GPL24076.txt")
    expr, meta = cp.parse_series_matrix(mp)
    lines = gp.read_text(encoding="utf-8", errors="replace").split("\n")
    b = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
    e = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    ann = pd.read_csv(io.StringIO("\n".join(lines[b + 1:e])), sep="\t", dtype=str, keep_default_na=False).set_index("ID")
    ga = ann["gene_assignment"]
    sym = ga.str.split(" // ").str[1].where(ga != "---").str.upper()
    sym = sym[~sym.fillna("LOC").str.startswith("LOC")]
    return expr, meta, sym


def main(n_null=1000):
    expr, meta, sym = load()
    expr.index = expr.index.astype(str)
    system = meta["title"].str.extract(r"^(flask|bioreactor)")[0].to_numpy()
    day = meta["title"].str.extract(r"day(\d)")[0].to_numpy()
    X = expr.T.to_numpy(dtype=float)
    X = cp.center_within(X, day)
    Z = (X - X.mean(axis=0)) / np.where(X.std(axis=0) > 1e-8, X.std(axis=0), 1.0)
    d = Z[system == "bioreactor"].mean(axis=0) - Z[system == "flask"].mean(axis=0)
    g = pd.DataFrame({"symbol": sym.reindex(expr.index).to_numpy(), "delta": d}).dropna().groupby("symbol").delta.mean()
    hall = cp.read_gmt(RAW / cp.SOURCES["hallmark"]["path"])
    rng = np.random.default_rng(SEED)
    vals = g.to_numpy()
    rows = []
    for name, genes in hall.items():
        m = [x for x in set(genes) if x in g.index]
        if len(m) < 10:
            continue
        s = float(g[m].mean())
        null = np.array([vals[rng.choice(len(vals), len(m), replace=False)].mean() for _ in range(n_null)])
        rows.append({"pathway": name, "n_genes": len(m), "delta": s, "p_vs_random_sets": float((1 + (np.abs(null) >= abs(s)).sum()) / (n_null + 1)),
                     "top_genes": "; ".join(f"{k} ({v:+.2f})" for k, v in g[m].sort_values(ascending=s < 0).head(6).items())})
    out = pd.DataFrame(rows)
    out["abs"] = out.delta.abs()
    out = out.sort_values("abs", ascending=False).drop(columns="abs").reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    out.round(4).to_csv(RESULTS / "scale_format_pathways.csv", index=False)
    # share of expression variance explained by sampling day and by culture system (median over probesets)
    Xr = expr.T.to_numpy(dtype=float)
    tot = ((Xr - Xr.mean(axis=0)) ** 2).sum(axis=0)
    ss_day = sum((day == k).sum() * (Xr[day == k].mean(axis=0) - Xr.mean(axis=0)) ** 2 for k in np.unique(day))
    ss_sys = sum((system == k).sum() * (Xr[system == k].mean(axis=0) - Xr.mean(axis=0)) ** 2 for k in np.unique(system))
    ok = tot > 1e-8
    var_day, var_sys = float(np.median(ss_day[ok] / tot[ok])), float(np.median(ss_sys[ok] / tot[ok]))
    hyp = out[out.pathway == "HALLMARK_HYPOXIA"].iloc[0]
    info = {"dataset": "GSE104787", "n_arrays": int(len(meta)), "n_flask_cultures": 3, "n_bioreactor_cultures": 3,
            "genes": int(len(g)), "hallmark_sets_ranked": int(len(out)), "hypoxia_rank": int(hyp["rank"]),
            "hypoxia_delta": float(hyp.delta), "hypoxia_p_vs_random_sets": float(hyp.p_vs_random_sets),
            "n_sets_p_lt_0.05": int((out.p_vs_random_sets < 0.05).sum()), "top5": out.pathway.head(5).tolist(),
            "median_variance_share_day": round(var_day, 3), "median_variance_share_system": round(var_sys, 3)}
    base = dict(phase="P5", box="Scale risk", evidence_type="blinded recovery (culture format, not scale)",
                modality="transcriptomics (GSE104787, Affymetrix CHO Gene array)", endpoint="bioreactor vs shake flask contrast",
                holdout_scheme="not applicable: descriptive contrast, no model fitted; 6 culture units", encoder="MSigDB Hallmark",
                model="delta (bioreactor minus flask), day-centred", n=int(len(meta)), n_unit="arrays (3 flasks, 3 bioreactors, 3 days)",
                source_data="GSE104787", ci_level="")
    info["metric_rows"] = [
        {**base, "metric": f"rank of HALLMARK_HYPOXIA among {len(out)} Hallmark sets (by |delta|)", "value": int(hyp["rank"]),
         "p_value": float(hyp.p_vs_random_sets), "note": f"pre-registered target from Gao et al. 2016; delta = {hyp.delta:+.3f} (positive = higher in bioreactor)"}]
    info["metric_rows"].append({**base, "encoder": "none (probeset level)", "model": "variance decomposition",
                                "metric": "median share of expression variance explained by culture system", "value": round(var_sys, 3),
                                "note": f"sampling day explains {var_day:.3f}; the deposit's summary reports large differences over time and much smaller ones between systems"})
    (RESULTS / "scale_format.json").write_text(json.dumps(info, indent=2))
    print(json.dumps({k: v for k, v in info.items() if k != "metric_rows"}, indent=1))
    print(out[["rank", "pathway", "n_genes", "delta", "p_vs_random_sets"]].head(15).to_string(index=False))
    return out, info


if __name__ == "__main__":
    main()
