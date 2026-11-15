"""P2: knockout ranking with blinded literature recovery.

Stage 1 (`freeze`): build the rankings. No ground truth is read.
Stage 2 (`evaluate`): annotate with literature ground truth and score.
Pre-registration is in DECISIONS.md; the freeze commit precedes ground-truth collection.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import cho_poc as cp  # noqa: E402

RESULTS = ROOT / "results"
ICHO_URL = "https://bigg.ucsd.edu/static/models/iCHOv1.json"
ICHO_PATH = ROOT / "raw" / "icho" / "iCHOv1.json"
RAW_KEGG = ROOT / "raw" / "kegg"
ICHO_SHA256 = "0b786b81c08142be3a88f33b366d4e0f05d8689baad877e5375a9fa843a5f87c"


def gene_delta():
    """Route A. Gene-level delta on GSE30321 (transcriptomics), within experiment and raw."""
    D = cp.prepare()
    meta = D["meta"]
    y = meta[cp.ENDPOINT].to_numpy(dtype=float)
    lab = ~np.isnan(y) & ~meta["qc_label_excluded"].to_numpy()
    X, y = D["X"][lab], y[lab]
    exps = meta["experiment_id"].to_numpy()[lab]
    sym = np.array([s if isinstance(s, str) else "" for s in D["symbols"]])

    def delta(Xm, ym):
        Z = (Xm - Xm.mean(axis=0)) / np.where(Xm.std(axis=0) > 1e-8, Xm.std(axis=0), 1.0)
        lo, hi = np.quantile(ym, [0.25, 0.75])
        return Z[ym >= hi].mean(axis=0) - Z[ym <= lo].mean(axis=0)

    multi = pd.Series(exps).map(pd.Series(exps).value_counts()).to_numpy() >= 4
    d_within = delta(cp.center_within(X[multi], exps[multi]), cp.center_within(y[multi][:, None], exps[multi]).ravel())
    d_raw = delta(X, y)
    df = pd.DataFrame({"symbol": sym, "probe": D["probes"], "delta_within_experiment": d_within, "delta_raw": d_raw})
    df = df[df.symbol != ""]
    g = df.groupby("symbol").agg(n_probesets=("probe", "size"), delta_within_experiment=("delta_within_experiment", "mean"),
                                 delta_raw=("delta_raw", "mean")).reset_index()
    g = g.sort_values("delta_within_experiment").reset_index(drop=True)
    g.insert(0, "knockout_rank", np.arange(1, len(g) + 1))            # 1 = most negative delta
    g["overexpression_rank"] = len(g) - g["knockout_rank"] + 1          # 1 = most positive delta
    g["modality"] = "transcriptomics (GSE30321)"
    info = {"n_genes": int(len(g)), "n_arrays_within": int(multi.sum()),
            "n_experiments_within": int(len(np.unique(exps[multi])))}
    return g, info


def fba_deletions(processes=2):
    """Route B. iCHO1766 single-gene deletions: growth and IgG-synthesis capacity."""
    import cobra
    from cobra.flux_analysis import single_gene_deletion

    if not ICHO_PATH.exists():
        import requests

        ICHO_PATH.parent.mkdir(parents=True, exist_ok=True)
        ICHO_PATH.write_bytes(requests.get(ICHO_URL, timeout=300).content)
    sha = cp.sha256_of(ICHO_PATH)
    if sha != ICHO_SHA256:
        print(f"  WARNING iCHOv1.json SHA256 {sha} differs from the pinned {ICHO_SHA256}")
    model = cobra.io.load_json_model(str(ICHO_PATH))
    model.solver = "glpk"
    wt_growth = model.slim_optimize()
    growth = single_gene_deletion(model, processes=processes)
    pm = model.copy()
    pm.reactions.DM_igg_g.bounds = (0.0, 1000.0)
    pm.reactions.BIOMASS_cho.lower_bound = 0.2 * wt_growth
    pm.objective = "DM_igg_g"
    wt_prod = pm.slim_optimize()
    prod = single_gene_deletion(pm, processes=processes)

    def tidy(res, wt, col):
        ids = [next(iter(s)) for s in res["ids"]]
        v = res["growth"].to_numpy(dtype=float)
        v = np.where(np.isnan(v) | (res["status"].to_numpy() != "optimal"), 0.0, v)
        return pd.DataFrame({"gene_id": ids, col: np.clip(v / wt, 0, None)})

    out = tidy(growth, wt_growth, "growth_ratio").merge(tidy(prod, wt_prod, "igg_capacity_ratio"), on="gene_id")
    names = {g.id: g.name for g in model.genes}
    nrx = {g.id: len(g.reactions) for g in model.genes}
    subs = {g.id: "; ".join(sorted({r.subsystem for r in g.reactions if r.subsystem})[:4]) for g in model.genes}
    sym = out.gene_id.map(names).str.upper()
    # model gene names that are bare LOC identifiers get the KEGG symbol where KEGG has one
    kegg = RAW_KEGG / "cge_genes.txt"
    if not kegg.exists():
        import requests

        kegg.parent.mkdir(parents=True, exist_ok=True)
        kegg.write_bytes(requests.get("https://rest.kegg.jp/list/cge", timeout=300).content)
    ks = {}
    for line in kegg.read_text().splitlines():
        p = line.split("\t")
        if ";" in p[-1]:
            ks[p[0].replace("cge:", "")] = p[-1].split(";")[0].split(",")[0].strip().upper()
    is_loc = sym.isna() | sym.astype(str).str.startswith("LOC")
    out["symbol"] = [ks.get(g, s) if loc and ks.get(g) else s for g, s, loc in zip(out.gene_id, sym, is_loc)]
    out["n_reactions"] = out.gene_id.map(nrx)
    out["subsystems"] = out.gene_id.map(subs)
    out["tolerated"] = (out.growth_ratio >= 0.9) & (out.igg_capacity_ratio >= 0.9)
    out["tradeoff_score"] = np.where((out.igg_capacity_ratio >= 0.99) & (out.growth_ratio >= 0.2),
                                     1 - out.growth_ratio, 0.0)
    out = out.sort_values(["tradeoff_score", "growth_ratio", "gene_id"], ascending=[False, True, True]).reset_index(drop=True)
    out["modality"] = "mechanistic model (iCHO1766 FBA), no omics data"
    out = out[["gene_id", "symbol", "growth_ratio", "igg_capacity_ratio", "tolerated", "tradeoff_score", "n_reactions",
               "subsystems", "modality"]]
    info = {"wt_growth_per_h": float(wt_growth), "wt_igg_flux": float(wt_prod), "n_genes": int(len(out)),
            "n_tolerated": int(out.tolerated.sum()), "n_lethal_growth_lt_0.01": int((out.growth_ratio < 0.01).sum()),
            "n_tradeoff_positive": int((out.tradeoff_score > 1e-6).sum())}
    return out, info


def freeze():
    g, ia = gene_delta()
    g.round(4).to_csv(RESULTS / "knockout_ranking.csv", index=False)
    print("Route A frozen:", ia)
    print(g.head(15)[["knockout_rank", "symbol", "delta_within_experiment", "delta_raw"]].to_string(index=False))
    f, ib = fba_deletions()
    f.round(5).to_csv(RESULTS / "knockout_fba.csv", index=False)
    print("Route B frozen:", ib)
    (RESULTS / "knockout_freeze_info.json").write_text(json.dumps({"route_a": ia, "route_b": ib}, indent=2))


if __name__ == "__main__":
    if "freeze" in sys.argv:
        freeze()
    elif "evaluate" in sys.argv:
        import knockout_eval

        knockout_eval.main()
