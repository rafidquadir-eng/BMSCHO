"""Collect every phase-2 metric into results/metrics_v2.csv.

One row per metric. Every row names its modality and its holdout scheme; a row
without both is not reportable (brief, non-negotiable 4).
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "results"
COLS = ["phase", "box", "evidence_type", "modality", "endpoint", "holdout_scheme", "encoder", "model", "metric", "value",
        "ci_lo", "ci_hi", "ci_level", "p_value", "n", "n_unit", "n_experiments", "n_holdout_groups", "source_data", "note"]
NA_HOLDOUT = "not applicable: descriptive contrast, no model fitted"


def rows_p1():
    info = json.loads((R / "yao_info.json").read_text())
    pr = pd.read_csv(R / "yao_pathway_ranking.csv")
    out = []
    base = dict(phase="P1", box="Metabolomics / delta branch", evidence_type="blinded recovery (concordance on published summary statistics)",
                modality="metabolomics (extracellular, untargeted LC-MS; Yao et al. 2021 Table S3)", endpoint="qP association",
                holdout_scheme=NA_HOLDOUT, encoder="KEGG reference metabolic pathways", model="delta statistic (mean Fisher z)",
                n=info["unique_kegg"], n_unit="metabolites (12 clones behind each correlation)", source_data="PMC8704136 supplement",
                ci_level="")
    for pid, nm in [("map00020", "citrate cycle"), ("map00250", "alanine, aspartate and glutamate metabolism")]:
        r = pr[pr.pathway == pid].iloc[0]
        out.append({**base, "metric": f"rank of {nm} among {len(pr)} pathways", "value": int(r.rank_delta), "p_value": float(r.p_perm),
                    "note": f"pre-registered; rank by authors' over-representation statistic = {int(r.rank_ora)}"})
    out.append({**base, "metric": "mean rank of the two target pathways", "value": info["mean_rank_targets"], "p_value": info["p_joint"],
                "note": f"pre-registered joint test; expected {info['expected_mean_rank_under_null']:.1f} under the null"})
    out.append({**base, "endpoint": "qP of library 2 (mAb C)", "metric": "cross-library Spearman / AUC", "value": "",
                "evidence_type": "predictive", "holdout_scheme": "train library 1, test library 2 (not run)", "n": 0, "n_unit": "clones",
                "note": "NOT RUN: per-clone metabolite data are not public (data availability: not applicable). Requires partner data."})
    return out


def rows_p2():
    e = json.loads((R / "knockout_eval.json").read_text())
    out = []
    a = dict(phase="P2", box="Ranked fixes / delta branch", evidence_type="blinded recovery", modality="transcriptomics (GSE30321 microarray)",
             endpoint="validated CHO knockout/knockdown targets (Schulz et al. 2026)", holdout_scheme=NA_HOLDOUT + "; ranking frozen before ground truth was loaded",
             encoder="gene-level delta, within experiment", model="rank by most negative delta", n_unit="validated targets on the array",
             source_data="GSE30321; PMC12818826 suppl. sheet B", ci_level=0.95)
    k = e["route_a_knockout_end_within_experiment"]["k100"]
    out.append({**a, "metric": "AUC (validated targets vs other array genes)", "value": k["auc"], "ci_lo": k["auc_null_95"][0], "ci_hi": k["auc_null_95"][1],
                "p_value": k["p_auc"], "n": k["n_positives"], "note": f"interval is the null 95% range, not a CI; {e['route_a_coverage']['ko_positives_total']} validated targets exist, {k['n_positives']} are on the array"})
    out.append({**a, "metric": "validated targets in top 100", "value": k["hits_in_top_k"], "p_value": k["p_hits"], "n": k["n_positives"],
                "note": f"expected {k['expected_hits']} by chance", "ci_level": ""})
    o = e["route_a_overexpression_end_within_experiment"]["k100"]
    out.append({**a, "endpoint": "validated CHO overexpression targets (Schulz et al. 2026)", "model": "rank by most positive delta",
                "metric": "AUC (validated targets vs other array genes)", "value": o["auc"], "ci_lo": o["auc_null_95"][0], "ci_hi": o["auc_null_95"][1],
                "p_value": o["p_auc"], "n": o["n_positives"], "n_unit": "validated targets on the array", "note": "interval is the null 95% range, not a CI"})
    b = e["route_b_tradeoff_score"]
    m = dict(a, modality="mechanistic model (iCHO1766 FBA); no omics data", encoder="single-gene deletion", model="trade-off score",
             source_data="BiGG iCHOv1; PMC12818826 suppl. sheet B", n_unit="validated targets in the model")
    out.append({**m, "metric": "AUC (validated targets vs other model genes)", "value": b["auc"], "ci_lo": b["auc_null_95"][0], "ci_hi": b["auc_null_95"][1],
                "p_value": b["p_auc"], "n": b["n_positives"], "note": "interval is the null 95% range, not a CI"})
    t = e["route_b_tolerance"]
    out.append({**m, "model": "tolerance call (growth and IgG capacity >= 0.9 of wild type)", "metric": "share of validated targets called tolerated",
                "value": round(t["share_tolerated_ko_positives"], 3), "n": b["n_positives"], "ci_level": "",
                "note": f"share for all model genes = {t['share_tolerated_all_model_genes']:.3f}; uninformative"})
    mz = e["route_b_marzluf"]
    out.append({**m, "endpoint": "growth-incompatible knockouts named by Marzluf et al. 2026", "model": "growth ratio after deletion",
                "metric": "AUC (named incompatible vs other screened candidates in the model)", "value": mz["auc_growth_ratio_low_predicts_incompatible"],
                "n": len(mz["named_incompatible_in_model"]), "n_unit": "named incompatible genes in the model", "ci_level": "",
                "source_data": "BiGG iCHOv1; PMC12902299", "note": "FBA calls every one of them tolerated; the 36 retained knockouts are not published by name"})
    return out


def rows_p3():
    m = pd.read_csv(R / "metrics_endpoints.csv")
    out = []
    for r in m.itertuples():
        out.append({"phase": "P3", "box": "Viability head" if r.endpoint == "Viability" else "Other culture-health endpoints",
                    "evidence_type": "predictive", "modality": r.modality, "endpoint": r.endpoint, "holdout_scheme": r.holdout_scheme,
                    "encoder": r.encoder, "model": r.model, "metric": r.metric, "value": r.value, "ci_lo": r.ci_lo, "ci_hi": r.ci_hi,
                    "ci_level": 0.95, "p_value": "", "n": r.n, "n_unit": "arrays", "n_experiments": r.n_experiments,
                    "n_holdout_groups": r.n_holdout_groups, "source_data": "GSE30321",
                    "note": (r.note if isinstance(r.note, str) else "") + ("; culture-level split, shown only next to the experiment-level one"
                                                                         if r.holdout_scheme == "culture_grouped_5fold" else "")})
    return out


def rows_p4():
    a = json.loads((R / "alignment.json").read_text())
    out = []
    for vname, blk in a["leg1_growth_phase_transcriptomics_vs_proteomics"].items():
        base = dict(phase="P4", box="Proteomics / pathway-constrained encoder", evidence_type="cross-modality alignment (no fused prediction)",
                    modality="transcriptomics (GSE30321) vs proteomics (PRIDE PXD053511); unpaired", endpoint="growth-phase contrast",
                    holdout_scheme=NA_HOLDOUT + "; each contrast computed inside one modality", model="Spearman of contrast vectors",
                    source_data="GSE30321; PXD053511", ci_level=0.95)
        g = blk["gene_level"]
        out.append({**base, "encoder": "none (gene level)", "metric": "Spearman", "value": g["spearman"], "p_value": g["p_perm"], "n": g["n_genes"],
                    "n_unit": "genes measured in both", "note": vname, "ci_level": ""})
        for basis in ["iCHO1766 subsystems", "KEGG pathways (Chinese hamster)"]:
            r = blk[basis]
            out.append({**base, "encoder": basis, "metric": "Spearman", "value": r["spearman"], "ci_lo": r["null_95"][0], "ci_hi": r["null_95"][1],
                        "p_value": r["p_perm"], "n": r["n_pathways"], "n_unit": "pathways",
                        "note": vname + "; interval is the null 95% range, not a CI"})
    blk = a["leg2_qp_transcriptomics_vs_metabolomics"]
    for basis in ["iCHO1766 subsystems", "KEGG pathways (Chinese hamster)"]:
        r = blk[basis]
        out.append({"phase": "P4", "box": "Metabolomics / pathway-constrained encoder", "evidence_type": "cross-modality alignment (no fused prediction)",
                    "modality": "transcriptomics (GSE30321) vs metabolomics (Yao 2021 Table S3); unpaired", "endpoint": "Qp association",
                    "holdout_scheme": NA_HOLDOUT + "; each contrast computed inside one modality", "encoder": basis, "model": "Spearman of contrast vectors",
                    "metric": "Spearman", "value": r["spearman"], "ci_lo": r["null_95"][0], "ci_hi": r["null_95"][1], "ci_level": 0.95,
                    "p_value": r["p_perm"], "n": r["n_pathways"], "n_unit": "pathways", "source_data": "GSE30321; PMC8704136",
                    "note": "pre-registered; interval is the null 95% range, not a CI"})
    return out


def rows_p5():
    f = R / "scale_format.json"
    if not f.exists():
        return []
    return json.loads(f.read_text()).get("metric_rows", [])


def main():
    rows = rows_p1() + rows_p2() + rows_p3() + rows_p4() + rows_p5()
    df = pd.DataFrame(rows).reindex(columns=COLS)
    assert df.modality.notna().all() and (df.modality != "").all(), "metric without a modality"
    assert df.holdout_scheme.notna().all() and (df.holdout_scheme != "").all(), "metric without a holdout scheme"
    df.to_csv(R / "metrics_v2.csv", index=False)
    print(f"metrics_v2.csv: {len(df)} rows; by phase {df.phase.value_counts().to_dict()}")
    return df


if __name__ == "__main__":
    main()
