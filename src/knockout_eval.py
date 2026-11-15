"""P2 stage 2: score the frozen knockout rankings against literature ground truth.

Ground truth:
  * Schulz et al. 2026 (mAbs; PMC12818826) supplementary sheet B_SingleTargetData:
    every engineered target with mechanism and maximum reported fold improvement.
  * Marzluf et al. 2026 (Front Bioeng Biotechnol; PMC12902299): only the knockouts
    the text names as growth-incompatible or performance-reducing. The identity of
    the 36 retained knockouts is not published by name (it is a bar in Figure 2G).
The rankings were frozen in commit de3fe41 before either file was opened.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
LIT = ROOT / "raw" / "lit"
SEED = 20261116
SCHULZ_SUPP = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12818826/supplementaryFiles"

# Target names as written in the review -> official symbols. Curated by hand, once.
ALIAS = {
    "BAK": ["BAK1"], "BLIMP": ["PRDM1"], "BLIMP1": ["PRDM1"], "BLIMP1Β": ["PRDM1"], "CMYC": ["MYC"],
    "BCKDHA/B": ["BCKDHA", "BCKDHB"], "BMPRIA": ["BMPR1A"], "BMPRII": ["BMPR2"], "BCL-2": ["BCL2"],
    "BCL-XL": ["BCL2L1"], "BIP": ["HSPA5"], "GRP94": ["HSP90B1"], "C12ORF35": ["RESF1"], "CASPASE-7": ["CASP7"],
    "CERT": ["CERT1"], "CHOP": ["DDIT3"], "COFILIN1": ["CFL1"], "CYPB": ["PPIB"], "ERGIC-53": ["LMAN1"],
    "ERGIC53": ["LMAN1"], "ERO1-LA": ["ERO1A"], "ERO1LA": ["ERO1A"], "ERP57": ["PDIA3"], "FAM60A": ["SINHCAF"],
    "GADD45": ["GADD45A"], "GFAT2": ["GFPT2"], "GS": ["GLUL"], "HSP27": ["HSPB1"], "IGF-1R": ["IGF1R"],
    "IL-17F": ["IL17F"], "LDH": ["LDHA"], "PDHK1": ["PDK1"], "MDHII": ["MDH2"], "MT-SP1": ["ST14"],
    "MCL-1": ["MCL1"], "NDPK-A": ["NME1"], "P115": ["USO1"], "PACESOL": ["FURIN"], "PDI": ["P4HB"],
    "PERK": ["EIF2AK3"], "PGC-1": ["PPARGC1A"], "PKM1": ["PKM"], "PPARΑ": ["PPARA"], "SNAP-23": ["SNAP23"],
    "SPR14": ["SRP14"], "SVIP1": ["SVIP"], "TAUT": ["SLC6A6"], "ALT1": ["GPT"], "XBP-1S": ["XBP1"],
    "XBP1S": ["XBP1"], "XBP-1": ["XBP1"], "XCT": ["SLC7A11"], "CATALASE": ["CAT"], "HQSOX1": ["QSOX1"],
    "HSQSOX1B": ["QSOX1"], "SURVIVIN": ["BIRC5"], "HSURVIVIN": ["BIRC5"], "P180": ["RRBP1"], "AMPK": ["PRKAA1"],
    "PRAK": ["MAPKAPK5"], "ATF6AC": ["ATF6"], "ATF6C": ["ATF6"], "ATF6Β": ["ATF6B"], "GLUT5": ["SLC2A5"],
    "PRDM1": ["PRDM1"], "(M)CDC25B": ["CDC25B"],
}
# Whole-string patterns that expand to several genes.
EXPAND = [
    (r"Neu 1, 2, 3", "NEU1, NEU2, NEU3"), (r"PDHK1, 2, 3", "PDK1, PDK2, PDK3"), (r"PKM1&2", "PKM"),
    (r"MPC1&2", "MPC1, MPC2"), (r"Erp27&57", "ERP27, PDIA3"), (r"DNMT1 & 3a", "DNMT1, DNMT3A"),
    (r"TET2 & 3", "TET2, TET3"), (r"SRP14, ?9, ?54", "SRP14, SRP9, SRP54"), (r"SR \(a ?& ?b\)", "SRPRA, SRPRB"),
    (r"Sec61 \(a, b, y\)", "SEC61A1, SEC61B, SEC61G"), (r"BCKDHa/b", "BCKDHA, BCKDHB"), (r"BAX/BAK", "BAX, BAK"),
    (r"KDEL receptor 1", "KDELR1"), (r"HSPA5/BIP", "HSPA5"), (r"BLIMP1/prdm1", "PRDM1"), (r"PRDM1/BLIMP1", "PRDM1"),
    (r"SLC2A5/ GLUT5", "SLC2A5"), (r"hSurvivin/BIRC5", "BIRC5"), (r"LDHa single allele", "LDHA"),
]
SKIP_NAME = re.compile(r"miR|let[-‐]7|lnc|mitosRNA|SINEUP|UCOE|Aneuploidy|Chromosome|directed evolution|30Kc|lvp1|^LC\b|"
                       r"mutat|mutant|const active|S132A|D816V|YAP5SA|delta|R98S|HEX", re.I)


def _fetch_schulz():
    import io
    import zipfile

    import requests

    dest = LIT / "PMC12818826" / "KMAB_A_2615475_SM0321.xlsx"
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        z = zipfile.ZipFile(io.BytesIO(requests.get(SCHULZ_SUPP, timeout=300).content))
        dest.write_bytes(z.read("KMAB_A_2615475_SM0321.xlsx"))
    return dest


def parse_schulz() -> pd.DataFrame:
    b = pd.read_excel(_fetch_schulz(), sheet_name="B_SingleTargetData")
    b.columns = [str(c).strip() for c in b.columns]
    rows = []
    for _, r in b.iterrows():
        name = str(r["Target engingeered"]).replace("\n", " ").strip()
        mech = str(r["Engineering mechanism"]).upper()
        fold = pd.to_numeric(r["Max. improvement fold"], errors="coerce")
        if SKIP_NAME.search(name):
            continue
        ko = ("KNOCKOUT" in mech) or ("KNOCKDOWN" in mech)
        oe = "OVEREXPRESSION" in mech
        if ("OEX" in mech) or ("KD &" in mech) or (ko and oe) or "OTHERS" in mech or "FUNCTIONAL VARIANT" in mech:
            cls = "mixed/other (excluded)"
        elif ko:
            cls = "KO/KD"
        elif oe:
            cls = "OE"
        else:
            cls = "mixed/other (excluded)"
        s = name
        for pat, rep in EXPAND:
            s = re.sub(pat, rep, s, flags=re.I)
        s = re.sub(r"\([^)]*\)", " ", s)
        toks = [t.strip().upper().replace("‑", "-").replace("‐", "-") for t in re.split(r"[,;&+]|\band\b", s) if t.strip()]
        genes = []
        for t in toks:
            t = t.strip()
            if t in ALIAS:
                genes += ALIAS[t]
            elif re.fullmatch(r"[A-Z0-9][A-Z0-9-]{1,14}", t):
                genes.append(t)
        combo = "COMBINATION" in mech or len(set(genes)) > 1
        for g in sorted(set(genes)):
            rows.append({"gene": g, "class": cls, "fold": float(fold) if pd.notna(fold) else np.nan, "combination": combo,
                         "target_as_written": name, "pathway": str(r.get("Pathway/ Mechanism", "")),
                         "reference": str(r.get("Reference/ DOI", "")).strip(), "year": r.get("Publication year", "")})
    return pd.DataFrame(rows)


def ground_truth():
    g = parse_schulz()
    g = g[g["class"].isin(["KO/KD", "OE"])]
    out = []
    for (gene, cls), d in g.groupby(["gene", "class"]):
        out.append({"gene": gene, "class": cls, "max_fold": d.fold.max(), "min_fold": d.fold.min(), "n_reports": len(d),
                    "raised_titer_or_qp": bool((d.fold > 1).any()), "single_target_evidence": bool((~d.combination).any()),
                    "targets_as_written": " || ".join(sorted(set(d.target_as_written)))[:200],
                    "references": "; ".join(sorted(set(d.reference)))[:300],
                    "source": "Schulz et al. 2026 (PMC12818826), suppl. sheet B"})
    gt = pd.DataFrame(out)
    # BOK appears only inside the BAX/BAK/BOK triple knockout; the literature reports that deleting
    # BOK adds nothing to BAX/BAK. It is treated as a named negative, not a positive.
    gt.loc[(gt.gene == "BOK") & (gt["class"] == "KO/KD"), "raised_titer_or_qp"] = False
    gt.loc[gt.gene == "BOK", "source"] = "Schulz et al. 2026 (triple KO only); no added effect of BOK per MacDonald et al."
    return gt


MARZLUF_INCOMPATIBLE = {  # named in the text of Marzluf et al. 2026 (PMC12902299)
    "GAPDH": "pool failed to recover", "ALDOA": "pool failed to recover", "EEF2": "pool failed to recover",
    "EEF1A1": "reverted to wild type (selective disadvantage)", "PGK1": "reverted to wild type (selective disadvantage)",
    "LDHA": "reverted to wild type (selective disadvantage)", "TINAGL1": "reverted to wild type (selective disadvantage)",
    "NME1": "reverted to wild type (selective disadvantage)", "PRDX1": "reduced growth or productivity",
    "CTSB": "reduced growth or productivity", "HSPA5": "reduced performance", "HSP90B1": "reduced performance",
    "GPC1": "reduced performance",
}
MARZLUF_OK = {"FN1": "improved titer (~2-fold)", "DAG1": "neutral/beneficial", "LGALS1": "neutral/beneficial"}
MARZLUF_CANDIDATES = ["ACLY", "ACTB", "AEBP1", "AKR1B1", "ALDOA", "ANXA2", "BGN", "CALR", "CLSTN1", "CLU", "CSPG4", "CTSA",
                      "CTSB", "CTSZ", "DAG1", "EEF1A1", "EEF2", "ENO3", "FABP4", "FASN", "FLNA", "FN1", "FSTL1", "FUCA1",
                      "GAPDH", "GPC1", "GSTM6", "GSTP1", "HIST1H1E", "HSP90AA1", "HSP90AB1", "HSP90B1", "HSPA5", "HSPA8",
                      "HSPG2", "HTRA1", "ITIH5", "LAMA5", "LAMB1", "LDHA", "LGALS1", "LGALS3BP", "LGMN", "LPL", "NCL",
                      "NID1", "NME1", "NUCB2", "PCOLCE", "PGK1", "PKM", "PLBD2", "PLOD1", "PPIA", "PPIB", "PRDX1", "PSAP",
                      "PXDN", "QSOX1", "SPP1", "THBS1", "TIMP1", "TINAGL1", "TKT", "TNXB", "TUBA1B", "TUBB1", "UBC",
                      "VIM", "YWHAZ"]
NAMED_IN_BRIEF = ["FUT8", "BAX", "BAK1", "BOK", "B4GALT1", "BCKDHA", "BCKDHB", "GLUL", "LDHA", "C1GALT1C1", "NEU1", "NEU2", "NEU3"]


def auc_of(score, pos):
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(pos.astype(int), score)) if 0 < pos.sum() < len(pos) else float("nan")


def enrichment(score, pos, k, n_perm=10000, seed=SEED):
    """score: higher = stronger candidate. Returns AUC, top-k hits, and permutation p-values."""
    rng = np.random.default_rng(seed)
    order = np.argsort(-score)
    top = np.zeros(len(score), bool)
    top[order[:k]] = True
    n_pos = int(pos.sum())
    hits = int((top & pos).sum())
    a = auc_of(score, pos)
    null_hits, null_auc = np.empty(n_perm), np.empty(n_perm)
    ranks = pd.Series(score).rank().to_numpy()
    for b in range(n_perm):
        idx = rng.choice(len(score), n_pos, replace=False)
        null_hits[b] = top[idx].sum()
        null_auc[b] = (ranks[idx].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * (len(score) - n_pos))
    return {"n_universe": int(len(score)), "n_positives": n_pos, "k": int(k), "hits_in_top_k": hits,
            "expected_hits": round(float(k * n_pos / len(score)), 2), "p_hits": float((1 + (null_hits >= hits).sum()) / (n_perm + 1)),
            "auc": round(a, 4), "p_auc": float((1 + (null_auc >= a).sum()) / (n_perm + 1)),
            "auc_null_95": [round(float(np.percentile(null_auc, 2.5)), 3), round(float(np.percentile(null_auc, 97.5)), 3)]}


def main():
    gt = ground_truth()
    gt.to_csv(RESULTS / "knockout_ground_truth.csv", index=False)
    A = pd.read_csv(RESULTS / "knockout_ranking.csv")
    B = pd.read_csv(RESULTS / "knockout_fba.csv", dtype={"gene_id": str})
    B["symbol"] = B["symbol"].astype(str).str.upper()
    ko_pos = set(gt[(gt["class"] == "KO/KD") & gt.raised_titer_or_qp].gene)
    ko_neg = set(gt[(gt["class"] == "KO/KD") & ~gt.raised_titer_or_qp].gene)
    ko_single = set(gt[(gt["class"] == "KO/KD") & gt.raised_titer_or_qp & gt.single_target_evidence].gene)
    oe_pos = set(gt[(gt["class"] == "OE") & gt.raised_titer_or_qp].gene)
    res = {"ground_truth": {"ko_kd_positive_genes": len(ko_pos), "ko_kd_positive_single_target": len(ko_single),
                            "ko_kd_no_benefit_genes": len(ko_neg), "oe_positive_genes": len(oe_pos)}}

    # ---------------- Route A ----------------
    sym = A.symbol.astype(str)
    a_ko, a_oe = sym.isin(ko_pos).to_numpy(), sym.isin(oe_pos).to_numpy()
    res["route_a_coverage"] = {"array_genes": int(len(A)), "ko_positives_on_array": int(a_ko.sum()), "ko_positives_total": len(ko_pos),
                               "oe_positives_on_array": int(a_oe.sum()), "oe_positives_total": len(oe_pos),
                               "ko_on_array": sorted(sym[a_ko]), "named_in_brief_on_array": [g for g in NAMED_IN_BRIEF if g in set(sym)],
                               "named_in_brief_not_on_array": [g for g in NAMED_IN_BRIEF if g not in set(sym)]}
    dw, dr = A.delta_within_experiment.to_numpy(), A.delta_raw.to_numpy()
    res["route_a_knockout_end_within_experiment"] = {k: enrichment(-dw, a_ko, k) for k in ("100", "50") for k in [int(k)]} if False else {
        "k100": enrichment(-dw, a_ko, 100), "k50": enrichment(-dw, a_ko, 50)}
    res["route_a_knockout_end_raw_delta"] = {"k100": enrichment(-dr, a_ko, 100)}
    res["route_a_knockout_end_single_target_only"] = {"k100": enrichment(-dw, sym.isin(ko_single).to_numpy(), 100)}
    res["route_a_overexpression_end_within_experiment"] = {"k100": enrichment(dw, a_oe, 100), "k50": enrichment(dw, a_oe, 50)}
    res["route_a_overexpression_end_raw_delta"] = {"k100": enrichment(dr, a_oe, 100)}
    cls = np.where(a_ko, "KO/KD raised titer or qP", np.where(sym.isin(ko_neg), "KO/KD no benefit",
                   np.where(a_oe, "overexpression raised titer or qP", "")))
    A["known_target"] = cls
    A["named_in_brief"] = sym.isin(NAMED_IN_BRIEF)
    A["marzluf_screen"] = sym.map({**MARZLUF_INCOMPATIBLE, **MARZLUF_OK}).fillna("")
    fold = gt.set_index(["gene", "class"]).max_fold
    A["max_reported_fold"] = [fold.get((s, "KO/KD"), fold.get((s, "OE"), np.nan)) if c else np.nan for s, c in zip(sym, cls)]
    A["ground_truth_source"] = np.where(cls != "", "Schulz et al. 2026 suppl. sheet B", "")
    A.to_csv(RESULTS / "knockout_ranking.csv", index=False)
    res["route_a_named_ranks"] = {r.symbol: {"knockout_rank": int(r.knockout_rank), "of": int(len(A)), "delta_within": float(r.delta_within_experiment),
                                             "class": r.known_target} for r in A[A.symbol.isin(ko_pos | ko_neg | set(NAMED_IN_BRIEF))].itertuples()}

    # ---------------- Route B ----------------
    bs = B.symbol
    b_ko = bs.isin(ko_pos).to_numpy()
    res["route_b_coverage"] = {"model_genes": int(len(B)), "ko_positives_in_model": int(b_ko.sum()), "ko_in_model": sorted(bs[b_ko]),
                               "named_in_brief_in_model": [g for g in NAMED_IN_BRIEF if g in set(bs)],
                               "named_in_brief_not_in_model": [g for g in NAMED_IN_BRIEF if g not in set(bs)]}
    res["route_b_tolerance"] = {"share_tolerated_ko_positives": float(B.tolerated[b_ko].mean()) if b_ko.sum() else None,
                                "share_tolerated_all_model_genes": float(B.tolerated.mean()),
                                "calls": {r.symbol: {"growth_ratio": float(r.growth_ratio), "igg_capacity_ratio": float(r.igg_capacity_ratio),
                                                     "tolerated": bool(r.tolerated), "tradeoff_score": float(r.tradeoff_score)}
                                          for r in B[bs.isin(ko_pos | ko_neg | set(NAMED_IN_BRIEF))].itertuples()}}
    res["route_b_tradeoff_score"] = enrichment(B.tradeoff_score.to_numpy(), b_ko, 50) if b_ko.sum() else None
    cand = B[bs.isin(MARZLUF_CANDIDATES)].copy()
    cand["named_incompatible"] = cand.symbol.isin(MARZLUF_INCOMPATIBLE)
    res["route_b_marzluf"] = {
        "candidates_total": len(MARZLUF_CANDIDATES), "candidates_in_model": int(cand.symbol.nunique()),
        "named_incompatible_in_model": sorted(set(cand.symbol[cand.named_incompatible])),
        "auc_growth_ratio_low_predicts_incompatible": round(auc_of(-cand.growth_ratio.to_numpy(), cand.named_incompatible.to_numpy()), 3)
        if 0 < cand.named_incompatible.sum() < len(cand) else None,
        "calls": {r.symbol: {"growth_ratio": float(r.growth_ratio), "tolerated": bool(r.tolerated),
                             "reported": MARZLUF_INCOMPATIBLE.get(r.symbol, MARZLUF_OK.get(r.symbol, "not individually reported"))}
                  for r in cand.sort_values("growth_ratio").itertuples()},
        "caveat": "the 36 retained knockouts are not published by name; unnamed candidates are a mix of retained and rejected"}
    B["known_target"] = np.where(b_ko, "KO/KD raised titer or qP", np.where(bs.isin(ko_neg), "KO/KD no benefit", ""))
    B["marzluf_screen"] = bs.map({**MARZLUF_INCOMPATIBLE, **MARZLUF_OK}).fillna("")
    B.to_csv(RESULTS / "knockout_fba.csv", index=False)
    (RESULTS / "knockout_eval.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
