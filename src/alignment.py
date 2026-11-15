"""P4: does a shared pathway basis align independent modalities?

No sample here has more than one modality. Nothing in this file is a fused or
multi-omic prediction. Each statistic is computed inside one modality; the test
is whether pathway-level contrast vectors from different modalities correlate.
Design pre-registered in DECISIONS.md (commit 83b8d9e) before computing.
"""
from __future__ import annotations

import json
import re
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
URLS = {
    "prot/AllQuantifiedProteinGroups.tsv": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/08/PXD053511/AllQuantifiedProteinGroups.tsv",
    "prot/uniprot_10029.tsv": "https://rest.uniprot.org/uniprotkb/stream?query=%28taxonomy_id%3A10029%29&fields=accession%2Cgene_primary%2Cgene_names%2Cxref_geneid%2Cprotein_name&format=tsv",
    "prot/uniprot_reviewed_mouse_human.tsv": "https://rest.uniprot.org/uniprotkb/stream?query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%20OR%20organism_id%3A10090%29&fields=accession%2Cgene_primary%2Cprotein_name%2Corganism_id&format=tsv",
    "kegg/cge_genes.txt": "https://rest.kegg.jp/list/cge",
    "kegg/cge_gene_pathway.txt": "https://rest.kegg.jp/link/pathway/cge",
    "kegg/cge_pathways.txt": "https://rest.kegg.jp/list/pathway/cge",
    "icho/iCHOv1.json": "https://bigg.ucsd.edu/static/models/iCHOv1.json",
}
CURRENCY = {"h", "h2o", "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph", "pi", "ppi", "co2", "o2", "coa", "fad", "fadh2",
            "h2o2", "nh4", "na1", "k", "hco3", "gtp", "gdp", "udp", "utp", "ctp", "cdp", "so4"}
SKIP_SUBSYSTEMS = re.compile(r"TRANSPORT|EXCHANGE|UNASSIGNED|BIOMASS|^$|MISCELLANEOUS", re.I)


def fetch():
    import yao

    for rel, url in URLS.items():
        yao._get(url, RAW / rel)
    yao.fetch()


# ----------------------------- gene-level statistics ----------------------------- #
def transcript_delta(endpoint, within=True, only_37=False):
    """Gene-level delta on GSE30321 (transcriptomics) for one endpoint."""
    D = cp.prepare()
    meta = D["meta"]
    y = meta[endpoint].to_numpy(dtype=float)
    ok = ~np.isnan(y) & ~meta["qc_label_excluded"].to_numpy()
    if only_37:
        ok &= (meta["temperature"] == 37).to_numpy()
    X, y, exps = D["X"][ok], y[ok], meta["experiment_id"].to_numpy()[ok]
    if within:
        multi = pd.Series(exps).map(pd.Series(exps).value_counts()).to_numpy() >= 4
        X, y, exps = X[multi], y[multi], exps[multi]
        X = cp.center_within(X, exps)
        y = cp.center_within(y[:, None], exps).ravel()
    Z = (X - X.mean(axis=0)) / np.where(X.std(axis=0) > 1e-8, X.std(axis=0), 1.0)
    lo, hi = np.quantile(y, [0.25, 0.75])
    d = Z[y >= hi].mean(axis=0) - Z[y <= lo].mean(axis=0)
    sym = pd.Series([s if isinstance(s, str) else None for s in D["symbols"]])
    g = pd.DataFrame({"symbol": sym, "stat": d}).dropna().groupby("symbol").stat.mean()
    return g, {"n_arrays": int(len(y)), "n_experiments": int(len(np.unique(exps)))}


def protein_names_map():
    def names(s):
        s = str(s)
        out = [re.split(r" \(", s)[0]] + re.findall(r"\(([^()]+)\)", s)
        return [n.strip().lower() for n in out if len(n.strip()) > 3 and not n.startswith("EC ")]

    m = {}
    for f in ["prot/uniprot_reviewed_mouse_human.tsv", "prot/uniprot_10029.tsv"]:  # hamster last, so it wins
        u = pd.read_csv(RAW / f, sep="\t", dtype=str)
        gcol = "Gene Names (primary)"
        for g, pn in zip(u[gcol], u["Protein names"]):
            if isinstance(g, str) and g:
                for n in names(pn):
                    m[n] = g.split(";")[0].strip().upper()
    return m


def proteomics_stat():
    """PXD053511: Welch t, day 4 (exponential) minus day 7 (stationary), per gene."""
    from scipy.stats import ttest_ind

    p = pd.read_csv(RAW / "prot/AllQuantifiedProteinGroups.tsv", sep="\t", index_col=False, low_memory=False)
    p = p[p["Protein Decoy/Contaminant/Target"] == "T"]
    d4 = [c for c in p.columns if c.startswith("Intensity_d4_")]
    d7 = [c for c in p.columns if c.startswith("Intensity_d7_")]
    I = p[d4 + d7].apply(pd.to_numeric, errors="coerce").replace(0, np.nan)
    L = np.log2(I)
    L = L - L.median(axis=0)
    keep = (L[d4].notna().sum(axis=1) >= 3) & (L[d7].notna().sum(axis=1) >= 3)
    t = ttest_ind(L.loc[keep, d4], L.loc[keep, d7], axis=1, equal_var=False, nan_policy="omit").statistic
    lfc = L.loc[keep, d4].mean(axis=1) - L.loc[keep, d7].mean(axis=1)
    q = p.loc[keep, ["Protein Accession", "Gene", "Protein Full Name"]].copy()
    q["t"], q["log2fc"] = np.asarray(t, dtype=float), lfc.to_numpy()
    nm = protein_names_map()
    locus = re.compile(r"^(I79_|H671_|CgPICR_|LOC\d)")

    def symbol(gene, pname):
        for g in str(gene).split("|"):
            g = g.strip()
            if g and g != "nan" and not locus.match(g):
                return g.upper(), "deposited gene column"
        n = str(pname).strip().lower()
        if n in nm:
            return nm[n], "protein name match"
        return None, "unmapped"

    sm = [symbol(g, n) for g, n in zip(q["Gene"], q["Protein Full Name"])]
    q["symbol"], q["mapping"] = [s[0] for s in sm], [s[1] for s in sm]
    info = {"protein_groups": int(len(p)), "quantified_3_of_4_per_group": int(keep.sum()),
            "mapping": q.mapping.value_counts().to_dict(), "mapped_fraction": float(q.symbol.notna().mean()),
            "unique_genes": int(q.symbol.nunique()), "n_up_exponential_t_gt_2": int((q.t > 2).sum()),
            "n_up_stationary_t_lt_-2": int((q.t < -2).sum())}
    g = q.dropna(subset=["symbol", "t"]).groupby("symbol").t.mean()
    q.to_csv(RESULTS / "proteomics_gene_stats.csv", index=False)
    return g, info


# ----------------------------- pathway bases ----------------------------- #
def kegg_symbols():
    m = {}
    for line in (RAW / "kegg/cge_genes.txt").read_text().splitlines():
        p = line.split("\t")
        if ";" in p[-1]:
            m[p[0].replace("cge:", "")] = p[-1].split(";")[0].split(",")[0].strip().upper()
    return m


def icho_basis():
    """iCHO1766 subsystems -> gene symbols and -> KEGG compounds."""
    m = json.loads((RAW / "icho/iCHOv1.json").read_text())
    ks = kegg_symbols()
    gsym = {}
    for g in m["genes"]:
        n = str(g.get("name") or "").upper()
        gsym[g["id"]] = ks.get(g["id"], n) if (not n or n.startswith("LOC")) else n
    met_kegg = {x["id"]: x.get("annotation", {}).get("kegg.compound", []) for x in m["metabolites"]}
    genes, cpds = {}, {}
    for r in m["reactions"]:
        s = r.get("subsystem", "")
        if SKIP_SUBSYSTEMS.search(s):
            continue
        for gid in re.findall(r"\d+", r.get("gene_reaction_rule", "")):
            if gsym.get(gid):
                genes.setdefault(s, set()).add(gsym[gid])
        for mid in r["metabolites"]:
            if mid.rsplit("_", 1)[0] in CURRENCY:
                continue
            for c in met_kegg.get(mid, []):
                cpds.setdefault(s, set()).add(c)
    return genes, cpds


def kegg_basis():
    ks = kegg_symbols()
    pname = {}
    for line in (RAW / "kegg/cge_pathways.txt").read_text().splitlines():
        p, n = line.split("\t")
        pname[p.replace("path:", "")] = n.replace(" - Cricetulus griseus (Chinese hamster)", "")
    genes = {}
    for line in (RAW / "kegg/cge_gene_pathway.txt").read_text().splitlines():
        g, p = line.split("\t")
        g, p = g.replace("cge:", ""), p.replace("path:", "")
        n = int(p[3:])
        if g in ks and not (1100 <= n <= 1299):  # drop global/overview maps
            genes.setdefault(pname.get(p, p), set()).add(ks[g])
    cpds = {}
    num2name = {k[3:]: v for k, v in pname.items()}
    for line in (RAW / "yao/kegg_compound_pathway.txt").read_text().splitlines():
        c, p = line.split("\t")
        num = p.replace("path:map", "")
        if num in num2name and not (1100 <= int(num) <= 1299):
            cpds.setdefault(num2name[num], set()).add(c.replace("cpd:", ""))
    return genes, cpds


# ----------------------------- the test ----------------------------- #
def zscore(s):
    return (s - s.mean()) / s.std()


def pathway_scores(stat: pd.Series, sets: dict, min_n: int):
    z = zscore(stat)
    idx = {k: i for i, k in enumerate(z.index)}
    out, members = {}, {}
    for name, feats in sets.items():
        ii = [idx[f] for f in feats if f in idx]
        if len(ii) >= min_n:
            out[name], members[name] = float(z.to_numpy()[ii].mean()), np.array(ii)
    return pd.Series(out), members, z


def align(stat_a, sets_a, min_a, stat_b, sets_b, min_b, n_perm=10000, seed=SEED):
    """Spearman between pathway vectors; null = permute feature labels of modality B."""
    from scipy.stats import spearmanr

    pa, _, _ = pathway_scores(stat_a, sets_a, min_a)
    pb, mb, zb = pathway_scores(stat_b, sets_b, min_b)
    common = sorted(set(pa.index) & set(pb.index))
    if len(common) < 5:
        return {"n_pathways": len(common), "spearman": None, "p_perm": None}, pd.DataFrame()
    a = pa[common].to_numpy()
    rho = float(spearmanr(a, pb[common].to_numpy()).correlation)
    rng = np.random.default_rng(seed)
    zv = zb.to_numpy()
    null = np.empty(n_perm)
    ra = pd.Series(a).rank().to_numpy()
    for b in range(n_perm):
        zp = rng.permutation(zv)
        pbp = np.array([zp[mb[c]].mean() for c in common])
        null[b] = np.corrcoef(ra, pd.Series(pbp).rank().to_numpy())[0, 1]
    tab = pd.DataFrame({"pathway": common, "score_a": a, "score_b": pb[common].to_numpy(),
                        "n_features_b": [len(mb[c]) for c in common]})
    return {"n_pathways": len(common), "spearman": round(rho, 4), "p_perm": float((1 + (null >= rho).sum()) / (n_perm + 1)),
            "null_95": [round(float(np.percentile(null, 2.5)), 3), round(float(np.percentile(null, 97.5)), 3)]}, tab


def gene_level(stat_a, stat_b, n_perm=10000, seed=SEED):
    from scipy.stats import spearmanr

    common = sorted(set(stat_a.index) & set(stat_b.index))
    a, b = stat_a[common].to_numpy(), stat_b[common].to_numpy()
    rho = float(spearmanr(a, b).correlation)
    rng = np.random.default_rng(seed)
    ra, rb = pd.Series(a).rank().to_numpy(), pd.Series(b).rank().to_numpy()
    null = np.array([np.corrcoef(ra, rng.permutation(rb))[0, 1] for _ in range(n_perm)])
    return {"n_genes": len(common), "spearman": round(rho, 4), "p_perm": float((1 + (null >= rho).sum()) / (n_perm + 1))}


def main():
    fetch()
    import yao

    res, tabs = {}, []
    ig, ic = icho_basis()
    kg, kc = kegg_basis()
    bases = {"iCHO1766 subsystems": (ig, ic), "KEGG pathways (Chinese hamster)": (kg, kc)}
    res["bases"] = {k: {"gene_sets": len(v[0]), "compound_sets": len(v[1])} for k, v in bases.items()}

    # ---- leg 1: growth phase, transcriptomics vs proteomics ----
    prot, pinfo = proteomics_stat()
    res["proteomics"] = pinfo
    variants = {"all temperatures, within experiment (pre-registered primary)": dict(within=True, only_37=False),
                "37 C arrays only, within experiment (sensitivity)": dict(within=True, only_37=True)}
    res["leg1_growth_phase_transcriptomics_vs_proteomics"] = {}
    for vname, kw in variants.items():
        tx, tinfo = transcript_delta("growth rate", **kw)
        block = {"transcriptomics": tinfo, "gene_level": gene_level(tx, prot)}
        for bname, (gs, _) in bases.items():
            r, tab = align(tx, gs, 5, prot, gs, 5)
            block[bname] = r
            if len(tab):
                tabs.append(tab.assign(leg="growth phase: transcriptomics (a) vs proteomics (b)", basis=bname, variant=vname))
        res["leg1_growth_phase_transcriptomics_vs_proteomics"][vname] = block

    # ---- leg 2: Qp association, transcriptomics vs metabolomics ----
    txq, qinfo = transcript_delta("qproductivity", within=False)
    s3 = yao.map_metabolites(yao.parse_table_s3(yao.fetch()), yao.kegg_tables()[0]).dropna(subset=["kegg"]).drop_duplicates("kegg")
    met = pd.Series(np.arctanh(np.clip(s3.r_qp.to_numpy(), -0.999, 0.999)), index=s3.kegg.to_numpy())
    res["leg2_qp_transcriptomics_vs_metabolomics"] = {"transcriptomics": qinfo, "metabolites": int(len(met))}
    for bname, (gs, cs) in bases.items():
        r, tab = align(txq, gs, 5, met, cs, 3)
        res["leg2_qp_transcriptomics_vs_metabolomics"][bname] = r
        if len(tab):
            tabs.append(tab.assign(leg="Qp association: transcriptomics (a) vs metabolomics (b)", basis=bname, variant="raw delta"))

    (RESULTS / "alignment.json").write_text(json.dumps(res, indent=2))
    pd.concat(tabs).round(4).to_csv(RESULTS / "alignment_pathway_scores.csv", index=False)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
