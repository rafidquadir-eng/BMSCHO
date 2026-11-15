"""P1: pathway-recovery test on Yao et al. 2021 (Metabolites 11:823) Table S3.

Only summary statistics are public (per-metabolite Pearson r with qP across 12
clones). No per-clone data exist publicly, so this is a concordance check of the
delta-style pathway ranking against the authors' own conclusion, pre-registered
in DECISIONS.md. Modality: extracellular metabolomics (untargeted LC-MS, day 7).
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw" / "yao"
RESULTS = ROOT / "results"
SEED = 20261116
SUPP_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8704136/supplementaryFiles"
TARGETS = {"map00250": "Alanine, aspartate and glutamate metabolism", "map00020": "Citrate cycle (TCA cycle)"}

# Names in Table S3 that do not match a KEGG synonym exactly. Resolved by hand.
MANUAL = {
    "alpha-D-Glucose": "C00267", "D-Gluconic acid": "C00257", "L-methionine-S-oxide": "C02989",
    "sn-Glycero-3-phospho-1-inositol": "C01225", "1D-myo-Inositol 1,3,4,5,6-pentakisphosphate": "C01284",
    "(S)-3-Methyl-2-oxopentanoic acid": "C00671", "D-myo-Inositol 1,2-cyclic phosphate": "C04299",
    "1-(5-Phospho-D-ribosyl)-5-amino-4-imidazolecarboxylate": "C04751",
    "7,8-Dihydroneopterin 3'-triphosphate": "C04895", "5-Acetylamino-6-formylamino-3-methyluracil": "C16365",
    "3-Methoxy-4-hydroxyphenylglycolaldehyde": "C05583", "2-Hydroxy-dATP": "C19969", "2-Hydroxy-dAMP": "C19970",
    "Pterine": "C00715",
}


def _get(url, dest: Path, binary=False):
    import requests

    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=300, headers={"User-Agent": "research-client"})
            r.raise_for_status()
            dest.write_bytes(r.content)
            return dest
        except Exception as exc:
            if attempt == 3:
                raise
            print(f"  retry {url}: {exc}")
            time.sleep(2 ** (attempt + 1))


def fetch():
    import io
    import zipfile

    z = _get(SUPP_URL, RAW / "supp.zip")
    pdf = RAW / "supplementary.pdf"
    if not pdf.exists():
        outer = zipfile.ZipFile(z)
        inner = zipfile.ZipFile(io.BytesIO(outer.read("metabolites-11-00823-s001.zip")))
        name = [n for n in inner.namelist() if n.endswith(".pdf")][0]
        pdf.write_bytes(inner.read(name))
    _get("https://rest.kegg.jp/list/compound", RAW / "kegg_compound.txt")
    _get("https://rest.kegg.jp/link/pathway/compound", RAW / "kegg_compound_pathway.txt")
    _get("https://rest.kegg.jp/list/pathway", RAW / "kegg_pathway.txt")
    return pdf


def parse_table_s3(pdf: Path) -> pd.DataFrame:
    from pypdf import PdfReader

    r = PdfReader(str(pdf))
    text = "\n".join((r.pages[i].extract_text() or "") for i in (6, 7, 8))
    text = text[text.index("p-value"):]
    toks = text.split()
    num = re.compile(r"^-?\d\.\d{3}$")
    rows, name = [], []
    i = 0
    while i < len(toks):
        if i + 3 < len(toks) and all(num.match(t) for t in toks[i:i + 4]):
            nm = " ".join(name)
            nm = re.sub(r"^(Correlation|coefficient|p-value|with|qP|growth|Annotated|metabolite|\s)+", "", nm).strip()
            nm = re.sub(r"-\s+", "-", nm)  # line-wrapped hyphenated names
            rows.append((nm, *map(float, toks[i:i + 4])))
            name, i = [], i + 4
        else:
            name.append(toks[i])
            i += 1
    df = pd.DataFrame(rows, columns=["metabolite", "r_qp", "p_qp", "r_growth", "p_growth"])
    df = df[df.metabolite.str.len() > 0].reset_index(drop=True)
    return df


def kegg_tables():
    syn = {}
    for line in (RAW / "kegg_compound.txt").read_text().splitlines():
        cid, names = line.split("\t", 1)
        cid = cid.replace("cpd:", "")
        for n in names.split("; "):
            syn.setdefault(n.strip().lower(), cid)
    c2p = {}
    for line in (RAW / "kegg_compound_pathway.txt").read_text().splitlines():
        c, p = line.split("\t")
        c, p = c.replace("cpd:", ""), p.replace("path:", "")
        n = int(p[3:])
        if p.startswith("map") and 10 <= n <= 999:  # metabolic maps only; excludes map011xx / map012xx overviews
            c2p.setdefault(c, set()).add(p)
    pname = {}
    for line in (RAW / "kegg_pathway.txt").read_text().splitlines():
        p, n = line.split("\t")
        pname[p.replace("path:", "")] = n
    return syn, c2p, pname


def map_metabolites(df, syn):
    ids, how = [], []
    for m in df.metabolite:
        if m in MANUAL:
            ids.append(MANUAL[m]); how.append("manual" if MANUAL[m] else "unresolved")
            continue
        cands = [m, m.replace("L-Cysteate", "L-Cysteic acid"), re.sub(r" acid$", "ate", m), re.sub(r"ate$", "ic acid", m)]
        hit = next((syn[c.lower()] for c in cands if c.lower() in syn), None)
        ids.append(hit); how.append("exact" if hit and m.lower() in syn else ("variant" if hit else "unresolved"))
    return df.assign(kegg=ids, match=how)


def pathway_rank(df, c2p, pname, n_perm=10000, min_members=3, seed=SEED):
    d = df.dropna(subset=["kegg"]).drop_duplicates("kegg").reset_index(drop=True)
    z = np.arctanh(np.clip(d.r_qp.to_numpy(), -0.999, 0.999))
    sig = (d.p_qp.to_numpy() < 0.05)
    paths = sorted({p for c in d.kegg for p in c2p.get(c, ())})
    M = np.array([[p in c2p.get(c, ()) for c in d.kegg] for p in paths])
    keep = M.sum(axis=1) >= min_members
    paths, M = [p for p, k in zip(paths, keep) if k], M[keep]

    def stat(zv):
        return np.array([zv[m].mean() - zv[~m].mean() for m in M])

    obs = stat(z)
    rng = np.random.default_rng(seed)
    perm = np.empty((n_perm, len(paths)))
    for b in range(n_perm):
        perm[b] = stat(rng.permutation(z))
    p_one = (1 + (perm >= obs).sum(axis=0)) / (n_perm + 1)
    from scipy.stats import hypergeom

    p_ora = np.array([hypergeom.sf(sig[m].sum() - 1, len(d), sig.sum(), m.sum()) for m in M])
    out = pd.DataFrame({"pathway": paths, "name": [pname.get(p, p) for p in paths], "n_measured": M.sum(axis=1),
                        "n_sig_qp": [int(sig[m].sum()) for m in M],
                        "mean_r_qp": [float(d.r_qp.to_numpy()[m].mean()) for m in M],
                        "delta_z": obs, "p_perm": p_one, "p_ora": p_ora,
                        "members": ["; ".join(d.metabolite[m]) for m in M]})
    out = out.sort_values("delta_z", ascending=False).reset_index(drop=True)
    out.insert(0, "rank_delta", np.arange(1, len(out) + 1))
    out["rank_ora"] = out.p_ora.rank(method="min").astype(int)
    out["target"] = out.pathway.isin(TARGETS)
    # joint test: mean rank of the two targets under the permutation null
    idx = [paths.index(t) for t in TARGETS if t in paths]
    joint = {}
    if len(idx) == len(TARGETS):
        def mean_rank(v):
            r = (-v).argsort().argsort() + 1
            return r[idx].mean()
        o = mean_rank(obs)
        null = np.array([mean_rank(perm[b]) for b in range(n_perm)])
        joint = {"mean_rank_targets": float(o), "p_joint": float((1 + (null <= o).sum()) / (n_perm + 1)),
                 "expected_mean_rank_under_null": float(null.mean())}
    info = {"metabolites_in_table_s3": int(len(df)), "mapped_to_kegg": int(df.kegg.notna().sum()),
            "unique_kegg": int(len(d)), "eligible_pathways": int(len(paths)), "n_sig_qp": int(sig.sum()),
            "frac_positive_r": float((d.r_qp > 0).mean()), **joint}
    return out, info, d


def main():
    pdf = fetch()
    df = parse_table_s3(pdf)
    syn, c2p, pname = kegg_tables()
    df = map_metabolites(df, syn)
    df.to_csv(RESULTS / "yao_metabolite_map.csv", index=False)
    out, info, d = pathway_rank(df, c2p, pname)
    num = ["mean_r_qp", "delta_z", "p_perm", "p_ora"]
    out[num] = out[num].round(4)
    out.to_csv(RESULTS / "yao_pathway_ranking.csv", index=False)
    (RESULTS / "yao_info.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))
    print(out[["rank_delta", "name", "n_measured", "n_sig_qp", "mean_r_qp", "delta_z", "p_perm", "rank_ora", "p_ora"]].to_string(index=False))
    # Post hoc sensitivity (NOT pre-registered): only pathways annotated in Chinese hamster.
    _get("https://rest.kegg.jp/list/pathway/cge", RAW / "kegg_pathway_cge.txt")
    cge = {line.split("\t")[0].replace("cge", "map") for line in (RAW / "kegg_pathway_cge.txt").read_text().splitlines()}
    out_h, info_h, _ = pathway_rank(df, {c: {p for p in ps if p in cge} for c, ps in c2p.items()}, pname)
    out_h[num] = out_h[num].round(4)
    out_h.to_csv(RESULTS / "yao_pathway_ranking_hamster_only.csv", index=False)
    (RESULTS / "yao_info_hamster_only.json").write_text(json.dumps(info_h, indent=2))
    return out, info


if __name__ == "__main__":
    if "--parse-only" in sys.argv:
        d = map_metabolites(parse_table_s3(fetch()), kegg_tables()[0])
        pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
        print(d.to_string()); print(d.match.value_counts())
    else:
        main()
