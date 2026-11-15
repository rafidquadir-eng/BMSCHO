#!/usr/bin/env python
"""CHO proof of concept: predict cell-specific productivity (Qp) from
transcriptomics on GSE30321 with leakage-safe grouped cross-validation.

One command reproduces everything:

    python cho_poc.py --all

Stages (each can be run alone):
    --self-test   offline synthetic checks of the machinery (~1 min)
    --download    fetch GEO / ArrayExpress / MSigDB files into raw/, write data/PROVENANCE.md
    --run         grouped-CV metrics, matched baseline, leakage comparison, delta decode
    --figures     slide figures from results/
    --pretrain    self-supervised pretraining learning curves (src/pretrain.py, needs torch)

Scope boundary: single modality (2011 microarray transcriptomics). Nothing here
validates multi-omic fusion, latent dynamics or scale-up prediction.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261116
ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

BASELINE_RMSE = 4.44  # Clarke et al. 2011 J Biotechnol 151:159, RMSE(CMV), pg/cell/day
ENDPOINT = "qproductivity"

SOURCES = {
    "series_matrix": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30321/matrix/GSE30321_series_matrix.txt.gz",
        "path": "geo/GSE30321_series_matrix.txt.gz",
        "sha256": "790c5ef6b48070b0a5af3ea9596466c78fc29d25b1ea20e0ecac37254b465cdf",
    },
    "platform": {
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13791&targ=self&form=text&view=data",
        "path": "gpl/GPL13791.txt",
        "sha256": "f102b8d336e68a88646306659526c7452cabdaefa2f25e0cef33cb961447208f",
    },
    "sdrf": {
        "url": "https://www.ebi.ac.uk/biostudies/files/E-GEOD-30321/E-GEOD-30321.sdrf.txt",
        "path": "ae/E-GEOD-30321.sdrf.txt",
        "sha256": "69ce5da3854675369ca2784a611ad9cd8b570a6dd3610a71bef87473ce0818ee",
    },
    "hallmark": {
        "url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
        "path": "genesets/h.all.v2023.2.Hs.symbols.gmt",
        "sha256": "8ff1f036f0988c99d42cf9e5ea24d355f10b8ece5e0d4f6edf26a0458f480780",
    },
    "go_bp": {
        "url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.symbols.gmt",
        "path": "genesets/c5.go.bp.v2023.2.Hs.symbols.gmt",
        "sha256": "b4b3b49a267c2239fad40d39a56f2146a57743992943915e80ebed34ec7ebaae",
    },
}

# Themes reported by the companion WGCNA paper (Clarke et al. 2011, J Biotechnol
# 155:350): cell cycle, protein secretion, vesicle transport. Fixed here, before
# looking at any delta-decode output, so the check is confirmatory.
THEME_SETS = {
    "HALLMARK_E2F_TARGETS": "cell cycle",
    "HALLMARK_G2M_CHECKPOINT": "cell cycle",
    "HALLMARK_MITOTIC_SPINDLE": "cell cycle",
    "HALLMARK_PROTEIN_SECRETION": "protein secretion / vesicle transport",
    "GOBP_CELL_CYCLE": "cell cycle",
    # GO:0009306 "protein secretion" is not a gene set in MSigDB 2023.2; its parent is used
    "GOBP_PROTEIN_LOCALIZATION_TO_EXTRACELLULAR_REGION": "protein secretion",
    "GOBP_VESICLE_MEDIATED_TRANSPORT": "vesicle transport",
    "GOBP_GOLGI_VESICLE_TRANSPORT": "vesicle transport",
}
GO_THEME_NAMES = [k for k in THEME_SETS if k.startswith("GOBP_")]


# --------------------------------------------------------------------------- #
# Data acquisition and parsing
# --------------------------------------------------------------------------- #
def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_all(force: bool = False) -> None:
    """Fetch every source file into raw/ and write data/PROVENANCE.md."""
    import requests

    DATA.mkdir(exist_ok=True)
    rows = []
    for name, src in SOURCES.items():
        dest = RAW / src["path"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        stamp = None
        if force or not dest.exists():
            for attempt in range(4):
                try:
                    r = requests.get(src["url"], timeout=300)
                    r.raise_for_status()
                    dest.write_bytes(r.content)
                    break
                except Exception as exc:  # network retry with backoff
                    if attempt == 3:
                        raise
                    print(f"  retry {name}: {exc}")
                    time.sleep(2 ** (attempt + 1))
            stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        else:
            stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(dest.stat().st_mtime))
        sha = sha256_of(dest)
        status = "matches pinned hash" if sha == src["sha256"] else "DIFFERS from pinned hash"
        if sha != src["sha256"]:
            # GEO regenerates series-matrix files from time to time; warn loudly, do not abort.
            print(f"  WARNING {name}: SHA256 {sha} differs from the pinned {src['sha256']}; "
                  "results may not reproduce bit-for-bit")
        rows.append((name, src["url"], src["path"], dest.stat().st_size, sha, stamp, status))
        print(f"  {name}: {dest.stat().st_size} bytes, {status}")
    write_provenance(rows)


def write_provenance(rows) -> None:
    expr, meta = load_dataset()
    gpl = parse_platform(RAW / SOURCES["platform"]["path"])
    hall = read_gmt(RAW / SOURCES["hallmark"]["path"])
    lines = [
        "# Data provenance",
        "",
        "Raw files live in `raw/` (git-ignored) and are re-fetched by `python cho_poc.py --download`.",
        "No raw data is committed. This file is regenerated on every download.",
        "",
        "| name | URL | local path (under raw/) | bytes | SHA256 | downloaded (UTC) | check |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, url, path, size, sha, stamp, status in rows:
        lines.append(f"| {name} | {url} | `{path}` | {size} | `{sha}` | {stamp} | {status} |")
    n_lab = int(meta[ENDPOINT].notna().sum())
    lines += [
        "",
        "## Contents",
        "",
        f"- Expression matrix: {expr.shape[0]} probesets x {expr.shape[1]} arrays "
        "(log2 RMA, as deposited by the submitters; no re-normalisation here).",
        f"- Platform GPL13791 (WyeHamster2a custom Affymetrix CHO array): {len(gpl)} probesets.",
        f"- Sample metadata: {meta.shape[0]} arrays, bioprocess variables parsed from "
        "`!Sample_characteristics_ch1`.",
        f"- Arrays with a Qp label (`qproductivity`): {n_lab} of {meta.shape[0]}.",
        f"- ArrayExpress SDRF (E-GEOD-30321): carries the same nine characteristics, no extra "
        "columns (no culture, cell-line or time-point identifier).",
        f"- MSigDB Hallmark v2023.2.Hs: {len(hall)} gene sets.",
        "- MSigDB C5 GO BP v2023.2.Hs: only the four pre-specified theme sets are used "
        f"({', '.join(GO_THEME_NAMES)}).",
        "- GEO `/suppl/` for GSE30321 does not exist (HTTP 404); per-sample supplementary "
        "file is `NONE` for all 295 arrays, so no CEL files are public.",
        "",
        "## Label and grouping provenance",
        "",
        "- Qp is taken verbatim from GEO sample characteristics (`qproductivity`, pg/cell/day).",
        "- There is no culture identifier in GEO, the SOFT file or the ArrayExpress SDRF. "
        "The culture and experiment IDs are derived (see `derive_groups` in `cho_poc.py` "
        "and `results/culture_map.csv`).",
    ]
    (DATA / "PROVENANCE.md").write_text("\n".join(lines) + "\n")


def parse_series_matrix(path_or_text) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse a GEO series matrix. Returns (expr: probes x samples, meta: samples x fields).

    Sample characteristics of the form "key: value" become one column per key.
    "NA" becomes NaN; numeric columns are converted.
    """
    if isinstance(path_or_text, Path) or "\n" not in str(path_or_text):
        p = Path(path_or_text)
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt", encoding="utf-8", errors="replace") as f:
            text = f.read()
    else:
        text = str(path_or_text)
    meta: dict[str, list] = {}
    table_lines: list[str] = []
    in_table = False
    for line in text.split("\n"):
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            in_table = False
            continue
        if in_table:
            table_lines.append(line)
            continue
        if not line.startswith("!Sample_"):
            continue
        parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
        key, vals = parts[0][len("!Sample_"):], parts[1:]
        if key == "characteristics_ch1":
            keys = {v.split(":", 1)[0].strip() for v in vals if ":" in v}
            if len(keys) != 1:
                raise ValueError(f"mixed characteristic keys on one line: {keys}")
            k = keys.pop()
            meta[k] = [v.split(":", 1)[1].strip() if ":" in v else None for v in vals]
        elif key in ("geo_accession", "title"):
            meta[key] = vals
    expr = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    m = pd.DataFrame(meta).set_index("geo_accession")
    for c in m.columns:
        if c == "title":
            continue
        m[c] = pd.to_numeric(m[c].replace({"NA": None, "": None}), errors="coerce")
    if list(expr.columns) != list(m.index):
        raise ValueError("expression columns and sample metadata are not aligned")
    return expr, m


def parse_platform(path: Path) -> pd.DataFrame:
    """GPL full-table text -> DataFrame indexed by probe ID with an upper-cased `symbol`."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").split("\n")
    b = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
    e = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    df = pd.read_csv(io.StringIO("\n".join(lines[b + 1:e])), sep="\t", keep_default_na=False, dtype=str)
    df = df.set_index("ID")
    sym = df["Gene Symbol"].str.strip()
    bad = sym.isin(["", "NA", "---"]) | sym.str.contains(r"[ /]") | sym.str.match(r"^[A-Za-z]{2}\.\d+$")
    df["symbol"] = sym.str.upper().where(~bad, None)
    return df


def read_gmt(path: Path, keep: list[str] | None = None) -> dict[str, list[str]]:
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if keep is not None and parts[0] not in keep:
                continue
            out[parts[0]] = [g.upper() for g in parts[2:] if g]
    return out


def load_dataset() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, meta = parse_series_matrix(RAW / SOURCES["series_matrix"]["path"])
    meta["idx"] = meta["title"].str.extract(r"_(\d+)$").astype(int)
    return expr, meta


# --------------------------------------------------------------------------- #
# Culture / experiment grouping (derived; GEO has no culture ID)
# --------------------------------------------------------------------------- #
def _r(a: int, b: int) -> list[int]:
    return list(range(a, b + 1))


# Each entry: (experiment_id, rule, sample indices [1-based index from the GEO
# title "..._<idx>"], ambiguous?). Rules:
#   seed   culture = distinct seeding density within the experiment
#   pos:k  arrays are deposited as consecutive time points of k replicate
#          cultures in fixed order; culture = position modulo k
#   single every array is its own culture (replicate cultures, one time point)
#   lump   all arrays form one group (cannot be resolved; merging is the
#          leakage-safe direction)
#   dict   explicit culture -> indices
# "ambiguous" marks blocks where the public metadata cannot settle the
# assignment. Every ambiguous block is unlabeled for Qp, so the ambiguity does
# not touch any reported predictive metric.
CULTURE_RULES = [
    ("E01", "seed", _r(1, 2), False),
    ("E02", "single", _r(3, 6), False),
    ("E03", "seed", _r(7, 15), False),
    ("E04", "seed", _r(16, 24), False),
    ("E05", "pos:4", _r(25, 32), False),
    ("E06", "seed", _r(33, 44), False),
    ("E07", "single", _r(45, 48), True),
    ("E08", "pos:4", _r(49, 56), False),
    ("E09", "single", _r(57, 59), False),
    ("E10", "single", _r(60, 63), False),
    ("E11", "seed", _r(64, 65), False),
    ("E12", {"a": [66, 69, 73], "b": [67, 70, 74], "c": [71, 75], "d": [68, 72, 76]}, _r(66, 76), False),
    ("E13", "lump", _r(77, 83), True),
    ("E13", "pos:3", _r(84, 94), True),
    ("E14", "pos:4", _r(95, 106), True),
    ("E15", "seed", _r(107, 109), False),
    ("E16", "pos:3", _r(110, 118), False),
    ("E17", "pos:3", _r(119, 127), False),
    ("E18", {"a": [128, 131, 133], "b": [129, 132, 134], "c": [130, 135]}, _r(128, 135), False),
    ("E19", "seed", _r(136, 138), False),
    ("E20", {"a": _r(139, 142), "b": _r(143, 146), "c": [147], "d": [148], "e": [149]}, _r(139, 149), True),
    ("E21", "seed", _r(150, 161), False),
    ("E23", "pos:3", _r(163, 171), False),
    ("E24", {"a": [172, 174, 176, 178, 180, 183], "b": [173, 175, 177, 179, 181, 184], "u": [182],
             "c": [209, 213, 217, 221]}, _r(172, 184) + [209, 213, 217, 221], True),
    ("E25", "seed", _r(185, 193), False),
    ("E26", "single", [194], True),
    ("E27", "seed", _r(195, 203), False),
    ("E28", "single", _r(204, 208), True),
    ("E30", "seed", _r(210, 212) + _r(214, 216) + _r(218, 220) + [222, 223], False),
    ("E31", "pos:3", _r(224, 232), False),
    ("E32", "seed", _r(233, 241), False),
    # E33: the deposited metadata and expression disagree for arrays 248-251
    # (see QC_EXCLUDE_IDX). Which array belongs to which of the three replicate
    # cultures is therefore uncertain, so the whole experiment is one group.
    ("E33", "lump", _r(242, 247) + _r(249, 251), True),
    ("E40", "single", [248], True),
    ("E34", "seed", _r(252, 265), False),
    ("E35", "seed", _r(266, 271) + _r(281, 286), False),
    ("E36", "seed", _r(272, 280) + _r(287, 289), False),
    ("E37", "single", [162, 290, 291], False),
    ("E38", "single", [292], True),
    ("E39", "single", _r(293, 295), False),
]


# Arrays whose expression contradicts their deposited metadata. Found from the
# replicate structure of the expression matrix alone (no model, no residuals):
# metadata says 248-250 are the 31 C triplicate of experiment E33 and 251 opens
# the next experiment (Qp 27.4); in the expression matrix 249-251 are the tight
# triplicate (r = 0.998-0.999) and 248 sits with a different experiment
# (r = 0.96 to its supposed replicates). The Qp labels of 248 and 251 are
# therefore not trusted. They stay in the unlabeled corpus.
QC_EXCLUDE_IDX = (248, 251)


def derive_groups(meta: pd.DataFrame) -> pd.DataFrame:
    """Derive culture_id and experiment_id for every array.

    GEO states "295 microarrays from 121 individual CHO cultures" but deposits no
    culture identifier. The arrays are deposited in contiguous experiment blocks
    in which replicate cultures recur across time points with the same seeding
    density, which is a per-culture constant. That structure is encoded in
    CULTURE_RULES. experiment_id is the coarser unit (one cell line / run).
    """
    by_idx = meta.reset_index().set_index("idx")
    recs = {}
    for exp, rule, idxs, amb in CULTURE_RULES:
        if isinstance(rule, dict):
            flat = sorted(i for v in rule.values() for i in v)
            if flat != sorted(idxs):
                raise ValueError(f"{exp}: explicit rule does not cover its indices")
            for lab, members in rule.items():
                for i in members:
                    recs[i] = (exp, f"{exp}_{lab}", amb)
        elif rule == "seed":
            for i in idxs:
                sd = by_idx.loc[i, "seeding density"]
                if pd.isna(sd):
                    raise ValueError(f"{exp}: seeding density missing for idx {i}")
                recs[i] = (exp, f"{exp}_sd{sd:g}", amb)
        elif rule.startswith("pos:"):
            k = int(rule.split(":")[1])
            off = min(idxs)
            for i in idxs:
                recs[i] = (exp, f"{exp}_p{off}_{(i - off) % k}", amb)
        elif rule == "single":
            for i in idxs:
                recs[i] = (exp, f"{exp}_s{i}", amb)
        elif rule == "lump":
            for i in idxs:
                recs[i] = (exp, f"{exp}_lump{min(idxs)}", amb)
        else:
            raise ValueError(rule)
    missing = set(by_idx.index) - set(recs)
    if missing:
        raise ValueError(f"arrays without a culture assignment: {sorted(missing)}")
    out = pd.DataFrame.from_dict(recs, orient="index", columns=["experiment_id", "culture_id", "ambiguous"])
    out.index.name = "idx"
    out = out.join(by_idx[["geo_accession"]]).reset_index().set_index("geo_accession")
    out["qc_label_excluded"] = out["idx"].isin(QC_EXCLUDE_IDX)
    return out.loc[meta.index]


def expression_clusters(X, experiments, r_cut=0.5, k=1000):
    """Group experiments whose mean expression profiles are similar (r > r_cut).

    Unsupervised (no Qp). Used only for the strictest split, where whole families
    of related experiments are held out together. Average linkage on 1 - r between
    experiment centroids of the k highest-variance probesets.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    X = np.asarray(X, dtype=float)
    experiments = np.asarray(experiments)
    v = np.sort(np.argsort(X.var(axis=0))[::-1][:k])
    Xc = (X - X.mean(axis=0))[:, v]
    ex = sorted(set(experiments))
    cent = np.array([Xc[experiments == e].mean(axis=0) for e in ex])
    Dm = 1.0 - np.corrcoef(cent)
    np.fill_diagonal(Dm, 0.0)
    lab = fcluster(linkage(squareform((Dm + Dm.T) / 2, checks=False), method="average"),
                   t=1.0 - r_cut, criterion="distance")
    m = {e: f"L{int(c):02d}" for e, c in zip(ex, lab)}
    return np.array([m[e] for e in experiments])


def matched_subset_mask(meta: pd.DataFrame) -> pd.Series:
    """Stationary-phase, temperature-shifted arrays with a Qp label.

    Clarke et al. 2011 give no sample list. Rule fixed before any model was run:
    temperature below 37 C (post-shift) AND specific growth rate <= 0.01 per hour
    (stationary) AND Qp measured. Arrays with missing growth rate are excluded.
    """
    ok = meta[ENDPOINT].notna()
    if "qc_label_excluded" in meta:
        ok &= ~meta["qc_label_excluded"]
    return (meta["temperature"] < 37) & (meta["growth rate"] <= 0.01) & ok


# --------------------------------------------------------------------------- #
# Encoders (all fitted inside the training fold)
# --------------------------------------------------------------------------- #
from sklearn.base import BaseEstimator, TransformerMixin  # noqa: E402


class TopVarianceSelector(BaseEstimator, TransformerMixin):
    """Raw-gene encoder: keep the k highest-variance probesets (training fold only)."""

    def __init__(self, k: int = 1000):
        self.k = k

    def fit(self, X, y=None):
        v = np.var(np.asarray(X), axis=0)
        k = min(self.k, X.shape[1])
        self.idx_ = np.sort(np.argsort(v)[::-1][:k])
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.idx_]


class PathwayScorer(BaseEstimator, TransformerMixin):
    """Pathway-constrained encoder.

    Each output is the mean of z-scored member probesets of one gene set. The gene
    set membership is prior knowledge (fixed mask); the z-scoring parameters are
    estimated on the training fold only.
    """

    def __init__(self, membership=None, names=None):
        self.membership = membership  # list of integer index arrays into the columns of X
        self.names = names

    def fit(self, X, y=None):
        X = np.asarray(X)
        self.mu_ = X.mean(axis=0)
        sd = X.std(axis=0)
        self.sd_ = np.where(sd > 1e-8, sd, 1.0)
        return self

    def transform(self, X):
        Z = (np.asarray(X) - self.mu_) / self.sd_
        return np.column_stack([Z[:, m].mean(axis=1) for m in self.membership])


class CoexpressionModules(BaseEstimator, TransformerMixin):
    """Data-driven fallback encoder: within-fold coexpression modules.

    Hierarchical clustering (average linkage, 1 - |r|) of the top-variance
    probesets on the training fold; each module is scored as the mean of its
    sign-aligned z-scored members. Fitted entirely inside the training fold.
    """

    def __init__(self, n_modules: int = 40, k: int = 1500):
        self.n_modules = n_modules
        self.k = k

    def fit(self, X, y=None):
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform

        X = np.asarray(X)
        v = X.var(axis=0)
        self.idx_ = np.sort(np.argsort(v)[::-1][: min(self.k, X.shape[1])])
        Xs = X[:, self.idx_]
        self.mu_ = Xs.mean(axis=0)
        sd = Xs.std(axis=0)
        self.sd_ = np.where(sd > 1e-8, sd, 1.0)
        Z = (Xs - self.mu_) / self.sd_
        C = np.corrcoef(Z, rowvar=False)
        D = 1.0 - np.abs(C)
        np.fill_diagonal(D, 0.0)
        D = np.clip((D + D.T) / 2.0, 0.0, None)
        lab = fcluster(linkage(squareform(D, checks=False), method="average"),
                       t=min(self.n_modules, Xs.shape[1]), criterion="maxclust")
        self.modules_, self.signs_ = [], []
        for m in np.unique(lab):
            members = np.where(lab == m)[0]
            ref = members[np.argmax(np.abs(C[np.ix_(members, members)]).sum(axis=0))]
            self.modules_.append(members)
            self.signs_.append(np.where(C[ref, members] >= 0, 1.0, -1.0))
        return self

    def transform(self, X):
        Z = (np.asarray(X)[:, self.idx_] - self.mu_) / self.sd_
        return np.column_stack([(Z[:, m] * s).mean(axis=1) for m, s in zip(self.modules_, self.signs_)])


def build_pathway_membership(probe_ids, symbols, gene_sets, min_genes=10):
    """Map gene sets onto array columns. Returns (membership, names, stats rows)."""
    sym = pd.Series(list(symbols), index=range(len(probe_ids)))
    by_symbol: dict[str, list[int]] = {}
    for i, s in sym.items():
        if isinstance(s, str) and s:
            by_symbol.setdefault(s, []).append(i)
    membership, names, stats = [], [], []
    for name, genes in gene_sets.items():
        hit = sorted({g for g in genes if g in by_symbol})
        cols = sorted({i for g in hit for i in by_symbol[g]})
        used = len(hit) >= min_genes
        stats.append({"gene_set": name, "genes_in_set": len(set(genes)), "genes_on_array": len(hit),
                      "probesets": len(cols), "used": used})
        if used:
            membership.append(np.array(cols))
            names.append(name)
    return membership, names, pd.DataFrame(stats)


# --------------------------------------------------------------------------- #
# Cross-validation with a hard no-leak assertion
# --------------------------------------------------------------------------- #
def assert_no_group_leak(groups, train_idx, test_idx) -> None:
    """Runtime guard: no group may appear on both sides of a split."""
    g = np.asarray(groups)
    overlap = set(g[train_idx]) & set(g[test_idx])
    assert not overlap, f"group leakage: {sorted(overlap)[:5]} appear in both train and test"
    assert not (set(train_idx) & set(test_idx)), "sample appears in both train and test"


def make_splits(groups, n_splits, seed, scheme="grouped"):
    """Outer splits.

    scheme="grouped"   GroupKFold on the supplied groups, shuffled with a fixed
                       seed; n_splits=None means leave-one-group-out.
    scheme="ungrouped" KFold on arrays. Only used by the leakage demonstration;
                       it deliberately ignores the culture structure.
    """
    from sklearn.model_selection import GroupKFold, KFold, LeaveOneGroupOut

    groups = np.asarray(groups)
    n = len(groups)
    dummy = np.zeros(n)
    if scheme == "ungrouped":
        return list(KFold(n_splits=n_splits, shuffle=True, random_state=seed).split(dummy))
    n_groups = len(np.unique(groups))
    if n_splits is None or n_splits >= n_groups:
        splits = list(LeaveOneGroupOut().split(dummy, groups=groups))
    else:
        splits = list(GroupKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(dummy, groups=groups))
    for tr, te in splits:
        assert_no_group_leak(groups, tr, te)
    return splits


def _inner_cv(groups_tr, scheme, seed, n_inner=3):
    from sklearn.model_selection import GroupKFold, KFold

    if scheme == "ungrouped":
        return list(KFold(n_splits=n_inner, shuffle=True, random_state=seed).split(np.zeros(len(groups_tr))))
    k = min(n_inner, len(np.unique(groups_tr)))
    sp = list(GroupKFold(n_splits=k, shuffle=True, random_state=seed).split(np.zeros(len(groups_tr)), groups=groups_tr))
    for tr, te in sp:
        assert_no_group_leak(groups_tr, tr, te)
    return sp


def make_encoder(kind, pathway=None):
    if kind == "raw":
        return TopVarianceSelector(k=1000)
    if kind == "pathway":
        return PathwayScorer(membership=pathway[0], names=pathway[1])
    if kind == "coexpr":
        return CoexpressionModules(n_modules=40, k=1500)
    raise ValueError(kind)


def make_model(name, encoder, seed):
    """Return (pipeline, param_grid, scoring). Every step is refit inside each fold."""
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
    from sklearn.linear_model import ElasticNet, LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    steps = [("enc", encoder), ("scale", StandardScaler())]
    if name == "pls":
        return (Pipeline(steps + [("m", PLSRegression(scale=False))]),
                {"m__n_components": [1, 2, 3, 4, 6, 8, 10]}, "neg_root_mean_squared_error")
    if name == "enet":
        return (Pipeline(steps + [("m", ElasticNet(max_iter=20000, random_state=seed))]),
                {"m__alpha": [0.03, 0.1, 0.3, 1.0, 3.0], "m__l1_ratio": [0.1, 0.5]},
                "neg_root_mean_squared_error")
    if name == "gbm":
        return (Pipeline(steps + [("m", GradientBoostingRegressor(
            learning_rate=0.05, subsample=0.8, max_features="sqrt", random_state=seed))]),
                {"m__n_estimators": [150, 400], "m__max_depth": [2, 3]}, "neg_root_mean_squared_error")
    if name == "logit":
        return (Pipeline(steps + [("m", LogisticRegression(max_iter=5000, random_state=seed))]),
                {"m__C": [0.001, 0.01, 0.1, 1.0]}, "roc_auc")
    if name == "gbc":
        return (Pipeline(steps + [("m", GradientBoostingClassifier(
            learning_rate=0.05, subsample=0.8, max_features="sqrt", random_state=seed))]),
                {"m__n_estimators": [150, 400], "m__max_depth": [2, 3]}, "roc_auc")
    raise ValueError(name)


CLASSIFIERS = {"logit", "gbc"}


def nested_cv_predict(X, y, groups, model_name, encoder_kind, pathway=None, n_splits=5,
                      scheme="grouped", seed=SEED, verbose=False):
    """Out-of-fold predictions with hyperparameters chosen inside each training fold.

    Returns a DataFrame with y, pred, null (training-fold mean of y, the
    no-information reference) and fold. For classifiers pred is P(class 1).
    """
    from sklearn.model_selection import GridSearchCV

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    groups = np.asarray(groups)
    splits = make_splits(groups, n_splits, seed, scheme=scheme)
    pred = np.full(len(y), np.nan)
    null = np.full(len(y), np.nan)
    fold = np.full(len(y), -1)
    is_clf = model_name in CLASSIFIERS
    for f, (tr, te) in enumerate(splits):
        if scheme == "grouped":
            assert_no_group_leak(groups, tr, te)  # non-negotiable: checked on every outer split
        pipe, grid, scoring = make_model(model_name, make_encoder(encoder_kind, pathway), seed)
        if model_name == "pls":
            grid = {"m__n_components": [c for c in grid["m__n_components"] if c < min(len(tr) * 2 // 3 - 1, 11)]}
        inner = _inner_cv(groups[tr], scheme, seed + f + 1)
        if is_clf:
            ytr = y[tr].astype(int)
            inner = [(a, b) for a, b in inner if len(np.unique(ytr[a])) == 2 and len(np.unique(ytr[b])) == 2]
            if len(np.unique(ytr)) < 2:
                pred[te] = float(ytr.mean())
                null[te] = float(ytr.mean())
                fold[te] = f
                continue
        if len(inner) >= 2:
            gs = GridSearchCV(pipe, grid, scoring=scoring, cv=inner, n_jobs=1, refit=True)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gs.fit(X[tr], y[tr].astype(int) if is_clf else y[tr])
            best = gs.best_estimator_
        else:  # not enough class balance for an inner search: fixed default
            best = pipe.fit(X[tr], y[tr].astype(int) if is_clf else y[tr])
        if is_clf:
            pred[te] = best.predict_proba(X[te])[:, 1]
        else:
            pred[te] = np.ravel(best.predict(X[te]))
        null[te] = y[tr].mean()
        fold[te] = f
        if verbose:
            print(f"    fold {f}: n_train={len(tr)} n_test={len(te)}")
    assert not np.isnan(pred).any()
    return pd.DataFrame({"y": y, "pred": pred, "null": null, "fold": fold, "group": groups})


# --------------------------------------------------------------------------- #
# Metrics and cluster bootstrap
# --------------------------------------------------------------------------- #
def rmse(y, p):
    return float(np.sqrt(np.mean((np.asarray(y) - np.asarray(p)) ** 2)))


def q2(y, p, null):
    """Cross-validated R2: 1 - SSE(model) / SSE(training-fold mean).

    Stricter than sklearn's r2_score on pooled predictions, which centres on the
    test labels' own mean.
    """
    y, p, null = map(np.asarray, (y, p, null))
    return float(1.0 - np.sum((y - p) ** 2) / np.sum((y - null) ** 2))


def auc(y, p):
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(np.asarray(y).astype(int), np.asarray(p)))


def spearman(y, p):
    from scipy.stats import spearmanr

    return float(spearmanr(y, p).correlation)


def cluster_bootstrap(df: pd.DataFrame, stat_fns: dict, cluster_col="group", n_boot=2000, seed=SEED):
    """Percentile CI by resampling whole clusters (cultures or experiments) with replacement.

    Resampling arrays would treat replicate arrays of one culture as independent
    and give intervals that are too narrow.
    """
    rng = np.random.default_rng(seed)
    clusters = df[cluster_col].to_numpy()
    uniq = np.unique(clusters)
    members = {u: np.where(clusters == u)[0] for u in uniq}
    point = {k: fn(df) for k, fn in stat_fns.items()}
    draws = {k: [] for k in stat_fns}
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([members[u] for u in pick])
        sub = df.iloc[idx]
        for k, fn in stat_fns.items():
            try:
                v = fn(sub)
            except ValueError:
                v = np.nan
            draws[k].append(v)
    out = {}
    for k in stat_fns:
        d = np.asarray(draws[k], dtype=float)
        d = d[np.isfinite(d)]
        out[k] = (point[k], float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)))
    return out


REG_STATS = {
    "rmse": lambda d: rmse(d["y"], d["pred"]),
    "r2_cv": lambda d: q2(d["y"], d["pred"], d["null"]),
    "spearman": lambda d: spearman(d["y"], d["pred"]),
    "rmse_null": lambda d: rmse(d["y"], d["null"]),
}
CLF_STATS = {"auc": lambda d: auc(d["y"], d["pred"])}


# --------------------------------------------------------------------------- #
# Delta decode
# --------------------------------------------------------------------------- #
def delta_decode(Z, y, names, clusters, n_boot=2000, seed=SEED, q=0.25):
    """mean(z | top-quartile outcome) - mean(z | bottom-quartile outcome), per feature.

    Z is standardised over the supplied samples so deltas are in SD units. The
    interval comes from a cluster bootstrap in which the quartile cut-offs are
    recomputed on every resample.
    """
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float)
    clusters = np.asarray(clusters)
    Zs = (Z - Z.mean(axis=0)) / np.where(Z.std(axis=0) > 1e-8, Z.std(axis=0), 1.0)

    def one(Zm, ym):
        lo, hi = np.quantile(ym, [q, 1 - q])
        top, bot = ym >= hi, ym <= lo
        if top.sum() < 2 or bot.sum() < 2:
            return np.full(Zm.shape[1], np.nan)
        return Zm[top].mean(axis=0) - Zm[bot].mean(axis=0)

    point = one(Zs, y)
    rng = np.random.default_rng(seed)
    uniq = np.unique(clusters)
    members = {u: np.where(clusters == u)[0] for u in uniq}
    boots = np.empty((n_boot, Zs.shape[1]))
    for b in range(n_boot):
        idx = np.concatenate([members[u] for u in rng.choice(uniq, size=len(uniq), replace=True)])
        boots[b] = one(Zs[idx], y[idx])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5], axis=0)
    out = pd.DataFrame({"feature": names, "delta": point, "ci_lo": lo, "ci_hi": hi})
    out["ci_excludes_zero"] = (out["ci_lo"] > 0) | (out["ci_hi"] < 0)
    return out


def center_within(values, clusters):
    """Subtract each cluster's mean (removes between-experiment differences)."""
    df = pd.DataFrame(np.asarray(values, dtype=float))
    return (df - df.groupby(np.asarray(clusters)).transform("mean")).to_numpy()


# --------------------------------------------------------------------------- #
# Self-test (offline, synthetic)
# --------------------------------------------------------------------------- #
def _synthetic(n_cultures=60, arrays_per=3, n_genes=400, seed=0, signal=1.0):
    """Synthetic data with the real data's nuisance: replicate arrays share a
    culture-specific expression signature that is unrelated to the outcome
    except through one planted pathway."""
    rng = np.random.default_rng(seed)
    n = n_cultures * arrays_per
    groups = np.repeat(np.arange(n_cultures), arrays_per)
    culture_sig = rng.normal(0, 1.0, (n_cultures, n_genes))
    y_culture = rng.normal(20, 8, n_cultures)
    X = culture_sig[groups] + rng.normal(0, 0.3, (n, n_genes))
    planted = np.arange(0, 20)
    X[:, planted] += signal * ((y_culture[groups] - 20) / 8)[:, None]
    y = y_culture[groups] + rng.normal(0, 0.5, n)
    return X, y, groups, planted


def self_test() -> int:
    t0 = time.time()
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name} {detail}")

    print("self-test (synthetic, offline)")
    # 1. series-matrix parser
    txt = "\n".join([
        '!Series_title\t"toy"',
        '!Sample_title\t"c_1"\t"c_2"\t"c_3"',
        '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"',
        '!Sample_characteristics_ch1\t"qproductivity: 4.5"\t"qproductivity: NA"\t"qproductivity: 7"',
        '!Sample_characteristics_ch1\t"temperature: 37"\t"temperature: 31"\t"temperature: 31"',
        "!series_matrix_table_begin",
        '"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"',
        '"p1_at"\t1.0\t2.0\t3.0',
        '"p2_at"\t4.0\t5.0\t6.0',
        "!series_matrix_table_end", ""])
    expr, meta = parse_series_matrix(txt)
    check("parser: matrix shape", expr.shape == (2, 3))
    check("parser: NA label -> NaN", np.isnan(meta.loc["GSM2", "qproductivity"]) and meta.loc["GSM3", "qproductivity"] == 7)
    check("parser: numeric characteristics", meta["temperature"].tolist() == [37, 31, 31])

    # 2. the no-leak assertion fires on a leaky split and passes on grouped splits
    X, y, groups, planted = _synthetic()
    try:
        assert_no_group_leak(groups, np.array([0, 1]), np.array([2, 3]))
        check("leak guard fires on a leaky split", False)
    except AssertionError:
        check("leak guard fires on a leaky split", True)
    sp = make_splits(groups, 5, SEED)
    check("grouped splits are group-disjoint", all(not (set(groups[a]) & set(groups[b])) for a, b in sp))
    check("grouped splits cover every sample once", sorted(np.concatenate([b for _, b in sp]).tolist()) == list(range(len(y))))

    # 3. grouped vs ungrouped CV: ungrouped must look better on replicate-structured data
    Xn, yn, gn, _ = _synthetic(signal=0.0, seed=1)  # no real signal at all
    g_df = nested_cv_predict(Xn, yn, gn, "pls", "raw", n_splits=5, scheme="grouped")
    u_df = nested_cv_predict(Xn, yn, gn, "pls", "raw", n_splits=5, scheme="ungrouped")
    rg, ru = rmse(g_df.y, g_df.pred), rmse(u_df.y, u_df.pred)
    check("null data: grouped CV shows no skill", q2(g_df.y, g_df.pred, g_df.null) < 0.15, f"(R2cv={q2(g_df.y, g_df.pred, g_df.null):.2f})")
    check("null data: ungrouped CV is inflated", ru < 0.8 * rg, f"(RMSE grouped={rg:.2f}, ungrouped={ru:.2f})")

    # 4. planted signal is recovered under grouped CV
    s_df = nested_cv_predict(X, y, groups, "pls", "raw", n_splits=5)
    check("planted signal: grouped R2cv > 0.5", q2(s_df.y, s_df.pred, s_df.null) > 0.5, f"(R2cv={q2(s_df.y, s_df.pred, s_df.null):.2f})")

    # 5. cluster bootstrap CI
    ci = cluster_bootstrap(s_df, {"rmse": REG_STATS["rmse"]}, n_boot=300)["rmse"]
    check("bootstrap CI brackets the point estimate", ci[1] <= ci[0] <= ci[2], f"({ci[1]:.2f} <= {ci[0]:.2f} <= {ci[2]:.2f})")
    iid = s_df.assign(group=np.arange(len(s_df)))
    ci_iid = cluster_bootstrap(iid, {"rmse": REG_STATS["rmse"]}, n_boot=300)["rmse"]
    check("cluster CI is wider than array-level CI", (ci[2] - ci[1]) > (ci_iid[2] - ci_iid[1]))

    # 6. pathway encoder + delta decode recover the planted pathway
    rng = np.random.default_rng(3)
    membership = [planted] + [rng.choice(np.arange(20, X.shape[1]), 20, replace=False) for _ in range(9)]
    names = ["PLANTED"] + [f"RANDOM_{i}" for i in range(9)]
    enc = PathwayScorer(membership=membership, names=names).fit(X)
    Zp = enc.transform(X)
    check("pathway encoder output shape", Zp.shape == (len(y), 10))
    dd = delta_decode(Zp, y, names, groups, n_boot=200).sort_values("delta", key=np.abs, ascending=False)
    check("delta decode ranks planted pathway first", dd.iloc[0]["feature"] == "PLANTED" and dd.iloc[0]["ci_excludes_zero"])
    check("delta decode: random sets mostly null", dd.iloc[1:]["ci_excludes_zero"].sum() <= 2)
    p_df = nested_cv_predict(X, y, groups, "pls", "pathway", pathway=(membership, names), n_splits=5)
    check("pathway encoder predicts under grouped CV", q2(p_df.y, p_df.pred, p_df.null) > 0.5)

    # 7. coexpression-module encoder is fold-safe and usable
    c_df = nested_cv_predict(X, y, groups, "pls", "coexpr", n_splits=5)
    check("coexpression encoder runs inside folds", np.isfinite(c_df.pred).all())

    # 8. classifier path
    hi, lo = np.quantile(y, 2 / 3), np.quantile(y, 1 / 3)
    keep = (y >= hi) | (y <= lo)
    k_df = nested_cv_predict(X[keep], (y[keep] >= hi).astype(float), groups[keep], "logit", "raw", n_splits=5)
    check("classifier AUC on planted signal > 0.8", auc(k_df.y, k_df.pred) > 0.8, f"(AUC={auc(k_df.y, k_df.pred):.2f})")

    # 9. determinism
    s2 = nested_cv_predict(X, y, groups, "pls", "raw", n_splits=5)
    check("fixed seed -> identical predictions", np.allclose(s_df.pred, s2.pred))

    print(f"self-test {'PASSED' if ok else 'FAILED'} in {time.time() - t0:.0f}s")
    return 0 if ok else 1


# --------------------------------------------------------------------------- #
# Main analysis
# --------------------------------------------------------------------------- #
def prepare():
    """Load everything the analysis needs. Returns a dict."""
    expr, meta = load_dataset()
    gpl = parse_platform(RAW / SOURCES["platform"]["path"])
    grp = derive_groups(meta)
    meta = meta.join(grp[["experiment_id", "culture_id", "ambiguous", "qc_label_excluded"]])
    probes = list(expr.index)
    symbols = gpl["symbol"].reindex(probes).tolist()
    hallmark = read_gmt(RAW / SOURCES["hallmark"]["path"])
    go = read_gmt(RAW / SOURCES["go_bp"]["path"], keep=GO_THEME_NAMES)
    membership, names, stats = build_pathway_membership(probes, symbols, hallmark, min_genes=10)
    X = expr.T.to_numpy(dtype=float)
    meta["expression_cluster"] = expression_clusters(X, meta["experiment_id"].to_numpy())
    return {"expr": expr, "meta": meta, "X": X, "probes": probes, "symbols": symbols,
            "hallmark": hallmark, "go": go, "pathway": (membership, names), "pathway_stats": stats}


def _row(analysis, subset, cv, encoder, n_features, model, df, stats, cluster_col, n_boot):
    ci = cluster_bootstrap(df, stats, cluster_col=cluster_col, n_boot=n_boot)
    rows = []
    for metric, (pt, lo, hi) in ci.items():
        rows.append({"analysis": analysis, "subset": subset, "cv_scheme": cv, "encoder": encoder,
                     "n_features": n_features, "model": model, "metric": metric,
                     "value": round(pt, 4), "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                     "n_arrays": len(df), "n_cultures": df["culture"].nunique(),
                     "n_experiments": df["experiment"].nunique(),
                     "bootstrap_unit": {"group": "cv group", "culture": "culture",
                                        "experiment": "experiment", "array": "array"}[cluster_col]})
    return rows


CV_SCHEMES = {
    # name: (grouping column, n_splits, scheme, bootstrap unit)
    "array_kfold_LEAKY": ("array", 5, "ungrouped", "array"),
    "culture_grouped_5fold": ("culture", 5, "grouped", "culture"),
    "leave_one_experiment_out": ("experiment", None, "grouped", "experiment"),
    "leave_one_expression_cluster_out": ("cluster", None, "grouped", "experiment"),
}


def run_cv(D, mask, target, model, encoder, cv_name, seed=SEED):
    meta = D["meta"]
    idx = np.where(mask)[0]
    col, n_splits, scheme, _ = CV_SCHEMES[cv_name]
    cultures = meta["culture_id"].to_numpy()[idx]
    experiments = meta["experiment_id"].to_numpy()[idx]
    groups = {"array": np.arange(len(idx)), "culture": cultures, "experiment": experiments,
              "cluster": meta["expression_cluster"].to_numpy()[idx]}[col]
    df = nested_cv_predict(D["X"][idx], target[idx], groups, model, encoder, pathway=D["pathway"],
                           n_splits=n_splits, scheme=scheme, seed=seed)
    df["culture"] = cultures
    df["experiment"] = experiments
    df["array"] = meta.index.to_numpy()[idx]
    if scheme == "grouped":  # belt and braces: culture-disjointness of every fold, re-checked on the output
        for f in np.unique(df["fold"]):
            te = df["fold"] == f
            assert not (set(df.loc[te, "culture"]) & set(df.loc[~te, "culture"])), "culture leak"
    return df


def n_features_of(D, encoder):
    return {"raw": min(1000, D["X"].shape[1]), "pathway": len(D["pathway"][1]), "coexpr": 40}[encoder]


def run_analysis(n_boot=2000, quick=False):
    RESULTS.mkdir(exist_ok=True)
    D = prepare()
    meta = D["meta"]
    y = meta[ENDPOINT].to_numpy(dtype=float)
    as_deposited = ~np.isnan(y)
    labeled = as_deposited & ~meta["qc_label_excluded"].to_numpy()
    matched = matched_subset_mask(meta).to_numpy()
    meta_out = meta[["idx", "experiment_id", "culture_id", "expression_cluster", "ambiguous", "qc_label_excluded",
                     "temperature", "growth rate", ENDPOINT]].copy()
    meta_out["labeled"] = labeled
    meta_out["matched_subset"] = matched
    meta_out.to_csv(RESULTS / "culture_map.csv")
    D["pathway_stats"].to_csv(RESULTS / "pathway_coverage.csv", index=False)

    summary = dataset_summary(D, labeled, matched)
    (RESULTS / "dataset_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    rows, preds = [], {}

    def do(analysis, subset_name, mask, target, model, encoder, cv_name):
        t = time.time()
        df = run_cv(D, mask, target, model, encoder, cv_name)
        stats = CLF_STATS if model in CLASSIFIERS else REG_STATS
        rows.extend(_row(analysis, subset_name, cv_name, encoder, n_features_of(D, encoder), model, df, stats,
                         CV_SCHEMES[cv_name][3], n_boot))
        key = f"{analysis}|{subset_name}|{cv_name}|{encoder}|{model}"
        preds[key] = df.assign(key=key)
        main = "auc" if model in CLASSIFIERS else "rmse"
        r = [x for x in rows[-len(stats):] if x["metric"] == main][0]
        print(f"  {key}: {main}={r['value']:.3f} [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]  ({time.time() - t:.0f}s)")
        pd.DataFrame(rows).to_csv(RESULTS / "metrics.csv", index=False)
        return df

    # ---- 1. matched baseline comparison (the only head-to-head against 4.44) ----
    print("matched subset (stationary, temperature-shifted)")
    for cv_name in ["culture_grouped_5fold", "leave_one_experiment_out", "array_kfold_LEAKY"]:
        do("matched_baseline", "matched_stationary_tempshift", matched, y, "pls", "raw", cv_name)

    # ---- 2. full labeled dataset, regression ----
    print("full labeled dataset: regression")
    reg_models = ["pls", "enet"] if quick else ["pls", "enet", "gbm"]
    encoders = ["raw", "pathway"] if quick else ["raw", "pathway", "coexpr"]
    for cv_name in ["culture_grouped_5fold", "leave_one_experiment_out"]:
        for enc in encoders:
            for model in reg_models:
                if enc == "coexpr" and model != "pls":
                    continue
                do("full_regression", "all_labeled", labeled, y, model, enc, cv_name)

    # ---- 3. full labeled dataset, classification (top vs bottom tercile) ----
    print("full labeled dataset: classification (top vs bottom Qp tercile)")
    lo, hi = np.nanquantile(y, [1 / 3, 2 / 3])
    cls_mask = labeled & ((y <= lo) | (y >= hi))
    ycls = (y >= hi).astype(float)
    clf_models = ["logit"] if quick else ["logit", "gbc"]
    for cv_name in ["culture_grouped_5fold", "leave_one_experiment_out"]:
        for enc in ["raw", "pathway"]:
            for model in clf_models:
                do("full_classification", "all_labeled_terciles", cls_mask, ycls, model, enc, cv_name)

    # ---- 4. leakage comparison: same model, four split units ----
    print("leakage comparison")
    for cv_name in ["array_kfold_LEAKY", "leave_one_expression_cluster_out"]:
        do("leakage", "all_labeled", labeled, y, "pls", "raw", cv_name)
        do("leakage", "all_labeled_terciles", cls_mask, ycls, "logit", "raw", cv_name)

    # ---- 4b. sensitivity: the two QC-excluded arrays put back, labels as deposited ----
    print("sensitivity: no QC exclusion")
    for cv_name in ["culture_grouped_5fold", "leave_one_experiment_out"]:
        do("sensitivity_labels_as_deposited", "all_labeled_incl_qc_flagged", as_deposited, y, "pls", "raw", cv_name)

    pd.concat(preds.values()).to_csv(RESULTS / "oof_predictions.csv", index=False)
    pd.DataFrame(rows).to_csv(RESULTS / "metrics.csv", index=False)

    # ---- 5. delta decode ----
    run_delta(D, labeled, n_boot=n_boot)
    return D


def dataset_summary(D, labeled, matched):
    meta = D["meta"]
    sym = pd.Series(D["symbols"])
    hall_genes = set().union(*[set(v) for v in D["hallmark"].values()])
    on_array = set(sym.dropna())
    st = D["pathway_stats"]
    y = meta[ENDPOINT]
    return {
        "arrays": int(len(meta)),
        "probesets": int(D["X"].shape[1]),
        "arrays_with_qp_as_deposited": int(meta[ENDPOINT].notna().sum()),
        "arrays_with_qp_after_qc": int(labeled.sum()),
        "qc_excluded_idx": list(QC_EXCLUDE_IDX),
        "expression_clusters_labeled": int(meta.loc[labeled, "expression_cluster"].nunique()),
        "expression_cluster_members": {k: sorted(set(v)) for k, v in
                                       meta.loc[labeled].groupby("expression_cluster")["experiment_id"]},
        "derived_cultures_total": int(meta["culture_id"].nunique()),
        "cultures_stated_by_geo": 121,
        "derived_cultures_unambiguous": int(meta.loc[~meta["ambiguous"], "culture_id"].nunique()),
        "derived_cultures_in_ambiguous_blocks": int(meta.loc[meta["ambiguous"], "culture_id"].nunique()),
        "arrays_in_ambiguous_blocks": int(meta["ambiguous"].sum()),
        "labeled_arrays_in_ambiguous_blocks": int((meta["ambiguous"].to_numpy() & labeled).sum()),
        "derived_experiments_total": int(meta["experiment_id"].nunique()),
        "labeled_cultures": int(meta.loc[labeled, "culture_id"].nunique()),
        "labeled_experiments": int(meta.loc[labeled, "experiment_id"].nunique()),
        "qp_mean": round(float(y[labeled].mean()), 2),
        "qp_sd": round(float(y[labeled].std()), 2),
        "qp_min": float(y[labeled].min()),
        "qp_max": float(y[labeled].max()),
        "qp_variance_between_experiments_fraction": round(float(
            1 - (y[labeled] - y[labeled].groupby(meta.loc[labeled, "experiment_id"]).transform("mean")).var(ddof=0)
            / y[labeled].var(ddof=0)), 3),
        "matched_arrays": int(matched.sum()),
        "matched_cultures": int(meta.loc[matched, "culture_id"].nunique()),
        "matched_experiments": int(meta.loc[matched, "experiment_id"].nunique()),
        "matched_qp_sd": round(float(y[matched].std()), 2),
        "probesets_with_symbol": int(sym.notna().sum()),
        "symbol_mapping_coverage": round(float(sym.notna().mean()), 3),
        "unique_symbols": int(sym.dropna().nunique()),
        "hallmark_genes_total": len(hall_genes),
        "hallmark_genes_on_array": len(hall_genes & on_array),
        "hallmark_gene_match_rate": round(len(hall_genes & on_array) / len(hall_genes), 3),
        "array_symbols_in_hallmark_rate": round(len(hall_genes & on_array) / len(on_array), 3),
        "hallmark_sets_used_min10_genes": int(st["used"].sum()),
        "hallmark_sets_dropped": st.loc[~st["used"], "gene_set"].tolist(),
        "feature_counts": {"raw": n_features_of(D, "raw"), "pathway": n_features_of(D, "pathway"),
                           "coexpr": n_features_of(D, "coexpr")},
    }


def run_delta(D, labeled, n_boot=2000):
    """Delta decode on the pathway encoder, with the adversarial within-experiment variant."""
    meta = D["meta"]
    idx = np.where(labeled)[0]
    y = meta[ENDPOINT].to_numpy(dtype=float)[idx]
    X = D["X"][idx]
    cultures = meta["culture_id"].to_numpy()[idx]
    exps = meta["experiment_id"].to_numpy()[idx]

    # Hallmark sets used by the encoder + the pre-specified GO theme sets
    sets = dict(D["hallmark"])
    sets.update(D["go"])
    membership, names, stats = build_pathway_membership(D["probes"], D["symbols"], sets, min_genes=10)
    Z = PathwayScorer(membership=membership, names=names).fit(X).transform(X)

    d_cult = delta_decode(Z, y, names, cultures, n_boot=n_boot)
    d_exp = delta_decode(Z, y, names, exps, n_boot=n_boot)
    # within-experiment: remove each experiment's mean from both z and Qp, so the
    # contrast cannot be driven by which cell line / run an array belongs to
    multi = pd.Series(exps).map(pd.Series(exps).value_counts()).to_numpy() >= 4
    Zw = center_within(Z[multi], exps[multi])
    yw = center_within(y[multi][:, None], exps[multi]).ravel()
    d_within = delta_decode(Zw, yw, names, exps[multi], n_boot=n_boot)

    # gene-level decode: rank member genes by their own top-vs-bottom-quartile delta
    Xz = (X - X.mean(axis=0)) / np.where(X.std(axis=0) > 1e-8, X.std(axis=0), 1.0)
    lo, hi = np.quantile(y, [0.25, 0.75])
    gene_delta = Xz[y >= hi].mean(axis=0) - Xz[y <= lo].mean(axis=0)
    sym = np.array([s if isinstance(s, str) else "" for s in D["symbols"]])

    st = stats.set_index("gene_set")
    out = d_cult.rename(columns={"feature": "pathway", "ci_lo": "ci_lo_culture_boot", "ci_hi": "ci_hi_culture_boot"})
    out = out.drop(columns="ci_excludes_zero")
    out["ci_lo_experiment_boot"] = d_exp["ci_lo"].to_numpy()
    out["ci_hi_experiment_boot"] = d_exp["ci_hi"].to_numpy()
    out["robust_to_experiment_bootstrap"] = d_exp["ci_excludes_zero"].to_numpy()
    out["delta_within_experiment"] = d_within["delta"].to_numpy()
    out["within_ci_lo"] = d_within["ci_lo"].to_numpy()
    out["within_ci_hi"] = d_within["ci_hi"].to_numpy()
    out["within_experiment_ci_excludes_zero"] = d_within["ci_excludes_zero"].to_numpy()
    out["source"] = ["MSigDB Hallmark" if n.startswith("HALLMARK_") else "GO BP (pre-specified theme)" for n in names]
    out["wgcna_theme"] = [THEME_SETS.get(n, "") for n in names]
    out["genes_on_array"] = [int(st.loc[n, "genes_on_array"]) for n in names]
    out["genes_in_set"] = [int(st.loc[n, "genes_in_set"]) for n in names]
    drivers = []
    for n, m, d in zip(names, membership, out["delta"]):
        gd = pd.Series(gene_delta[m], index=sym[m]).groupby(level=0).mean()
        gd = gd.sort_values(ascending=(d < 0))
        drivers.append("; ".join(f"{g} ({v:+.2f})" for g, v in gd.head(8).items()))
    out["top_driver_genes"] = drivers
    out["abs_delta"] = out["delta"].abs()
    out = out.sort_values("abs_delta", ascending=False).drop(columns="abs_delta")
    out.insert(0, "rank", range(1, len(out) + 1))
    num = out.select_dtypes("number").columns.difference(["rank", "genes_on_array", "genes_in_set"])
    out[num] = out[num].round(3)
    out.to_csv(RESULTS / "delta_pathways.csv", index=False)
    info = {"n_arrays": int(len(y)), "top_quartile_cut": float(hi), "bottom_quartile_cut": float(lo),
            "n_arrays_within_experiment_variant": int(multi.sum()),
            "n_experiments_within_experiment_variant": int(len(np.unique(exps[multi])))}
    (RESULTS / "delta_info.json").write_text(json.dumps(info, indent=2))
    print("delta decode: top 12 pathways")
    print(out[["rank", "pathway", "delta", "ci_lo_culture_boot", "ci_hi_culture_boot",
               "robust_to_experiment_bootstrap", "delta_within_experiment",
               "within_experiment_ci_excludes_zero", "wgcna_theme"]].head(12).to_string(index=False))
    print("delta decode: pre-specified WGCNA themes")
    print(out.loc[out["wgcna_theme"] != "", ["rank", "pathway", "delta", "ci_lo_culture_boot", "ci_hi_culture_boot",
                                              "robust_to_experiment_bootstrap", "delta_within_experiment",
                                              "within_experiment_ci_excludes_zero", "genes_on_array"]].to_string(index=False))
    return out


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--figures", action="store_true")
    ap.add_argument("--pretrain", action="store_true")
    ap.add_argument("--all", action="store_true", help="self-test, download, run, pretrain, figures")
    ap.add_argument("--quick", action="store_true", help="smaller model grid and fewer bootstraps (smoke run)")
    ap.add_argument("--n-boot", type=int, default=2000)
    a = ap.parse_args(argv)
    if not any([a.self_test, a.download, a.run, a.figures, a.pretrain, a.all]):
        ap.print_help()
        return 2
    if a.self_test or a.all:
        rc = self_test()
        if rc:
            return rc
    if a.download or a.all:
        print("download")
        download_all()
    if a.run or a.all:
        run_analysis(n_boot=200 if a.quick else a.n_boot, quick=a.quick)
    if a.pretrain or a.all:
        sys.path.insert(0, str(ROOT / "src"))
        try:
            import pretrain
        except ImportError as exc:  # torch missing: the core result must not depend on it
            print(f"pretraining skipped ({exc}); install torch to run it")
        else:
            pretrain.main(quick=a.quick)
    if a.figures or a.all:
        sys.path.insert(0, str(ROOT / "src"))
        import figures

        figures.make_all()
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "0")
    sys.exit(main())
