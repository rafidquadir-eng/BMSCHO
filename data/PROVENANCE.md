# Data provenance

Raw files live in `raw/` (git-ignored) and are re-fetched by `python cho_poc.py --download`.
No raw data is committed. This file is regenerated on every download.

| name | URL | local path (under raw/) | bytes | SHA256 | downloaded (UTC) | check |
|---|---|---|---|---|---|---|
| series_matrix | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30321/matrix/GSE30321_series_matrix.txt.gz | `geo/GSE30321_series_matrix.txt.gz` | 5723843 | `790c5ef6b48070b0a5af3ea9596466c78fc29d25b1ea20e0ecac37254b465cdf` | 2026-11-15T01:31:18Z | matches pinned hash |
| platform | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13791&targ=self&form=text&view=data | `gpl/GPL13791.txt` | 255331 | `f102b8d336e68a88646306659526c7452cabdaefa2f25e0cef33cb961447208f` | 2026-11-15T01:36:28Z | matches pinned hash |
| sdrf | https://www.ebi.ac.uk/biostudies/files/E-GEOD-30321/E-GEOD-30321.sdrf.txt | `ae/E-GEOD-30321.sdrf.txt` | 155039 | `69ce5da3854675369ca2784a611ad9cd8b570a6dd3610a71bef87473ce0818ee` | 2026-11-15T01:32:19Z | matches pinned hash |
| hallmark | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt | `genesets/h.all.v2023.2.Hs.symbols.gmt` | 48690 | `8ff1f036f0988c99d42cf9e5ea24d355f10b8ece5e0d4f6edf26a0458f480780` | 2026-11-15T01:34:51Z | matches pinned hash |
| go_bp | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.symbols.gmt | `genesets/c5.go.bp.v2023.2.Hs.symbols.gmt` | 4929531 | `b4b3b49a267c2239fad40d39a56f2146a57743992943915e80ebed34ec7ebaae` | 2026-11-15T01:40:22Z | matches pinned hash |

## Contents

- Expression matrix: 3592 probesets x 295 arrays (log2 RMA, as deposited by the submitters; no re-normalisation here).
- Platform GPL13791 (WyeHamster2a custom Affymetrix CHO array): 3592 probesets.
- Sample metadata: 295 arrays, bioprocess variables parsed from `!Sample_characteristics_ch1`.
- Arrays with a Qp label (`qproductivity`): 169 of 295.
- ArrayExpress SDRF (E-GEOD-30321): carries the same nine characteristics, no extra columns (no culture, cell-line or time-point identifier).
- MSigDB Hallmark v2023.2.Hs: 50 gene sets.
- MSigDB C5 GO BP v2023.2.Hs: only the four pre-specified theme sets are used (GOBP_CELL_CYCLE, GOBP_PROTEIN_LOCALIZATION_TO_EXTRACELLULAR_REGION, GOBP_VESICLE_MEDIATED_TRANSPORT, GOBP_GOLGI_VESICLE_TRANSPORT).
- GEO `/suppl/` for GSE30321 does not exist (HTTP 404); per-sample supplementary file is `NONE` for all 295 arrays, so no CEL files are public.

## Label and grouping provenance

- Qp is taken verbatim from GEO sample characteristics (`qproductivity`, pg/cell/day).
- There is no culture identifier in GEO, the SOFT file or the ArrayExpress SDRF. The culture and experiment IDs are derived (see `derive_groups` in `cho_poc.py` and `results/culture_map.csv`).
