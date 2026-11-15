# Decision log

Running log of every autonomous choice. Times are UTC. Newest entries at the bottom.

---

### 2026-11-15 01:30 — The repository was empty; `cho_poc.py` did not exist
- **Found:** `rafidquadir-eng/BMSCHO` had no commits, no branches, no `cho_poc.py`, no `main`.
- **Decision:** Write `cho_poc.py` from scratch with the capabilities the brief lists (series-matrix parser, grouped CV, bootstrap CIs, pathway encoder, Δ decode, within-fold coexpression modules, `--self-test`).
- **Alternatives:** Stop and ask (ruled out by the brief).
- **Consequence:** "Already tested" does not apply. The self-test is new and is the only test of this code.

### 2026-11-15 01:31 — Branch and PR base
- **Found:** No `main` exists, so a PR "to `main`" has no base. A push of an initial commit to `main` was blocked by the session's permission layer (this session may only push to its own working branch without further approval).
- **Decision:** Commit on `poc/gse30321` as instructed. Do not create `main` by another route. The PR is opened if a base branch can be created; otherwise this is reported as the one undelivered item.
- **Alternatives:** Push the work to `main` directly (not authorised); make `poc/gse30321` the default branch and skip the PR.

### 2026-11-15 01:31 — `.gitignore` deviates from the brief by one line
- **Decision:** `data/` became `data/*` plus `!data/PROVENANCE.md`.
- **Reasoning:** The brief asks for `data/` to be ignored and also for `data/PROVENANCE.md` as a deliverable. Git cannot re-include a file inside an ignored directory, so the directory pattern had to change. Everything else under `data/` is still ignored.

### 2026-11-15 01:32 — GPU instance not used; everything runs locally on CPU
- **Found:** This session runs in a cloud container. It has no `ssh` client and the key file lives on a laptop path that does not exist here. No GPU is present (`nvidia-smi` missing). 4 CPU cores, 15 GB RAM, Python 3.13.
- **Decision:** Run every phase locally, including §6. Nothing is staged on the instance.
- **Reasoning:** The brief says a broken instance must never cost the deadline. The expression matrix is 295 × 3,592, so the denoising autoencoder trains in seconds on CPU.
- **Consequence:** `nvidia-smi`, `df -h` and the Python version of the instance were not recorded.

### 2026-11-15 01:35 — Qp labels are in GEO; no fallback needed
- **Found:** `!Sample_characteristics_ch1` carries nine variables per array, including `qproductivity`. 169 of 295 arrays have a Qp value. GEO has no `/suppl/` directory for this series (HTTP 404) and every sample lists supplementary file `NONE`, so there are no CEL files. The ArrayExpress SDRF (E-GEOD-30321) holds the same nine columns and nothing more. Fallbacks 3–5 were not needed for Qp.
- **Decision:** Primary endpoint stays Qp. No pivot.

### 2026-11-15 01:36 — The array is far smaller than the brief assumed
- **Found:** GPL13791 (WyeHamster2a) has 3,592 probesets, not ~20,000. Values are log2 RMA as deposited.
- **Decision:** Use the deposited RMA values. No re-normalisation (impossible without CEL files).
- **Consequence:** Pathway coverage is thin (next entry). Any gene the array lacks cannot appear in the Δ decode.

### 2026-11-15 01:45 — Annotation and gene sets
- **Found:** 2,045 of 3,592 probesets (56.9%) map to a usable gene symbol (1,707 unique). Coverage is above the 40% threshold, so the pathway encoder stays in play. MSigDB Hallmark v2023.2 downloaded without registration from the Broad's public data server. 653 of 4,384 Hallmark genes (14.9%) are on the array; 38% of annotated array genes are in at least one Hallmark set.
- **Decision:** Match by upper-cased symbol. Keep Hallmark sets with at least 10 genes on the array: 44 of 50. Dropped: ANGIOGENESIS, APICAL_SURFACE, HEDGEHOG_SIGNALING, INTERFERON_ALPHA_RESPONSE, NOTCH_SIGNALING, WNT_BETA_CATENIN_SIGNALING.
- **Alternatives:** Keep all 50 (scores from 3–8 genes are noise); orthology mapping through Ensembl (no time, and symbols on this array are already human-style).
- **Also:** Symbols that are UniGene IDs (`Rn.147997`) or free text are treated as unmapped.

### 2026-11-15 01:50 — WGCNA themes fixed before running the Δ decode
- **Decision:** The three themes of the companion paper are mapped to gene sets now, before any decode output exists: cell cycle → HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT, HALLMARK_MITOTIC_SPINDLE, GOBP_CELL_CYCLE; protein secretion → HALLMARK_PROTEIN_SECRETION, GOBP_PROTEIN_LOCALIZATION_TO_EXTRACELLULAR_REGION; vesicle transport → GOBP_VESICLE_MEDIATED_TRANSPORT, GOBP_GOLGI_VESICLE_TRANSPORT.
- **Reasoning:** Choosing the sets after seeing the ranking would let me find the themes by construction. Hallmark has no vesicle-transport set, so four GO BP sets were added for the decode only. They are not features of the predictive encoder.
- **Note:** GO "protein secretion" (GO:0009306) is not a set in MSigDB 2023.2; its parent term is used.

### 2026-11-15 02:00 — Culture ID does not exist publicly; it is derived
- **Found:** No culture, cell-line, time-point or growth-phase identifier in the series matrix, the SOFT family file or the SDRF.
- **Decision:** Derive two IDs from the deposit structure (rules are in `CULTURE_RULES` in `cho_poc.py`; output in `results/culture_map.csv`):
  - `culture_id` — arrays are deposited in contiguous experiment blocks in which replicate cultures recur across time points with identical seeding density. Seeding density is constant per culture, so within a block it identifies the culture. Where replicates share a seeding density, position in the repeating replicate order is used. Where one time point was deposited, each array is its own culture.
  - `experiment_id` — the block itself (one cell line / run as deposited).
- **Result:** 122 derived cultures against 121 stated by GEO; 38 experiment blocks. Blocks where the metadata cannot settle the assignment are flagged `ambiguous`; apart from E33 (below) they contain no Qp-labeled array, so they touch no predictive metric. Where unsure, arrays are merged into one group, which is the leakage-safe direction.
- **Check against expression (label-free):** for 93% of arrays the most-correlated other array is in the same derived experiment; the exceptions are arrays with no metadata at all, plus the E33 case below.
- **Limit:** I cannot prove the map equals the authors' 121 cultures. The count agrees to within one.

### 2026-11-15 02:05 — Two split units are reported, not one
- **Found:** 79% of the variance of Qp lies between experiment blocks. Replicate flasks of one experiment have near-identical Qp and near-identical expression (r ≈ 0.998).
- **Decision:** Culture-grouped 5-fold CV is run exactly as the brief specifies and satisfies non-negotiable 1. Leave-one-experiment-out is reported next to it everywhere. A third, stricter split holds out whole clusters of experiments with similar expression profiles (9 clusters; unsupervised, average linkage, r > 0.5 between experiment centroids).
- **Reasoning:** Culture-grouped CV still places sister flasks of the same experiment in train and test. A reviewer who has built these models will ask what happens on an unseen cell line. Experiment-level splits are coarser than culture-level ones, so every one of them also satisfies "no culture in both train and test".
- **Alternatives:** Report culture-grouped only (meets the letter of the brief, invites the obvious objection).

### 2026-11-15 02:10 — Definition of the matched subset
- **Found:** The paper's 70 samples are not listed anywhere public and growth phase is not annotated.
- **Decision (fixed before any model was run):** temperature < 37 °C AND specific growth rate ≤ 0.01 h⁻¹ AND Qp measured. Arrays with missing growth rate are excluded. This gives 75 arrays, 38 cultures, 12 experiments.
- **Alternatives considered and their sizes:** threshold 0.005 → 66; 0.0075 → 72; 0.0125 → 80; all temperature-shifted arrays with Qp → 95. The threshold was not tuned on model output.
- **Limit:** This is "a subset matching the paper's description", not the paper's samples.

### 2026-11-15 02:20 — Two arrays removed from the labeled set after a label-free QC check
- **Found:** Deposited metadata say arrays 248–250 are the 31 °C triplicate of experiment E33 (Qp ≈ 4) and 251 begins the next experiment (Qp 27.4). In the expression matrix 249–251 are the tight triplicate (r = 0.998–0.999) and 248 sits with a different experiment (r = 0.96 with its supposed replicates). Metadata and expression are offset by one position here. No other labeled array shows this.
- **Decision:** Arrays 248 and 251 keep their expression in the unlabeled corpus but lose their Qp label (167 labeled arrays remain). E33 is treated as a single group because the culture assignment inside it is no longer certain. 251 is assigned to E33, 248 to its own block.
- **Reasoning:** The check uses replicate structure only, not model residuals, so it cannot favour any model. Training on a label that belongs to a different sample is wrong regardless of its effect.
- **Sensitivity:** the primary model is also run with all 169 labels exactly as deposited (`sensitivity_labels_as_deposited` in `results/metrics.csv`).
- **Alternatives:** Shift the labels to realign them (speculative); keep everything as deposited (trains on a label I believe is wrong).

### 2026-11-15 02:25 — Statistical conventions
- **R²:** reported as cross-validated R² = 1 − SSE(model) / SSE(training-fold mean). This is stricter than `r2_score` on pooled predictions, which centres on the test labels.
- **Confidence intervals:** percentile cluster bootstrap (2,000 resamples) of the out-of-fold predictions. The resampling unit is the split unit: cultures for culture-grouped CV, experiments for leave-one-experiment-out. Resampling arrays would treat replicates as independent. These intervals do not include the variance from refitting the model.
- **Classification target:** top vs bottom tercile of Qp over the labeled arrays; the middle tercile is dropped. Cut-offs are a label definition and are computed once on all labeled arrays.
- **Hyperparameters:** chosen by a 3-fold inner CV inside each training fold, grouped by the same unit as the outer split. Feature selection, scaling and encoder fitting are pipeline steps refit in every fold.
- **Encoders:** raw = 1,000 highest-variance probesets (variance computed in the training fold); pathway = 44 Hallmark scores (mean of z-scored member probesets, z-scoring fitted in the training fold); coexpr = 40 within-fold coexpression modules.
- **Primary model, named in advance:** PLS on the raw encoder for regression (mirrors the published method), logistic regression on the raw encoder for classification. The figures show these, whichever model scores best.
