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
- **Decision:** Commit on `poc/gse30321` as instructed. Do not create `main` by another route.
- **Alternatives:** Push the work to `main` directly (not authorised); make `poc/gse30321` the default branch and skip the PR.
- **Resolution (2026-11-15 05:25):** The PR is opened against a placeholder base, `claude/cool-curie-wthg7f`, which holds only the first commit (`.gitignore`) and is the one branch this session is cleared to push besides the working branch. The PR therefore shows the full diff, figures and ledger. **It is not a PR to `main`.** To finish: create `main` from commit `25d650a` and change the PR's base branch to `main` (GitHub: Edit → base), then delete the placeholder branch.

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

### 2026-11-15 02:40 — The published baseline is not beaten; reported as is
- **Found:** On the matched subset, culture-grouped nested-CV PLS gives RMSE 6.6 pg/cell/day, above the published 4.44. Held-out-experiment RMSE is 11.1.
- **Decision:** Report it plainly. No search over models, feature counts or subset definitions to get under 4.44. PLS on the raw encoder was named as the head-to-head model before the run; elastic net on the same subset is listed as a secondary row only.
- **Reasoning:** Any number reached by trying variants until one wins is not a head-to-head. The published figure also comes from a protocol I cannot reproduce exactly (different sample list, jackknife gene selection, cross-model validation whose handling of replicate cultures the abstract does not state).

### 2026-11-15 02:45 — Adversarial checks added to the main run
- **Process-covariate baseline:** Qp predicted from temperature, growth rate, cell density and viability alone (fixed small GBM, no expression). Answers "does the array add anything over the batch record?". Titre is excluded because Qp is computed from it.
- **Random-gene-set null for the Δ decode:** each pathway's |Δ| is compared with 1,000 random probe sets of the same size. Cell-line identity makes arbitrary gene sets differ between high- and low-Qp arrays, so a bootstrap interval that excludes zero is weak evidence on its own.
- **Within-experiment Δ:** the decode repeated after removing each experiment's mean from both the pathway scores and Qp, so it cannot be driven by which cell line an array comes from.
- **Repeated partitions for the leakage figure:** the array-level and culture-level bars are means over 10 random partitions, with a cluster bootstrap over cultures. A single partition moved RMSE by about ±0.7, which is larger than the array-vs-culture difference.

### 2026-11-15 03:00 — Leakage figure shows four split units
- **Decision:** `fig2_leakage.png` compares arrays, cultures, experiments and families of related experiments, not just grouped vs ungrouped.
- **Reasoning:** On these data, grouping by culture barely changes the estimate, because the leak runs through sister flasks of the same experiment. A two-bar figure would have shown "no leakage", which is the wrong lesson. The real step is at the experiment level.

### 2026-11-15 03:10 — §6 adapted to the labeled sample that exists
- **Found:** Only 62 labeled cultures exist, so train sizes of 20–100 cultures are not possible.
- **Decision:** Hold out 15 cultures; train on 10/20/30/40/47 cultures; 5 seeds. A second panel holds out 5 whole experiments and trains on 10/20/30/40 cultures. Everything on CPU.
- **In-fold corpus:** every array not belonging to a held-out group, labeled or not, labels never used. Asserted in code.
- **Stretch (pooling other public CHO data) skipped:** of the related accessions in the brief, only GSE37251 is CHO expression (30 arrays, a different array design, GPL15441). GSE296860, GSE179451 and E-MTAB-8060 are unrelated organisms or studies. Thirty arrays on another platform would add 10% to the corpus and a platform batch effect; not worth the harmonisation work.

### 2026-11-15 04:55 — Pretraining protocol revised once, after seeing the first run
- **First protocol (v1):** encoder with no normalisation on `z`, fine-tuned with weight decay 1e-4 and no dropout. Result: pretraining made test RMSE worse at every training size (for example 17.9 vs 12.4 at 10 cultures), and the transductive arm was worse than predicting the training mean. Files kept in `results/superseded/`.
- **Diagnosis (training-set quantities only):** both networks reach zero training error within about 40 epochs, and the pretrained encoder emits activations about 40 times larger than a fresh one (SD 4.4 vs 0.11). v1 was measuring an unstable fine-tune, not the quality of the representation.
- **Revision (v2), applied identically to every arm:** LayerNorm on `z`; input dropout 0.2 and weight decay 1e-2 during supervised training; one extra arm that freezes the in-fold encoder and fits ridge regression on `z` (a linear probe).
- **Commitment made before v2 results exist:** v2 is final. Whatever it shows is what gets reported, and v1 is disclosed next to it. No third protocol.
- **Risk I am accepting:** changing a protocol after seeing test results is a forking path. Disclosing both runs is the mitigation.

### 2026-11-15 05:20 — Pretraining result (v2, final): no reliable gain; reported as a negative
- **Result:** in-fold pretrained minus scratch = −0.61 RMSE with cultures held out (17 of 25 runs better; per-seed means −0.21, −0.01, +0.05, −0.64, −2.21) and −0.29 with experiments held out (11 of 20). PLS beats both networks at every size on held-out cultures.
- **Decision:** `fig3_learning_curve.png` is titled "no consistent gain". The title is computed from the paired differences in `src/figures.py`, not typed by hand.
- **Not done:** no further tuning of the autoencoder, no extra seeds chosen after the fact, no switch of the headline arm to the frozen-encoder probe (which did slightly better on held-out experiments).
- **Engineering note:** the first v2 launch stalled because strong weight decay pushed weights into denormal floats, which are about 100× slower on CPU. `torch.set_flush_denormal(True)` fixed it. Values printed before and after the fix were identical.

### 2026-11-15 05:25 — Figure choices
- `fig1_scatter.png` has two panels so that the ±4.44 band is never read against the wrong samples: left is the matched subset (the head-to-head), right is all labeled arrays, labelled as a different experiment.
- `fig1_roc.png` shows three curves (held-out cultures, experiments, families) because a single culture-grouped curve at AUC 0.97 would overstate what the model can do on a new cell line.
- Palette: three categorical colours validated for colour-vision deficiency, plus grey for the leaky split and the PLS reference. All type ≥ 14 pt.

### 2026-11-15 05:30 — Reproduction check
- **Decision:** before shipping, clone the branch into an empty directory and run `python cho_poc.py --all` there, then compare its results with the committed ones.
- **Outcome (05:47):** the run completed in about 25 minutes on 4 CPU cores (self-test, download, analysis, pretraining, figures). `metrics.csv`, `delta_pathways.csv`, `culture_map.csv`, `oof_predictions.csv` and `pretrain_runs.csv` matched the committed files with a maximum absolute numeric difference of 0.0. All five downloads matched their pinned SHA256. Only `data/PROVENANCE.md` differed, in its download timestamps.
- **Limit:** same machine and same package versions. Bit-for-bit agreement on other hardware or BLAS builds is not guaranteed, especially for the torch runs.

### 2026-11-15 05:28 — What is and is not committed
- **Committed:** code, figures, and result tables. `results/culture_map.csv` and `results/oof_predictions.csv` include each array's Qp, temperature and growth rate, copied from public GEO metadata, because a reviewer needs them to check the predictions and the grouping.
- **Not committed:** the expression matrix, the platform file, gene sets, any credential, key or host address. `raw/` and `data/*` (except `PROVENANCE.md`) are git-ignored; `git grep` for the instance address and key name returns nothing.

---

# Phase 2 — multi-modal extension (branch `poc/multimodal`)

Times in this section are read from the system clock (UTC) when the entry is written.

### 2026-11-15 10:04 UTC — P1: what Yao et al. 2021 actually published
- **Found:** Full text and supplement retrieved from Europe PMC (PMC8704136; MDPI and CORE returned HTTP 403). The data availability statement reads "Not applicable". The supplement is one PDF. It contains **no per-clone metabolite table** for either library.
- **What is public:** Table S3 — Pearson correlation of each annotated metabolite with qP and with peak VCD across the 12 clones of library 1 (day 7). Table S1 — the 34 targeted metabolites (names and MRM transitions only, no values). Figures: a 97-row heatmap at about 2 pixels per row, and scatter plots for two metabolites (aspartate, cystine) in both libraries.
- **Decision:** (a) Run the pathway-recovery test on Table S3, as a summary-statistic test. (b) **Do not run cross-library prediction.** It needs per-clone values for more than the two metabolites the authors chose to plot after seeing their own result; a "prediction" built from those two panels would be circular. (c) Do not digitise the heatmap: the resolution cannot support per-clone values and I will not present guessed numbers as data.
- **Consequence:** The cross-library box is marked "requires partner data". BMS holds these data.

### 2026-11-15 10:04 UTC — P1 PRE-REGISTRATION (written before any pathway statistic is computed)
- **Honest status of "blind":** I have read the paper's conclusion (alanine/aspartate/glutamate metabolism and the TCA cycle). The method below is fixed before computing, and it uses the authors' per-metabolite statistics. This is therefore a **concordance check of the Δ-style pathway ranking against the authors' own pathway analysis on the same summary data**, not an independent discovery and not a prediction.
- **Input:** every row of Table S3 (significant or not): Pearson r between day-7 metabolite level and qP.
- **Metabolite → KEGG compound:** case-insensitive exact match of the name against KEGG compound synonyms; unmatched names resolved by hand and listed in `results/yao_metabolite_map.csv`. Unresolvable names are dropped and counted.
- **Pathways:** KEGG reference metabolic pathways (map00010–map00999), from the KEGG compound–pathway links. Global and overview maps (map011xx, map012xx) excluded. Eligible if ≥ 3 Table S3 metabolites are members.
- **Statistic (Δ analogue for summary data):** mean Fisher-z of r(qP) over pathway members minus the mean over all other measured metabolites. Pathways ranked by this value, descending.
- **Null:** 10,000 permutations of metabolite labels. One-sided p per pathway.
- **Pre-specified targets:** map00250 (alanine, aspartate and glutamate metabolism) and map00020 (citrate cycle). Reported: rank of each among eligible pathways, per-pathway p, and a joint p for "mean rank of the two targets at least this good".
- **Secondary, also fixed now:** the same ranking with the over-representation statistic the authors used (members among metabolites with p(qP) < 0.05, hypergeometric), to separate "the data say so" from "our statistic says so".
- **Reported whatever the outcome.**
