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

### 2026-11-15 10:06 UTC — P1 RESULT (pre-registered test, reported as it came out)
- **Input:** 92 metabolites parsed from Table S3 (the paper states 97 unique features; 92 rows are in the table). All 92 mapped to KEGG compounds. 53 KEGG metabolic pathways had ≥ 3 measured members.
- **Citrate cycle:** rank **6 of 53** by the Δ statistic (permutation p = 0.020); rank **1 of 53** by the authors' own over-representation statistic.
- **Alanine, aspartate and glutamate metabolism:** rank **15 of 53** (p = 0.050); rank 7 by over-representation.
- **Joint test:** mean rank of the two targets 10.5 against 27.2 expected by chance, **p = 0.061**. This does not reach 0.05.
- **Reading:** the TCA cycle is recovered; alanine/aspartate/glutamate metabolism is in the top third, not at the top. The two targets do not "come out on top" under the pre-registered statistic. Ranks 1–5 went to small pathways (3–4 measured metabolites each: terpenoid backbone, butanoate, C5-branched dibasic acid, nicotinate, glutathione).
- **Why the ranking is flat:** 97% of the metabolites correlate positively with qP, so every pathway has a high mean r (0.38–0.85). A competitive statistic has little contrast to work with.
- **Post hoc sensitivity (not pre-registered):** restricting to the 38 pathways annotated in Chinese hamster gives TCA rank 6, alanine/aspartate/glutamate rank 13, joint p = 0.10. The pre-registered list included KEGG reference maps that do not exist in mammals (for example carbapenem biosynthesis); that was a flaw in my pre-registration and I am not using the sensitivity run to rescue the result.
- **Not done:** no cross-library prediction (see the entry above).

### 2026-11-15 10:06 UTC — P3 design: other endpoints on GSE30321
- **Endpoints:** cell viability, growth rate, cell density, q-lactate, q-ammonium, as deposited in GEO. Modality: transcriptomics only.
- **Splits:** leave-one-experiment-out is the headline; leave-one-expression-family-out is the strict check; culture-grouped 5-fold is shown only next to them. The leak assertion runs on every split.
- **Why the family split matters more here than for Qp:** many arrays labeled for these endpoints sit in experiment blocks whose culture assignment I flagged as ambiguous in phase 1. Experiment blocks are firmer than culture IDs, and the family grouping is derived from expression alone, so it does not depend on my reading of the metadata.
- **The two arrays with offset metadata (248, 251) are excluded from every endpoint**, since all of their deposited variables are suspect, not just Qp.
- **Extra reference:** a temperature-only linear model under the same split. Growth rate and viability move with the temperature shift, so "expression predicts growth rate" has to beat "the culture was at 31 °C".
- **Pipeline otherwise unchanged:** raw (1,000 probesets) and pathway (44 Hallmark) encoders; PLS, elastic net, gradient boosting; logistic regression for top-vs-bottom tercile; nested hyperparameters; 2,000-resample cluster bootstrap on the split unit.

### 2026-11-15 10:08 UTC — P2 PRE-REGISTRATION (written before any ground-truth list is loaded)
- **What I already know and cannot un-know:** the brief names FUT8, BAX, BAK, BOK, B4GALT1, BCKDHA/B, GLUL, LDHA, COSMC and NEU1/2/3. None of them was used to design the scores below. The Schulz et al. review and the Marzluf et al. screen have **not** been opened at the time of this entry.
- **Route A, data-driven (modality: transcriptomics, GSE30321):** gene-level Δ = mean standardised expression in the top Qp quartile minus the bottom quartile, computed **within experiment** (each experiment's mean removed from expression and from Qp; experiments with ≥ 4 labeled arrays). A knockout acts inside one cell line, so the within-experiment contrast is the relevant one. Probesets are averaged per gene symbol. Knockout candidates = most negative Δ. Overexpression candidates = most positive Δ. The raw (between-experiment) Δ is reported as a secondary column.
- **Route B, mechanistic (iCHO1766 from BiGG, COBRApy FBA, default medium bounds):** for every model gene, single-gene deletion gives g = maximum growth relative to wild type, and p = maximum IgG synthesis flux relative to wild type with growth held at ≥ 20% of wild-type maximum.
  - *Stated limit, before running:* an FBA optimum cannot increase when a gene is deleted, so FBA cannot predict benefit directly. Route B therefore yields (B1) a **tolerance call** — tolerated if g ≥ 0.9 and p ≥ 0.9 — and (B2) an exploratory **trade-off score** = 1 − g for deletions with p ≥ 0.99 and g ≥ 0.2 (growth capacity cut, product capacity intact), 0 otherwise.
- **Frozen artefacts:** `results/knockout_ranking.csv` (Route A, all array genes, ranked) and `results/knockout_fba.csv` (Route B, all 1,766 model genes) are committed before the ground truth is fetched. The commit hash is the timestamp.
- **Ground truth, to be collected afterwards:** (i) genes whose knockout, knockdown or inhibition raised titer or qP in CHO; (ii) genes whose overexpression did; (iii) reported neutral or harmful knockouts (for example BOK); (iv) the Marzluf screen's compatible vs incompatible knockouts. Source priority: the Schulz et al. systematic review if its target table can be obtained; otherwise open-access reviews, with the source of every gene recorded.
- **Metrics, fixed now:**
  - Route A, knockout end: among genes on the array, AUC of −Δ for set (i) against all other array genes; number of set-(i) genes in the 100 most negative; permutation p from 10,000 random gene sets of the same size. k = 50 reported as secondary.
  - Route A, overexpression end: the same with +Δ and set (ii).
  - Route B: share of set (i) called tolerated against the share for all model genes; for the Marzluf screen, AUC of g for compatible vs incompatible knockouts; AUC of the trade-off score for set (i).
  - Known negatives: their ranks are reported individually.
  - **Coverage is reported first.** The array carries about 1,700 genes and the model 1,766, so many validated targets will be untestable. Untestable is reported as untestable, not as a miss.
- **Expectation, stated in advance:** weak. Most famous knockouts act on glycosylation or apoptosis, which neither a Qp expression contrast nor steady-state FBA is built to see.

### 2026-11-15 10:16 UTC — P2 freeze and ground truth
- **Freeze:** commit `de3fe41` (2026-11-15 10:10 UTC) holds both rankings. The ground-truth files were downloaded after it (Europe PMC, 10:11 UTC).
- **Ground truth obtained:** Schulz et al. 2026, *mAbs* (PMC12818826), supplementary sheet `B_SingleTargetData`: 466 rows, each an engineered target with mechanism and maximum fold improvement. After removing miRNAs, lncRNAs, engineered variants and non-gene interventions, and resolving aliases by hand (`ALIAS` in `src/knockout_eval.py`): **60 genes whose knockout or knockdown raised titer or qP** (40 with single-target evidence), 10 knockout/knockdown genes with no benefit, **153 genes whose overexpression raised titer or qP**.
- **Marzluf et al. 2026 (PMC12902299):** the 36 retained knockouts are **not published by name** (they are one bar in Figure 2G). Only the genes the text names are usable: 8 growth-incompatible, 5 performance-reducing, 3 neutral or beneficial. This is a much weaker labeled set than the brief assumed.
- **BOK:** appears in the review only inside the BAX/BAK/BOK triple knockout. Treated as a named negative, not a positive, per the brief.
- **Rows with mixed mechanisms** (knockout plus overexpression in one construct) are excluded from both sets because the direction per gene is not machine-readable.

### 2026-11-15 10:16 UTC — P2 RESULT: knockout recovery is at chance
- **Coverage first:** of the 60 validated knockout/knockdown targets, **14 are on the array** and **17 are in iCHO1766**. BAX, BAK1, BOK, BCKDHA/B, COSMC, NEU1 and NEU3 are not on the array; BAX, BAK1, BOK, COSMC and NEU3 are not in the metabolic model. BOK is untestable by either route.
- **Route A (transcriptomics Δ, within experiment), knockout end:** AUC **0.50** (permutation p = 0.50; null 95% range 0.35–0.65). **2 of 14** validated targets in the 100 most negative genes against 0.8 expected (p = 0.19). The two are GLUL (rank 29 of 1,707) and DAG1 (rank 30). With k = 50 the same two hits give p = 0.06. FUT8 ranks 1,197, LDHA 182, MYC 1,700.
- **Route A, overexpression end:** AUC **0.39** (p = 0.99); 2 of 37 in the top 100 against 2.2 expected. If anything the sign is the wrong way round.
- **Route A with the raw (between-experiment) Δ:** AUC 0.39, 0 of 14 in the top 100.
- **Route B (iCHO1766 FBA):** uninformative. 93% of all model genes are called tolerated, and 94% of the validated targets. The trade-off score has AUC 0.49 and no validated target in its top 50. FBA calls GAPDH, ALDOA, PGK1 and LDHA knockouts tolerated (growth ratio 1.0) where the Marzluf screen found them lethal or selected against; the model's isozyme rules make single deletions harmless. One mechanistically sensible call: FUT8 deletion drops the model's IgG synthesis capacity to zero, because the model's IgG carries a core-fucosylated glycan.
- **Reading:** neither route recovers validated knockout targets above chance. With 14 testable positives the test could only have detected a strong effect, so this is "no evidence of recovery", not "evidence of no ability". The two top-30 hits are not enough to claim anything.
- **Not done:** no re-scoring, no alternative statistics tried after seeing this.

### 2026-11-15 10:16 UTC — P4 PRE-REGISTRATION: cross-modality alignment on a pathway basis (written before any cross-modality statistic is computed)
- **Claim under test:** independent modalities, projected onto a common pathway basis, give concordant pathway-level contrasts for the same biological contrast.
- **Leg 1 — growth phase, transcriptomics vs proteomics (the main test).**
  - *Proteomics:* PRIDE PXD053511 (CHO-K1GS, label-free, day 4 exponential n = 4 vs day 7 stationary n = 4, 37 °C throughout). Deposited `AllQuantifiedProteinGroups.tsv`. Proteins quantified in ≥ 3 of 4 samples in each group; log2 intensity, median-centred per sample; statistic = Welch t (day 4 − day 7). Positive = higher in exponential phase.
  - *Transcriptomics:* GSE30321. Gene-level Δ for **growth rate, within experiment**: top quartile minus bottom quartile of experiment-centred growth rate (experiments with ≥ 4 arrays). Positive = higher when the culture grows faster.
  - *Known mismatch, stated now:* in GSE30321 slow growth mostly coincides with the temperature shift to 31 °C; the proteomics study has no temperature shift. A 37 °C-only variant is reported as a sensitivity run. It has only 44 arrays in 5 experiments with a usable growth range, which is why it is not primary.
  - *Protein → gene:* the deposited gene column where it is a real symbol; otherwise exact match of the protein name to UniProt names (hamster, then mouse and human reviewed entries). Mapping rate reported.
- **Leg 2 — Qp association, transcriptomics vs metabolomics.** Yao Table S3 (r with qP across 12 clones) against the GSE30321 gene Δ for Qp. The **raw, between-experiment Δ** is used because Yao's contrast is between clones. Expected to be weak: 97% of Yao's correlations are positive.
- **Pathway bases:** (a) **iCHO1766 subsystems** (primary, as the architecture specifies): genes through the model's gene–reaction rules; metabolites through the reactions they take part in, currency metabolites excluded. (b) **KEGG pathways annotated for Chinese hamster** (secondary; includes non-metabolic processes the metabolic model cannot represent). A pathway is used if it has ≥ 5 measured genes in each gene-level modality, or ≥ 3 metabolites for metabolomics.
- **Pathway score:** mean over member features of the modality's statistic, after converting that statistic to a z-score across all features of the modality.
- **Test statistic:** Spearman correlation between the two modalities' pathway vectors. **P-value:** 10,000 permutations of feature labels within one modality (protein labels for leg 1, metabolite labels for leg 2), pathway scores recomputed each time; one-sided for a positive correlation.
- **Also reported:** the gene-level Spearman correlation for genes measured in both modalities in leg 1. If genes already agree, the pathway basis is not what creates the agreement, and I will say so.
- **Not run, with reason:** the hamster proteome atlas (tissues and cell lines, no contrast shared with any other dataset here); Metabolomics Workbench (one Chinese hamster study, ST001177, not examined for lack of time); any three-way fusion (no sample has more than one modality).
- **A null result is reported as a null result.**

### 2026-11-15 10:18 UTC — P4 RESULT: pathway-level alignment is not demonstrated by the pre-registered test
- **Proteomics processing:** 6,959 protein groups; 4,207 quantified in ≥ 3 of 4 samples per group; 87% mapped to a gene symbol (1,488 from the deposited gene column, 2,178 by protein-name match; UniProt has retired most of the accessions the deposit used); 3,303 unique genes. 557 genes are shared with the array.
- **Leg 1, pre-registered primary (all temperatures, within experiment; 205 arrays, 22 experiments):**
  - iCHO1766 subsystems: Spearman **0.20** over only **18** testable subsystems, p = **0.24**.
  - KEGG hamster pathways: Spearman **0.16** over 266 pathways, p = **0.13**.
  - Neither is significant.
- **Gene level, no pathway basis:** Spearman **0.18** over 557 shared genes, p < 0.0001. Transcript and protein contrasts do agree, weakly, gene by gene. The pathway projection did not strengthen that agreement, so the pathway basis cannot be credited with it.
- **Leg 1, sensitivity (37 °C arrays only; 69 arrays, 11 experiments):** KEGG Spearman **0.34**, p < 0.0001; iCHO1766 Spearman 0.36, p = 0.08. This variant matches the proteomics design better (no temperature shift) and I said so before running it, but it was registered as a sensitivity run and is reported as one.
- **Leg 2 (Qp association, transcriptomics vs Yao metabolomics):** Spearman −0.12 (iCHO1766, 10 subsystems) and −0.08 (KEGG, 40 pathways); p = 0.63 and 0.66. Null, as expected from the lack of contrast in Yao's table.
- **Reading:** the primary test is a null. There is a real but weak transcript–protein agreement for the growth-phase contrast, visible at gene level and, once the temperature confound is removed, at KEGG pathway level. The iCHO1766 basis is too sparse on this array to test: 18 subsystems.
- **Why iCHO1766 is sparse here:** it covers metabolism only (1,766 genes) and the array carries 1,707 genes, so few subsystems reach 5 genes in both modalities.

### 2026-11-15 10:21 UTC — P5: what exists for scale risk, and a pre-registration
- **GSE30321 has no scale or format variable.** The deposit carries nine per-sample variables (growth rate, temperature, q-ammonia, q-lactate, viability, seeding density, cell density, titre, Qp). The summary says samples span "a wide range of process scales and formats", but that information was not deposited. A small-scale-to-large-scale degradation test cannot be run on it.
- **Gao et al. 2016** (*Biotechnol J* 11:1190; hypoxia on scale-up to 5,000 L; authors include BMS staff) is closed access with no deposited data that I could find. Not usable.
- **Found instead:** GSE104787 — one IgG-producing CHO-K1 line grown from the same pre-culture in 3 uncontrolled shake flasks and 3 controlled bioreactors, arrays on days 2, 4 and 5 (18 arrays, genome-wide Affymetrix CHO array, 16,485 gene symbols). This is a **culture-format contrast, not a scale contrast**, and it has 6 culture units, so no predictive metric is possible.
- **PRE-REGISTRATION (before computing):**
  - Δ per gene = mean z-scored expression in bioreactor arrays minus flask arrays, after removing each day's mean. Hallmark pathway score = mean Δ of member genes (≥ 10 genes on the array). Pathways ranked by |score|. Null: 1,000 random gene sets of the same size.
  - **Target 1:** HALLMARK_HYPOXIA. I know from its title that Gao et al. report hypoxia at 5,000 L. Whether oxygen biology also separates a flask from a benchtop bioreactor is a different question; a hit would be suggestive only, and a miss says nothing about Gao et al.
  - **Target 2:** whatever the dataset's own publication concluded. I have not read it. The Hallmark ranking is committed first, then the abstract is fetched and compared.
  - The scale-risk box stays "requires partner data" whatever this shows.

### 2026-11-15 10:22 UTC — P5 RESULT (frozen in commit `0785634` before the dataset's publication was looked up)
- **Hypoxia:** HALLMARK_HYPOXIA ranks **17 of 50** in the bioreactor-vs-flask contrast (Δ = +0.14, p = 0.012 against random gene sets). It is not near the top. Sixteen of 50 sets pass p < 0.05, so passing means little here.
- **What does lead:** proliferation programmes are higher in the bioreactor (E2F targets, MYC targets, G2M checkpoint, DNA repair) together with the unfolded protein response; protein secretion is lower.
- **Target 2:** I could not find a publication for this deposit, and its GEO summary names no pathway. The one checkable statement in it — large differences over time, much smaller ones between systems — holds: sampling day explains a median 24.6% of expression variance, culture system 2.1%.
- **Reading:** no support here for hypoxia as the main difference between a shake flask and a benchtop bioreactor, which is not what Gao et al. claimed anyway (they compared 5,000 L with small scale). The scale-risk box stays "requires partner data".
- **Lead not followed for lack of time:** Alsayyari et al. 2018 (*Biotechnol Bioeng*), transcriptomes of a 10 L process and its ambr 15 scale-down. That is a true scale contrast and may be on GEO.

### 2026-11-15 10:22 UTC — P6: glycans, stopped inside the timebox
- **Searched:** GEO (Chinese hamster series mentioning glycan, glycosylation or glycoform: 3 hits, none with glycan profiles) and Europe PMC (open-access CHO papers with glycan profiles and omics on the same clones: glycoengineering panels and vaccine glycoprofiling, no clone-level omics-plus-glycan matrix).
- **Decision:** the glycan head is marked "requires partner data". The FUT8 and B4GALT1 knockout literature gives mechanism, not a prediction set. One mechanistic note from P2: deleting FUT8 in iCHO1766 removes the model's capacity to make its (core-fucosylated) IgG, which is the model encoding known biology, not a prediction.

### 2026-11-15 16:50 UTC — P3 RESULT: other endpoints (transcriptomics, GSE30321)
- **Viability:** high vs low AUC 0.84 (0.70–0.94) with experiments held out, 0.77 (0.60–0.90) with families held out. Regression has no skill: R² 0.09 (−0.84–0.26). Temperature alone: R² 0.01.
- **Growth rate:** AUC 0.98, R² 0.76 (0.68–0.83). Temperature alone already gives R² 0.62.
- **Viable cell density:** AUC 0.82 (0.68–0.96), R² 0.08. Temperature alone 0.17.
- **q-lactate:** AUC 0.97, R² 0.62. Temperature alone 0.61. Expression adds nothing measurable.
- **q-ammonium:** AUC 0.93, R² 0.42. Temperature alone 0.35.
- **Leak size:** culture-grouped R² for viability is 0.77 and for density 0.77, against 0.09 and 0.08 with experiments held out. The culture-level split overstates these endpoints even more than it did Qp.
- **Decision:** the viability box is marked positive on the strength of the high-vs-low AUC only, and the box says in the same breath that regression has no skill.
- **Added after seeing the regression table:** a temperature-only reference for the AUCs as well, since "temperature alone explains the R²" raised the same question for the classification numbers. It is a reference row, not a model choice.

### 2026-11-15 16:50 UTC — fig0 design choices
- **Three statuses, not two.** "Tested: positive evidence", "tested: weak or null", "requires partner data". A two-state figure (metric vs grey) would have forced null results into the same colour as successes.
- **A box is green only if** its headline metric has an interval that excludes chance under experiment-level holdout. That gives two green boxes plus the transcriptomics input that feeds them.
- **Status is carried by a text tag and the border style** as well as colour, so it survives a greyscale print.
- **Every number on the figure is read from a results file** by `src/figures_v2.py`.
- **Dark background** to match the deck the diagram came from; phase-1 figures stay light.
- **"iMEMORY initialisation: not tested"** is written into the culture-state box. The brief's diagram names it; nothing public lets me test it.

### 2026-11-15 16:50 UTC — Reproducibility of the frozen P2 files
- Re-running `src/knockout.py freeze` and `src/knockout_eval.py` from the committed scripts reproduces the frozen ranking file exactly (all 1,707 rows), the FBA ratios with zero numeric difference, and an identical evaluation JSON.
- One fix made along the way: the gene-symbol fill-in for model genes named only by a LOC identifier had been applied by hand before the freeze; it now lives in `src/knockout.py`. Symbols are identifiers, not scores, and the regenerated table matches the frozen one on every symbol.

### 2026-11-15 16:50 UTC — What was cut or not attempted
- Cross-library prediction on the Yao data (no per-clone data).
- The hamster proteome atlas and Metabolomics Workbench study ST001177 (no shared contrast / no time).
- The Alsayyari et al. 2018 scale-down transcriptomes (a real scale contrast; found late).
- Any attempt to improve the knockout or alignment results after seeing them.
- iCHO1766 as a predictive encoder: only 18 subsystems reach five genes on this array.
- **Temperature-only AUC, now computed:** viability 0.44 (0.25–0.69), growth rate 0.90, viable cell density 0.75, q-lactate 0.94, q-ammonium 0.80, against 0.84, 0.98, 0.82, 0.97 and 0.93 from expression. Viability and q-ammonium are the endpoints where expression clearly adds to temperature. The rerun reproduced every earlier number.

### 2026-11-15 17:24 UTC — Phase-2 reproduction check
- **What was run:** a clean clone of `poc/multimodal`, then `python cho_poc.py --download && python src/run_v2.py`. About 20 minutes on 4 CPU cores.
- **Outcome:** `metrics_v2.csv`, `metrics_endpoints.csv`, `knockout_ranking.csv`, `yao_pathway_ranking.csv`, `alignment_pathway_scores.csv`, `scale_format_pathways.csv`, `knockout_ground_truth.csv` and all four result JSON files matched the committed versions with zero numeric difference. Phase-1 `metrics.csv` was untouched.
- **One difference found and fixed:** `knockout_fba.csv` had the same values in a different row order, because LP solutions carry float noise near 1e-12 and tied genes sorted differently. Ratios are now rounded to six decimals before sorting; two consecutive runs give byte-identical files, and every value still equals the frozen commit `de3fe41`.
- **Limit:** same machine, same package versions. The public sources (Europe PMC, PRIDE, KEGG, BiGG, UniProt, GEO) are fetched live and only the phase-1 files and iCHO1766 have pinned hashes; KEGG and UniProt content can change.
