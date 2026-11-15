# CHO proof of concept: Qp from transcriptomics on GSE30321

Predicts cell-specific productivity (Qp) of CHO production cultures from 2011-era microarray
expression, with leakage-safe cross-validation, and compares against the published baseline
(Clarke et al. 2011, J Biotechnol 151:159).

**Read `REPORT.md` first.** It holds the results, the claims ledger and the list of things
this work does not show. `DECISIONS.md` logs every choice made along the way.

## Reproduce everything

```bash
pip install -r requirements.txt && python cho_poc.py --all
```

That one command runs, in order: the offline self-test, the downloads (GEO, ArrayExpress,
MSigDB; about 11 MB into `raw/`), the full analysis, the pretraining experiment and the figures.
Expect about 45 minutes on 4 CPU cores; no GPU is needed. Seeds are fixed.

Individual stages:

| command | what it does | time |
|---|---|---|
| `python cho_poc.py --self-test` | synthetic checks of the parser, the no-leak guard, bootstrap, encoders, Δ decode | 10 s |
| `python cho_poc.py --download` | fetches sources, checks SHA256 against pinned values, writes `data/PROVENANCE.md` | 20 s |
| `python cho_poc.py --run` | grouped-CV metrics, matched baseline, leakage comparison, Δ decode | 25 min |
| `python cho_poc.py --pretrain` | self-supervised pretraining learning curves (needs `torch`) | 20 min |
| `python cho_poc.py --figures` | slide figures from `results/` | 10 s |

Add `--quick` for a smoke run with a smaller model grid.

## Outputs

| file | content |
|---|---|
| `results/metrics.csv` | every reported metric with its 95% CI, split scheme, encoder, model and sample counts |
| `results/delta_pathways.csv` | Δ decode: pathways ranked by top-quartile minus bottom-quartile Qp, with driver genes |
| `results/culture_map.csv` | derived culture and experiment ID for each of the 295 arrays |
| `results/oof_predictions.csv` | out-of-fold predictions behind every metric |
| `results/pretrain_*.csv` | pretraining runs, per-arm summary, paired differences |
| `figures/fig1_roc.png`, `fig1_scatter.png` | the two slide figures |
| `figures/fig2_leakage.png` | the same model under four hold-out units |
| `figures/fig3_learning_curve.png` | pretrained vs from scratch |
| `data/PROVENANCE.md` | URLs, SHA256, timestamps, row and column counts |

## Layout

- `cho_poc.py` — parser, culture derivation, encoders, nested grouped CV, bootstrap, Δ decode, CLI
- `src/pretrain.py` — masked-feature denoising autoencoder experiment
- `src/figures.py` — plotting only; no model is fitted there
- `raw/`, `data/` — downloads; git-ignored except `data/PROVENANCE.md`

## Guarantees enforced in code

- Every grouped split passes `assert_no_group_leak`; the assertion stays in the committed code.
- Scaling, feature selection, encoders and hyperparameters are fitted inside the training fold.
- The in-fold pretraining corpus is asserted to contain no array from a held-out group.

## Scope

Single modality: transcriptomics on a 3,592-probeset custom array. No proteomics, metabolomics,
glycan data or paired small/large-scale runs. Nothing here tests multi-omic fusion, latent
dynamics or scale-up prediction.
