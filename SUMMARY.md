# Summary of results and findings

The short version. Full detail and the claims ledger are in `REPORT.md`; every choice made along the way is in `DECISIONS.md`.

RMSE is in pg/cell/day. Brackets are 95% bootstrap intervals over cultures or experiments.

## Headline results

| | Held-out cultures | Held-out experiments |
|---|---|---|
| Matched subset vs the published 4.44 (75 arrays) | 6.6 (5.3–7.9) | 11.1 (7.1–15.2) |
| All labeled arrays, RMSE (167 arrays) | 5.8 (4.8–6.7) | 10.8 (8.1–13.3) |
| All labeled arrays, R² | 0.83 | 0.45 |
| Top vs bottom Qp tercile, AUC | 0.97 (0.93–1.00) | 0.92 (0.80–0.99) |

A model that predicts the training mean scores about 14. Process variables alone (temperature, growth rate, cell density, viability) score 13.5.

## Findings

1. **The published baseline is not beaten.** On the matched subset (stationary, temperature-shifted), culture-grouped PLS scores 6.6, and its interval excludes 4.44. No variants were tried to close the gap, so it stays a single pre-specified test.
2. **Expression predicts Qp well for a held-out culture.** But sister flasks from the same experiment are in the training set.
3. **Most of that is recognising the experiment, not predicting a new one.** Holding out whole experiments nearly doubles the error. Holding out families of related experiments gives AUC 0.71 (0.44–0.89), an interval that includes chance.
4. **Grouping by culture does not remove the leak in this dataset.** Splitting by array gives 5.7 and by culture 6.0, the same answer. The real step is at the experiment level, because replicate flasks of one experiment are near-copies.
5. **The pathway encoder loses to raw genes**: 10.1 vs 5.8 on held-out cultures. The array covers only 15% of Hallmark genes, so this says little about a pathway encoder on a full transcriptome.
6. **The Δ decode recovers one of the three WGCNA themes cleanly.** Golgi vesicle transport ranks 5 of 48 (p = 0.010 against random gene sets) and holds within experiments. Cell cycle is suggestive (p = 0.05–0.07). Protein secretion is not recovered. The sign is negative (high-Qp arrays score lower), and the contrast is confounded with cell line and temperature.
7. **Self-supervised pretraining gives no reliable gain.** Pretraining inside the fold beats training from scratch by 0.6 with cultures held out and 0.3 with experiments held out. Both are within seed-to-seed noise, and plain PLS beats both networks. The protocol was revised once after an unstable first run; both runs are disclosed.

## Data issues found

- **No culture ID is public.** It was derived from deposit order and seeding density: 122 cultures against the 121 GEO states.
- **Two arrays carry labels that contradict their expression** (arrays 248 and 251). They were removed from the labeled set, leaving 167 of 169. A sensitivity run with them put back is included.
- **The array is small**: 3,592 probesets, not about 20,000.
- **The labeled set is small.** Only 62 cultures in 19 experiments have Qp. The learning curves therefore use 10–47 training cultures.

## Deviations from the brief

- **The repository was empty.** There was no `cho_poc.py` and no `main` branch, so the pipeline and its self-test were written from scratch.
- **No GPU instance was used.** This session had no SSH client and no access to the key, so everything ran on CPU.
- **The PR does not target `main`.** It targets a placeholder branch holding only the `.gitignore` commit. To fix: create `main` from commit `25d650a`, then change the PR's base to `main`.

## Reproducibility

```bash
pip install -r requirements.txt && python cho_poc.py --all
```

The run takes about 25 to 45 minutes on 4 CPU cores. A clean clone reproduced every committed result file exactly.

## What this does not demonstrate

- **Not tested:** multi-omic fusion, latent dynamics or scale-up prediction. The data to test them are not in this dataset.
- **No prediction for a new cell line.** Strict splits give an error of about 11 pg/cell/day and intervals that reach chance.
- **No improvement on the 2011 published model.**
