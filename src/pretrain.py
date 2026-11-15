"""Self-supervised pretraining experiment (brief section 6).

Question: does masked-feature denoising pretraining on unlabeled CHO expression
improve small-sample Qp prediction, compared with the same network trained from
scratch?

No external weights. Five arms, same test arrays in every arm:

    scratch            encoder + linear head, random init, labeled training arrays only
    pretrained_infold  encoder initialised from a denoising autoencoder trained on
                       every array that is NOT in a held-out group (labeled or not,
                       labels never used), then fine-tuned. This is the headline arm.
    pretrained_infold_frozen_ridge
                       the same in-fold encoder, frozen, with ridge regression on z
                       (a linear probe: tests the representation without fine-tuning)
    pretrained_transductive
                       encoder pretrained on all 295 arrays, including the held-out
                       arrays' expression (never their labels). Upper bound only.
    pls                PLS on the 1,000 highest-variance probesets, the classical
                       reference, components chosen inside the training set

Two hold-out units, reported separately:
    culture     15 cultures held out (the unit the brief specifies)
    experiment  5 whole experiment blocks held out (stricter)

Everything that is tuned (number of fine-tuning epochs, PLS components) is tuned
on a grouped split inside the training set. Hyperparameters of the autoencoder
are fixed in advance and were not tuned.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import cho_poc as cp  # noqa: E402

RESULTS = ROOT / "results"

# fixed in advance, not tuned
HIDDEN = 256
Z_DIM = 32
MASK_FRAC = 0.25
PRE_EPOCHS = 200
PRE_LR = 1e-3
FT_LR = 5e-4
WEIGHT_DECAY = 1e-4       # autoencoder
FT_WEIGHT_DECAY = 1e-2    # supervised training, both arms
INPUT_DROPOUT = 0.2       # supervised training, both arms
FT_MAX_EPOCHS = 200
FT_CHECK_EVERY = 20
BATCH = 32
PRE_BATCH = 64
N_SEEDS = 5
TRAIN_SIZES = {"culture": [10, 20, 30, 40, 47], "experiment": [10, 20, 30, 40]}
N_TEST = {"culture": 15, "experiment": 5}


def _torch():
    import torch

    torch.set_num_threads(4)
    return torch


def build_encoder(n_in, torch):
    nn = torch.nn
    # LayerNorm on z keeps the latent on one scale whether the encoder was pretrained or not;
    # without it a pretrained encoder hands the head activations ~40x larger than a fresh one.
    return nn.Sequential(nn.Linear(n_in, HIDDEN), nn.GELU(), nn.Linear(HIDDEN, Z_DIM), nn.LayerNorm(Z_DIM))


def pretrain_dae(Xs, seed, epochs=PRE_EPOCHS):
    """Masked-feature denoising autoencoder. Xs is standardised expression (numpy).

    A random 25% of features per sample is replaced by 0 (the feature mean after
    standardisation); the loss is the reconstruction error on the masked entries.
    Returns the encoder state dict and the final masked-reconstruction loss.
    """
    torch = _torch()
    torch.manual_seed(seed)
    g = torch.Generator().manual_seed(seed)
    X = torch.tensor(Xs, dtype=torch.float32)
    enc = build_encoder(X.shape[1], torch)
    dec = torch.nn.Sequential(torch.nn.Linear(Z_DIM, HIDDEN), torch.nn.GELU(), torch.nn.Linear(HIDDEN, X.shape[1]))
    opt = torch.optim.AdamW(list(enc.parameters()) + list(dec.parameters()), lr=PRE_LR, weight_decay=WEIGHT_DECAY)
    n = X.shape[0]
    last = float("nan")
    for _ in range(epochs):
        perm = torch.randperm(n, generator=g)
        tot, cnt = 0.0, 0
        for i in range(0, n, PRE_BATCH):
            xb = X[perm[i:i + PRE_BATCH]]
            mask = torch.rand(xb.shape, generator=g) < MASK_FRAC
            out = dec(enc(xb.masked_fill(mask, 0.0)))
            loss = ((out - xb) ** 2)[mask].mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item() * len(xb)
            cnt += len(xb)
        last = tot / cnt
    return {k: v.clone() for k, v in enc.state_dict().items()}, last


def _fit_supervised(Xtr, ytr, seed, init_state, epochs, eval_sets=None):
    """Train encoder + linear head for `epochs`. Optionally record RMSE on eval_sets every FT_CHECK_EVERY epochs."""
    torch = _torch()
    torch.manual_seed(seed)
    g = torch.Generator().manual_seed(seed)
    X = torch.tensor(Xtr, dtype=torch.float32)
    y = torch.tensor(ytr, dtype=torch.float32)
    enc = build_encoder(X.shape[1], torch)
    if init_state is not None:
        enc.load_state_dict(init_state)
    head = torch.nn.Linear(Z_DIM, 1)
    model = torch.nn.Sequential(torch.nn.Dropout(INPUT_DROPOUT), enc, head)
    opt = torch.optim.AdamW(model.parameters(), lr=FT_LR, weight_decay=FT_WEIGHT_DECAY)
    n = X.shape[0]
    curve = {}
    for ep in range(1, epochs + 1):
        model.train()
        perm = torch.randperm(n, generator=g)
        for i in range(0, n, BATCH):
            b = perm[i:i + BATCH]
            loss = ((model(X[b]).squeeze(-1) - y[b]) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
        if eval_sets is not None and ep % FT_CHECK_EVERY == 0:
            model.eval()
            with torch.no_grad():
                curve[ep] = [float(((model(torch.tensor(a, dtype=torch.float32)).squeeze(-1)
                                     - torch.tensor(b, dtype=torch.float32)) ** 2).mean().sqrt())
                             for a, b in eval_sets]
    return model, curve


def fit_predict_nn(Xtr, ytr, gtr, Xte, seed, init_state):
    """Choose the number of epochs on a grouped split inside the training set,
    refit on the whole training set, predict the test arrays."""
    torch = _torch()
    mu, sd = ytr.mean(), ytr.std() + 1e-8
    yz = (ytr - mu) / sd
    rng = np.random.default_rng(seed)
    ug = rng.permutation(np.unique(gtr))
    n_val = max(2, len(ug) // 4)
    val = np.isin(gtr, ug[:n_val])
    cp.assert_no_group_leak(gtr, np.where(~val)[0], np.where(val)[0])
    _, curve = _fit_supervised(Xtr[~val], yz[~val], seed, init_state, FT_MAX_EPOCHS, eval_sets=[(Xtr[val], yz[val])])
    best_ep = min(curve, key=lambda e: curve[e][0])
    model, _ = _fit_supervised(Xtr, yz, seed, init_state, best_ep)
    model.eval()
    with torch.no_grad():
        return model(torch.tensor(Xte, dtype=torch.float32)).squeeze(-1).numpy() * sd + mu, best_ep


def fit_predict_frozen(Xtr, ytr, gtr, Xte, seed, init_state):
    """Linear probe: frozen pretrained encoder, ridge regression on z (alpha chosen on grouped inner folds)."""
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import GridSearchCV
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    torch = _torch()
    enc = build_encoder(Xtr.shape[1], torch)
    enc.load_state_dict(init_state)
    enc.eval()
    with torch.no_grad():
        ztr = enc(torch.tensor(Xtr, dtype=torch.float32)).numpy()
        zte = enc(torch.tensor(Xte, dtype=torch.float32)).numpy()
    gs = GridSearchCV(Pipeline([("s", StandardScaler()), ("m", Ridge())]), {"m__alpha": [0.01, 0.1, 1, 10, 100, 1000]},
                      scoring="neg_root_mean_squared_error", cv=cp._inner_cv(gtr, "grouped", seed))
    return gs.fit(ztr, ytr).predict(zte)


def fit_predict_pls(Xtr, ytr, gtr, Xte, seed):
    from sklearn.model_selection import GridSearchCV

    pipe, grid, scoring = cp.make_model("pls", cp.make_encoder("raw"), seed)
    grid = {"m__n_components": [c for c in grid["m__n_components"] if c < min(len(ytr) * 2 // 3 - 1, 11)]}
    inner = cp._inner_cv(gtr, "grouped", seed)
    return np.ravel(GridSearchCV(pipe, grid, scoring=scoring, cv=inner).fit(Xtr, ytr).predict(Xte))


def standardise(X, ref):
    mu = ref.mean(axis=0)
    sd = ref.std(axis=0)
    return (X - mu) / np.where(sd > 1e-8, sd, 1.0)


def run_unit(D, unit, n_seeds, sizes, log):
    meta = D["meta"]
    X = D["X"]
    y = meta[cp.ENDPOINT].to_numpy(dtype=float)
    labeled = ~np.isnan(y) & ~meta["qc_label_excluded"].to_numpy()
    culture = meta["culture_id"].to_numpy()
    holdout_group = culture if unit == "culture" else meta["experiment_id"].to_numpy()
    rows = []
    for seed in range(n_seeds):
        s = cp.SEED + seed
        rng = np.random.default_rng(s)
        test_groups = rng.permutation(np.unique(holdout_group[labeled]))[: N_TEST[unit]]
        in_test_group = np.isin(holdout_group, test_groups)
        te = np.where(labeled & in_test_group)[0]
        pool = np.where(labeled & ~in_test_group)[0]
        pool_cultures = rng.permutation(np.unique(culture[pool]))
        # in-fold pretraining corpus: every array outside the held-out groups, labeled or not
        corpus = np.where(~in_test_group)[0]
        assert not (set(holdout_group[corpus]) & set(test_groups)), "held-out group leaked into pretraining corpus"
        assert not (set(corpus) & set(te)), "held-out array leaked into pretraining corpus"
        t0 = time.time()
        state_in, loss_in = pretrain_dae(standardise(X[corpus], X[corpus]), s)
        state_tr, loss_tr = pretrain_dae(standardise(X, X), s)
        log(f"  [{unit}] seed {seed}: test arrays={len(te)}, corpus in-fold={len(corpus)} / transductive={len(X)}, "
            f"DAE masked-MSE in-fold={loss_in:.3f} transductive={loss_tr:.3f} ({time.time() - t0:.0f}s)")
        for n_cult in sizes:
            if n_cult > len(pool_cultures):
                continue
            tr = pool[np.isin(culture[pool], pool_cultures[:n_cult])]
            cp.assert_no_group_leak(holdout_group, tr, te)
            cp.assert_no_group_leak(culture, tr, te)
            ytr, yte, gtr = y[tr], y[te], culture[tr]
            preds = {}
            # scratch: sees only the labeled training arrays, so it is standardised on them
            preds["scratch"], ep_s = fit_predict_nn(standardise(X[tr], X[tr]), ytr, gtr,
                                                    standardise(X[te], X[tr]), s, None)
            preds["pretrained_infold"], ep_i = fit_predict_nn(standardise(X[tr], X[corpus]), ytr, gtr,
                                                              standardise(X[te], X[corpus]), s, state_in)
            preds["pretrained_infold_frozen_ridge"] = fit_predict_frozen(
                standardise(X[tr], X[corpus]), ytr, gtr, standardise(X[te], X[corpus]), s, state_in)
            preds["pretrained_transductive"], ep_t = fit_predict_nn(standardise(X[tr], X), ytr, gtr,
                                                                    standardise(X[te], X), s, state_tr)
            preds["pls"] = fit_predict_pls(X[tr], ytr, gtr, X[te], s)
            preds["train_mean"] = np.full(len(te), ytr.mean())
            for arm, p in preds.items():
                rows.append({"holdout_unit": unit, "seed": seed, "n_train_cultures": n_cult, "n_train_arrays": len(tr),
                             "n_test_arrays": len(te), "arm": arm, "rmse": cp.rmse(yte, p)})
            log(f"    n={n_cult:>2} cultures ({len(tr):>3} arrays): " + "  ".join(
                f"{a}={cp.rmse(yte, p):.2f}" for a, p in preds.items()) + f"  epochs s/i/t={ep_s}/{ep_i}/{ep_t}")
    return rows


def main(quick: bool = False):
    try:
        torch = _torch()
    except ImportError as exc:
        print(f"pretraining skipped: {exc}")
        return None
    torch.use_deterministic_algorithms(True)
    print(f"pretraining experiment (torch {torch.__version__}, CPU)")
    D = cp.prepare()
    n_seeds = 2 if quick else N_SEEDS
    rows = []
    for unit in ["culture", "experiment"]:
        sizes = TRAIN_SIZES[unit][:2] if quick else TRAIN_SIZES[unit]
        rows += run_unit(D, unit, n_seeds, sizes, print)
        pd.DataFrame(rows).to_csv(RESULTS / "pretrain_runs.csv", index=False)
    runs = pd.DataFrame(rows)
    summ = (runs.groupby(["holdout_unit", "n_train_cultures", "arm"])["rmse"]
            .agg(mean="mean", sd="std", n_seeds="count").reset_index())
    # paired difference pretrained_infold - scratch (same seed, same train and test arrays)
    wide = runs.pivot_table(index=["holdout_unit", "n_train_cultures", "seed"], columns="arm", values="rmse")
    diff = (wide["pretrained_infold"] - wide["scratch"]).groupby(level=[0, 1]).agg(["mean", "std", "count"])
    diff["wins"] = (wide["pretrained_infold"] < wide["scratch"]).groupby(level=[0, 1]).sum()
    diff = diff.reset_index().rename(columns={"mean": "paired_diff_mean", "std": "paired_diff_sd", "count": "n_seeds",
                                              "wins": "seeds_where_pretrained_better"})
    summ.round(4).to_csv(RESULTS / "pretrain_summary.csv", index=False)
    diff.round(4).to_csv(RESULTS / "pretrain_paired_diff.csv", index=False)
    print(summ.pivot_table(index=["holdout_unit", "n_train_cultures"], columns="arm", values="mean").round(2).to_string())
    print("paired difference, pretrained_infold minus scratch (negative = pretraining helps):")
    print(diff.round(2).to_string(index=False))
    return runs


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
