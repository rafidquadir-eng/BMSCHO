"""Slide figures. Reads results/, writes figures/. No model fitting happens here."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import cho_poc as cp  # noqa: E402

RESULTS, FIGURES = ROOT / "results", ROOT / "figures"

# Colorblind-safe categorical slots (validated all-pairs for three series) + neutrals
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, NEUTRAL = "#0b0b0b", "#52514e", "#e4e3df", "#8a8985"
SURFACE = "#ffffff"

plt.rcParams.update({
    "font.size": 16, "axes.titlesize": 19, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "legend.fontsize": 15, "figure.titlesize": 21, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.titleweight": "bold",
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 1.0, "axes.axisbelow": True,
    "legend.frameon": False, "axes.titlelocation": "left", "axes.titlepad": 14,
})


def _metric(m, analysis, subset, cv, encoder, model, metric):
    r = m[(m.analysis == analysis) & (m.subset == subset) & (m.cv_scheme == cv) & (m.encoder == encoder)
          & (m.model == model) & (m.metric == metric)]
    if len(r) != 1:
        raise KeyError((analysis, subset, cv, encoder, model, metric, len(r)))
    r = r.iloc[0]
    return float(r.value), float(r.ci_lo), float(r.ci_hi)


def _preds(p, key):
    d = p[p.key == key]
    if d.empty:
        raise KeyError(key)
    return d


def fig_roc(m, p):
    from sklearn.metrics import roc_curve

    fig, ax = plt.subplots(figsize=(14.5, 8.2))
    series = [
        ("full_classification", "culture_grouped_5fold", "Held-out cultures", BLUE, "-"),
        ("full_classification", "leave_one_experiment_out", "Held-out experiments", ORANGE, "-"),
        ("leakage", "leave_one_expression_cluster_out", "Held-out families of\nrelated experiments", AQUA, "-"),
    ]
    for analysis, cv, label, color, ls in series:
        d = _preds(p, f"{analysis}|all_labeled_terciles|{cv}|raw|logit")
        fpr, tpr, _ = roc_curve(d.y.astype(int), d.pred)
        a, lo, hi = _metric(m, analysis, "all_labeled_terciles", cv, "raw", "logit", "auc")
        ax.plot(fpr, tpr, color=color, lw=3, ls=ls, label=f"{label}\nAUC {a:.2f}  (95% CI {lo:.2f}–{hi:.2f})")
    ax.plot([0, 1], [0, 1], color=NEUTRAL, lw=1.5, ls="--")
    ax.text(0.66, 0.60, "chance", color=INK2, rotation=45, fontsize=15, rotation_mode="anchor")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_aspect("equal")
    n = len(_preds(p, "full_classification|all_labeled_terciles|culture_grouped_5fold|raw|logit"))
    ax.set_title("High vs low Qp from expression:\nAUC falls as the held-out unit gets stricter")
    leg = ax.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), title=f"Top vs bottom Qp tercile\nn = {n} arrays",
                    title_fontsize=15, labelspacing=1.3, handlelength=1.6, fontsize=16)
    leg._legend_box.align = "left"
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_roc.png", dpi=200)
    plt.close(fig)


def fig_scatter(m, p):
    fig, axes = plt.subplots(1, 2, figsize=(17, 8.4), sharex=True, sharey=True)
    panels = [
        ("matched_baseline", "matched_stationary_tempshift", "Matched subset: the head-to-head",
         "stationary, temperature-shifted arrays"),
        ("full_regression", "all_labeled", "All labeled arrays: a different experiment",
         "not comparable with the published 4.44"),
    ]
    lim = (-3, 60)
    for ax, (analysis, subset, title, sub) in zip(axes, panels):
        d = _preds(p, f"{analysis}|{subset}|culture_grouped_5fold|raw|pls")
        r, lo, hi = _metric(m, analysis, subset, "culture_grouped_5fold", "raw", "pls", "rmse")
        q, qlo, qhi = _metric(m, analysis, subset, "culture_grouped_5fold", "raw", "pls", "r2_cv")
        xs = np.array(lim)
        ax.fill_between(xs, xs - cp.BASELINE_RMSE, xs + cp.BASELINE_RMSE, color=NEUTRAL, alpha=0.22, lw=0,
                        label=f"±{cp.BASELINE_RMSE} pg/cell/day\n(published RMSE, Clarke 2011)")
        ax.plot(xs, xs, color=INK2, lw=1.5, label="1:1")
        ax.scatter(d.y, d.pred, s=70, color=BLUE, edgecolor=SURFACE, linewidth=1.2, alpha=0.9, zorder=3)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.set_title(f"{title}\n", fontsize=18)
        ax.text(0.0, 1.015, sub, transform=ax.transAxes, color=INK2, fontsize=15)
        ax.set_xlabel("Measured Qp (pg/cell/day)")
        txt = (f"n = {len(d)} arrays, {d.culture.nunique()} cultures\n"
               f"RMSE {r:.1f}  (95% CI {lo:.1f}–{hi:.1f})\n"
               f"R² {q:.2f}  (95% CI {qlo:.2f}–{qhi:.2f})")
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", ha="left", fontsize=15,
                bbox=dict(boxstyle="round,pad=0.4", fc=SURFACE, ec=GRID))
    axes[0].set_ylabel("Predicted Qp, held-out cultures (pg/cell/day)")
    axes[1].legend(loc="lower right", handlelength=1.4)
    r0 = _metric(m, "matched_baseline", "matched_stationary_tempshift", "culture_grouped_5fold", "raw", "pls", "rmse")
    verdict = "above" if r0[0] > cp.BASELINE_RMSE else "below"
    fig.suptitle(f"Culture-grouped CV: RMSE {r0[0]:.1f} on the matched subset, {verdict} the published "
                 f"{cp.BASELINE_RMSE}", x=0.012, ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURES / "fig1_scatter.png", dpi=200)
    plt.close(fig)


def fig_leakage(m, p):
    schemes = [
        ("leakage", "array_kfold_LEAKY", "Arrays\n(leaky)", NEUTRAL),
        (None, "culture_grouped_5fold", "Cultures", BLUE),
        (None, "leave_one_experiment_out", "Experiments", BLUE),
        ("leakage", "leave_one_expression_cluster_out", "Families of\nrelated\nexperiments", BLUE),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(17, 8.2))
    specs = [
        (axes[0], "full_regression", "all_labeled", "pls", "rmse", "Regression error (RMSE, pg/cell/day)", "{:.1f}"),
        (axes[1], "full_classification", "all_labeled_terciles", "logit", "auc", "High vs low Qp (AUC)", "{:.2f}"),
    ]
    for ax, base, subset, model, metric, title, fmt in specs:
        vals = []
        for analysis, cv, label, color in schemes:
            # the two randomly partitioned schemes use the 10-partition estimate; the
            # leave-one-group-out schemes have a single possible partition
            random_partition = cv in ("array_kfold_LEAKY", "culture_grouped_5fold")
            vals.append(_metric(m, "leakage_repeated_partitions" if random_partition else (analysis or base),
                                subset, cv, "raw", model, metric))
        x = np.arange(len(schemes))
        pts = np.array([v[0] for v in vals])
        err = np.array([[v[0] - v[1] for v in vals], [v[2] - v[0] for v in vals]])
        ax.bar(x, pts, width=0.62, color=[s[3] for s in schemes], edgecolor=SURFACE, linewidth=2, zorder=2)
        ax.errorbar(x, pts, yerr=err, fmt="none", ecolor=INK, elinewidth=2, capsize=6, capthick=2, zorder=3)
        for xi, v in zip(x, vals):
            # value inside the bar, clear of the CI whisker and the reference line
            ax.text(xi, 0.7 if metric == "rmse" else 0.05, fmt.format(v[0]), ha="center", va="bottom",
                    fontsize=21, fontweight="bold", color=SURFACE)
        ax.set_xticks(x)
        ax.set_xticklabels([s[2] for s in schemes])
        ax.set_title(title, fontsize=18)
        ax.grid(axis="x", visible=False)
        if metric == "rmse":
            null = _metric(m, base, subset, "leave_one_experiment_out", "raw", model, "rmse_null")[0]
            ax.axhline(null, color=INK2, ls="--", lw=1.5)
            ax.text(-0.42, null + 0.25, f"no-skill reference: predict the training mean ({null:.1f})", color=INK2,
                    fontsize=14, va="bottom")
            ax.set_ylim(0, max(null, max(v[2] for v in vals)) + 2.6)
            ax.set_ylabel("lower is better")
        else:
            ax.axhline(0.5, color=INK2, ls="--", lw=1.5)
            ax.text(-0.42, 0.51, "chance (0.5)", color=INK2, fontsize=14, va="bottom")
            ax.set_ylim(0, 1.12)
            ax.set_ylabel("higher is better")
        ax.set_xlabel("What is held out together")
    fig.suptitle("Same model, same data: what you hold out decides the answer", x=0.012, ha="left",
                 fontweight="bold")
    fig.text(0.012, 0.905, "Splitting by culture is not enough here: sister flasks of one experiment land in "
             "both train and test. Bars: point estimate, 95% bootstrap CI.", color=INK2, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(FIGURES / "fig2_leakage.png", dpi=200)
    plt.close(fig)


def fig_learning_curve():
    path = RESULTS / "pretrain_summary.csv"
    if not path.exists():
        print("  fig3 skipped: no pretraining results")
        return
    s = pd.read_csv(path)
    arms = [
        ("scratch", "From scratch", ORANGE, "-", "o"),
        ("pretrained_infold", "Pretrained inside the fold", BLUE, "-", "s"),
        ("pretrained_transductive", "Pretrained on all arrays\n(transductive: sees test expression)", AQUA, "--", "^"),
        ("pls", "PLS reference", NEUTRAL, "-", "D"),
    ]
    units = [("culture", "15 cultures held out"), ("experiment", "5 whole experiments held out")]
    fig, axes = plt.subplots(1, 2, figsize=(17, 8.4), sharey=True)
    for ax, (unit, title) in zip(axes, units):
        d = s[s.holdout_unit == unit]
        offs = np.linspace(-0.6, 0.6, len(arms))
        for (arm, label, color, ls, mk), off in zip(arms, offs):
            a = d[d.arm == arm].sort_values("n_train_cultures")
            ax.errorbar(a.n_train_cultures + off, a["mean"], yerr=a["sd"], color=color, ls=ls, lw=2.6, marker=mk,
                        markersize=9, markeredgecolor=SURFACE, markeredgewidth=1.2, capsize=4, elinewidth=1.6,
                        label=label)
        null = d[d.arm == "train_mean"]["mean"].mean()
        ax.axhline(null, color=INK2, ls=":", lw=2, label="Predict the training mean (no skill)")
        ax.set_title(title, fontsize=18)
        ax.set_xlabel("Labeled cultures used for training")
        ax.set_xticks(sorted(d.n_train_cultures.unique()))
    axes[0].set_ylabel("Test RMSE (pg/cell/day), mean ± SD over 5 seeds")
    axes[0].set_ylim(bottom=0)
    axes[1].legend(loc="lower left", handlelength=2.2, labelspacing=0.7)
    verdict = pretrain_verdict()
    fig.suptitle(verdict, x=0.012, ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURES / "fig3_learning_curve.png", dpi=200)
    plt.close(fig)


def pretrain_verdict() -> str:
    """Headline for fig3, computed from the paired differences so it cannot drift from the data."""
    d = pd.read_csv(RESULTS / "pretrain_paired_diff.csv")
    n = int(d.n_seeds.max())
    helped = int((d.seeds_where_pretrained_better == d.n_seeds).sum())
    hurt = int((d.seeds_where_pretrained_better == 0).sum())
    mean_diff = d.paired_diff_mean.mean()
    if helped == len(d):
        return f"Self-supervised pretraining lowers test error at every training size (all {n} seeds)"
    if hurt == len(d):
        return f"Self-supervised pretraining raises test error at every training size (all {n} seeds)"
    if abs(mean_diff) < 0.5 * d.paired_diff_sd.mean():
        return "Self-supervised pretraining gives no consistent gain over training from scratch"
    word = "lowers" if mean_diff < 0 else "raises"
    return (f"Self-supervised pretraining {word} test error on average ({mean_diff:+.1f} pg/cell/day), "
            "not consistently across seeds")


def make_all():
    FIGURES.mkdir(exist_ok=True)
    m = pd.read_csv(RESULTS / "metrics.csv")
    p = pd.read_csv(RESULTS / "oof_predictions.csv")
    fig_roc(m, p)
    fig_scatter(m, p)
    fig_leakage(m, p)
    fig_learning_curve()
    print("figures written to", FIGURES)


if __name__ == "__main__":
    make_all()
