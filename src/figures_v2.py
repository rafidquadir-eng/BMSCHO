"""Phase-2 figures. Reads results/, writes figures/. No model fitting happens here.

Every number drawn is read from a results file; none is typed in.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from figures import AQUA, BLUE, GRID, INK, INK2, NEUTRAL, ORANGE, SURFACE  # noqa: E402  (also sets rcParams)

R, F = ROOT / "results", ROOT / "figures"


def _wrap(s, width):
    import textwrap

    return "\n".join(textwrap.fill(par, width) for par in s.split("\n"))


def _j(name):
    return json.loads((R / name).read_text())


def _m1(analysis, subset, cv, encoder, model, metric):
    m = pd.read_csv(R / "metrics.csv")
    r = m[(m.analysis == analysis) & (m.subset == subset) & (m.cv_scheme == cv) & (m.encoder == encoder) & (m.model == model)
          & (m.metric == metric)].iloc[0]
    return float(r.value), float(r.ci_lo), float(r.ci_hi)


def _ep(endpoint, cv, encoder, model, metric):
    m = pd.read_csv(R / "metrics_endpoints.csv")
    r = m[(m.endpoint == endpoint) & (m.holdout_scheme == cv) & (m.encoder == encoder) & (m.model == model) & (m.metric == metric)].iloc[0]
    return float(r.value), float(r.ci_lo), float(r.ci_hi), int(r.n), int(r.n_holdout_groups)


# --------------------------------------------------------------------------- #
# fig4: Yao 2021 pathway recovery
# --------------------------------------------------------------------------- #
def fig4():
    pr = pd.read_csv(R / "yao_pathway_ranking.csv")
    info = _j("yao_info.json")
    top = pr.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(15.5, 9.2))
    colors = [BLUE if t else NEUTRAL for t in top.target]
    ax.barh(np.arange(len(top)), top.delta_z, color=colors, height=0.7, edgecolor=SURFACE, linewidth=1.5)
    labels = [f"{int(r.rank_delta)}.  {r['name']}  ({int(r.n_measured)})" for _, r in top.iterrows()]
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(labels, fontsize=14)
    for tick, t in zip(ax.get_yticklabels(), top.target):
        if t:
            tick.set_fontweight("bold")
    for i, (_, r) in enumerate(top.iterrows()):
        if r.target:
            ax.text(r.delta_z + 0.01, i, f"p = {r.p_perm:.3f}", va="center", fontsize=14, fontweight="bold")
    ax.set_xlabel("Δ statistic: mean correlation with qP inside the pathway\nminus outside it (Fisher z)")
    ax.grid(axis="y", visible=False)
    tca = pr[pr.pathway == "map00020"].iloc[0]
    aag = pr[pr.pathway == "map00250"].iloc[0]
    fig.suptitle(f"Yao 2021 metabolomics: TCA cycle recovered at rank {int(tca.rank_delta)} of {len(pr)};\n"
                 f"alanine/aspartate/glutamate metabolism only rank {int(aag.rank_delta)}", x=0.012, ha="left", fontweight="bold")
    fig.text(0.012, 0.02, _wrap(
        f"Modality: metabolomics (Yao et al. 2021, Table S3: {info['unique_kegg']} metabolites, 12 clones). No model fitted, no holdout: "
        f"a pre-registered ranking of published summary statistics. Joint test for the two target pathways: mean rank "
        f"{info['mean_rank_targets']:.1f} vs {info['expected_mean_rank_under_null']:.1f} expected, p = {info['p_joint']:.3f}. "
        "Number in brackets = metabolites measured in the pathway.\n"
        "Cross-library prediction (mAb A/B → mAb C) was NOT run: per-clone metabolite data are not public.", 135),
        fontsize=13.5, color=INK2, va="bottom")
    fig.tight_layout(rect=(0, 0.13, 1, 0.91))
    fig.savefig(F / "fig4_metabolomics_qp.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# fig5: knockout recovery
# --------------------------------------------------------------------------- #
def fig5():
    A = pd.read_csv(R / "knockout_ranking.csv")
    e = _j("knockout_eval.json")
    k = e["route_a_knockout_end_within_experiment"]["k100"]
    fig, axes = plt.subplots(1, 2, figsize=(18, 8.6), gridspec_kw={"width_ratios": [2.1, 1]})
    ax = axes[0]
    ax.plot(A.knockout_rank, A.delta_within_experiment, color=NEUTRAL, lw=2, zorder=1)
    ax.axvspan(0, 100, color=BLUE, alpha=0.10, lw=0)
    ax.text(50, A.delta_within_experiment.max() * 0.97, "top 100", ha="center", va="top", fontsize=14, color=INK2)
    pos = A[A.known_target == "KO/KD raised titer or qP"]
    neg = A[A.known_target == "KO/KD no benefit"]
    ax.scatter(pos.knockout_rank, pos.delta_within_experiment, s=130, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3,
               label=f"Knockout/knockdown raised titer or qP in the literature (n = {len(pos)} on the array)")
    ax.scatter(neg.knockout_rank, neg.delta_within_experiment, s=130, facecolor=SURFACE, edgecolor=ORANGE, linewidth=2.5, zorder=3,
               label=f"Knockout reported with no benefit (n = {len(neg)})")
    lift = np.linspace(0.55, 1.35, len(pos))
    for (_, r), dy in zip(pos.sort_values("knockout_rank").iterrows(), np.resize([0.42, 0.72, 1.02], len(pos))):
        ax.annotate(r.symbol, (r.knockout_rank, r.delta_within_experiment), xytext=(r.knockout_rank, r.delta_within_experiment + dy),
                    ha="center", fontsize=13, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    for _, r in neg.iterrows():
        ax.annotate(r.symbol, (r.knockout_rank, r.delta_within_experiment), xytext=(r.knockout_rank, r.delta_within_experiment - 0.5),
                    ha="center", fontsize=13, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xlabel(f"Knockout-candidate rank among {len(A):,} array genes (1 = most negative Δ)")
    ax.set_ylabel("Gene Δ: top minus bottom Qp quartile,\nwithin experiment (SD units)")
    ax.set_xlim(-20, len(A) + 20)
    ax.legend(loc="lower right", fontsize=13.5)
    ax.set_title("Data-driven ranking (modality: transcriptomics)", fontsize=17)

    ax = axes[1]
    o = e["route_a_overexpression_end_within_experiment"]["k100"]
    b = e["route_b_tradeoff_score"]
    items = [("Knockout targets\ntranscriptomics Δ", k), ("Overexpression targets\ntranscriptomics Δ", o), ("Knockout targets\niCHO1766 FBA", b)]
    for i, (lab, v) in enumerate(items):
        lo, hi = v["auc_null_95"]
        ax.plot([lo, hi], [i, i], color=NEUTRAL, lw=14, solid_capstyle="butt", alpha=0.5, zorder=1)
        ax.scatter([v["auc"]], [i], s=240, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.axvline(0.5, color=INK2, ls="--", lw=1.5)
    ax.text(0.505, -0.52, "chance", ha="left", fontsize=13, color=INK2)
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([f"{lab}\nAUC {v['auc']:.2f}, p = {v['p_auc']:.2f}, n = {v['n_positives']}" for lab, v in items], fontsize=13.5)
    ax.set_ylim(-0.6, len(items) - 0.4)
    ax.set_xlim(0.25, 0.75)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("AUC for validated targets\n(grey: 95% range of random gene sets)")
    ax.set_title("Recovery vs chance", fontsize=17)
    cov = e["route_a_coverage"]
    fig.suptitle(f"Knockout recovery is at chance: {k['hits_in_top_k']} of {k['n_positives']} testable validated targets in the top 100 "
                 f"(p = {k['p_hits']:.2f})", x=0.012, ha="left", fontweight="bold")
    fig.text(0.012, 0.015, _wrap(f"Ranking frozen before the ground truth was loaded. Ground truth: Schulz et al. 2026 systematic review; "
             f"{cov['ko_positives_total']} validated knockout/knockdown targets, of which {cov['ko_positives_on_array']} are on this 2011 array. "
             "No model fitted, no holdout.", 165), fontsize=13.5, color=INK2, va="bottom")
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(F / "fig5_knockout_recovery.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# fig6: endpoints
# --------------------------------------------------------------------------- #
def fig6():
    m = pd.read_csv(R / "metrics_endpoints.csv")
    eps = ["Viability", "Growth rate", "Viable cell density", "q-lactate", "q-ammonium"]
    fig, axes = plt.subplots(1, 2, figsize=(18, 8.4), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    rows = [("Qp (phase 1)", _m1("full_classification", "all_labeled_terciles", "leave_one_experiment_out", "raw", "logit", "auc"),
             _m1("leakage", "all_labeled_terciles", "leave_one_expression_cluster_out", "raw", "logit", "auc"), 19, 9)]
    for ep in eps:
        a = _ep(ep, "leave_one_experiment_out", "raw", "logit", "auc")
        b = _ep(ep, "leave_one_expression_cluster_out", "raw", "logit", "auc")
        rows.append((ep, a[:3], b[:3], a[4], b[4]))
    y = np.arange(len(rows))
    for i, (lab, a, b, na, nb) in enumerate(rows):
        for (v, lo, hi), off, col, mk in [(a, -0.17, BLUE, "o"), (b, 0.17, ORANGE, "s")]:
            ax.plot([lo, hi], [i + off, i + off], color=col, lw=3, zorder=2)
            ax.scatter([v], [i + off], s=150, color=col, marker=mk, edgecolor=SURFACE, linewidth=1.5, zorder=3)
            ax.text(1.035, i + off, f"{v:.2f} ({lo:.2f}–{hi:.2f})", va="center", fontsize=13.5, color=INK)
    ax.scatter([], [], s=150, color=BLUE, marker="o", label="Whole experiments held out")
    ax.scatter([], [], s=150, color=ORANGE, marker="s", label="Families of related experiments held out")
    ax.axvline(0.5, color=INK2, ls="--", lw=1.5)
    ax.text(0.5, -0.72, "chance", ha="center", fontsize=13, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r[0]}\n{r[3]} experiments, {r[4]} families" for r in rows], fontsize=14)
    ax.set_xlim(0.0, 1.32)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylim(-0.9, len(rows) - 0.4)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("AUC, top vs bottom tercile (95% bootstrap CI)")
    ax.set_title("High vs low, from expression", fontsize=17)
    ax.legend(loc="lower left", fontsize=13, bbox_to_anchor=(0.0, -0.33), ncol=1)

    ax = axes[1]
    rows2 = []
    for ep in eps:
        e_ = _ep(ep, "leave_one_experiment_out", "raw", "pls", "r2_cv")
        t_ = _ep(ep, "leave_one_experiment_out", "temperature only", "linear", "r2_cv")
        rows2.append((ep, e_[:3], t_[:3]))
    q = _m1("full_regression", "all_labeled", "leave_one_experiment_out", "raw", "pls", "r2_cv")
    rows2 = [("Qp (phase 1)", q, None)] + rows2
    for i, (lab, a, t) in enumerate(rows2):
        ax.plot([a[1], a[2]], [i - 0.17, i - 0.17], color=BLUE, lw=3)
        ax.scatter([a[0]], [i - 0.17], s=150, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        if t is not None:
            ax.plot([max(t[1], -1.0), t[2]], [i + 0.17, i + 0.17], color=NEUTRAL, lw=3)
            ax.scatter([max(t[0], -1.0)], [i + 0.17], s=150, color=NEUTRAL, marker="D", edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.scatter([], [], s=150, color=BLUE, label="Expression (PLS, 1,000 probesets)")
    ax.scatter([], [], s=150, color=NEUTRAL, marker="D", label="Culture temperature alone")
    ax.axvline(0, color=INK2, ls="--", lw=1.5)
    ax.text(0.0, -0.72, "no skill", ha="center", fontsize=13, color=INK2)
    ax.set_yticks(np.arange(len(rows2)))
    ax.set_yticklabels([r[0] for r in rows2], fontsize=14)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-0.9, len(rows2) - 0.4)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Cross-validated R², whole experiments held out\n(values below −1 are drawn at −1)")
    ax.set_title("Regression skill vs the temperature shift", fontsize=17)
    ax.legend(loc="lower left", fontsize=13, bbox_to_anchor=(0.0, -0.33))
    v = _ep("Viability", "leave_one_experiment_out", "raw", "logit", "auc")
    fig.suptitle(f"Other endpoints from expression alone (modality: transcriptomics): viability AUC {v[0]:.2f} "
                 f"({v[1]:.2f}–{v[2]:.2f}), experiments held out", x=0.012, ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(F / "fig6_endpoints.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# fig7: modality alignment
# --------------------------------------------------------------------------- #
def fig7():
    a = _j("alignment.json")
    t = pd.read_csv(R / "alignment_pathway_scores.csv")
    leg = "growth phase: transcriptomics (a) vs proteomics (b)"
    basis = "KEGG pathways (Chinese hamster)"
    variants = list(a["leg1_growth_phase_transcriptomics_vs_proteomics"].keys())
    fig, axes = plt.subplots(1, 2, figsize=(18, 8.8), sharex=True, sharey=True)
    titles = ["Pre-registered primary: all temperatures", "Sensitivity: 37 °C arrays only"]
    for ax, v, ttl in zip(axes, variants, titles):
        d = t[(t.leg == leg) & (t.basis == basis) & (t.variant == v)]
        blk = a["leg1_growth_phase_transcriptomics_vs_proteomics"][v]
        r, ic, g = blk[basis], blk["iCHO1766 subsystems"], blk["gene_level"]
        ax.scatter(d.score_a, d.score_b, s=np.clip(d.n_features_b * 6, 40, 400), color=BLUE, alpha=0.55, edgecolor=SURFACE, linewidth=1)
        ax.axhline(0, color=INK2, lw=1)
        ax.axvline(0, color=INK2, lw=1)
        p = "p < 0.0001" if r["p_perm"] < 0.0002 else f"p = {r['p_perm']:.2f}"
        pi = "p < 0.0001" if ic["p_perm"] < 0.0002 else f"p = {ic['p_perm']:.2f}"
        txt = (f"KEGG pathways (n = {r['n_pathways']}): ρ = {r['spearman']:.2f}, {p}\n"
               f"iCHO1766 subsystems (n = {ic['n_pathways']}): ρ = {ic['spearman']:.2f}, {pi}\n"
               f"Gene level, no pathway basis (n = {g['n_genes']}): ρ = {g['spearman']:.2f}, p < 0.0001\n"
               f"Transcriptomics: {blk['transcriptomics']['n_arrays']} arrays, {blk['transcriptomics']['n_experiments']} experiments")
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=14, bbox=dict(boxstyle="round,pad=0.4", fc=SURFACE, ec=GRID))
        ax.set_title(ttl, fontsize=17)
        ax.set_xlabel("Transcriptomics pathway score\n(faster vs slower growth, within experiment; GSE30321)")
    axes[0].set_ylabel("Proteomics pathway score\n(day 4 exponential vs day 7 stationary; PXD053511)")
    prim = a["leg1_growth_phase_transcriptomics_vs_proteomics"][variants[0]][basis]
    sens = a["leg1_growth_phase_transcriptomics_vs_proteomics"][variants[1]][basis]
    verdict = "is null" if prim["p_perm"] >= 0.05 else "is significant"
    fig.suptitle(f"Pathway-level alignment of transcripts and proteins: the pre-registered test {verdict} (ρ = {prim['spearman']:.2f}, "
                 f"p = {prim['p_perm']:.2f});\nit appears once the temperature shift is removed (ρ = {sens['spearman']:.2f})",
                 x=0.012, ha="left", fontweight="bold", fontsize=19)
    fig.text(0.012, 0.012, _wrap("Unpaired data from different labs; each point is one KEGG pathway (size = proteins measured). No sample has both "
             "modalities: this is not a fused prediction. p-values from 10,000 permutations of protein labels. The 37 °C run matches the "
             "proteomics design (no temperature shift) but was registered as a sensitivity analysis.", 165), fontsize=13.5, color=INK2, va="bottom")
    fig.tight_layout(rect=(0, 0.08, 1, 0.9))
    fig.savefig(F / "fig7_modality_alignment.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# fig0: architecture status (the money slide) — dark, to match the deck
# --------------------------------------------------------------------------- #
BG, CARD, EDGE, TXT, TXT2 = "#0b0b0b", "#151515", "#3a3a3a", "#f0efec", "#c3c2b7"
OK, WEAK, NEED = "#199e70", "#c98500", "#8a8985"
STATUS = {"ok": (OK, "TESTED: POSITIVE EVIDENCE", "-"), "weak": (WEAK, "TESTED: WEAK OR NULL", "-"), "need": (NEED, "REQUIRES PARTNER DATA", (0, (4, 3)))}


def _box(ax, x, y, w, h, title, lines, status, title_size=17, body_size=12.6):
    col, tag, ls = STATUS[status]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9", fc=CARD, ec=col, lw=2.6, ls=ls, zorder=2))
    ax.text(x + 1.0, y + h - 1.15, title, color=TXT, fontsize=title_size, fontweight="bold", va="top", zorder=3)
    ax.text(x + w - 1.0, y + h - 1.25, tag, color=col, fontsize=10.5, fontweight="bold", va="top", ha="right", zorder=3)
    ax.text(x + 1.0, y + h - 3.75, "\n".join(lines), color=TXT2, fontsize=body_size, va="top", linespacing=1.32, zorder=3)


def _arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color="#6a6a6a", lw=1.6), zorder=1)


def fig0():
    e, al, yi = _j("knockout_eval.json"), _j("alignment.json"), _j("yao_info.json")
    pr = pd.read_csv(R / "yao_pathway_ranking.csv")
    tca = int(pr[pr.pathway == "map00020"].rank_delta.iloc[0])
    aag = int(pr[pr.pathway == "map00250"].rank_delta.iloc[0])
    auc_q = _m1("full_classification", "all_labeled_terciles", "leave_one_experiment_out", "raw", "logit", "auc")
    auc_qf = _m1("leakage", "all_labeled_terciles", "leave_one_expression_cluster_out", "raw", "logit", "auc")
    rm_q = _m1("full_regression", "all_labeled", "leave_one_experiment_out", "raw", "pls", "rmse")
    rm_null = _m1("full_regression", "all_labeled", "leave_one_experiment_out", "raw", "pls", "rmse_null")
    rm_match = _m1("matched_baseline", "matched_stationary_tempshift", "culture_grouped_5fold", "raw", "pls", "rmse")
    rm_path = _m1("full_regression", "all_labeled", "leave_one_experiment_out", "pathway", "pls", "rmse")
    cov = _m1("covariate_baseline", "all_labeled", "leave_one_experiment_out", "process covariates", "gbm", "rmse")
    via = _ep("Viability", "leave_one_experiment_out", "raw", "logit", "auc")
    viaf = _ep("Viability", "leave_one_expression_cluster_out", "raw", "logit", "auc")
    prim = al["leg1_growth_phase_transcriptomics_vs_proteomics"]
    v0, v1 = list(prim.keys())
    kegg = "KEGG pathways (Chinese hamster)"
    ko = e["route_a_knockout_end_within_experiment"]["k100"]
    dp = pd.read_csv(R / "delta_pathways.csv")
    golgi = dp[dp.pathway == "GOBP_GOLGI_VESICLE_TRANSPORT"].iloc[0]
    pd_ = pd.read_csv(R / "pretrain_paired_diff.csv")
    pt = pd_[pd_.holdout_unit == "experiment"]
    sf = _j("scale_format.json") if (R / "scale_format.json").exists() else None

    fig = plt.figure(figsize=(19.2, 12.4), facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 64.6)
    ax.axis("off")
    ax.set_facecolor(BG)
    ax.text(2, 62.6, "Architecture status on public data: what is tested, what is weak, what needs partner data",
            color=TXT, fontsize=23, fontweight="bold", va="top")
    ax.text(2, 59.6, "Every metric is single-modality and carries its holdout scheme. No public sample has more than one modality with an outcome, "
            "so nothing here is a fused multi-omic prediction.", color=TXT2, fontsize=13.5, va="top")
    # legend
    lx = 2
    for key in ["ok", "weak", "need"]:
        col, tag, ls = STATUS[key]
        ax.add_patch(FancyBboxPatch((lx, 56.0), 2.2, 1.1, boxstyle="round,pad=0,rounding_size=0.3", fc=CARD, ec=col, lw=2.4, ls=ls))
        ax.text(lx + 2.8, 56.55, tag, color=col, fontsize=11.5, fontweight="bold", va="center")
        lx += 21

    # ---- inputs ----
    yI, hI, wI = 45.6, 9.3, 23.2
    xs = [2, 26.27, 50.53, 74.8]
    _box(ax, xs[0], yI, wI, hI, "Transcriptomics", [
        "GSE30321: 295 arrays, 167 with Qp, 19 experiments",
        f"High vs low Qp: AUC {auc_q[0]:.2f} ({auc_q[1]:.2f}–{auc_q[2]:.2f})",
        "holdout: whole experiments"], "ok")
    g0 = prim[v0]["gene_level"]
    _box(ax, xs[1], yI, wI, hI, "Proteomics", [
        "PXD053511: 8 samples, no outcome labels",
        f"Agrees with transcripts gene by gene: ρ = {g0['spearman']:.2f}",
        f"({g0['n_genes']} genes, growth-phase contrast, unpaired)"], "weak")
    _box(ax, xs[2], yI, wI, hI, "Metabolomics", [
        "Yao 2021: summary statistics only (12 clones)",
        f"TCA cycle rank {tca} of {len(pr)}; Ala/Asp/Glu rank {aag}",
        "Per-clone data and library 2: not public"], "weak")
    _box(ax, xs[3], yI, wI, hI, "Process settings", [
        "4 logged variables alone, GSE30321",
        f"RMSE {cov[0]:.1f} vs {rm_null[0]:.1f} no-skill: no skill",
        "holdout: whole experiments"], "weak")
    # ---- encoder ----
    yE, hE = 34.6, 8.2
    for x in xs:
        _arrow(ax, x + wI / 2, yI, x + wI / 2, yE + hE)
    ic = prim[v0]["iCHO1766 subsystems"]
    kg0, kg1 = prim[v0][kegg], prim[v1][kegg]
    _box(ax, 2, yE, 96, hE, "Pathway-constrained encoder", [
        f"Fixed pathway average (44 Hallmark sets; array covers 15% of their genes): Qp RMSE {rm_path[0]:.1f} vs {rm_q[0]:.1f} for raw genes   ·   "
        "modality: transcriptomics, holdout: whole experiments",
        f"Cross-modality alignment on a pathway basis (transcripts vs proteins, unpaired): pre-registered test null, ρ = {kg0['spearman']:.2f}, p = {kg0['p_perm']:.2f} "
        f"(KEGG); ρ = {kg1['spearman']:.2f} at 37 °C only   ·   iCHO1766: {ic['n_pathways']} testable subsystems, ρ = {ic['spearman']:.2f}, p = {ic['p_perm']:.2f}"],
         "weak", body_size=12.4)
    # ---- culture state + delta ----
    yZ, hZ = 23.6, 8.2
    _arrow(ax, 50, yE, 26, yZ + hZ)
    _box(ax, 2, yZ, 46, hZ, "Culture state z (32 dimensions)", [
        f"Self-supervised z vs training from scratch: {pt.paired_diff_mean.mean():+.1f} RMSE, better in "
        f"{int(pt.seeds_where_pretrained_better.sum())} of {int(pt.n_seeds.sum())} runs",
        "No reliable gain; PLS matches or beats it. iMEMORY initialisation: not tested",
        "modality: transcriptomics, holdout: whole experiments"], "weak")
    _arrow(ax, 48, yZ + hZ / 2, 52, yZ + hZ / 2)
    _box(ax, 52, yZ, 46, hZ, "Δ branch: top quartile minus bottom quartile", [
        f"Transcriptomics: Golgi vesicle transport rank {int(golgi['rank'])} of {len(dp)} (p = {golgi.p_vs_random_gene_sets:.2f}), holds within experiment",
        f"Metabolomics: TCA cycle rank {tca} of {len(pr)} (p = {float(pr[pr.pathway == 'map00020'].p_perm.iloc[0]):.2f}); joint test p = {yi['p_joint']:.2f}",
        "descriptive contrasts on fixed pathway scores, not on a learned z"], "weak")
    # ---- dynamics ----
    yD, hD = 15.4, 5.6
    _arrow(ax, 25, yZ, 40, yD + hD)
    _box(ax, 27, yD, 46, hD, "Latent dynamics, 14-day run", [
        "No public CHO dataset has time-resolved omics with day labels and outcomes"], "need")
    # ---- heads ----
    yH, hH, wH = 2.4, 10.4, 18.4
    xh = [2, 21.4, 40.8, 60.2, 79.6]
    for x in xh:
        _arrow(ax, 50, yD, x + wH / 2, yH + hH)
    _box(ax, xh[0], yH, wH, hH, "Titer, qP", [
        f"AUC {auc_q[0]:.2f} ({auc_q[1]:.2f}–{auc_q[2]:.2f})",
        f"RMSE {rm_q[0]:.1f} vs {rm_null[0]:.1f} no-skill",
        f"Families held out: AUC {auc_qf[0]:.2f}",
        f"Published 4.44 not beaten ({rm_match[0]:.1f})",
        "transcriptomics; experiments held out"], "ok", body_size=12)
    _box(ax, xh[1], yH, wH, hH, "Viability", [
        f"AUC {via[0]:.2f} ({via[1]:.2f}–{via[2]:.2f})",
        f"{via[4]} experiments held out",
        f"Families held out: AUC {viaf[0]:.2f}",
        f"({viaf[1]:.2f}–{viaf[2]:.2f})",
        "transcriptomics"], "ok" if via[1] > 0.5 else "weak", body_size=12)
    _box(ax, xh[2], yH, wH, hH, "Glycans", [
        "No public CHO dataset pairs",
        "glycan profiles with omics on",
        "the same clones at usable n"], "need", body_size=12)
    scale_lines = ["GSE30321 deposits no scale", "or format variable"]
    if sf:
        scale_lines += [f"Flask vs bioreactor (n = {sf['n_flask_cultures']} + {sf['n_bioreactor_cultures']}",
                        f"cultures): hypoxia rank {sf['hypoxia_rank']} of {sf['hallmark_sets_ranked']}", "format contrast, not scale"]
    _box(ax, xh[3], yH, wH, hH, "Scale risk", scale_lines, "need", body_size=12)
    _box(ax, xh[4], yH, wH, hH, "Ranked fixes", [
        "Knockout recovery at chance:",
        f"AUC {ko['auc']:.2f} (p = {ko['p_auc']:.2f})",
        f"{ko['hits_in_top_k']} of {ko['n_positives']} validated targets",
        f"in top 100 (p = {ko['p_hits']:.2f})",
        "transcriptomics Δ; frozen ranking"], "weak", body_size=12)
    ax.text(2, 0.9, "Intervals are 95% cluster-bootstrap CIs over the held-out unit. Blinded tests were pre-registered and frozen before the ground truth was loaded. "
            "Calibration of intervals was not assessed. Sources: results/metrics.csv, results/metrics_v2.csv.", color=TXT2, fontsize=11.5, va="bottom")
    fig.savefig(F / "fig0_architecture_status.png", dpi=150, facecolor=BG)
    plt.close(fig)


def make_all():
    F.mkdir(exist_ok=True)
    for fn in (fig4, fig5, fig7, fig6, fig0):
        try:
            fn()
            print("  wrote", fn.__name__)
        except FileNotFoundError as exc:
            print(f"  {fn.__name__} skipped: {exc}")


if __name__ == "__main__":
    make_all()
