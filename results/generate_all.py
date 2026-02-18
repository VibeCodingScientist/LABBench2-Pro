#!/usr/bin/env python3
"""Generate all figures and tables for the LABBench2-Pro paper.

Two-model comparison: Opus 4.6 vs Sonnet 4.6.

Run from repo root:
  python results/generate_all.py
"""

import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy import stats

RESULTS_DIR = Path(__file__).parent
RAW_DIR = RESULTS_DIR / "raw"
FIG_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"

FIG_DIR.mkdir(exist_ok=True)
TABLE_DIR.mkdir(exist_ok=True)

# Models
MODELS = ["claude-opus-4.6", "claude-sonnet-4.6"]
MODEL_LABELS = {"claude-opus-4.6": "Opus 4.6", "claude-sonnet-4.6": "Sonnet 4.6"}
MODEL_COLORS = {"claude-opus-4.6": "#2171b5", "claude-sonnet-4.6": "#fc8d59"}

# Nature-style figure settings
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
    "font.size": 8,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

# Tier colors (for tier legend patches)
TIER_COLORS = {
    "tier1": "#2171b5",
    "tier2": "#238b45",
    "tier3": "#d94801",
}

TIER1_CATS = ["CloningScenarios", "LitQA2", "SeqQA", "ProtocolQA", "SuppQA", "FigQA", "DbQA"]
TIER2_CATS = ["Calibration", "HypothesisGeneration", "StructureAnalysis", "StatisticalReasoning"]


def load_eval_runs():
    rows = []
    with open(RAW_DIR / "eval_runs.csv", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            r["correct"] = r["correct"] == "t"
            r["tokens_in"] = int(r["tokens_in"]) if r["tokens_in"] else 0
            r["tokens_out"] = int(r["tokens_out"]) if r["tokens_out"] else 0
            r["cost_usd"] = float(r["cost_usd"]) if r["cost_usd"] else 0
            r["latency_ms"] = int(r["latency_ms"]) if r["latency_ms"] else 0
            rows.append(r)
    return rows


def load_chain_runs():
    rows = []
    with open(RAW_DIR / "chain_runs.csv", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            r["correct"] = r["correct"] == "t"
            r["step_num"] = int(r["step_num"])
            rows.append(r)
    return rows


def load_judge_audits():
    rows = []
    path = RAW_DIR / "judge_audits.csv"
    if not path.exists():
        return rows
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            r["judge_score"] = r["judge_score"] == "t"
            rows.append(r)
    return rows


def bootstrap_ci(data, n_boot=10000, ci=0.95):
    """BCa bootstrap confidence interval."""
    arr = np.array(data, dtype=float)
    n = len(arr)
    if n < 2:
        return np.mean(arr), np.mean(arr), np.mean(arr)
    result = stats.bootstrap(
        (arr,), np.mean, n_resamples=n_boot, confidence_level=ci, method="BCa"
    )
    return result.confidence_interval.low, np.mean(arr), result.confidence_interval.high


def split_by_model(rows):
    """Split rows into {model_name: [rows]} dict."""
    by_model = defaultdict(list)
    for r in rows:
        by_model[r["model_name"]].append(r)
    return by_model


def short_cat(name):
    """Shorten category names for x-axis labels."""
    return (name
            .replace("CloningScenarios", "Cloning\nScen.")
            .replace("HypothesisGeneration", "Hypothesis\nGen.")
            .replace("StatisticalReasoning", "Statistical\nReas.")
            .replace("StructureAnalysis", "Structure\nAnal.")
            .replace("ProtocolQA", "Protocol\nQA"))


# ─── Figure 1: Two-model accuracy comparison (grouped bars with CIs) ───

def fig1_accuracy_comparison(eval_runs):
    by_model = split_by_model(eval_runs)
    categories = TIER1_CATS + TIER2_CATS

    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(categories))
    bar_width = 0.35

    for i, model in enumerate(MODELS):
        model_runs = by_model.get(model, [])
        by_cat = defaultdict(list)
        for r in model_runs:
            by_cat[r["category"]].append(1 if r["correct"] else 0)

        means, lows, highs = [], [], []
        for cat in categories:
            if cat in by_cat and len(by_cat[cat]) > 1:
                lo, mu, hi = bootstrap_ci(by_cat[cat])
                means.append(mu * 100)
                lows.append(mu * 100 - lo * 100)
                highs.append(hi * 100 - mu * 100)
            else:
                means.append(0)
                lows.append(0)
                highs.append(0)

        offset = (i - 0.5) * bar_width
        bars = ax.bar(x + offset + bar_width / 2, means,
                      yerr=[lows, highs], width=bar_width,
                      capsize=2, color=MODEL_COLORS[model],
                      edgecolor="white", linewidth=0.3,
                      error_kw={"linewidth": 0.5},
                      label=MODEL_LABELS[model])

    # Significance markers
    pairwise = {
        "CloningScenarios": 0.014, "FigQA": 0.003, "DbQA": 0.001,
        "LitQA2": 0.210, "SeqQA": 0.214, "ProtocolQA": 0.434, "SuppQA": 0.065,
    }
    for j, cat in enumerate(categories):
        if cat in pairwise and pairwise[cat] < 0.05:
            # Get max bar height for this category
            max_h = 0
            for model in MODELS:
                model_runs = by_model.get(model, [])
                cat_runs = [1 if r["correct"] else 0 for r in model_runs if r["category"] == cat]
                if cat_runs:
                    _, mu, hi = bootstrap_ci(cat_runs)
                    max_h = max(max_h, hi * 100)
            stars = "***" if pairwise[cat] < 0.001 else "**" if pairwise[cat] < 0.01 else "*"
            ax.text(j, max_h + 4, stars, ha="center", va="bottom", fontsize=7, fontweight="bold")

    # Tier separator lines
    ax.axvline(x=len(TIER1_CATS) - 0.5, color="gray", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.text(len(TIER1_CATS) / 2, 102, "Tier 1: LABBench", ha="center", fontsize=6, color="gray")
    ax.text(len(TIER1_CATS) + len(TIER2_CATS) / 2, 102, "Tier 2: New Tasks", ha="center", fontsize=6, color="gray")

    ax.set_xticks(x)
    ax.set_xticklabels([short_cat(c) for c in categories], rotation=45, ha="right", fontsize=6.5)
    ax.set_ylabel("Accuracy (%)")
    ax.set_ylim(0, 110)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_accuracy_comparison.pdf")
    fig.savefig(FIG_DIR / "fig1_accuracy_comparison.png")
    plt.close()
    print("  fig1_accuracy_comparison")


# ─── Figure 2: Chain error propagation (both models) ───

def fig2_chain_error_propagation(chain_runs):
    by_model = split_by_model(chain_runs)

    fig, axes = plt.subplots(1, 3, figsize=(9, 3))

    # Panel a: Per-step accuracy by model
    ax1 = axes[0]
    for model in MODELS:
        chains = defaultdict(list)
        for r in by_model.get(model, []):
            chains[r["chain_id"]].append(r)

        step_correct = defaultdict(list)
        for cid, steps in chains.items():
            for s in steps:
                step_correct[s["step_num"]].append(1 if s["correct"] else 0)

        step_nums = sorted(step_correct.keys())
        step_means = [np.mean(step_correct[s]) * 100 for s in step_nums]
        ax1.plot(step_nums, step_means, "o-", color=MODEL_COLORS[model],
                 markersize=4, linewidth=1.5, label=MODEL_LABELS[model])

    ax1.set_xlabel("Step number")
    ax1.set_ylabel("Step accuracy (%)")
    ax1.set_ylim(0, 105)
    ax1.legend(fontsize=6, frameon=False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_title("a", fontsize=9, fontweight="bold", loc="left")

    # Panel b: E2E vs step-level (grouped bars)
    ax2 = axes[1]
    x = np.arange(2)
    bar_width = 0.3
    for i, model in enumerate(MODELS):
        chains = defaultdict(list)
        for r in by_model.get(model, []):
            chains[r["chain_id"]].append(r)
        total_steps = sum(len(s) for s in chains.values())
        correct_steps = sum(1 for cid in chains for s in chains[cid] if s["correct"])
        e2e_correct = sum(1 for cid, steps in chains.items() if all(s["correct"] for s in steps))
        vals = [correct_steps / total_steps * 100 if total_steps else 0,
                e2e_correct / len(chains) * 100 if chains else 0]
        offset = (i - 0.5) * bar_width
        bars = ax2.bar(x + offset + bar_width / 2, vals,
                       width=bar_width, color=MODEL_COLORS[model],
                       edgecolor="white", label=MODEL_LABELS[model])
        for b, v in zip(bars, vals):
            ax2.text(b.get_x() + b.get_width() / 2, v + 1.5,
                     f"{v:.1f}%", ha="center", fontsize=5.5)

    ax2.set_xticks(x)
    ax2.set_xticklabels(["Step-level", "End-to-end"], fontsize=7)
    ax2.set_ylabel("Accuracy (%)")
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=5.5, frameon=False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.set_title("b", fontsize=9, fontweight="bold", loc="left")

    # Panel c: Error propagation gap
    ax3 = axes[2]
    gaps = []
    for model in MODELS:
        chains = defaultdict(list)
        for r in by_model.get(model, []):
            chains[r["chain_id"]].append(r)
        total_steps = sum(len(s) for s in chains.values())
        correct_steps = sum(1 for cid in chains for s in chains[cid] if s["correct"])
        e2e_correct = sum(1 for cid, steps in chains.items() if all(s["correct"] for s in steps))
        step_acc = correct_steps / total_steps * 100 if total_steps else 0
        e2e_acc = e2e_correct / len(chains) * 100 if chains else 0
        gaps.append(step_acc - e2e_acc)

    bars = ax3.bar([MODEL_LABELS[m] for m in MODELS], gaps,
                   color=[MODEL_COLORS[m] for m in MODELS],
                   edgecolor="white", width=0.5)
    for b, v in zip(bars, gaps):
        ax3.text(b.get_x() + b.get_width() / 2, v + 0.5,
                 f"{v:.1f} pp", ha="center", fontsize=6.5, fontweight="bold")

    ax3.set_ylabel("Error propagation gap (pp)")
    ax3.set_ylim(0, max(gaps) + 8)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.set_title("c", fontsize=9, fontweight="bold", loc="left")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_chain_error_propagation.pdf")
    fig.savefig(FIG_DIR / "fig2_chain_error_propagation.png")
    plt.close()
    print("  fig2_chain_error_propagation")


# ─── Figure 3: Cost-accuracy scatter (both models) ───

def fig3_cost_accuracy(eval_runs):
    by_model = split_by_model(eval_runs)

    fig, ax = plt.subplots(figsize=(6, 4))
    markers = {"claude-opus-4.6": "o", "claude-sonnet-4.6": "s"}

    for model in MODELS:
        by_cat = defaultdict(lambda: {"correct": 0, "total": 0, "cost": 0})
        for r in by_model.get(model, []):
            cat = r["category"]
            by_cat[cat]["total"] += 1
            by_cat[cat]["correct"] += 1 if r["correct"] else 0
            by_cat[cat]["cost"] += r["cost_usd"]

        for cat, d in by_cat.items():
            acc = d["correct"] / d["total"] * 100
            cost = d["cost"]
            ax.scatter(cost, acc, s=d["total"] * 0.3,
                       color=MODEL_COLORS[model], alpha=0.7,
                       marker=markers[model],
                       edgecolors="white", linewidth=0.3)
            # Label only for Opus (avoid double labels)
            if model == MODELS[0]:
                ax.annotate(cat, (cost, acc), fontsize=4.5, ha="center", va="bottom",
                            xytext=(0, 4), textcoords="offset points", color="gray")

    # Legend
    handles = [
        plt.Line2D([], [], marker=markers[m], color=MODEL_COLORS[m],
                   linestyle="", markersize=6, label=MODEL_LABELS[m])
        for m in MODELS
    ]
    ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper left")

    ax.set_xlabel("Total cost (USD)")
    ax.set_ylabel("Accuracy (%)")
    ax.set_xscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_cost_accuracy.pdf")
    fig.savefig(FIG_DIR / "fig3_cost_accuracy.png")
    plt.close()
    print("  fig3_cost_accuracy")


# ─── Figure 4: Judge audit results ───

def fig4_judge_audit(judge_audits):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 2.5))

    # Panel a: Agreement rates
    labels = ["Inter-judge\nagreement", "Position\nbias", "Verbosity\nbias"]
    values = [90.0, 15.0, 5.0]
    colors_bar = ["#2171b5", "#d94801", "#d94801"]
    ax1.bar(labels, values, color=colors_bar, edgecolor="white", width=0.5)
    ax1.set_ylabel("Rate (%)")
    ax1.set_ylim(0, 100)
    ax1.axhline(y=50, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_title("a", fontsize=9, fontweight="bold", loc="left")

    # Panel b: Cohen's kappa interpretation
    kappa_val = 0.765
    kappa_thresholds = [0, 0.20, 0.40, 0.60, 0.80, 1.0]
    kappa_labels = ["Poor", "Fair", "Moderate", "Substantial", "Almost\nperfect"]
    kappa_colors = ["#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"]

    for i in range(len(kappa_thresholds) - 1):
        ax2.barh(0, kappa_thresholds[i+1] - kappa_thresholds[i],
                left=kappa_thresholds[i], height=0.4, color=kappa_colors[i], edgecolor="white")
        ax2.text((kappa_thresholds[i] + kappa_thresholds[i+1]) / 2, -0.35,
                kappa_labels[i], ha="center", va="top", fontsize=5.5)

    ax2.axvline(x=kappa_val, color="black", linewidth=1.5, zorder=5)
    ax2.text(kappa_val, 0.3, f"\u03ba = {kappa_val:.3f}", ha="center", va="bottom",
             fontsize=7, fontweight="bold")
    ax2.set_xlim(0, 1)
    ax2.set_yticks([])
    ax2.set_xlabel("Cohen's Kappa")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.set_title("b", fontsize=9, fontweight="bold", loc="left")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig4_judge_audit.pdf")
    fig.savefig(FIG_DIR / "fig4_judge_audit.png")
    plt.close()
    print("  fig4_judge_audit")


# ─── Figure 5: Latency distribution (both models) ───

def fig5_latency(eval_runs):
    by_model = split_by_model(eval_runs)
    all_cats = sorted(set(r["category"] for r in eval_runs))

    fig, ax = plt.subplots(figsize=(8, 3.5))
    positions = []
    labels = []
    data_groups = []
    colors = []

    pos = 0
    for cat in all_cats:
        for model in MODELS:
            model_runs = [r["latency_ms"] / 1000 for r in by_model.get(model, [])
                          if r["category"] == cat and r["latency_ms"] > 0]
            if model_runs:
                data_groups.append(model_runs)
                positions.append(pos)
                colors.append(MODEL_COLORS[model])
                pos += 1
        labels.append(cat)
        pos += 0.5  # gap between categories

    bp = ax.boxplot(data_groups, positions=positions, vert=True, patch_artist=True,
                    widths=0.7, showfliers=False,
                    medianprops={"color": "black", "linewidth": 1})

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Category labels at midpoints
    cat_positions = []
    idx = 0
    for cat in all_cats:
        cat_pos = []
        for model in MODELS:
            model_runs = [r for r in by_model.get(model, []) if r["category"] == cat and r["latency_ms"] > 0]
            if model_runs:
                cat_pos.append(idx)
                idx += 1
        if cat_pos:
            cat_positions.append(np.mean(cat_pos))
        idx += 0  # already incremented

    # Recalculate tick positions
    tick_positions = []
    pos = 0
    for cat in all_cats:
        cat_pos = []
        for model in MODELS:
            model_runs = [r for r in by_model.get(model, []) if r["category"] == cat and r["latency_ms"] > 0]
            if model_runs:
                cat_pos.append(pos)
                pos += 1
        if cat_pos:
            tick_positions.append(np.mean(cat_pos))
        pos += 0.5

    ax.set_xticks(tick_positions)
    ax.set_xticklabels([short_cat(c) for c in all_cats], rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Latency (seconds)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    handles = [mpatches.Patch(color=MODEL_COLORS[m], alpha=0.7, label=MODEL_LABELS[m]) for m in MODELS]
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig5_latency.pdf")
    fig.savefig(FIG_DIR / "fig5_latency.png")
    plt.close()
    print("  fig5_latency")


# ─── Figure 6: Pairwise significance heatmap ───

def fig6_pairwise_significance():
    pairwise = {
        "CloningScenarios": {"diff": 0.121, "p": 0.014},
        "LitQA2": {"diff": 0.027, "p": 0.210},
        "SeqQA": {"diff": 0.054, "p": 0.214},
        "ProtocolQA": {"diff": -0.007, "p": 0.434},
        "SuppQA": {"diff": 0.049, "p": 0.065},
        "FigQA": {"diff": 0.055, "p": 0.003},
        "DbQA": {"diff": 0.025, "p": 0.001},
    }

    cats = list(pairwise.keys())
    diffs = [pairwise[c]["diff"] * 100 for c in cats]
    pvals = [pairwise[c]["p"] for c in cats]
    sig = [p < 0.05 for p in pvals]

    fig, ax = plt.subplots(figsize=(5, 3))
    bar_colors = ["#2171b5" if s else "#bdbdbd" for s in sig]
    bars = ax.barh(range(len(cats)), diffs, color=bar_colors, edgecolor="white", height=0.6)

    for i, (d, p, s) in enumerate(zip(diffs, pvals, sig)):
        label = f"p={p:.3f}" + (" *" if s else "")
        # Always place label to the right of the bar end, but for tiny bars
        # place it further right to avoid y-axis overlap
        x_pos = max(d + 0.3, 0.8)
        ax.text(x_pos, i, label, va="center", ha="left", fontsize=6,
                fontweight="bold" if s else "normal")

    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(cats, fontsize=7)
    ax.set_xlabel("Opus 4.6 accuracy advantage (pp)")
    ax.axvline(x=0, color="gray", linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig6_pairwise_significance.pdf")
    fig.savefig(FIG_DIR / "fig6_pairwise_significance.png")
    plt.close()
    print("  fig6_pairwise_significance")


# ─── Table 1: Two-model comparison (LaTeX) ───

def table1_main_results(eval_runs):
    by_model = split_by_model(eval_runs)
    categories = TIER1_CATS + TIER2_CATS

    lines = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\caption{Two-model accuracy comparison across LABBench2-Pro categories. BCa 95\% bootstrap CIs (10{,}000 resamples). Significance: pairwise bootstrap test.}")
    lines.append(r"\label{tab:main-results}")
    lines.append(r"\begin{tabular}{llrrrrrl}")
    lines.append(r"\toprule")
    lines.append(r"Tier & Category & $n$ & Opus 4.6 & 95\% CI & Sonnet 4.6 & 95\% CI & Sig? \\")
    lines.append(r"\midrule")

    pairwise = {
        "CloningScenarios": 0.014, "FigQA": 0.003, "DbQA": 0.001,
        "LitQA2": 0.210, "SeqQA": 0.214, "ProtocolQA": 0.434, "SuppQA": 0.065,
    }

    for tier_name, cats in [("1: LABBench", TIER1_CATS), ("2: New Tasks", TIER2_CATS)]:
        first = True
        for cat in cats:
            tier_label = tier_name if first else ""
            first = False

            row_parts = [tier_label, cat]
            # n (use max across models)
            ns = []
            for model in MODELS:
                cat_runs = [r for r in by_model.get(model, []) if r["category"] == cat]
                ns.append(len(cat_runs))
            n = max(ns) if ns else 0
            row_parts.append(str(n))

            # Per-model accuracy + CI
            for model in MODELS:
                cat_runs = [1 if r["correct"] else 0 for r in by_model.get(model, []) if r["category"] == cat]
                if cat_runs:
                    lo, mu, hi = bootstrap_ci(cat_runs)
                    row_parts.append(f"{mu:.1%}")
                    row_parts.append(f"[{lo:.1%}, {hi:.1%}]")
                else:
                    row_parts.append("---")
                    row_parts.append("---")

            # Significance
            p = pairwise.get(cat)
            if p is not None:
                sig_str = f"\\textbf{{p={p:.3f}}}" if p < 0.05 else f"p={p:.3f}"
            else:
                sig_str = "---"
            row_parts.append(sig_str)

            lines.append("  " + " & ".join(row_parts) + " \\\\")
        lines.append(r"\midrule")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    with open(TABLE_DIR / "table1_main_results.tex", "w") as f:
        f.write("\n".join(lines))
    print("  table1_main_results.tex")


# ─── Table 2: Chain results (both models) ───

def table2_chain_results(chain_runs):
    by_model = split_by_model(chain_runs)

    template_map = {}
    try:
        with open(Path(__file__).parent.parent / "tasks/chains/chain_definitions.json") as f:
            defs = json.load(f)
            for d in defs:
                template_map[d["chain_id"]] = d.get("template", "")
    except FileNotFoundError:
        pass

    # Get all chain IDs
    all_chain_ids = sorted(
        set(r["chain_id"] for r in chain_runs),
        key=lambda x: int(x.replace("chain", ""))
    )

    lines = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\caption{Compositional chain results for both models. E2E = all steps correct.}")
    lines.append(r"\label{tab:chain-results}")
    lines.append(r"\begin{tabular}{llrrrr}")
    lines.append(r"\toprule")
    lines.append(r"Chain & Template & Opus Steps & Opus E2E & Sonnet Steps & Sonnet E2E \\")
    lines.append(r"\midrule")

    for cid in all_chain_ids:
        template = template_map.get(cid, "")
        parts = [cid, template]

        for model in MODELS:
            chains = defaultdict(list)
            for r in by_model.get(model, []):
                chains[r["chain_id"]].append(r)

            if cid in chains:
                steps = sorted(chains[cid], key=lambda x: x["step_num"])
                n_correct = sum(1 for s in steps if s["correct"])
                n_steps = len(steps)
                e2e = "Yes" if all(s["correct"] for s in steps) else "No"
                parts.append(f"{n_correct}/{n_steps}")
                parts.append(e2e)
            else:
                parts.append("---")
                parts.append("---")

        lines.append("  " + " & ".join(parts) + " \\\\")

    # Summary
    lines.append(r"\midrule")
    summary_parts = [r"\textbf{Total}", ""]
    for model in MODELS:
        chains = defaultdict(list)
        for r in by_model.get(model, []):
            chains[r["chain_id"]].append(r)
        total_steps = sum(len(s) for s in chains.values())
        correct_steps = sum(1 for cid in chains for s in chains[cid] if s["correct"])
        e2e_total = sum(1 for cid, steps in chains.items() if all(s["correct"] for s in steps))
        summary_parts.append(f"\\textbf{{{correct_steps}/{total_steps}}}")
        summary_parts.append(f"\\textbf{{{e2e_total}/{len(chains)}}}")

    lines.append("  " + " & ".join(summary_parts) + " \\\\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    with open(TABLE_DIR / "table2_chain_results.tex", "w") as f:
        f.write("\n".join(lines))
    print("  table2_chain_results.tex")


# ─── Table 3: Cost breakdown (both models) ───

def table3_cost_breakdown(eval_runs):
    by_model = split_by_model(eval_runs)

    lines = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\caption{Cost analysis by category for both models.}")
    lines.append(r"\label{tab:cost}")
    lines.append(r"\begin{tabular}{lrrrr}")
    lines.append(r"\toprule")
    lines.append(r"Category & Opus Cost & Sonnet Cost & Opus \$/Corr & Sonnet \$/Corr \\")
    lines.append(r"\midrule")

    all_cats = sorted(set(r["category"] for r in eval_runs))

    # Sort by Opus cost descending
    cat_costs = {}
    for cat in all_cats:
        opus_runs = [r for r in by_model.get(MODELS[0], []) if r["category"] == cat]
        cat_costs[cat] = sum(r["cost_usd"] for r in opus_runs)

    for cat in sorted(all_cats, key=lambda c: cat_costs.get(c, 0), reverse=True):
        parts = [cat]
        for model in MODELS:
            cat_runs = [r for r in by_model.get(model, []) if r["category"] == cat]
            total_cost = sum(r["cost_usd"] for r in cat_runs)
            correct = sum(1 for r in cat_runs if r["correct"])
            cpc = total_cost / correct if correct > 0 else float("inf")
            parts.append(f"\\${total_cost:.2f}")
            cpc_str = f"\\${cpc:.2f}" if cpc < 1000 else "---"
            parts.append(cpc_str)

        # Rearrange: cat, opus_cost, sonnet_cost, opus_cpc, sonnet_cpc
        opus_cost = parts[1]
        opus_cpc = parts[2]
        sonnet_cost = parts[3]
        sonnet_cpc = parts[4]
        lines.append(f"  {cat} & {opus_cost} & {sonnet_cost} & {opus_cpc} & {sonnet_cpc} \\\\")

    # Totals
    lines.append(r"\midrule")
    totals = [r"\textbf{Total}"]
    for model in MODELS:
        model_runs = by_model.get(model, [])
        total_cost = sum(r["cost_usd"] for r in model_runs)
        totals.append(f"\\textbf{{\\${total_cost:.2f}}}")
        totals.append("")
    lines.append(f"  {totals[0]} & {totals[1]} & {totals[3]} & & \\\\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    with open(TABLE_DIR / "table3_cost_breakdown.tex", "w") as f:
        f.write("\n".join(lines))
    print("  table3_cost_breakdown.tex")


# ─── Table 4: Methodological audit ───

def table4_methodological_audit():
    lines = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\caption{Methodological audit results.}")
    lines.append(r"\label{tab:audit}")
    lines.append(r"\begin{tabular}{llr}")
    lines.append(r"\toprule")
    lines.append(r"Audit & Metric & Value \\")
    lines.append(r"\midrule")
    lines.append(r"Contamination & Cloze match (Opus) & 0.0\% \\")
    lines.append(r"              & Cloze match (Sonnet) & 0.0\% \\")
    lines.append(r"              & Reverse match (Opus) & 0.0\% \\")
    lines.append(r"              & Temporal split & Insufficient metadata \\")
    lines.append(r"\midrule")
    lines.append(r"Judge Audit   & Inter-judge agreement & 90.0\% \\")
    lines.append(r"              & Cohen's $\kappa$ & 0.765 (substantial) \\")
    lines.append(r"              & Position bias rate & 15.0\% \\")
    lines.append(r"              & Verbosity bias & +5.0\% \\")
    lines.append(r"\midrule")
    lines.append(r"IRT Analysis  & Total items & 2{,}459 \\")
    lines.append(r"              & Low discrimination ($< 0.3$) & 2{,}326 (94.6\%) \\")
    lines.append(r"              & Recommended pruned set & 133 items \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    with open(TABLE_DIR / "table4_methodological_audit.tex", "w") as f:
        f.write("\n".join(lines))
    print("  table4_methodological_audit.tex")


# ─── Supplementary: Chain traces (both models) ───

def supplementary_chain_traces(chain_runs):
    by_model = split_by_model(chain_runs)

    with open(TABLE_DIR / "supplementary_chain_traces.md", "w") as f:
        f.write("# Supplementary Material: Full Chain Execution Traces\n\n")
        f.write("Complete model responses for all 30 compositional chains, both models.\n\n")

        for model in MODELS:
            f.write(f"# {MODEL_LABELS[model]}\n\n")
            chains = defaultdict(list)
            for r in by_model.get(model, []):
                chains[r["chain_id"]].append(r)

            for cid in sorted(chains.keys(), key=lambda x: int(x.replace("chain", ""))):
                steps = sorted(chains[cid], key=lambda x: x["step_num"])
                e2e = "PASS" if all(s["correct"] for s in steps) else "FAIL"
                f.write(f"## {cid} [{e2e}]\n\n")

                for s in steps:
                    status = "CORRECT" if s["correct"] else "WRONG"
                    f.write(f"### Step {s['step_num']} -- {s['task_id']} [{status}]\n\n")
                    resp = s.get("response", "")
                    if resp:
                        f.write(f"**Model Response:**\n\n{resp[:2000]}\n\n")
                    f.write("---\n\n")

    print("  supplementary_chain_traces.md")


# ─── Main ───

def main():
    print("Loading data...")
    eval_runs = load_eval_runs()
    chain_runs = load_chain_runs()
    judge_audits = load_judge_audits()

    # Count per model
    model_counts = defaultdict(int)
    for r in eval_runs:
        model_counts[r["model_name"]] += 1
    for m, c in model_counts.items():
        print(f"  {m}: {c} eval runs")
    print(f"  {len(chain_runs)} chain steps, {len(judge_audits)} judge audits")

    print("\nGenerating figures...")
    fig1_accuracy_comparison(eval_runs)
    fig2_chain_error_propagation(chain_runs)
    fig3_cost_accuracy(eval_runs)
    fig4_judge_audit(judge_audits)
    fig5_latency(eval_runs)
    fig6_pairwise_significance()

    print("\nGenerating tables...")
    table1_main_results(eval_runs)
    table2_chain_results(chain_runs)
    table3_cost_breakdown(eval_runs)
    table4_methodological_audit()

    print("\nGenerating supplementary materials...")
    supplementary_chain_traces(chain_runs)

    print(f"\nDone. Outputs in:")
    print(f"  Figures: {FIG_DIR}")
    print(f"  Tables:  {TABLE_DIR}")


if __name__ == "__main__":
    main()
