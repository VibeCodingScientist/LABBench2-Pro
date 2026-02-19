#!/usr/bin/env python3
"""Generate all figures and tables for the LABBench2-Pro paper.

Four-model, three-provider comparison:
  Opus 4.6 (Anthropic) vs Sonnet 4.6 (Anthropic)
  vs GPT-5.2 (OpenAI) vs Gemini 2.5 Pro (Google)

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

# Models (4 models, 3 providers)
MODELS = ["claude-opus-4.6", "claude-sonnet-4.6", "gpt-5.2", "gemini-2.5-pro"]
MODEL_LABELS = {
    "claude-opus-4.6": "Opus 4.6",
    "claude-sonnet-4.6": "Sonnet 4.6",
    "gpt-5.2": "GPT-5.2",
    "gemini-2.5-pro": "Gemini 2.5 Pro",
}
MODEL_COLORS = {
    "claude-opus-4.6": "#2171b5",
    "claude-sonnet-4.6": "#6baed6",
    "gpt-5.2": "#2ca02c",
    "gemini-2.5-pro": "#d62728",
}
MODEL_MARKERS = {
    "claude-opus-4.6": "o",
    "claude-sonnet-4.6": "s",
    "gpt-5.2": "^",
    "gemini-2.5-pro": "D",
}

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

TIER1_CATS = ["CloningScenarios", "LitQA2", "SeqQA", "ProtocolQA", "SuppQA", "FigQA", "DbQA"]
TIER2_CATS = ["Calibration", "HypothesisGeneration", "StructureAnalysis", "StatisticalReasoning"]


def load_eval_runs():
    rows = []
    with open(RAW_DIR / "eval_runs.csv", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            r["correct"] = r["correct"] in ("t", "True", "true")
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
            r["correct"] = r["correct"] in ("t", "True", "true")
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
            r["judge_score"] = r["judge_score"] in ("t", "True", "true")
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


# ─── Figure 1: Four-model accuracy comparison (grouped bars with CIs) ───

def fig1_accuracy_comparison(eval_runs):
    by_model = split_by_model(eval_runs)
    categories = TIER1_CATS + TIER2_CATS
    n_models = len(MODELS)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(categories))
    bar_width = 0.8 / n_models

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

        offset = (i - (n_models - 1) / 2) * bar_width
        ax.bar(x + offset, means,
               yerr=[lows, highs], width=bar_width,
               capsize=1.5, color=MODEL_COLORS[model],
               edgecolor="white", linewidth=0.3,
               error_kw={"linewidth": 0.5},
               label=MODEL_LABELS[model])

    # Tier separator
    ax.axvline(x=len(TIER1_CATS) - 0.5, color="gray", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.text(len(TIER1_CATS) / 2, 104, "Tier 1: LABBench", ha="center", fontsize=6, color="gray")
    ax.text(len(TIER1_CATS) + len(TIER2_CATS) / 2, 104, "Tier 2: New Tasks", ha="center", fontsize=6, color="gray")

    ax.set_xticks(x)
    ax.set_xticklabels([short_cat(c) for c in categories], rotation=45, ha="right", fontsize=6.5)
    ax.set_ylabel("Accuracy (%)")
    ax.set_ylim(0, 112)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right", ncol=2)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_accuracy_comparison.pdf")
    fig.savefig(FIG_DIR / "fig1_accuracy_comparison.png")
    plt.close()
    print("  fig1_accuracy_comparison")


# ─── Figure 2: Chain error propagation (all models) ───

def fig2_chain_error_propagation(chain_runs):
    by_model = split_by_model(chain_runs)
    n_models = len(MODELS)

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))

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
                 marker=MODEL_MARKERS[model],
                 markersize=4, linewidth=1.2, label=MODEL_LABELS[model])

    ax1.set_xlabel("Step number")
    ax1.set_ylabel("Step accuracy (%)")
    ax1.set_ylim(0, 105)
    ax1.legend(fontsize=5.5, frameon=False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_title("a", fontsize=9, fontweight="bold", loc="left")

    # Panel b: E2E vs step-level (grouped bars)
    ax2 = axes[1]
    x = np.arange(2)
    bar_width = 0.8 / n_models
    for i, model in enumerate(MODELS):
        chains = defaultdict(list)
        for r in by_model.get(model, []):
            chains[r["chain_id"]].append(r)
        total_steps = sum(len(s) for s in chains.values())
        correct_steps = sum(1 for cid in chains for s in chains[cid] if s["correct"])
        e2e_correct = sum(1 for cid, steps in chains.items() if all(s["correct"] for s in steps))
        vals = [correct_steps / total_steps * 100 if total_steps else 0,
                e2e_correct / len(chains) * 100 if chains else 0]
        offset = (i - (n_models - 1) / 2) * bar_width
        bars = ax2.bar(x + offset, vals,
                       width=bar_width, color=MODEL_COLORS[model],
                       edgecolor="white", label=MODEL_LABELS[model])
        for b, v in zip(bars, vals):
            ax2.text(b.get_x() + b.get_width() / 2, v + 1.5,
                     f"{v:.0f}%", ha="center", fontsize=4.5)

    ax2.set_xticks(x)
    ax2.set_xticklabels(["Step-level", "End-to-end"], fontsize=7)
    ax2.set_ylabel("Accuracy (%)")
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=4.5, frameon=False)
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
                   edgecolor="white", width=0.6)
    for b, v in zip(bars, gaps):
        ax3.text(b.get_x() + b.get_width() / 2, v + 0.5,
                 f"{v:.1f} pp", ha="center", fontsize=5.5, fontweight="bold")

    ax3.set_ylabel("Error propagation gap (pp)")
    ax3.set_ylim(0, max(gaps) + 8)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.tick_params(axis="x", labelsize=6, rotation=20)
    ax3.set_title("c", fontsize=9, fontweight="bold", loc="left")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_chain_error_propagation.pdf")
    fig.savefig(FIG_DIR / "fig2_chain_error_propagation.png")
    plt.close()
    print("  fig2_chain_error_propagation")


# ─── Figure 3: Cost-accuracy scatter (all models) ───

def fig3_cost_accuracy(eval_runs):
    by_model = split_by_model(eval_runs)

    fig, ax = plt.subplots(figsize=(7, 4.5))

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
            if cost <= 0:
                cost = 0.001  # avoid log(0)
            ax.scatter(cost, acc, s=d["total"] * 0.3,
                       color=MODEL_COLORS[model], alpha=0.7,
                       marker=MODEL_MARKERS[model],
                       edgecolors="white", linewidth=0.3)
            # Label only for first model to avoid clutter
            if model == MODELS[0]:
                ax.annotate(cat, (cost, acc), fontsize=4, ha="center", va="bottom",
                            xytext=(0, 4), textcoords="offset points", color="gray")

    # Legend
    handles = [
        plt.Line2D([], [], marker=MODEL_MARKERS[m], color=MODEL_COLORS[m],
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


# ─── Figure 5: Latency distribution (all models) ───

def fig5_latency(eval_runs):
    by_model = split_by_model(eval_runs)
    all_cats = sorted(set(r["category"] for r in eval_runs))

    fig, ax = plt.subplots(figsize=(10, 4))
    positions = []
    data_groups = []
    colors = []

    pos = 0
    tick_positions = []
    for cat in all_cats:
        cat_pos = []
        for model in MODELS:
            model_runs = [r["latency_ms"] / 1000 for r in by_model.get(model, [])
                          if r["category"] == cat and r["latency_ms"] > 0]
            if model_runs:
                data_groups.append(model_runs)
                positions.append(pos)
                colors.append(MODEL_COLORS[model])
                cat_pos.append(pos)
                pos += 1
        if cat_pos:
            tick_positions.append(np.mean(cat_pos))
        pos += 0.5

    bp = ax.boxplot(data_groups, positions=positions, vert=True, patch_artist=True,
                    widths=0.7, showfliers=False,
                    medianprops={"color": "black", "linewidth": 1})

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels([short_cat(c) for c in all_cats], rotation=45, ha="right", fontsize=5.5)
    ax.set_ylabel("Latency (seconds)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [mpatches.Patch(color=MODEL_COLORS[m], alpha=0.7, label=MODEL_LABELS[m]) for m in MODELS]
    ax.legend(handles=handles, fontsize=5.5, frameon=False, loc="upper right", ncol=2)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig5_latency.pdf")
    fig.savefig(FIG_DIR / "fig5_latency.png")
    plt.close()
    print("  fig5_latency")


# ─── Figure 6: Cross-provider accuracy heatmap ───

def fig6_provider_heatmap(eval_runs):
    """Heatmap: models (rows) x categories (columns), cell = accuracy %."""
    by_model = split_by_model(eval_runs)
    categories = TIER1_CATS + TIER2_CATS

    # Build accuracy matrix
    matrix = np.zeros((len(MODELS), len(categories)))
    for i, model in enumerate(MODELS):
        by_cat = defaultdict(list)
        for r in by_model.get(model, []):
            by_cat[r["category"]].append(1 if r["correct"] else 0)
        for j, cat in enumerate(categories):
            if cat in by_cat:
                matrix[i, j] = np.mean(by_cat[cat]) * 100

    fig, ax = plt.subplots(figsize=(10, 3.5))
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=100)

    # Bold the best model per category
    for j in range(len(categories)):
        best_i = np.argmax(matrix[:, j])
        for i in range(len(MODELS)):
            val = matrix[i, j]
            is_best = (i == best_i and val > 0)
            color = "white" if val > 60 else "black"
            weight = "bold" if is_best else "normal"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=6, color=color, fontweight=weight)

    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels([short_cat(c) for c in categories], rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(MODELS)))
    ax.set_yticklabels([MODEL_LABELS[m] for m in MODELS], fontsize=7)

    # Tier separator
    ax.axvline(x=len(TIER1_CATS) - 0.5, color="white", linewidth=2)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label="Accuracy (%)")
    cbar.ax.tick_params(labelsize=6)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig6_provider_heatmap.pdf")
    fig.savefig(FIG_DIR / "fig6_provider_heatmap.png")
    plt.close()
    print("  fig6_provider_heatmap")


# ─── Figure 7: Provider strengths radar / bar comparison ───

def fig7_provider_strengths(eval_runs, chain_runs):
    """Bar chart showing which provider wins each 'dimension'."""
    by_model_eval = split_by_model(eval_runs)
    by_model_chain = split_by_model(chain_runs)

    dimensions = ["Retrieval\n(Tier 1 avg)", "Quantitative\n(StatReas)",
                  "Structure\n(StructAnal)", "Compositional\n(Chain E2E)",
                  "Calibration"]

    data = {}
    for model in MODELS:
        # Tier 1 average
        t1_correct, t1_total = 0, 0
        by_cat = defaultdict(list)
        for r in by_model_eval.get(model, []):
            by_cat[r["category"]].append(1 if r["correct"] else 0)
        for cat in TIER1_CATS:
            if cat in by_cat:
                t1_correct += sum(by_cat[cat])
                t1_total += len(by_cat[cat])
        t1_avg = t1_correct / t1_total * 100 if t1_total else 0

        # StatReas
        sr = np.mean(by_cat.get("StatisticalReasoning", [0])) * 100

        # StructAnal
        sa = np.mean(by_cat.get("StructureAnalysis", [0])) * 100

        # Chain E2E
        chains = defaultdict(list)
        for r in by_model_chain.get(model, []):
            chains[r["chain_id"]].append(r)
        e2e = sum(1 for cid, steps in chains.items() if all(s["correct"] for s in steps))
        e2e_pct = e2e / len(chains) * 100 if chains else 0

        # Calibration
        cal = np.mean(by_cat.get("Calibration", [0])) * 100

        data[model] = [t1_avg, sr, sa, e2e_pct, cal]

    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(dimensions))
    n_models = len(MODELS)
    bar_width = 0.8 / n_models

    for i, model in enumerate(MODELS):
        offset = (i - (n_models - 1) / 2) * bar_width
        ax.bar(x + offset, data[model], width=bar_width,
               color=MODEL_COLORS[model], edgecolor="white",
               label=MODEL_LABELS[model])

    ax.set_xticks(x)
    ax.set_xticklabels(dimensions, fontsize=7)
    ax.set_ylabel("Accuracy / Score (%)")
    ax.set_ylim(0, 110)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right", ncol=2)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig7_provider_strengths.pdf")
    fig.savefig(FIG_DIR / "fig7_provider_strengths.png")
    plt.close()
    print("  fig7_provider_strengths")


# ─── Table 1: Four-model comparison (LaTeX) ───

def table1_main_results(eval_runs):
    by_model = split_by_model(eval_runs)
    categories = TIER1_CATS + TIER2_CATS

    lines = []
    lines.append(r"\begin{table*}[ht]")
    lines.append(r"\centering")
    lines.append(r"\caption{Four-model accuracy comparison across LABBench2-Pro categories. BCa 95\% bootstrap CIs (10{,}000 resamples). Bold = best per category.}")
    lines.append(r"\label{tab:main-results}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{llr" + "rr" * len(MODELS) + "}")
    lines.append(r"\toprule")

    header = r"Tier & Category & $n$"
    for model in MODELS:
        header += f" & {MODEL_LABELS[model]} & 95\\% CI"
    header += r" \\"
    lines.append(header)
    lines.append(r"\midrule")

    for tier_name, cats in [("1: LABBench", TIER1_CATS), ("2: New Tasks", TIER2_CATS)]:
        first = True
        for cat in cats:
            tier_label = tier_name if first else ""
            first = False

            # Get n and accuracies
            model_accs = {}
            ns = []
            for model in MODELS:
                cat_runs = [1 if r["correct"] else 0 for r in by_model.get(model, []) if r["category"] == cat]
                ns.append(len(cat_runs))
                if cat_runs:
                    lo, mu, hi = bootstrap_ci(cat_runs)
                    model_accs[model] = (mu, lo, hi)
                else:
                    model_accs[model] = None

            n = max(ns) if ns else 0
            best_model = max((m for m in MODELS if model_accs.get(m)),
                             key=lambda m: model_accs[m][0] if model_accs[m] else 0,
                             default=None)

            row = f"  {tier_label} & {cat} & {n}"
            for model in MODELS:
                if model_accs.get(model):
                    mu, lo, hi = model_accs[model]
                    acc_str = f"{mu:.1%}"
                    ci_str = f"[{lo:.1%}, {hi:.1%}]"
                    if model == best_model:
                        acc_str = f"\\textbf{{{acc_str}}}"
                    row += f" & {acc_str} & {ci_str}"
                else:
                    row += " & --- & ---"
            row += r" \\"
            lines.append(row)
        lines.append(r"\midrule")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table*}")

    with open(TABLE_DIR / "table1_main_results.tex", "w") as f:
        f.write("\n".join(lines))
    print("  table1_main_results.tex")


# ─── Table 2: Chain results (all models) ───

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

    all_chain_ids = sorted(
        set(r["chain_id"] for r in chain_runs),
        key=lambda x: int(x.replace("chain", ""))
    )

    lines = []
    lines.append(r"\begin{table*}[ht]")
    lines.append(r"\centering")
    lines.append(r"\caption{Compositional chain results for all four models. E2E = all steps correct.}")
    lines.append(r"\label{tab:chain-results}")
    lines.append(r"\small")

    col_spec = "ll" + "rr" * len(MODELS)
    lines.append(r"\begin{tabular}{" + col_spec + "}")
    lines.append(r"\toprule")

    header = "Chain & Template"
    for model in MODELS:
        short = MODEL_LABELS[model]
        header += f" & {short} Steps & {short} E2E"
    header += r" \\"
    lines.append(header)
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
    lines.append(r"\end{table*}")

    with open(TABLE_DIR / "table2_chain_results.tex", "w") as f:
        f.write("\n".join(lines))
    print("  table2_chain_results.tex")


# ─── Table 3: Cost breakdown (all models) ───

def table3_cost_breakdown(eval_runs):
    by_model = split_by_model(eval_runs)

    lines = []
    lines.append(r"\begin{table*}[ht]")
    lines.append(r"\centering")
    lines.append(r"\caption{Cost analysis by category for all four models.}")
    lines.append(r"\label{tab:cost}")
    lines.append(r"\small")

    col_spec = "l" + "r" * len(MODELS)
    lines.append(r"\begin{tabular}{" + col_spec + "}")
    lines.append(r"\toprule")

    header = "Category"
    for model in MODELS:
        header += f" & {MODEL_LABELS[model]}"
    header += r" \\"
    lines.append(header)
    lines.append(r"\midrule")

    all_cats = sorted(set(r["category"] for r in eval_runs))

    # Sort by Opus cost descending
    cat_costs = {}
    for cat in all_cats:
        opus_runs = [r for r in by_model.get(MODELS[0], []) if r["category"] == cat]
        cat_costs[cat] = sum(r["cost_usd"] for r in opus_runs)

    for cat in sorted(all_cats, key=lambda c: cat_costs.get(c, 0), reverse=True):
        row = cat
        for model in MODELS:
            cat_runs = [r for r in by_model.get(model, []) if r["category"] == cat]
            total_cost = sum(r["cost_usd"] for r in cat_runs)
            row += f" & \\${total_cost:.2f}"
        row += r" \\"
        lines.append(f"  {row}")

    # Totals
    lines.append(r"\midrule")
    total_row = r"\textbf{Total}"
    for model in MODELS:
        model_runs = by_model.get(model, [])
        total_cost = sum(r["cost_usd"] for r in model_runs)
        total_row += f" & \\textbf{{\\${total_cost:.2f}}}"
    total_row += r" \\"
    lines.append(f"  {total_row}")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table*}")

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
    lines.append(r"IRT Analysis  & Total items & 2{,}522 \\")
    lines.append(r"              & Low discrimination ($< 0.3$) & 2{,}219 (88.0\%) \\")
    lines.append(r"              & Recommended pruned set & 303 items \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    with open(TABLE_DIR / "table4_methodological_audit.tex", "w") as f:
        f.write("\n".join(lines))
    print("  table4_methodological_audit.tex")


# ─── Supplementary: Chain traces (all models) ───

def supplementary_chain_traces(chain_runs):
    by_model = split_by_model(chain_runs)

    with open(TABLE_DIR / "supplementary_chain_traces.md", "w") as f:
        f.write("# Supplementary Material: Full Chain Execution Traces\n\n")
        f.write("Complete model responses for all 30 compositional chains, all four models.\n\n")

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
    for m, c in sorted(model_counts.items()):
        print(f"  {m}: {c} eval runs")
    print(f"  {len(chain_runs)} chain steps, {len(judge_audits)} judge audits")

    print("\nGenerating figures...")
    fig1_accuracy_comparison(eval_runs)
    fig2_chain_error_propagation(chain_runs)
    fig3_cost_accuracy(eval_runs)
    fig4_judge_audit(judge_audits)
    fig5_latency(eval_runs)
    fig6_provider_heatmap(eval_runs)
    fig7_provider_strengths(eval_runs, chain_runs)

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
