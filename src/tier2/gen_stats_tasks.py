"""Generate statistical reasoning tasks from gene expression data.

Usage: python -m src.tier2.gen_stats_tasks --output-dir tasks/stats_reasoning --count 200
"""

import argparse
import json
import os
import random

import numpy as np
from scipy import stats


def generate_test_selection_task(task_id: int) -> dict:
    """Which statistical test is appropriate for this design?"""
    designs = [
        {
            "design": "Two groups (treatment vs control), continuous outcome, normal distribution, equal variance",
            "answer": "Independent samples t-test",
            "n1": random.randint(10, 50),
            "n2": random.randint(10, 50),
        },
        {
            "design": "Two groups, continuous outcome, non-normal distribution",
            "answer": "Mann-Whitney U test (Wilcoxon rank-sum)",
            "n1": random.randint(10, 50),
            "n2": random.randint(10, 50),
        },
        {
            "design": "Same subjects measured before and after treatment, continuous outcome, normal",
            "answer": "Paired t-test",
            "n1": random.randint(10, 50),
            "n2": None,
        },
        {
            "design": "Three or more groups, continuous outcome, normal distribution",
            "answer": "One-way ANOVA",
            "n1": random.randint(10, 30),
            "n2": random.randint(3, 5),
        },
        {
            "design": "Two categorical variables, testing independence",
            "answer": "Chi-squared test of independence",
            "n1": random.randint(50, 200),
            "n2": None,
        },
        {
            "design": "Two groups, small sample size (<5 per cell), categorical outcome",
            "answer": "Fisher's exact test",
            "n1": random.randint(3, 8),
            "n2": random.randint(3, 8),
        },
        {
            "design": "Correlation between two continuous variables, normal distribution",
            "answer": "Pearson correlation",
            "n1": random.randint(20, 100),
            "n2": None,
        },
        {
            "design": "Correlation between two continuous variables, non-normal or ordinal",
            "answer": "Spearman rank correlation",
            "n1": random.randint(20, 100),
            "n2": None,
        },
    ]

    d = random.choice(designs)
    return {
        "id": f"pro_stats_{task_id:04d}",
        "category": "StatisticalReasoning",
        "question": f"A researcher has the following experimental design: {d['design']}. "
                     f"Sample sizes: n={d['n1']}" + (f", groups={d['n2']}" if d['n2'] else "") + ". "
                     f"Which statistical test is most appropriate?",
        "ideal": d["answer"],
        "verification": "programmatic",
        "verification_fn": "exact_match",
        "source": "pro",
        "meta": {"template": "test_selection", "difficulty": "easy"},
    }


def generate_pvalue_task(task_id: int) -> dict:
    """Compute the p-value for differential expression."""
    n1, n2 = random.randint(8, 30), random.randint(8, 30)
    # Generate two groups with known effect
    effect = random.choice([0.5, 1.0, 1.5, 2.0, 3.0])
    np.random.seed(task_id)
    group1 = np.random.normal(loc=5.0, scale=1.5, size=n1)
    group2 = np.random.normal(loc=5.0 + effect, scale=1.5, size=n2)

    t_stat, p_val = stats.ttest_ind(group1, group2)
    gene_name = f"GENE{random.randint(100, 999)}"

    # Format data as CSV-like string
    g1_str = ", ".join(f"{v:.2f}" for v in group1)
    g2_str = ", ".join(f"{v:.2f}" for v in group2)

    return {
        "id": f"pro_stats_{task_id:04d}",
        "category": "StatisticalReasoning",
        "question": f"Gene {gene_name} expression values:\n"
                     f"Control (n={n1}): {g1_str}\n"
                     f"Treatment (n={n2}): {g2_str}\n\n"
                     f"Compute the p-value using an independent samples t-test. "
                     f"Report the p-value to 4 decimal places.",
        "ideal": f"{p_val:.4f}",
        "verification": "programmatic",
        "verification_fn": "numeric_tolerance",
        "source": "pro",
        "meta": {"template": "p_value", "difficulty": "medium", "tolerance": 0.01},
    }


def generate_sig_count_task(task_id: int) -> dict:
    """How many genes are significant at FDR < 0.05?"""
    n_genes = random.choice([20, 50, 100])
    n_samples_per_group = random.randint(5, 15)

    np.random.seed(task_id)
    # Some genes are DE, most are not
    n_de = random.randint(2, n_genes // 3)
    p_values = []
    gene_names = []

    for i in range(n_genes):
        gene_names.append(f"Gene_{i+1}")
        if i < n_de:
            # True DE: low p-value
            p_values.append(np.random.uniform(0.0001, 0.04))
        else:
            p_values.append(np.random.uniform(0.01, 0.99))

    # BH correction
    p_arr = np.array(p_values)
    n = len(p_arr)
    sorted_idx = np.argsort(p_arr)
    sorted_p = p_arr[sorted_idx]
    bh_critical = np.arange(1, n + 1) / n * 0.05
    significant = sorted_p <= bh_critical
    # Find the largest k where p(k) <= k/n * alpha
    if significant.any():
        max_k = np.max(np.where(significant)[0])
        sig_count = max_k + 1
    else:
        sig_count = 0

    # Format as table
    table_lines = ["Gene | p-value"]
    table_lines.append("--- | ---")
    for name, pv in zip(gene_names, p_values):
        table_lines.append(f"{name} | {pv:.4f}")
    table = "\n".join(table_lines)

    return {
        "id": f"pro_stats_{task_id:04d}",
        "category": "StatisticalReasoning",
        "question": f"Given the following p-values from differential expression analysis "
                     f"of {n_genes} genes:\n\n{table}\n\n"
                     f"How many genes are significant after Benjamini-Hochberg FDR correction "
                     f"at alpha = 0.05? Report the exact count.",
        "ideal": str(sig_count),
        "verification": "programmatic",
        "verification_fn": "integer_match",
        "source": "pro",
        "meta": {"template": "sig_count", "difficulty": "hard", "n_genes": n_genes},
    }


def main():
    parser = argparse.ArgumentParser(description="Generate statistical reasoning tasks")
    parser.add_argument("--output-dir", default="tasks/stats_reasoning")
    parser.add_argument("--count", type=int, default=200)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    generators = [generate_test_selection_task, generate_pvalue_task, generate_sig_count_task]
    tasks = []

    for i in range(args.count):
        gen = generators[i % len(generators)]
        task = gen(i + 1)
        tasks.append(task)

        filepath = os.path.join(args.output_dir, f"{task['id']}.json")
        with open(filepath, "w") as f:
            json.dump(task, f, indent=2)

    print(f"Generated {len(tasks)} tasks in {args.output_dir}")
    by_template = {}
    for t in tasks:
        tpl = t["meta"]["template"]
        by_template[tpl] = by_template.get(tpl, 0) + 1
    for tpl, count in by_template.items():
        print(f"  {tpl}: {count}")


if __name__ == "__main__":
    main()
