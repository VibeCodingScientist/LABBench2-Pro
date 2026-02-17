"""Generate hypothesis generation tasks from PubMed Central.

Usage: python -m src.tier2.gen_hypothesis --output-dir tasks/hypothesis --count 100
"""

import argparse
import json
import os
import time

from Bio import Entrez


HYPOTHESIS_RUBRIC = """Evaluate the proposed hypotheses on these criteria:
1. Falsifiable: Each hypothesis must be testable with a concrete experiment
2. Connected: Hypotheses must logically follow from the presented finding
3. Novel: Hypotheses should suggest non-obvious extensions of the work
4. Specific: Hypotheses should name specific genes, pathways, or mechanisms

Score: 0-4 points (one per criterion met). A score of 3+ is "correct"."""


def fetch_papers(email: str, count: int = 200) -> list[dict]:
    """Fetch papers with unexpected findings from PubMed Central."""
    Entrez.email = email

    # Search for papers with interesting findings
    query = '"unexpected finding" OR "surprising result" OR "contrary to expectations"'
    handle = Entrez.esearch(db="pmc", term=query, retmax=count, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        return []

    papers = []
    # Fetch in batches of 10
    for batch_start in range(0, len(ids), 10):
        batch_ids = ids[batch_start:batch_start + 10]
        try:
            handle = Entrez.efetch(db="pmc", id=",".join(batch_ids), rettype="xml")
            # Parse just the abstracts for simplicity
            from Bio import Medline
            # Use abstract fetching instead
            handle.close()
        except Exception:
            pass
        time.sleep(0.5)  # Respect NCBI rate limits

    # Fallback: fetch abstracts via pubmed
    handle = Entrez.esearch(db="pubmed", term=query, retmax=count, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    pubmed_ids = record.get("IdList", [])
    for batch_start in range(0, min(len(pubmed_ids), count), 10):
        batch = pubmed_ids[batch_start:batch_start + 10]
        try:
            handle = Entrez.efetch(db="pubmed", id=",".join(batch), rettype="abstract", retmode="text")
            text = handle.read()
            handle.close()

            # Split into individual abstracts
            abstracts = text.strip().split("\n\n\n")
            for abstract in abstracts:
                if len(abstract) > 100:
                    papers.append({"text": abstract.strip(), "pmid": batch[0] if batch else "unknown"})
        except Exception:
            pass
        time.sleep(0.5)

    return papers[:count]


def generate_from_abstract(task_id: int, paper: dict) -> dict:
    """Generate a hypothesis task from a paper abstract."""
    text = paper["text"]
    # Truncate to a reasonable length
    if len(text) > 2000:
        text = text[:2000] + "..."

    return {
        "id": f"pro_hyp_{task_id:04d}",
        "category": "HypothesisGeneration",
        "question": f"Based on the following research finding, propose 3 testable hypotheses "
                     f"that could extend this work. Each hypothesis should be falsifiable, "
                     f"connected to the finding, novel, and specific.\n\n"
                     f"Finding:\n{text}",
        "ideal": "Evaluated by rubric: falsifiable, connected, novel, specific (3+/4 = correct)",
        "verification": "llm-judge",
        "source": "pro",
        "meta": {
            "template": "hypothesis",
            "pmid": paper.get("pmid", "unknown"),
            "rubric": HYPOTHESIS_RUBRIC,
            "difficulty": "hard",
        },
    }


def generate_placeholder_tasks(count: int) -> list[dict]:
    """Generate placeholder tasks when PubMed is unavailable."""
    findings = [
        "CRISPR-Cas9 editing of the FOXP2 gene in organoids led to unexpected enhancement of synaptic plasticity markers.",
        "Single-cell RNA sequencing revealed a novel population of immune cells in the tumor microenvironment that express both T-cell and macrophage markers.",
        "Gut microbiome transplant from centenarians to germ-free mice resulted in increased lifespan and reduced inflammatory markers.",
        "Protein folding simulations suggest that intrinsically disordered regions of p53 adopt transient alpha-helical structures under oxidative stress.",
        "Metagenomic analysis of deep-sea hydrothermal vent samples identified bacteria using a previously unknown carbon fixation pathway.",
        "Long-read sequencing of the human Y chromosome revealed extensive structural variation between populations that correlates with fertility metrics.",
        "Optogenetic stimulation of astrocytes, rather than neurons, was sufficient to restore memory formation in an Alzheimer's disease mouse model.",
        "ATAC-seq data from developing embryos shows that enhancer accessibility patterns predict cell fate 48 hours before transcriptional changes.",
        "A screen of 10,000 natural compounds identified a moss-derived molecule that selectively inhibits cancer stem cell self-renewal.",
        "Spatial transcriptomics of the aging brain reveals that glial cell heterogeneity increases while neuronal diversity decreases.",
    ]

    tasks = []
    for i in range(count):
        finding = findings[i % len(findings)]
        tasks.append({
            "id": f"pro_hyp_{i+1:04d}",
            "category": "HypothesisGeneration",
            "question": f"Based on the following research finding, propose 3 testable hypotheses "
                         f"that could extend this work. Each hypothesis should be falsifiable, "
                         f"connected to the finding, novel, and specific.\n\n"
                         f"Finding: {finding}",
            "ideal": "Evaluated by rubric: falsifiable, connected, novel, specific (3+/4 = correct)",
            "verification": "llm-judge",
            "source": "pro",
            "meta": {"template": "hypothesis", "rubric": HYPOTHESIS_RUBRIC, "difficulty": "hard"},
        })
    return tasks


def main():
    parser = argparse.ArgumentParser(description="Generate hypothesis tasks")
    parser.add_argument("--output-dir", default="tasks/hypothesis")
    parser.add_argument("--count", type=int, default=100)
    parser.add_argument("--email", default="labbench2pro@example.com", help="Email for NCBI Entrez")
    parser.add_argument("--use-placeholders", action="store_true", help="Skip PubMed, use placeholder findings")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.use_placeholders:
        tasks = generate_placeholder_tasks(args.count)
    else:
        print(f"Fetching papers from PubMed (email={args.email})...")
        try:
            papers = fetch_papers(args.email, args.count * 2)
            print(f"Fetched {len(papers)} papers")
        except Exception as e:
            print(f"PubMed fetch failed ({e}), using placeholders...")
            tasks = generate_placeholder_tasks(args.count)
            papers = []

        if papers:
            tasks = []
            for i, paper in enumerate(papers[:args.count]):
                tasks.append(generate_from_abstract(i + 1, paper))

    for task in tasks:
        filepath = os.path.join(args.output_dir, f"{task['id']}.json")
        with open(filepath, "w") as f:
            json.dump(task, f, indent=2)

    print(f"Generated {len(tasks)} hypothesis tasks in {args.output_dir}")


if __name__ == "__main__":
    main()
