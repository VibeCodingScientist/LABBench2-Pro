"""Generate hypothesis generation tasks from PubMed abstracts.

Fetches real papers via NCBI Entrez API. Requires NCBI_EMAIL env var.

Usage: python -m src.tier2.gen_hypothesis --output-dir tasks/hypothesis --count 100
"""

import argparse
import json
import os
import time

from Bio import Entrez, Medline

from src.config import NCBI_EMAIL

HYPOTHESIS_RUBRIC = """Evaluate the proposed hypotheses on these criteria:
1. Falsifiable: Each hypothesis must be testable with a concrete experiment
2. Connected: Hypotheses must logically follow from the presented finding
3. Novel: Hypotheses should suggest non-obvious extensions of the work
4. Specific: Hypotheses should name specific genes, pathways, or mechanisms

Score: 0-4 points (one per criterion met). A score of 3+ is "correct"."""

# Diverse biology search queries to get varied papers
SEARCH_QUERIES = [
    '"unexpected finding" AND biology[MeSH]',
    '"surprising result" AND molecular[Title]',
    '"novel mechanism" AND cell biology[MeSH]',
    '"contrary to expectations" AND genetics[MeSH]',
    '"previously unknown" AND biochemistry[MeSH]',
    '"first evidence" AND neuroscience[MeSH]',
    '"unexpected role" AND immunology[MeSH]',
    '"novel pathway" AND cancer[MeSH]',
    '"surprising discovery" AND microbiology[MeSH]',
    '"previously uncharacterized" AND genomics[MeSH]',
]


def fetch_abstracts(email: str, count: int) -> list[dict]:
    """Fetch paper abstracts from PubMed via Entrez."""
    Entrez.email = email
    papers = []
    per_query = max(count // len(SEARCH_QUERIES), 15)

    for query in SEARCH_QUERIES:
        if len(papers) >= count:
            break

        try:
            # Search
            handle = Entrez.esearch(db="pubmed", term=query, retmax=per_query, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            if not ids:
                continue

            time.sleep(0.4)

            # Fetch as Medline format (easy to parse)
            handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
            records = list(Medline.parse(handle))
            handle.close()

            for rec in records:
                abstract = rec.get("AB", "")
                title = rec.get("TI", "")
                pmid = rec.get("PMID", "unknown")
                authors = rec.get("AU", [])
                journal = rec.get("JT", "")
                year = rec.get("DP", "").split()[0] if rec.get("DP") else ""

                if len(abstract) < 200:
                    continue

                papers.append({
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "authors": authors[:3],
                    "journal": journal,
                    "year": year,
                })

            time.sleep(0.4)

        except Exception as e:
            print(f"  Warning: query '{query[:40]}...' failed: {e}")
            time.sleep(1)
            continue

    return papers[:count]


def generate_task(task_id: int, paper: dict) -> dict:
    """Generate a hypothesis task from a paper abstract."""
    # Build a rich context from the abstract
    context = f"Title: {paper['title']}\n"
    if paper.get("journal"):
        context += f"Journal: {paper['journal']}"
        if paper.get("year"):
            context += f" ({paper['year']})"
        context += "\n"
    context += f"\nAbstract:\n{paper['abstract']}"

    if len(context) > 2500:
        context = context[:2500] + "..."

    return {
        "id": f"pro_hyp_{task_id:04d}",
        "category": "HypothesisGeneration",
        "question": (
            "Based on the following research paper, propose 3 testable hypotheses "
            "that could extend this work. Each hypothesis should be:\n"
            "- Falsifiable (testable with a concrete experiment)\n"
            "- Connected (logically following from the findings)\n"
            "- Novel (non-obvious extensions)\n"
            "- Specific (naming genes, pathways, or mechanisms)\n\n"
            f"{context}"
        ),
        "ideal": "Evaluated by rubric: falsifiable, connected, novel, specific (3+/4 = correct)",
        "verification": "llm-judge",
        "source": "pro",
        "meta": {
            "template": "hypothesis",
            "pmid": paper["pmid"],
            "title": paper["title"],
            "rubric": HYPOTHESIS_RUBRIC,
            "difficulty": "hard",
        },
    }


def main():
    parser = argparse.ArgumentParser(description="Generate hypothesis tasks from PubMed")
    parser.add_argument("--output-dir", default="tasks/hypothesis")
    parser.add_argument("--count", type=int, default=100)
    args = parser.parse_args()

    email = NCBI_EMAIL
    if not email:
        print("NCBI_EMAIL not set. Using default (set it in .env for better rate limits).")
        email = "labbench2pro@research.org"

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Fetching ~{args.count} abstracts from PubMed (email={email})...")
    papers = fetch_abstracts(email, args.count * 2)  # fetch extra in case some are too short
    print(f"Fetched {len(papers)} usable abstracts")

    if not papers:
        print("ERROR: Could not fetch any papers from PubMed. Check network connectivity.", flush=True)
        raise SystemExit(1)

    tasks = []
    for i, paper in enumerate(papers[:args.count]):
        tasks.append(generate_task(i + 1, paper))

    for task in tasks:
        filepath = os.path.join(args.output_dir, f"{task['id']}.json")
        with open(filepath, "w") as f:
            json.dump(task, f, indent=2)

    print(f"Generated {len(tasks)} hypothesis tasks in {args.output_dir}")


if __name__ == "__main__":
    main()
