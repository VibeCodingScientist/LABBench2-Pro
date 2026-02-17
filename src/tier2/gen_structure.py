"""Generate PDB structure + gel image tasks via BioPython + matplotlib.

Downloads real PDB files, parses them for facts, generates questions with
deterministic ground truth. No hardcoded data.

Usage: python -m src.tier2.gen_structure --output-dir tasks/structures --count 150
"""

import argparse
import io
import json
import os
import random
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from Bio.PDB import PDBList, PDBParser, DSSP

# Well-known PDB IDs spanning diverse biology — all guaranteed to exist in RCSB
PDB_IDS = [
    "1MBO",  # Myoglobin
    "4HHB",  # Hemoglobin
    "1EMA",  # GFP
    "1TIM",  # Triosephosphate isomerase (TIM barrel)
    "2PTC",  # Trypsin-BPTI complex
    "1CRN",  # Crambin (ultra-high resolution)
    "1UBQ",  # Ubiquitin
    "1A2Y",  # Insulin
    "1LYZ",  # Lysozyme
    "3CLN",  # Calmodulin
    "1ATP",  # cAMP-dependent protein kinase
    "1GZM",  # Alpha-hemolysin
    "1BRS",  # Barnase-barstar complex
    "1A3N",  # Oxy T-state hemoglobin
    "2AW4",  # Aquaporin
    "1FME",  # Villin headpiece
    "1L2Y",  # Trp-cage
    "1IGT",  # IgG1 antibody
    "3PQR",  # GFP superfolder
    "1AKE",  # Adenylate kinase
    "1BNA",  # B-DNA dodecamer
    "1THR",  # Thrombin
    "1HEL",  # Hen egg-white lysozyme
    "1PPE",  # Porcine pancreatic elastase
    "1RN1",  # Ribonuclease A
    "1HOE",  # Alpha-amylase inhibitor
    "2LZM",  # T4 lysozyme
    "1OVA",  # Ovalbumin
    "3HHR",  # Human growth hormone receptor
    "1GCN",  # Glucagon
    "1BEB",  # Beta-lactoglobulin
    "4INS",  # Insulin (porcine)
    "1AZI",  # Azurin
    "1AB1",  # Alpha-bungarotoxin
    "1CON",  # Concanavalin A
    "2SOD",  # Superoxide dismutase
    "1CA2",  # Carbonic anhydrase II
    "1REI",  # Bence-Jones protein
    "3LZM",  # T4 lysozyme mutant
    "1CHO",  # Alpha-chymotrypsin
    "5CHA",  # Chymotrypsinogen A
    "7RSA",  # Ribonuclease A (high res)
    "1ACX",  # Actinidin
    "2RNT",  # RNase T1
    "1CCR",  # Cytochrome c
    "1MBN",  # Metmyoglobin
    "2CAB",  # Carbonic anhydrase B
    "1EST",  # Elastase
    "1PHT",  # Phosphotransferase
    "1STP",  # Streptavidin
]


def download_and_parse_pdb(pdb_id: str, tmp_dir: str) -> dict | None:
    """Download a PDB file and extract structural facts."""
    pdbl = PDBList()
    parser = PDBParser(QUIET=True)

    try:
        filename = pdbl.retrieve_pdb_file(pdb_id, pdir=tmp_dir, file_format="pdb")
        structure = parser.get_structure(pdb_id, filename)
    except Exception as e:
        print(f"  Warning: could not fetch {pdb_id}: {e}")
        return None

    model = structure[0]

    # Count chains and residues
    chains = list(model.get_chains())
    n_chains = len(chains)
    residues = [r for r in model.get_residues() if r.id[0] == " "]  # standard residues only
    n_residues = len(residues)

    # Header info
    header = structure.header
    name = header.get("name", pdb_id).strip().title() or pdb_id
    method = header.get("structure_method", "unknown").strip()
    resolution = header.get("resolution")
    organism = header.get("source", {})
    if isinstance(organism, dict):
        organism = organism.get("1", {}).get("organism_scientific", "unknown")
    elif isinstance(organism, list) and organism:
        organism = organism[0].get("organism_scientific", "unknown") if isinstance(organism[0], dict) else "unknown"
    else:
        organism = "unknown"

    return {
        "pdb_id": pdb_id,
        "name": name if name != pdb_id else pdb_id,
        "chains": n_chains,
        "residues": n_residues,
        "resolution": round(resolution, 2) if resolution else None,
        "method": method if method != "unknown" else "X-ray diffraction",
        "organism": organism,
    }


QUESTION_TEMPLATES = [
    ("chains", "How many polypeptide chains does the protein structure {pdb_id} ({name}) contain? Answer with just the number.",
     lambda f: str(f["chains"]), "integer_match"),
    ("residues", "How many standard amino acid residues are in PDB structure {pdb_id} ({name})? Answer with just the number.",
     lambda f: str(f["residues"]), "integer_match"),
    ("method", "What experimental method was used to determine the structure of {pdb_id} ({name})?",
     lambda f: f["method"], "exact_match"),
    ("organism", "From which organism was the protein in PDB entry {pdb_id} ({name}) derived?",
     lambda f: f["organism"], "exact_match"),
]

RESOLUTION_TEMPLATE = (
    "resolution",
    "What is the resolution (in Angstroms) of PDB structure {pdb_id} ({name})? Answer with just the number.",
    lambda f: str(f["resolution"]),
    "numeric_tolerance",
)


def generate_pdb_tasks(facts: list[dict]) -> list[dict]:
    """Generate questions from parsed PDB facts."""
    tasks = []
    task_id = 0

    for fact in facts:
        # Cycle through question types for each structure
        for q_type, q_template, answer_fn, verify_fn in QUESTION_TEMPLATES:
            task_id += 1
            tasks.append({
                "id": f"pro_struct_{task_id:04d}",
                "category": "StructureAnalysis",
                "question": q_template.format(**fact),
                "ideal": answer_fn(fact),
                "verification": "programmatic",
                "verification_fn": verify_fn,
                "source": "pro",
                "meta": {"template": "pdb", "pdb_id": fact["pdb_id"], "question_type": q_type},
            })

        # Resolution question (only if resolution is known)
        if fact.get("resolution"):
            task_id += 1
            q_type, q_template, answer_fn, verify_fn = RESOLUTION_TEMPLATE
            tasks.append({
                "id": f"pro_struct_{task_id:04d}",
                "category": "StructureAnalysis",
                "question": q_template.format(**fact),
                "ideal": answer_fn(fact),
                "verification": "programmatic",
                "verification_fn": verify_fn,
                "source": "pro",
                "meta": {"template": "pdb", "pdb_id": fact["pdb_id"], "question_type": q_type, "tolerance": 0.1},
            })

    return tasks


# --- Gel Image Tasks ---

def generate_gel_image(bands: list[tuple[float, float]], task_id: int) -> bytes:
    """Generate a synthetic gel electrophoresis image."""
    fig, ax = plt.subplots(figsize=(3, 6))
    ax.set_facecolor("#1a1a2e")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 250)
    ax.invert_yaxis()
    ax.set_ylabel("Molecular Weight (kDa)")
    ax.set_title("Lane 1", fontsize=10)
    ax.set_xticks([])

    for pos_kda, intensity in bands:
        y = pos_kda
        height = 3 + intensity * 2
        alpha = 0.3 + intensity * 0.5
        rect = plt.Rectangle((0.25, y - height / 2), 0.5, height,
                              color="white", alpha=min(alpha, 1.0))
        ax.add_patch(rect)

    ladder = [10, 15, 25, 35, 55, 70, 100, 130, 250]
    for mw in ladder:
        ax.axhline(y=mw, color="gray", linewidth=0.3, alpha=0.3)
        ax.text(-0.05, mw, f"{mw}", fontsize=6, ha="right", va="center", color="gray")

    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor="#0e0e1a")
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def generate_gel_tasks(count: int, start_id: int, output_dir: str) -> list[dict]:
    """Generate gel image tasks with known band patterns."""
    tasks = []
    for i in range(count):
        task_id = start_id + i
        np.random.seed(task_id + 1000)
        random.seed(task_id + 1000)

        n_bands = random.randint(1, 5)
        possible_sizes = [15, 25, 35, 42, 55, 70, 100, 130]
        chosen = random.sample(possible_sizes, min(n_bands, len(possible_sizes)))
        bands = [(size, random.uniform(0.3, 1.0)) for size in sorted(chosen)]

        img_data = generate_gel_image(bands, task_id)
        img_path = os.path.join(output_dir, f"gel_{task_id:04d}.png")
        with open(img_path, "wb") as f:
            f.write(img_data)

        band_sizes = [int(b[0]) for b in bands]
        largest = max(band_sizes)

        tasks.append({
            "id": f"pro_struct_{task_id:04d}",
            "category": "StructureAnalysis",
            "question": (
                "Examine this gel electrophoresis image. "
                "How many distinct bands are visible, and what is the approximate "
                "molecular weight (in kDa) of the largest band?"
            ),
            "ideal": f"{len(bands)} bands, largest at approximately {largest} kDa",
            "verification": "llm-judge",
            "source": "pro",
            "meta": {
                "template": "gel", "image_path": img_path,
                "band_sizes": band_sizes, "n_bands": len(bands),
            },
        })
    return tasks


def main():
    parser = argparse.ArgumentParser(description="Generate structure analysis tasks")
    parser.add_argument("--output-dir", default="tasks/structures")
    parser.add_argument("--pdb-count", type=int, default=50, help="Number of PDB structures to fetch")
    parser.add_argument("--gel-count", type=int, default=60, help="Number of gel images to generate")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- PDB Tasks: download and parse real structures ---
    pdb_ids = PDB_IDS[:args.pdb_count]
    print(f"Downloading {len(pdb_ids)} PDB structures...")

    with tempfile.TemporaryDirectory() as tmp_dir:
        facts = []
        for pdb_id in pdb_ids:
            fact = download_and_parse_pdb(pdb_id, tmp_dir)
            if fact:
                facts.append(fact)
                print(f"  {pdb_id}: {fact['name']} — {fact['chains']} chains, {fact['residues']} residues")

    print(f"Successfully parsed {len(facts)} structures")
    pdb_tasks = generate_pdb_tasks(facts)

    # --- Gel Tasks ---
    print(f"Generating {args.gel_count} gel image tasks...")
    gel_tasks = generate_gel_tasks(args.gel_count, start_id=len(pdb_tasks) + 1, output_dir=args.output_dir)

    # --- Write all tasks ---
    all_tasks = pdb_tasks + gel_tasks
    for task in all_tasks:
        filepath = os.path.join(args.output_dir, f"{task['id']}.json")
        with open(filepath, "w") as f:
            json.dump(task, f, indent=2)

    print(f"\nGenerated {len(all_tasks)} tasks ({len(pdb_tasks)} PDB, {len(gel_tasks)} gel) in {args.output_dir}")


if __name__ == "__main__":
    main()
