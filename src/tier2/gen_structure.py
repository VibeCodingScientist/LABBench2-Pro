"""Generate PDB structure + gel image tasks.

Usage: python -m src.tier2.gen_structure --output-dir tasks/structures --count 150
"""

import argparse
import io
import json
import os
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# --- PDB Tasks ---

# Pre-defined PDB facts (avoids needing network access during generation)
PDB_FACTS = [
    {"pdb_id": "1MBO", "name": "Myoglobin", "chains": 1, "residues": 153, "resolution": 2.0,
     "organism": "Physeter catodon", "method": "X-ray diffraction",
     "secondary": {"helix": 0.75, "sheet": 0.0, "coil": 0.25}},
    {"pdb_id": "1HHO", "name": "Hemoglobin", "chains": 4, "residues": 574, "resolution": 2.1,
     "organism": "Homo sapiens", "method": "X-ray diffraction",
     "secondary": {"helix": 0.72, "sheet": 0.0, "coil": 0.28}},
    {"pdb_id": "1EMA", "name": "Green fluorescent protein", "chains": 1, "residues": 238, "resolution": 1.9,
     "organism": "Aequorea victoria", "method": "X-ray diffraction",
     "secondary": {"helix": 0.05, "sheet": 0.45, "coil": 0.50}},
    {"pdb_id": "1BNA", "name": "B-DNA dodecamer", "chains": 2, "residues": 24, "resolution": 1.9,
     "organism": "Synthetic", "method": "X-ray diffraction",
     "secondary": {"helix": 0.0, "sheet": 0.0, "coil": 1.0}},
    {"pdb_id": "3J3Y", "name": "Ribosome 80S", "chains": 80, "residues": 13000, "resolution": 3.0,
     "organism": "Saccharomyces cerevisiae", "method": "Cryo-EM",
     "secondary": {"helix": 0.15, "sheet": 0.20, "coil": 0.65}},
    {"pdb_id": "1TIM", "name": "Triosephosphate isomerase", "chains": 2, "residues": 496, "resolution": 2.5,
     "organism": "Gallus gallus", "method": "X-ray diffraction",
     "secondary": {"helix": 0.40, "sheet": 0.25, "coil": 0.35}},
    {"pdb_id": "4HHB", "name": "Deoxyhemoglobin", "chains": 4, "residues": 574, "resolution": 1.74,
     "organism": "Homo sapiens", "method": "X-ray diffraction",
     "secondary": {"helix": 0.70, "sheet": 0.0, "coil": 0.30}},
    {"pdb_id": "1GFL", "name": "GFP S65T mutant", "chains": 1, "residues": 238, "resolution": 1.9,
     "organism": "Aequorea victoria", "method": "X-ray diffraction",
     "secondary": {"helix": 0.06, "sheet": 0.44, "coil": 0.50}},
    {"pdb_id": "2PTC", "name": "Trypsin", "chains": 2, "residues": 223, "resolution": 1.55,
     "organism": "Bos taurus", "method": "X-ray diffraction",
     "secondary": {"helix": 0.10, "sheet": 0.35, "coil": 0.55}},
    {"pdb_id": "1CRN", "name": "Crambin", "chains": 1, "residues": 46, "resolution": 0.83,
     "organism": "Crambe hispanica", "method": "X-ray diffraction",
     "secondary": {"helix": 0.40, "sheet": 0.15, "coil": 0.45}},
]


STRUCTURE_QUESTIONS = [
    ("chains", "How many polypeptide chains does {name} ({pdb_id}) contain?", lambda f: str(f["chains"])),
    ("residues", "How many amino acid residues are in {name} ({pdb_id})?", lambda f: str(f["residues"])),
    ("method", "What experimental method was used to determine the structure of {name} ({pdb_id})?", lambda f: f["method"]),
    ("organism", "From which organism was {name} ({pdb_id}) derived?", lambda f: f["organism"]),
    ("resolution", "What is the resolution (in Angstroms) of the {name} ({pdb_id}) structure?", lambda f: str(f["resolution"])),
    ("helix_frac", "What fraction of {name} ({pdb_id}) is alpha-helix? Report to 2 decimal places.", lambda f: f"{f['secondary']['helix']:.2f}"),
]


def generate_pdb_task(task_id: int) -> dict:
    fact = PDB_FACTS[task_id % len(PDB_FACTS)]
    q_type, q_template, answer_fn = STRUCTURE_QUESTIONS[task_id % len(STRUCTURE_QUESTIONS)]
    question = q_template.format(**fact)
    ideal = answer_fn(fact)

    verification_fn = "exact_match"
    if q_type in ("resolution", "helix_frac"):
        verification_fn = "numeric_tolerance"

    return {
        "id": f"pro_struct_{task_id:04d}",
        "category": "StructureAnalysis",
        "question": question,
        "ideal": ideal,
        "verification": "programmatic",
        "verification_fn": verification_fn,
        "source": "pro",
        "meta": {"template": "pdb", "pdb_id": fact["pdb_id"], "question_type": q_type,
                 "tolerance": 0.1 if verification_fn == "numeric_tolerance" else None},
    }


# --- Gel Image Tasks ---

def generate_gel_image(bands: list[tuple[float, float]], task_id: int) -> bytes:
    """Generate a synthetic gel electrophoresis image.

    bands: list of (position_kda, intensity) tuples
    """
    fig, ax = plt.subplots(figsize=(3, 6))
    ax.set_facecolor("#1a1a2e")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 250)
    ax.invert_yaxis()
    ax.set_ylabel("Molecular Weight (kDa)")
    ax.set_title(f"Lane 1", fontsize=10)
    ax.set_xticks([])

    # Draw bands
    for pos_kda, intensity in bands:
        y = pos_kda
        width = 0.5
        height = 3 + intensity * 2
        alpha = 0.3 + intensity * 0.5
        rect = plt.Rectangle((0.25, y - height/2), width, height,
                              color="white", alpha=min(alpha, 1.0))
        ax.add_patch(rect)

    # Ladder markers
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


def generate_gel_task(task_id: int) -> dict:
    """Generate a gel image task with known band pattern."""
    np.random.seed(task_id + 1000)

    n_bands = random.randint(1, 5)
    bands = []
    for _ in range(n_bands):
        pos = random.choice([15, 25, 35, 42, 55, 70, 100, 130])
        intensity = random.uniform(0.3, 1.0)
        bands.append((pos, intensity))
    bands.sort(key=lambda b: b[0])

    # Generate image
    img_data = generate_gel_image(bands, task_id)
    img_path = f"tasks/structures/gel_{task_id:04d}.png"
    os.makedirs("tasks/structures", exist_ok=True)
    with open(img_path, "wb") as f:
        f.write(img_data)

    band_sizes = [int(b[0]) for b in bands]
    largest = max(band_sizes)

    return {
        "id": f"pro_struct_{task_id:04d}",
        "category": "StructureAnalysis",
        "question": f"Examine this gel electrophoresis image. "
                     f"How many distinct bands are visible, and what is the approximate "
                     f"molecular weight (in kDa) of the largest band?",
        "ideal": f"{len(bands)} bands, largest at approximately {largest} kDa",
        "verification": "llm-judge",
        "source": "pro",
        "meta": {"template": "gel", "image_path": img_path, "band_sizes": band_sizes,
                 "n_bands": len(bands)},
    }


def main():
    parser = argparse.ArgumentParser(description="Generate structure analysis tasks")
    parser.add_argument("--output-dir", default="tasks/structures")
    parser.add_argument("--count", type=int, default=150)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    tasks = []

    # Split: 60% PDB, 40% gel
    n_pdb = int(args.count * 0.6)
    n_gel = args.count - n_pdb

    for i in range(n_pdb):
        task = generate_pdb_task(i + 1)
        tasks.append(task)

    for i in range(n_gel):
        task = generate_gel_task(n_pdb + i + 1)
        tasks.append(task)

    for task in tasks:
        filepath = os.path.join(args.output_dir, f"{task['id']}.json")
        with open(filepath, "w") as f:
            json.dump(task, f, indent=2)

    print(f"Generated {len(tasks)} tasks ({n_pdb} PDB, {n_gel} gel) in {args.output_dir}")


if __name__ == "__main__":
    main()
