#!/usr/bin/env python
"""Create a deterministic mini-GCTX and run the complete four-disease screen."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from twaslib.four_disease_screen import run_four_disease_screen  # noqa: E402
from twaslib.screen_figures import render_figure_suite  # noqa: E402


DISEASES = ["AD", "LSC", "atopic_dermatitis", "IPF"]
CONTEXTS = {
    "AD": {"tissue": "Brain_Cortex", "primary": ["NEU", "NPC"], "secondary": ["A549"]},
    "LSC": {"tissue": "Skin_Not_Sun_Exposed_Suprapubic", "primary": ["MCH58"], "secondary": ["A549"]},
    "atopic_dermatitis": {"tissue": "Skin_Not_Sun_Exposed_Suprapubic", "primary": ["MCH58"], "secondary": ["A549"]},
    "IPF": {"tissue": "Lung", "primary": ["HCC515"], "secondary": ["A549"]},
}


def build_demo(root):
    """Write deterministic synthetic inputs that exercise production interfaces."""
    inputs = root / "inputs"
    inputs.mkdir()
    rng = np.random.default_rng(20261112)
    gene_count = 60
    gene_ids = [str(1000 + i) for i in range(gene_count)]
    symbols = [f"GENE{i:03}" for i in range(gene_count)]
    orthogonal, _ = np.linalg.qr(rng.normal(size=(gene_count, len(DISEASES))))
    disease_z = {disease: orthogonal[:, i] / np.max(np.abs(orthogonal[:, i])) * 5
                 for i, disease in enumerate(DISEASES)}
    for disease in DISEASES:
        context = CONTEXTS[disease]
        pd.DataFrame({
            "gene_name": symbols, "tissue": context["tissue"],
            "zscore": disease_z[disease], "pvalue": np.linspace(.0001, .8, gene_count),
        }).to_csv(inputs / f"{disease}_twas.tsv", sep="\t", index=False)

    cells = ["NEU", "NPC", "MCH58", "HCC515", "A549"]
    compounds = []
    for disease in DISEASES:
        compounds.extend([(f"REV_{disease}_{rep}", disease, "reverse") for rep in [1, 2]])
        compounds.append((f"MIM_{disease}", disease, "mimic"))
    compounds.extend([(f"NULL_{i}", None, "null") for i in range(1, 5)])
    signatures, vectors = [], []
    for compound_index, (compound, target, mode) in enumerate(compounds):
        for cell_index, cell in enumerate(cells):
            for condition, (dose, time, noise_scale) in enumerate([(1, 6, .08), (10, 24, .14)], 1):
                sid = f"{compound}_{cell}_{condition}"
                if mode == "reverse": base = -disease_z[target]
                elif mode == "mimic": base = disease_z[target]
                else: base = rng.normal(0, 1, gene_count)
                cell_scale = 1 - .035 * cell_index
                noise = rng.normal(0, noise_scale, gene_count)
                vectors.append(base * cell_scale + noise)
                signatures.append({
                    "sig_id": sid, "pert_id": compound, "pert_iname": compound.lower(),
                    "pert_type": "trt_cp", "cell_id": cell, "pert_dose": dose,
                    "pert_dose_unit": "uM", "pert_time": time, "pert_time_unit": "h",
                })
    with h5py.File(inputs / "mini_level5.gctx", "w") as handle:
        handle["/0/DATA/0/matrix"] = np.asarray(vectors, dtype=np.float32).T
        handle["/0/META/ROW/id"] = np.asarray(gene_ids, dtype="S")
        handle["/0/META/COL/id"] = np.asarray([item["sig_id"] for item in signatures], dtype="S")
    pd.DataFrame(signatures).to_csv(inputs / "sig_info.tsv", sep="\t", index=False)
    pd.DataFrame({"pr_gene_id": gene_ids, "pr_gene_symbol": symbols,
                  "pr_is_lm": [1 if i < 40 else 0 for i in range(gene_count)]}).to_csv(
        inputs / "gene_info.tsv", sep="\t", index=False)
    pd.DataFrame({"pert_id": [item[0] for item in compounds],
                  "pert_iname": [item[0].lower() for item in compounds]}).to_csv(
        inputs / "pert_info.tsv", sep="\t", index=False)
    config = {
        "schema_version": 1, "input_kind": "synthetic_test", "output": "results",
        "analysis": {"z_clip": 5, "minimum_genes": 40, "permutations": 999, "seed": 20261112, "alpha": .05},
        "cmap": {"gctx": "inputs/mini_level5.gctx", "sig_info": "inputs/sig_info.tsv",
                 "gene_info": "inputs/gene_info.tsv", "pert_info": "inputs/pert_info.tsv", "cells": cells},
        "diseases": [{"disease": disease, "twas": f"inputs/{disease}_twas.tsv",
                      "tissue": CONTEXTS[disease]["tissue"], "primary_cells": CONTEXTS[disease]["primary"],
                      "secondary_cells": CONTEXTS[disease]["secondary"]} for disease in DISEASES],
    }
    config_path = root / "demo_config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True))
    return config_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="new demo directory; existing paths are refused")
    parser.add_argument("--formats", nargs="+", choices=["png", "svg", "pdf"], default=["png", "svg", "pdf"])
    args = parser.parse_args()
    root = Path(args.output).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=False)
    config_path = build_demo(root)
    batch = run_four_disease_screen(config_path)
    render_figure_suite(batch, formats=tuple(args.formats))
    print(f"DEMO / synthetic complete: {batch}")
    print("These outputs validate software behavior only; they are not biological candidate claims.")


if __name__ == "__main__":
    main()
