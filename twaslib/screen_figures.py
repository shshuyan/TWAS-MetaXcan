"""Publication-style figures rendered only from four-disease tidy outputs."""
from pathlib import Path
import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

from .screen_io import sha256


DISEASES = ["AD", "LSC", "atopic_dermatitis", "IPF"]
LABELS = {"AD": "Alzheimer’s disease", "LSC": "Lichen simplex chronicus",
          "atopic_dermatitis": "Atopic dermatitis", "IPF": "Idiopathic pulmonary fibrosis"}
COLORS = {"AD": "#0072B2", "LSC": "#009E73", "atopic_dermatitis": "#E69F00", "IPF": "#CC79A7"}
FIGURES = [
    "figure_01_pipeline", "figure_02_twas_signatures", "figure_03_reversal_heatmap",
    "figure_04_permutation_calibration", "figure_05_prioritization",
    "figure_06_condition_consistency", "figure_07_cross_disease_candidates",
    "figure_08_context_sensitivity",
]


def _style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
        "axes.labelsize": 9, "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "svg.fonttype": "none",
        "pdf.fonttype": 42, "axes.linewidth": .8,
    })


def _panel(ax, label):
    ax.text(-.11, 1.06, label, transform=ax.transAxes, weight="bold", fontsize=12, va="top")


def _empty(ax, disease=None):
    title = LABELS.get(disease, disease) if disease else "No analysable data"
    ax.text(.5, .55, title, ha="center", va="center", weight="bold", transform=ax.transAxes)
    ax.text(.5, .42, "No candidates", ha="center", va="center", color="#666666", transform=ax.transAxes)
    ax.set_axis_off()


def _load(batch):
    summary = json.loads((batch / "summary.json").read_text())
    frames = {}
    for disease in DISEASES:
        root = batch / disease
        frames[disease] = {
            "twas": pd.read_csv(root / "aligned_twas_genes.tsv.gz", sep="\t"),
            "primary": pd.read_csv(root / "all_primary_tests.tsv", sep="\t"),
            "secondary": pd.read_csv(root / "secondary_contexts.tsv", sep="\t"),
            "candidates": pd.read_csv(root / "candidate_compounds.tsv", sep="\t"),
            "evidence": pd.read_csv(root / "signature_evidence.tsv.gz", sep="\t"),
            "summary": json.loads((root / "summary.json").read_text()),
        }
    return summary, frames


def _pipeline(frames):
    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    steps = ["Frozen public\ninputs", "GWAS → TWAS\nZ-score", "Level 5 MODZ\nblock reader",
             "Weighted cosine\n+ permutations", "Disease-wide\nBH-FDR", "4 independent\ncandidate lists"]
    xs = np.linspace(.08, .92, len(steps))
    for index, (x, text) in enumerate(zip(xs, steps)):
        box = FancyBboxPatch((x - .072, .42), .144, .24, boxstyle="round,pad=.02,rounding_size=.025",
                             fc="#F3F6F8", ec="#3A506B", lw=1.2)
        ax.add_patch(box); ax.text(x, .54, text, ha="center", va="center", weight="bold", fontsize=8.5)
        if index < len(steps) - 1:
            ax.annotate("", xy=(xs[index + 1] - .08, .54), xytext=(x + .08, .54),
                        arrowprops={"arrowstyle": "-|>", "color": "#3A506B", "lw": 1.1})
    ax.text(.5, .18, "Primary cells define candidates  •  Secondary cells remain sensitivity evidence  •  Complete SHA-256 lineage",
            ha="center", color="#4D4D4D")
    data = pd.DataFrame({"step": np.arange(1, len(steps) + 1), "label": [s.replace("\n", " ") for s in steps]})
    return fig, data, "Preregistered analysis and provenance flow; four disease correction families remain separate."


def _twas(frames):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.1), sharex=True, sharey=True)
    rows = []
    for label, disease, ax in zip("ABCD", DISEASES, axes.ravel()):
        frame = frames[disease]["twas"]
        values = frame.zscore.to_numpy(float) if len(frame) else np.array([])
        if len(values):
            ax.hist(values, bins=min(30, max(5, len(values))), color=COLORS[disease], alpha=.85, edgecolor="white")
            ax.axvline(0, color="#333333", lw=.8)
            ax.set_title(LABELS[disease]); ax.set_xlabel("S-PrediXcan Z"); ax.set_ylabel("Genes")
            rows.extend({"disease": disease, "gene": gene, "zscore": z} for gene, z in zip(frame.gene_name, values))
        else: _empty(ax, disease)
        _panel(ax, label)
    fig.suptitle("Disease-specific genetic expression signatures", weight="bold")
    fig.tight_layout()
    return fig, pd.DataFrame(rows, columns=["disease", "gene", "zscore"]), "Tissue-specific TWAS Z distributions after CMap gene intersection; no cross-disease pooling."


def _heatmap(frames):
    rows = pd.concat([frames[d]["primary"] for d in DISEASES], ignore_index=True)
    rows["context"] = rows.disease + " · " + rows.cell
    table = rows.pivot_table(index="context", columns="drug", values="score", aggfunc="median")
    fig, ax = plt.subplots(figsize=(max(6.5, .55 * max(1, table.shape[1])), max(3.8, .42 * max(1, table.shape[0]))))
    if table.empty: _empty(ax)
    else:
        values = table.to_numpy(float)
        limit = max(.2, np.nanmax(np.abs(values)))
        image = ax.imshow(values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
        ax.set_xticks(range(len(table.columns)), table.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(table.index)), table.index)
        ax.set_xlabel("Compound"); ax.set_ylabel("Disease · primary cell")
        fig.colorbar(image, ax=ax, label="Weighted cosine (negative = reversal)", fraction=.025, pad=.02)
    ax.set_title("Primary-context reversal landscape", weight="bold"); _panel(ax, "A"); fig.tight_layout()
    return fig, rows, "Compound × primary-cell weighted cosine scores. Negative values indicate reversal."


def _permutation(frames):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.1), sharex=True, sharey=True)
    rows = []
    for label, disease, ax in zip("ABCD", DISEASES, axes.ravel()):
        frame = frames[disease]["primary"]
        if len(frame):
            y = -np.log10(frame.empirical_p.clip(lower=np.finfo(float).tiny))
            ax.scatter(frame.score, y, s=28, color=COLORS[disease], alpha=.8, edgecolor="white", linewidth=.4)
            floor = frames[disease]["summary"].get("minimum_empirical_p", 1)
            ax.axhline(-np.log10(floor), color="#777777", ls="--", lw=.8, label="Permutation resolution")
            ax.axvline(0, color="#333333", lw=.8); ax.set_title(LABELS[disease])
            ax.set_xlabel("Observed median score"); ax.set_ylabel("−log₁₀ empirical P")
            rows.extend(frame.assign(neg_log10_empirical_p=y).to_dict("records"))
        else: _empty(ax, disease)
        _panel(ax, label)
    fig.suptitle("Observed reversal and permutation calibration", weight="bold"); fig.tight_layout()
    return fig, pd.DataFrame(rows), "Observed scores versus matched signature-label permutation P values; dashed lines show finite permutation resolution."


def _prioritization(frames):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.1), sharex=True, sharey=True)
    rows = []
    for label, disease, ax in zip("ABCD", DISEASES, axes.ravel()):
        frame = frames[disease]["primary"].copy()
        if len(frame):
            frame["neg_log10_q"] = -np.log10(frame.q_bh_primary_family.clip(lower=np.finfo(float).tiny))
            colors = np.where(frame.candidate, "#D55E00", "#A6A6A6")
            ax.scatter(frame.score, frame.neg_log10_q, c=colors, s=34, alpha=.9, edgecolor="white", linewidth=.4)
            ax.axvline(0, color="#333333", lw=.8); ax.axhline(-np.log10(.05), color="#333333", ls="--", lw=.8)
            ax.set_title(LABELS[disease]); ax.set_xlabel("Weighted cosine"); ax.set_ylabel("−log₁₀ global q")
            rows.extend(frame.to_dict("records"))
        else: _empty(ax, disease)
        _panel(ax, label)
    fig.suptitle("Primary-family compound prioritisation", weight="bold"); fig.tight_layout()
    return fig, pd.DataFrame(rows), "Reversal score and within-disease global BH-FDR. Orange points pass the compound × primary-cell criterion."


def _consistency(frames):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.1), sharey=True)
    rows = []
    rng = np.random.default_rng(812)
    for label, disease, ax in zip("ABCD", DISEASES, axes.ravel()):
        evidence = frames[disease]["evidence"].copy()
        primary_cells = set(frames[disease]["summary"]["primary_cells"])
        evidence = evidence[evidence.cell_id.isin(primary_cells)]
        if len(evidence):
            ranked = evidence.groupby("pert_iname").score.median().abs().sort_values(ascending=False).head(10).index
            evidence = evidence[evidence.pert_iname.isin(ranked)].copy()
            drugs = sorted(evidence.pert_iname.unique())
            for x, drug in enumerate(drugs):
                values = evidence.loc[evidence.pert_iname == drug, "score"].to_numpy(float)
                jitter = rng.normal(0, .035, len(values))
                ax.scatter(np.full(len(values), x) + jitter, values, s=24, color=COLORS[disease], alpha=.75)
                ax.plot([x - .18, x + .18], [np.median(values)] * 2, color="#222222", lw=1.5)
            ax.axhline(0, color="#555555", lw=.8); ax.set_xticks(range(len(drugs)), drugs, rotation=35, ha="right")
            ax.set_title(LABELS[disease]); ax.set_ylabel("Signature-level weighted cosine")
            rows.extend(evidence.to_dict("records"))
        else: _empty(ax, disease)
        _panel(ax, label)
    fig.suptitle("Dose/time signature consistency in primary cells", weight="bold"); fig.tight_layout()
    return fig, pd.DataFrame(rows), "Individual Level 5 signatures with median bars; favourable single conditions are not selected post hoc."


def _candidates(frames):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.1), sharex=True)
    rows = []
    for label, disease, ax in zip("ABCD", DISEASES, axes.ravel()):
        frame = frames[disease]["candidates"].copy()
        if len(frame):
            frame = frame.sort_values("median_primary_score").head(10)
            ax.barh(frame.drug, frame.median_primary_score, color=COLORS[disease], alpha=.85)
            ax.axvline(0, color="#333333", lw=.8); ax.set_title(LABELS[disease])
            rows.extend(frame.assign(disease=disease).to_dict("records"))
        else: _empty(ax, disease)
        _panel(ax, label)
    fig.suptitle("Independent disease candidate lists (no pooled rank)", weight="bold")
    fig.supxlabel("Median primary-cell score", y=.015)
    fig.tight_layout(rect=(0, .035, 1, 1))
    return fig, pd.DataFrame(rows), "Top compounds are displayed within each disease only; absence of a bar is reported as No candidates."


def _sensitivity(frames):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.1), sharex=True)
    rows = []
    for label, disease, ax in zip("ABCD", DISEASES, axes.ravel()):
        primary = frames[disease]["primary"].groupby(["compound_id", "drug"], as_index=False).score.median().rename(columns={"score": "primary_score"})
        secondary = frames[disease]["secondary"].groupby(["compound_id", "drug"], as_index=False).score.median().rename(columns={"score": "secondary_score"})
        merged = primary.merge(secondary, on=["compound_id", "drug"], how="outer")
        if len(merged):
            merged["display_priority"] = merged[["primary_score", "secondary_score"]].abs().max(axis=1)
            merged = merged.sort_values(["display_priority", "compound_id"], ascending=[False, True]).head(10).reset_index(drop=True)
            y = np.arange(len(merged))
            for index, row in merged.iterrows():
                if pd.notna(row.primary_score) and pd.notna(row.secondary_score):
                    ax.plot([row.primary_score, row.secondary_score], [index, index], color="#B0B0B0", lw=1)
            ax.scatter(merged.primary_score, y, color=COLORS[disease], label="Primary", s=34, zorder=3)
            ax.scatter(merged.secondary_score, y, facecolor="white", edgecolor="#444444", label="Secondary", s=34, zorder=3)
            ax.axvline(0, color="#333333", lw=.8); ax.set_yticks(y, merged.drug); ax.set_title(LABELS[disease]); ax.set_xlabel("Median score")
            if disease == "AD": ax.legend(frameon=False, fontsize=8)
            rows.extend(merged.assign(disease=disease).to_dict("records"))
        else: _empty(ax, disease)
        _panel(ax, label)
    fig.suptitle("Primary versus secondary cell-context sensitivity", weight="bold"); fig.tight_layout()
    return fig, pd.DataFrame(rows), "Primary and secondary cell medians are displayed separately; secondary evidence cannot promote a candidate."


BUILDERS = [_pipeline, _twas, _heatmap, _permutation, _prioritization, _consistency, _candidates, _sensitivity]


def render_figure_suite(batch_output, formats=("png", "svg", "pdf")):
    """Render eight figures and their exact plotting tables from tidy outputs."""
    _style()
    batch = Path(batch_output).resolve()
    summary, frames = _load(batch)
    demo = summary.get("input_kind") == "synthetic_test"
    invalid = set(formats) - {"png", "svg", "pdf"}
    if invalid or not formats:
        raise ValueError("formats must be a nonempty subset of png, svg, pdf")
    figure_dir = batch / "figures"
    figure_dir.mkdir(exist_ok=True)
    runtime_path = batch / "runtime_versions.json"
    if runtime_path.is_file():
        runtime = json.loads(runtime_path.read_text())
        runtime["matplotlib"] = matplotlib.__version__
        runtime_path.write_text(json.dumps(runtime, indent=2, sort_keys=True) + "\n")
    input_manifest_path = batch / "input_manifest.json"
    if input_manifest_path.is_file():
        inputs = json.loads(input_manifest_path.read_text())
        role = "implementation:screen_figures.py"
        inputs = [item for item in inputs if item.get("role") != role]
        implementation = Path(__file__).resolve()
        inputs.append({"role": role, "file": str(implementation), "bytes": implementation.stat().st_size,
                       "sha256": sha256(implementation)})
        input_manifest_path.write_text(json.dumps(inputs, indent=2, sort_keys=True) + "\n")
    generated, captions = [], []
    for number, (stem, builder) in enumerate(zip(FIGURES, BUILDERS), 1):
        fig, data, caption = builder(frames)
        if demo:
            fig.text(.995, .006, "DEMO / synthetic", ha="right", va="bottom", color="#B2182B", fontsize=8, weight="bold")
        data.to_csv(figure_dir / f"{stem}.data.tsv", sep="\t", index=False)
        for extension in formats:
            path = figure_dir / f"{stem}.{extension}"
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", message="'mode' parameter is deprecated.*", category=DeprecationWarning,
                )
                fig.savefig(path, bbox_inches="tight", facecolor="white")
            generated.append(path)
        plt.close(fig)
        disclosure = " DEMO / synthetic; software validation only." if demo else ""
        captions.append(f"## Figure {number}\n\n{caption} Statistics are read from the audited result tables; BH-FDR is disease-specific.{disclosure}\n")
    (figure_dir / "figure_captions.md").write_text("\n".join(captions))
    # Figures are part of the final audit surface, so refresh the batch manifest.
    manifest = []
    for path in sorted(batch.rglob("*")):
        if path.is_file() and path != batch / "output_manifest.json":
            manifest.append({"role": "output", "file": str(path.relative_to(batch)),
                             "bytes": path.stat().st_size, "sha256": sha256(path)})
    (batch / "output_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return generated
