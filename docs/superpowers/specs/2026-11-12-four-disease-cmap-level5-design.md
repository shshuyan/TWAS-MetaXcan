# Four-disease TWAS × CMap Level 5 design

## Intent

Extend the existing MetaXcan pipeline into a reproducible compound-repurposing workflow for four separately analysed indications: Alzheimer’s disease (AD), lichen simplex chronicus / neurodermatitis (LSC), atopic dermatitis (ADerm), and idiopathic pulmonary fibrosis (IPF). Each disease must produce its own ranked compound list. The workflow is for research prioritisation, not clinical recommendation.

The immediate deliverable is a deterministic minimum closed loop that proves the software path end to end with a small GCTX fixture and synthetic TWAS inputs. Demo outputs and figures must be visibly labelled `DEMO / synthetic` and cannot be interpreted as biological findings. The same interfaces must accept the full GSE92742 Level 5 matrix and real TWAS outputs without changing the analysis code.

## Goals

- Register public GWAS and CMap inputs with URLs, freeze dates, checksums, ancestry/build metadata, sample sizes, phenotype definitions and overlap notes.
- Read the official GSE92742 Level 5 GCTX layout in bounded memory and join `sig_info`, `pert_info`, and `gene_info` metadata.
- Use a preregistered continuous disease–compound reversal score as the single primary analysis.
- Aggregate signatures within compound × cell while retaining dose/time/signature evidence.
- Calculate empirical permutation P values and disease-wide BH FDR across all compound × primary-cell tests.
- Keep secondary cells and thresholded Up/Down methods as sensitivity analyses only.
- Produce four independent candidate tables, complete audit manifests, a Chinese-first report, and publication-quality SVG/PDF/PNG figures.
- Provide a deterministic four-disease demo command and automated regression tests.

## Non-goals

- The demo will not download or commit the 21.3 GB Level 5 matrix.
- Synthetic candidates will not be presented as real drugs or wet-lab-ready results.
- This change will not claim TWAS causality, patient differential expression, clinical efficacy, or independence between correlated signatures.
- It will not merge the four disease signatures or create one pooled candidate ranking.
- It will not silently relax thresholds when LSC or another disease has insufficient power.

## Constraints and frozen scientific decisions

- Primary score: weighted cosine similarity between clipped tissue-specific S-PrediXcan Z and CMap Level 5 MODZ. Negative values indicate reversal.
- TWAS clipping: symmetric absolute cap from config, preregistered at `|Z| = 5` for the main analysis.
- Gene weights: `min(|TWAS Z|, 5) / 5`; all finite genes shared by TWAS and CMap are used. S-MultiXcan remains annotation and does not determine direction.
- Compound aggregation: median score over all eligible signatures for each compound × cell; no best-dose selection.
- Null: deterministic signature-label permutations matched hierarchically by cell/dose/time, then cell/time, then cell when an exact stratum has no external signature. The observed signature count per compound × cell is preserved and the current compound is excluded whenever another matched signature exists. This makes the formal 473,647-signature screen bounded-memory; gene-label permutation remains a small-matrix sensitivity function. Empirical lower-tail P is `(1 + count(null <= observed)) / (B + 1)`.
- Multiplicity: one BH family per disease across all tested compound × primary-cell hypotheses. Secondary cells receive descriptive/sensitivity statistics and cannot promote a primary failure.
- Formal run default: `B = 10000`, `alpha = 0.05`; `0.05 <= q < 0.10` is exploratory. The demo may use fewer permutations and must record that reduced resolution.
- Candidate: negative primary score and global primary-family `q < 0.05`. No candidates is a valid result.
- Four primary contexts: AD uses neuronal cells (`NEU`, with `NPC` as a prespecified additional primary context if available); LSC and ADerm use `MCH58`; IPF uses `HCC515` as an imperfect lung proxy. Proxy limitations are explicit in every report. Other configured cells are sensitivity contexts.
- Every run records resolved config, input/output SHA-256, software/runtime versions, random seed, thresholds, tested family size, skipped rows and reasons.
- Figures use a colour-blind-safe palette, vector text, explicit units/denominators, panel labels, and no significance stars without the exact test and correction in the caption.

## Public data registry

The registry supports, at minimum:

- AD primary: Bellenguez 2022 / GWAS Catalog `GCST90027158`; replication: FinnGen R12 `G6_ALZHEIMER`.
- LSC primary: FinnGen R12 `L12_LICHENSIMPLEX`; `L12_LICHSIMPANDPRURIGO` is phenotype sensitivity, not an automatic replacement.
- ADerm primary: GWAS Catalog `GCST90244787` European-ancestry analysis; FinnGen R12 `L12_ATOPIC` is a replication/sensitivity dataset after overlap review.
- IPF primary: clinically curated European IPF meta-analysis frozen in the registry after accession/licence verification; FinnGen R12 `IPF` is replication.
- CMap: `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz` plus official `sig_info`, `pert_info`, `gene_info`, README and checksum files.

Downloads are resumable, checksum-verified and stored outside Git. A manifest can point at local scratch or a Google Drive archive, but active GCTX computation should use local scratch rather than a sync-mounted path. Downloading a resource never implies permission to redistribute it.

## Architecture

1. `dataset_registry` validates frozen metadata and creates/downloads a checksum-verifiable data snapshot.
2. `cmap_level5` opens GCTX metadata, resolves matrix orientation, selects compound signatures/cells, maps Entrez IDs to symbols, and streams signature blocks.
3. `continuous_screen` aligns a disease TWAS vector, streams condition scores, aggregates compound × cell, runs matched signature-label permutations, and applies primary-family FDR.
4. `four_disease_screen` runs diseases independently, writes one output tree per disease, and creates a cross-disease comparison containing no pooled ranking.
5. `screen_figures` creates a fixed figure suite from tidy output tables. Captions state whether inputs are demo or real.
6. A demo builder writes deterministic mini GCTX/metadata/TWAS fixtures and invokes the exact production APIs.

## Output contract

Each disease directory contains `candidate_compounds.tsv`, `all_primary_tests.tsv`, `secondary_contexts.tsv`, `signature_evidence.tsv.gz`, `aligned_twas_genes.tsv.gz`, `permutation_summary.tsv`, `resolved_config.yaml`, `input_manifest.json`, `output_manifest.json`, `summary.json`, `report_zh.md`, and figure files. The batch root contains four disease links/tables, a cross-disease descriptive matrix, a provenance graph, and an index report.

The minimum figure suite is:

1. preregistered pipeline/provenance schematic;
2. disease TWAS Z distribution and measured-gene coverage;
3. compound × primary-cell reversal heatmap;
4. observed-versus-permutation null panels;
5. score-versus-global-FDR prioritisation plot;
6. dose/time/signature consistency plot;
7. four-disease top-candidate comparison without pooled ranking;
8. sensitivity-analysis forest plot.

## Approaches considered

- **A — extend the existing config-driven pipeline (selected):** reuses CLI, result conventions and audit hooks while replacing the primary drug-screen statistic with Level 5 continuous scoring.
- **B — independent package:** cleaner isolation but duplicates configuration, provenance and pipeline integration.
- **C — notebook-led analysis:** fastest for exploratory plots but inadequate for regression testing and formal audit.

## Acceptance criteria

- A test-created GCTX using official group names is read correctly in both possible matrix orientations.
- Four demo diseases complete independently with fixed seed and unchanged reruns are byte-stable except explicitly excluded timestamps.
- Candidate/FDR calculations are validated against hand-computable fixtures, including no-candidate and insufficient-gene cases.
- Secondary-cell significance never promotes a primary-context failure.
- The full unit/integration suite passes in the declared environment.
- The demo command emits all contracted tables, manifests, Chinese reports and figure formats; rendered PNGs are visually inspected.
- Git contains code/config/docs/small fixtures only, never real large GWAS/CMap data or generated run directories.

## Builder self-check

No subagent review was used because this session is not authorised to delegate. The design was checked against the existing code. The principal implementation risks are GCTX orientation/metadata variance, memory use, multiple-testing family leakage, deterministic output, and accidental interpretation of synthetic figures; each is assigned an explicit acceptance test above.
