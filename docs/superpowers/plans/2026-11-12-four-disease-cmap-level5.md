# Four-disease TWAS × CMap Level 5 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and run a deterministic four-disease TWAS × CMap Level 5 minimum closed loop with formal statistics, auditable lineage and publication-quality figures.

**Architecture:** Add bounded-memory GCTX ingestion and a continuous weighted-cosine screening engine alongside the legacy screen, then orchestrate four independent disease runs and render a fixed reporting suite. A synthetic mini-GCTX demo exercises the same production APIs while real public inputs remain external, checksum-registered files.

**Tech Stack:** Python 3.9, NumPy, pandas, SciPy, h5py, PyYAML, matplotlib, unittest.

**Spec:** `docs/superpowers/specs/2026-11-12-four-disease-cmap-level5-design.md`

## Global Constraints

- Four disease candidate lists remain independent; no pooled ranking.
- Primary score is weighted cosine on clipped tissue-specific S-PrediXcan Z and Level 5 MODZ; negative means reversal.
- Primary-family BH is computed within disease across compound × primary-cell tests only.
- Secondary contexts cannot promote a primary failure.
- Formal defaults are `B=10000`, `alpha=0.05`, `|Z|<=5`; demo reductions are labelled and recorded.
- Synthetic outputs are visibly labelled and never described as biological findings.
- Real large inputs and generated runs stay out of Git.

## Review Focus

- GCTX matrix orientation or byte-string metadata differs from the fixture: reader must detect or fail with a precise error.
- TWAS/CMap gene intersection is empty, duplicated or too small: retain an auditable insufficient-data result without threshold relaxation.
- A secondary-only hit is strong: it must never appear in the candidate table.
- Permutation resolution cannot support the requested FDR: report the bound rather than manufacture significance.
- Re-running with the same seed/config must preserve statistical tables and figure data deterministically.

---

### Task 1: Frozen dataset registry and reproducible downloader

**Files:**
- Create: `config/datasets/four_disease_public.yaml`
- Create: `twaslib/dataset_registry.py`
- Modify: `twas.py`
- Test: `tests/test_dataset_registry.py`

**Interfaces:**
- Produces: `load_registry(path) -> dict`, `download_snapshot(registry, destination, selected_ids=None) -> Path` and CLI `twas.py datasets`.
- Consumes: public URLs, optional expected bytes/checksum, freeze metadata.

- [ ] Write tests for schema validation, resumable local-URL download, SHA-256 enforcement and manifest creation; run and observe the expected missing-module/command failures.
- [ ] Implement the minimal registry/downloader and CLI; never overwrite a mismatched completed file.
- [ ] Run task tests and the full suite in the declared MetaXcan environment.
- [ ] Commit `feat: add frozen public dataset registry`.

### Task 2: Level 5 GCTX reader and continuous reversal statistics

**Files:**
- Create: `twaslib/cmap_level5.py`
- Create: `twaslib/continuous_screen.py`
- Modify: `twaslib/screen_stats.py`
- Test: `tests/test_cmap_level5.py`
- Test: `tests/test_continuous_screen.py`

**Interfaces:**
- Produces: `Level5Reference`, `weighted_cosine_scores`, `matched_signature_calibration`, sensitivity `permutation_calibration`, and `screen_disease(config, reference, output) -> Path`.
- Consumes: GCTX, signature/gene/compound metadata, complete tissue-specific S-PrediXcan table.

- [ ] Write failing tests for both GCTX orientations, metadata filtering, hand-computed weighted cosine, median aggregation, empirical P, primary-only global BH, secondary non-promotion and insufficient overlap.
- [ ] Implement bounded-block GCTX reads and continuous scoring with deterministic permutations.
- [ ] Run task tests and the full suite.
- [ ] Commit `feat: screen continuous TWAS signatures against CMap Level 5`.

### Task 3: Four-disease independent orchestration and audit lineage

**Files:**
- Create: `twaslib/four_disease_screen.py`
- Create: `config/drug_screens/four_disease_level5.yaml`
- Modify: `twas.py`
- Test: `tests/test_four_disease_screen.py`

**Interfaces:**
- Produces: `run_four_disease_screen(config_path, output=None) -> Path` and CLI `twas.py four-disease-screen`.
- Consumes: four disease-specific TWAS specs and one Level 5 reference spec.

- [ ] Write failing tests proving four isolated correction families, exact disease set, per-disease candidate tables, no pooled candidate rank, hashes and rerun protection.
- [ ] Implement orchestration, resolved configuration, input/output manifests and Chinese-first reports.
- [ ] Run task tests and the full suite.
- [ ] Commit `feat: orchestrate four independent disease screens`.

### Task 4: Publication-quality figure suite

**Files:**
- Create: `twaslib/screen_figures.py`
- Test: `tests/test_screen_figures.py`

**Interfaces:**
- Produces: `render_figure_suite(batch_output, formats=("png", "svg", "pdf")) -> list[Path]`.
- Consumes: tidy outputs from Task 3 only, so figures cannot silently recalculate primary statistics.

- [ ] Write failing tests for all eight figure identities, all three formats, demo watermark/caption and missing-data panels.
- [ ] Implement a shared journal style, colour-blind-safe palette, explicit panel labels and deterministic plot data exports.
- [ ] Render the suite and visually inspect PNGs at publication dimensions.
- [ ] Run task tests and the full suite.
- [ ] Commit `feat: add publication figure suite`.

### Task 5: Deterministic minimum closed loop, documentation and verification

**Files:**
- Create: `scripts/run_minimal_four_disease_demo.py`
- Create: `tests/test_minimal_demo.py`
- Modify: `envs_metaxcan.yaml`
- Modify: `.gitignore`
- Modify: `README.md`
- Modify: `docs/twas_drug_screen.md`

**Interfaces:**
- Produces: one command that builds mini GCTX/TWAS inputs outside tracked paths, runs all four diseases and renders the full output contract.
- Consumes: Tasks 1–4 production APIs.

- [ ] Write the failing end-to-end test first and verify the missing command/output failure.
- [ ] Implement the fixture builder/demo runner, dependency declaration and user documentation.
- [ ] Run the demo in a fresh temporary output; inspect candidate tables, manifests, reports and PNG renders.
- [ ] Run the full test suite, `git diff --check`, CLI help and a clean-tree/generated-data audit.
- [ ] Commit `test: add four-disease minimum closed loop`.

### Task 6: Whole-branch review and publication

**Files:** all changed files.

**Interfaces:**
- Consumes: completed tasks and their verification evidence.
- Produces: reviewed branch `codex/four-disease-cmap-level5` pushed to `origin`.

- [ ] Build the review package and perform the required fresh whole-branch review; classify findings before edits.
- [ ] Fix Critical/Important findings with RED→GREEN tests and rerun the full suite.
- [ ] Re-run the demo and final verification commands from a clean output location.
- [ ] Push the verified branch to `origin` and report the commit SHA, branch URL, demo output and explicit scientific limitations.
