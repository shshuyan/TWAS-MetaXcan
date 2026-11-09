# Replacement implementation validation — 2026-11-09

## Change in workflow

The previous fixed-shortlist evidence-recheck commit (`e5dafff`) was reverted
with `be08604`. It remains recoverable in Git history; the new tree does not
expose its `recheck` command or `candidate_recheck` pipeline step.

The replacement adds `drug-screen` / `drug_screen` after TWAS summarization:
full signed single-tissue TWAS + joint associations → supplied perturbation
library → candidate compounds. Neither GEO nor an upstream drug shortlist is
required. Shared rank reading, calibration and provenance helpers were extracted
into standalone modules rather than importing the superseded workflow.

## Executed checks

- Fifteen tests passed via `python -m unittest discover -s tests -p 'test_*.py' -v`:
  twelve discovery/input/CLI/pipeline checks and three statistical regression
  checks. Synthetic fixtures distinguish reversing and mimicking drugs, opposing
  cells, median versus best-condition selection, single-tissue direction versus
  joint z_mean, conflict handling, global multiplicity and insufficient queries.
- The standalone CLI completed a synthetic smoke run over four LINCS compound/cell
  groups and three CMap compound/cell groups. Only the known reversing primary
  compound was shortlisted; the mimicking compound was not.
- An independent SciPy BH calculation agreed with the seven reported adjusted
  P-values. All twelve input/implementation hashes and thirteen output hashes
  matched the completed smoke output.
- The available real neural LINCS reference passed whole-file validation:
  9,220 signatures, 3,645 unique compound IDs, 6,269 compound/cell groups;
  4,096 NEU signatures and 5,124 NPC signatures. This is reference compatibility
  checking, not a completed real TWAS-to-drug screen.
- `git diff --check` passed. Local configs, reference matrices, credentials and
  synthetic run artifacts are excluded from the source submission.

## Biological limits

No complete real AD S-PrediXcan/S-MultiXcan result tables were available locally
for this test. The replacement was therefore verified with explicitly labelled
synthetic TWAS inputs; it does not supply a new AD candidate-drug conclusion.
The actual GWAS-to-TWAS run and reference/environment preparation remain required
before a real discovery run. The workflow does not implement LD-corrected
heritability enrichment, colocalization or causal drug-efficacy validation.
