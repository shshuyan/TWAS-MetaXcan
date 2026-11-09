# AD candidate recheck: local analysis plan, 2026-11-08

This is an optional downstream virtual experiment in TWAS-MetaXcan. It evaluates a
frozen upstream drug list; it does not search for new drugs or automatically approve
wet-lab work. No LLM is used. Inputs are local and hash-recorded. No commits/pushes
are part of this change.

## Questions fixed before the real-data run

1. Do the neural perturbation response genes for each frozen drug overrepresent
   genes implicated by the specified AD genetic study?
2. Do the drug perturbations oppose the expression signatures of both the discovery
   and validation GEO cohorts under the same scoring/calibration method?
3. Where the frozen drug exists in the separate legacy CMap matrix, does that
   reference support the same direction? Missing legacy coverage is not failure.

## Frozen real-data inputs

- Twelve candidates from the existing AD agent run `ad-study-01`, discovery
  GSE132903. Candidate selection is not recomputed or optimized in this experiment.
- GSE132903 (97 AD/98 controls, middle temporal gyrus; group+sex) discovery DE.
- GSE118553 temporal donor-corrected DE (45 AD/25 controls; group+sex). Same-donor
  repeated arrays were previously averaged on the log2 scale. Cross-study donor
  independence and full age/postmortem/cell-composition adjustment remain unresolved.
- Bellenguez 2022 / PMID35379992 / GCST90027158: 89 Catalog associations at
  P<=5e-8, mapped to 112 unique labels in the existing preserved source export.
  The local run uses *locus-mapped* genes, not full TWAS or causal mapping. An
  explicitly selected alternative mode accepts complete S-MultiXcan or a single
  preselected tissue table, using tested genes as background and Bonferroni hits.
- Existing neural LINCS top250 up/down JSONL, exact perturbagen ID matching;
  NEU primary and NPC supplementary. GSE92742 gene annotation defines an explicit
  assayed/inferred gene universe; tail genes absent from it are excluded from the
  test and counted in the genetic outputs. This cross-release annotation is a
  limitation, not a claim of complete modern LINCS coverage.
- Legacy CMap build02 full ranks (lower means up), with the previously verified
  rank orientation and Entrez-symbol mapping. Exact names or explicit manually
  listed aliases only. No fuzzy identity or structure equivalence is claimed.

## Genetic enrichment

For each candidate, take the union of its NEU tail genes across available dose/time
conditions. For legacy CMap take the union of top/bottom250 genes across available
reference cell signatures. Test the resulting two sets separately against genetic
hits with one-sided Fisher. Background is the respective assay universe, intersected
with all tested TWAS genes in TWAS mode. Correct all covered candidate x reference
tests together by BH, alpha .05. Preserve overlap genes, 2x2 counts, expected overlap,
fold enrichment, odds ratio, raw P and adjusted q. Direct target overlap is a
separate annotation, not substituted for perturbation response enrichment.

This is exploratory gene-set overrepresentation: genes at the same GWAS locus are
not independent, mapping favors some genes, and perturbation ascertainment varies.
It is NOT an LD-aware SNP/heredity enrichment and does not identify causal drug
targets. It supplies no genetic risk-allele-to-expression direction.

## Reversal

Preserve the full DE tables. Per cohort/reference use FDR<.05 genes with positive and
negative logFC, ranked by P, up to100 each; require at least10 in each direction.
Background is genes both assayed and measured in the cohort. Do not relax thresholds.

LINCS score per signature is (same-direction overlap - opposite-direction overlap)
divided by query size. Aggregate by median per candidate/cell, preserving every
condition. CMap uses signed unweighted two-set KS, zero when enrichment signs agree,
then median across available reference cells. Neither is official CLUE tau.

Per cohort run10,000 random queries without replacement with fixed up/down sizes,
seed20261108. Empirical P=(count(null<=observed)+1)/10001. LINCS BH covers the frozen
candidate/cell groups; legacy CMap BH covers all1,309 matrix compounds. Negative
score and q<.05 is a supported check. This is a gene-query null, not a patient-label
test or clinical probability; coexpression is not preserved.

NEU support in every declared cohort plus NEU genetic enrichment defines two
co-occurring exploratory signals. Legacy CMap evidence is reported alongside and
never silently filled from LINCS. No combined P or efficacy score is produced.
Uncovered, insufficient query, tested-but-not-supported and supported are distinct.

## Dependence, reproducibility and limits

The same LINCS reference was used upstream: this is a consistency check, NOT new
independent drug evidence. Multiple conditions are not independent biological
replications. All candidate hits/negatives/missing records remain in output.
The real run does not execute the expensive upstream TWAS or retrieve new binding,
brain-exposure or clinical evidence. Tests using synthetic data validate statistics,
signs, failure states and integration, not biological success.

Each run uses a fresh output directory, writes hashes and parameters, and only
writes summary.json after success. The pipeline hook reuses a completion marker
only if config, input hashes and code hashes are unchanged. A changed run requires
a new output path; old outputs are never deleted automatically.
