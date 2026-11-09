# Local AD integration test — 2026-11-08

Base revision: `2958be2214e0127bf7499c80c37887660fe26a3a`.
This note records the local test performed before submission. The existing
AD-agent checkout was not modified. Large inputs are referenced locally, not
duplicated into the Git tree or included with the code submission.

## Executed

1. Eleven unit/integration tests passed, covering Fisher background counts, BH,
   KS sign conventions, same-sign zero scores, calibration with more than two
   groups, explicit genetic modes, insufficient queries, missing coverage,
   duplicate TWAS genes, identical-DE pseudo-replication, CLI, pipeline hook and
   changed-input rejection. Biological assertions use synthetic fixtures only.
2. Actual CLI: `twas.py recheck config/rechecks/local_ad.yaml` using the existing
   AD agent Python environment. Twelve frozen candidates, two real GEO DE inputs,
   real legacy CMap full ranks and real neural LINCS signatures. Each cohort used
   10,000 random queries for each reference calibration.
3. Actual pipeline integration:
   `twas.py run work/recheck-inputs/local_AD_hook.yaml --steps candidate_recheck`.
   The hook correctly reused the identical completed real output after checking
   all source hashes and resolved settings. It did not execute upstream TWAS.
4. `tests/verify_recheck_output.py` independently checked 14 source hashes, 15
   covered Fisher/BH genetic tests, 72 candidate/cohort/reference records, 380
   condition records, score medians and 96 structured virtual-experiment records.
   `git diff --check` passed.

The first real attempt found a calibration group-count bug before any result or
completion marker was written. It was fixed, regression-tested with at least three
groups, and rerun in a fresh directory. The incomplete attempt remains explicitly
labelled; it is not included in biological conclusions.

## Biological results

Completed output:
`results/local_rechecks/ad-20261108-02/`.

- Genetic mode is **gwas_mapped**, not a completed full-GWAS TWAS. Bellenguez 2022
  mapped loci supply 62 genes in the NEU annotation background and 60 in legacy
  CMap's mapped background. No candidate enrichment survived BH q<.05.
- BX-795 had the smallest NEU enrichment P, .008556, but q=.128346 after all15
  covered primary/reference tests. Its18 overlapping mapped genes are traceable
  in `genetic_enrichment.json`; overlap alone is not statistical support.
- None of the12 candidates passed the primary NEU median-reversal calibration in
  both GSE132903 and donor-corrected GSE118553. There are no co-occurring primary
  genetic-plus-reversal hits under the declared method.
- Legacy CMap covers vorinostat, capsaicin and explicitly aliased trichostatin A.
  Vorinostat passed for GSE132903 (score=-.1735, q≈.0119), but not GSE118553
  (score=-.1165, q≈.3179). The other two did not pass both cohorts. The nine absent
  compounds are marked missing, not biologically negative.
- Prior upstream screen hits are not contradicted by a numerical mismatch alone:
  this recheck uses per-reference top100 up/down disease queries and median over
  all available candidate/cell conditions, plus matched random-query calibration.
  The upstream screen's best condition and earlier FDR-based signatures answer a
  different statistical question. The same LINCS source is reused, so even a
  successful recheck would not establish independent drug efficacy.

## Next evidence needed

Run the actual full-GWAS→TWAS pipeline with matched models/reference, then set the
recheck mode to `twas` with the full tested S-MultiXcan table. Follow significant
gene associations with locus-aware fine-mapping/colocalization and direction
checks before drug-target interpretation. Improve assay annotation harmonization,
donor/age/postmortem/cell-composition control, and independent neural perturbation
validation. The current local extension implements and tests the evidence join;
it does not claim these missing biological steps have been completed.
