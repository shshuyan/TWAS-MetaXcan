# TWAS-led AD drug discovery

## Placement and scope

The new `drug_screen` step follows `summarize`:

GWAS → standardize/align → S-PrediXcan → S-MultiXcan → signed TWAS queries →
whole supplied LINCS/CMap library → candidate drugs.

This replaces the earlier fixed-shortlist recheck extension; that version is
recoverable from Git history but is not an active command or pipeline step.
`drug-screen` neither requires that shortlist nor runs GEO differential expression.
The new command reads completed TWAS results; it does not itself execute the upstream
GWAS/TWAS tools, download reference data, retrieve target literature or approve wet lab.

## Required inputs

1. Full S-PrediXcan expression results, including nonsignificant tests: `gene_name`,
   `tissue`, `pvalue`, `zscore`. The existing `spredixcan_eqtl_all.tsv.gz` is suitable.
2. Full S-MultiXcan results: `gene_name`, `pvalue`, one row per gene. The existing
   `smultixcan_eqtl_all.tsv.gz` is suitable. No signed joint effect is assumed.
3. At least one perturbation reference:
   - LINCS curated compound JSONL: `sig_id`, `pert_id`, `pert_name`, `cell_id`,
     `up`, `down`, optional `dose`, `time`, `signature_url`. If `pert_type` is
     supplied, only `trt_cp` is accepted. Declare `compound_only: true` and an
     explicit cell list. Include the actual assay/inference gene annotation TSV;
     the union of signature tails is not a valid measurement background.
   - Optional legacy CMap build02 full rank HDF5 and Entrez-symbol mapping, using
     the legacy HDF5 contract (`assay`, `rownames`, `colnames` in compound__cell__... form); lower ranks mean drug upregulation. Every
     available CMap compound/cell is screened. Arbitrary GCTX/logFC formats are
     not silently accepted as this reference format.

Ambiguous duplicate genes, missing tissues/Z-scores, unsupported splicing models,
contradictory signature directions and an uncovered primary cell fail explicitly.
Invalid/missing gene-association rows are retained with `valid_test: false` and
excluded from queries; nonmissing `status` must be blank, `ok` or `success`.
`complete_tested_results: true` is the user's input contract, not proof that no
rows were removed upstream. Study, build, ancestry, allele alignment and expression
prediction model quality must be checked in the upstream pipeline. These inputs
must come from the same GWAS/trait. Existing SNP coverage/performance columns are
preserved, not converted into an unreported quality filter.

## Prespecified direction and gene selection

- Explicitly confirm `trait_direction: higher_is_AD_risk` after checking GWAS
  beta/Z orientation. This is not inferred from a filename or allele label.
- Set a list of relevant tissues and one primary tissue **before** examining hits.
  S-PrediXcan Z > 0 means predicted expression is positively associated with AD
  risk, Z < 0 means negatively associated. Neither is measured case-control DE.
- By default, select single-tissue BH q < .05 across all prespecified valid
  gene-tissue tests **and** S-MultiXcan Bonferroni P < .05/N across all valid joint
  tests. Compute both corrections before assay intersection. These are correlated
  views of the same GWAS, not independent replication.
- `require_joint_significance: false` is an explicit alternate analysis using
  single-tissue BH selection while retaining joint results as annotation. It is
  never activated automatically to obtain candidates. A second analysis needs a
  fresh output directory and should be reported as sensitivity analysis.
- Significant opposite directions among the selected tissues are marked. Default
  `exclude_direction_conflicts: true` excludes those genes from queries; false
  keeps separate tissue directions, never an averaged sign. This only assesses
  conflicts in the configured tissues, not every possible tissue.
- Preserve all eligible genes by default (`max_genes_per_direction: null`). An
  optional explicit cap is applied downstream by single-tissue P, with gene-name
  tie breaks. At least 10 genes per direction per reference are required by the
  template. Fewer yields `insufficient_query`, not relaxed P/FDR or filled signs.

Official field definitions:
https://predictdb.hakyimlab.org/post/2022/03/08/metaxcan-output-file-formats/

## Reversal and candidate definition

For LINCS, score each condition as (same-direction overlap − opposite-direction
overlap) / query size. For CMap, use the existing signed two-query unweighted KS.
Negative is reversal. These scores are not CLUE tau and are not interchangeable.
Aggregate by median **within each compound × cell**, retaining every dose/time
condition; do not pool a favorable cell with an unfavorable cell or select its
best condition after seeing the result.

Calibrate against 10,000 random queries with the same positive/negative lengths.
Background is all eligible tested TWAS genes in that tissue intersected with the
reference measurement universe (and the joint-tested universe when requiring
joint significance), not just significant hits. Empirical P is
`(1 + count(null <= observed)) / (B + 1)`. One BH family covers all tested
source × tissue × drug × cell groups. Tests with insufficient queries have no P
and are retained as missing evidence. Finite B limits the smallest resolvable P.
The report explicitly flags insufficient resolution for an isolated hit: with M
tests, its attainable BH q is approximately M/(B+1) when all other P-values are
uninformative. A 10,000-query run across thousands of groups can therefore miss
an isolated signal even if its score is extreme. More permutations can be
prespecified for a new run; the program does not manufacture smaller P-values or
switch to an uncalibrated tail extrapolation.

The candidate list uses one prespecified primary source/cell/tissue with negative
median score and global q < alpha. Every other context is retained in the full
ranking and candidate annotations. A secondary-context hit is not silently
promoted when the primary analysis fails. CMap/LINCS chemical identities are not
merged by fuzzy names; no cross-library agreement or disagreement is asserted
without verified mapping. CMap-only screening is supported if explicitly primary.

The random-gene null does not preserve LD, gene coexpression, model ascertainment
or perturbation-tail biases. Thus q-values are exploratory calibration, not an
LD-corrected heritability test or clinical efficacy probabilities. TWAS may reflect
shared LD rather than causal expression mediation. Follow candidates with
colocalization/fine-mapping, target-direction evidence, independent neural
perturbations, exposure and toxicity. GEO can be a later orthogonal expression
check; it is not obligatory for this discovery step.

## Commands and pipeline hook

Copy `config/drug_screens/ad_template.yaml` to `local_ad.yaml`, set available
reference paths and full TWAS result paths. Relative paths resolve against the
YAML directory; `${root}` means this repository. `${results}`, `${work}` and
`${trait_id}` are supplied only by the trait hook, so standalone configs must
replace them with actual paths.

```bash
python twas.py drug-screen config/drug_screens/local_ad.yaml
python twas.py drug-screen config/drug_screens/local_ad.yaml --output /path/to/new-run
python twas.py run config/traits/AD_Bellenguez2022.yaml --steps drug_screen
```

The last command needs this trait setting:

```yaml
drug_screen:
  enabled: true
  config: ${root}/config/drug_screens/local_ad.yaml
```

For an end-to-end configured run omit `--steps`; upstream TWAS runs first. The
default AD configuration does not force an unavailable local perturbation resource
on other users. No LLM call or API key is required. The workflow outputs
`candidate_drugs.tsv/json`, `drug_screen.tsv/json`, all annotated TWAS input rows,
actual query/background genes, per-condition evidence, a Chinese report and
structured virtual experiments. `input_kind: synthetic_test` marks software
fixtures; `real_twas` is a declaration, not automatic biological validation.

Fresh runs never overwrite an existing directory. Pipeline reuse checks resolved
configuration, input/implementation hashes and output hashes. Changed, incomplete,
tampered or forced runs need a new directory. A completed run with zero queries
is explicitly `query_status: insufficient_query`, not negative drug biology.

## Software validation

`python -m unittest discover -s tests -p 'test_*.py' -v` covers known synthetic
reversing and mimicking compounds, discovery without candidates/GEO, cell isolation,
global correction, reversing single-tissue signs while joint z_mean is fixed,
direction conflicts, joint-gate semantics, all-gene backgrounds, insufficient
queries, CLI and pipeline resume integrity. These fixtures validate implementation,
not a completed real AD GWAS-to-drug discovery experiment.

## Four-disease Level 5 primary analysis

The newer four-disease command is separate from the legacy AD tail-overlap/KS screen above:

```bash
python scripts/run_minimal_four_disease_demo.py --output demo_output/four_disease
python twas.py four-disease-screen config/drug_screens/four_disease_level5.yaml
```

It reads GSE92742 Level 5 MODZ in bounded signature blocks and uses all finite TWAS genes shared with
CMap. TWAS Z is symmetrically clipped at 5 and weighted by `min(|Z|, 5) / 5` in a weighted cosine.
Negative is reversal. Signatures are aggregated by the median within compound × cell before matched
signature-label permutations (cell/dose/time hierarchy with replicate count preserved). This avoids
rescoring the 473,647 × 12,328 matrix for every null draw. One BH family per disease covers every compound × primary-cell test. A
candidate compound must have at least one primary-cell test with global q < .05 and no positive score
in any tested primary cell; it appears once in the disease candidate table. Secondary cells remain
separate sensitivity evidence.

The formal configuration fixes 10,000 permutations; finite permutation resolution is reported. The
minimum demo deliberately uses 999 and synthetic data. It proves I/O, statistics, multiplicity,
lineage and figure generation, but it does not validate any drug–disease relationship. Real runs must
retain the original compressed CMap archive/checksums, compute from a decompressed local-scratch GCTX,
freeze all GWAS accessions and sample-overlap decisions, and retain all generated manifests.
