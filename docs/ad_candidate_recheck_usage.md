# Local AD candidate recheck

This extension adds `twas.py recheck <YAML>` and the opt-in final pipeline step
`candidate_recheck`. It accepts the AD agent's existing candidate JSON directly.
It does not add a new upstream drug screen to MetaXcan and does not execute TWAS
when invoked as an isolated recheck. Standard MetaXcan environment dependencies
already include pandas, NumPy, SciPy, h5py and PyYAML; no LLM API or paid service is
required. The local verification used the AD agent's existing Python environment.

## Candidate and reference contracts

- Candidate JSON: list with `pert_id`, `pert_name`; or TSV with `candidate_id`,
  `drug`. Optional semicolon-separated `targets` or the AD agent's JSON
  `target_context` can annotate direct targets. IDs must be unique.
- LINCS JSONL: `sig_id`, `pert_id`, `cell_id`, `up`, `down`; optional `dose`, `time`,
  `signature_url`. Same pert_id, not approximate name, is required. Duplicate
  signature IDs and contradictory up/down sets are rejected.
- LINCS universe TSV: one `gene` column (override `gene_column` if needed), from an
  actual assay/inference annotation, never the union of observed hit tails alone.
- Optional CMap HDF5: `assay`, `rownames` (gene identifiers), `colnames` in
  `compound__cell__...` form. This implementation explicitly accepts lower-is-up
  legacy ranks, not arbitrary logFC/GCTX formats. Mapping TSV has `cmap_gene_id`
  and `symbol`; ambiguous duplicate symbols are excluded from queries. Gene ranks
  are preserved over the complete reference gene universe.
- Disease DE: `gene`, `log_fc`, `p_value`, `fdr_bh` (limma aliases accepted).
  Require at least two different files, including discovery and validation roles.
  Identical DE content with different filenames is rejected as replication.

## Genetic input modes (no silent fallback)

For completed TWAS, use:

```yaml
genetics:
  mode: twas
  complete_tested_results: true
  path: ${results}/smultixcan_eqtl_all.tsv.gz
  gene_column: gene_name
  pvalue_column: pvalue
  alpha: 0.05
```

Provide **all tested genes**, including nonsignificant ones. There must be exactly
one result per nonmissing gene symbol, e.g. full S-MultiXcan or one preselected
S-PrediXcan tissue. A significant-only table or pooled repeated tissue rows cannot
define the required background. Bonferroni uses the full valid input before
intersection with each assay. TWAS P-values supply no direction to reversal.

For an explicitly limited local GWAS mapping test, use:

```yaml
genetics:
  mode: gwas_mapped
  path: /path/to/single-study-locus-mapped-genes.tsv
  study: Bellenguez2022 PMID35379992 GCST90027158
  gene_column: mapped_genes
  pvalue_column: pvalue
  threshold: 5.0e-8
```

The table has one variant association per row and semicolon-separated mapped genes.
Do not mix studies or count linked SNPs as independent replication. Its background
is each declared assay gene universe. This mode is exploratory overrepresentation
of locus-mapped genes, not an LD-aware genetic/heritability enrichment, a completed
TWAS, causal mapping or risk-allele expression direction. It must stay labelled
`gwas_mapped` in reports.

## Execution and output

All relative paths resolve against the YAML directory; `${root}` is the local
TWAS-MetaXcan checkout. `${results}`, `${work}`, `${trait_id}` are only supplied by
the trait pipeline hook. For standalone use, specify actual input paths.

```bash
python twas.py recheck config/rechecks/local_ad.yaml
python twas.py recheck config/rechecks/local_ad.yaml --output /path/to/new-run
python twas.py run config/traits/your_AD_trait.yaml --steps candidate_recheck
python -m unittest discover -s tests -p 'test_*.py' -v
```

Outputs include Chinese `report_zh.md`, candidate evidence matrix, genetic Fisher
tables, per-cohort reversal, per-signature conditions/genes, actual disease queries,
known-target annotations and input SHA-256 manifest. Full CMap results are retained
for all reference compounds used in multiple-testing correction. `summary.json`
is written only after success. Failures leave no completion marker and are not
treated as negatives. Existing directories are never overwritten; the pipeline
can reuse unchanged complete output, but changed inputs/code need a new directory.

The combined record uses primary NEU genetic enrichment and primary NEU reversal
in every supplied cohort. Legacy CMap and NPC results remain explicit secondary
columns, including inconsistencies and coverage gaps. No P-values from correlated
sources are combined; `wet_lab_ready` is always false. Reused LINCS cannot become
new independent evidence merely because a second score was calculated.

## Deliberate boundaries

Not implemented: SNP-level LD-aware enrichment, colocalization/fine-mapping,
protein target direction, structure identity verification, BBB exposure, safety,
binding simulation or automatic wet-lab approval. These require separate evidence.
The gene-set null also does not model coexpression or the upstream selection of
candidates from the same reference, so all local biological conclusions remain
exploratory.
