# TWAS-MetaXcan: Alzheimer's disease TWAS with S-PrediXcan / S-MultiXcan

A transcriptome-wide association study (TWAS) pipeline for **Alzheimer's disease (AD)**, following the
MetaXcan tutorial workflow (harmonize GWAS → summary imputation → S-PrediXcan → S-MultiXcan →
gene-level results) with GTEx v8 MASHR expression (eQTL) and splicing (sQTL) models. Built on
[MetaXcan](https://github.com/hakyimlab/MetaXcan) and
[summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation).

The pipeline is config-driven: each GWAS is one YAML file in `config/traits/`. Every input layout is
first converted to **GWAS-VCF** ([Lyon et al. 2021](https://doi.org/10.1186/s13059-020-02248-0), the
MRC-IEU OpenGWAS format) and handled with `bcftools`, so other diseases (e.g. chemotherapy-induced
peripheral neuropathy) only need a new config file.

## AD GWAS datasets

| Config | Study | Cases / controls | Build | Source |
|---|---|---|---|---|
| **`AD_Bellenguez2022.yaml`** (primary) | Bellenguez et al. 2022, *Nat Genet* (EADB stage 1) | 85,934 (39,106 clinical + 46,828 UKB proxy) / 401,577 | hg38 (GWAS Catalog harmonised) | [GCST90027158](https://www.ebi.ac.uk/gwas/studies/GCST90027158) |
| `AD_Kunkle2019.yaml` | Kunkle et al. 2019, *Nat Genet* (IGAP stage 1), clinically diagnosed late-onset AD | 21,982 / 41,944 | hg19 | [GCST007511](https://www.ebi.ac.uk/gwas/studies/GCST007511) |
| `AD_FinnGen_R12.yaml` | FinnGen R12 endpoint `G6_ALZHEIMER` (replication) | 13,964 / 486,384 | hg38 | [FinnGen](https://www.finngen.fi/en/access_results) |
| `IPF_FinnGen_R12.yaml` | Idiopathic pulmonary fibrosis, FinnGen R12 `IPF` | 2,700 / 494,583 | hg38 | FinnGen |
| `TEMPLATE.yaml` | Template for new traits | | | |

All three AD datasets are European ancestry, matching the European LD reference and GTEx models.
Bellenguez 2022 includes UK Biobank "proxy" cases (family history of dementia); Kunkle 2019 is the
largest clinically-diagnosed-only set, useful as a sensitivity analysis.

## Run

```bash
python twas.py run config/traits/AD_Bellenguez2022.yaml        # download → ... → gene-level results
bash run_detached.sh config/traits/AD_Bellenguez2022.yaml      # same, in the background (survives logout)
python twas.py run config/traits/AD_*.yaml                     # all three AD datasets in turn
python twas.py compare config/traits/AD_Bellenguez2022.yaml config/traits/AD_Kunkle2019.yaml config/traits/AD_FinnGen_R12.yaml
```

Runtime and resources: the full run (with summary imputation, 220 jobs) takes about 10–12 h with
`njobs: 2` and needs ~16 GB RAM (each imputation job peaks at ~6–7 GB). Setting `skip_imputation: true`
in the trait config runs S-PrediXcan directly on the harmonized GWAS (well under 1 h for 49 tissues).

## Workflow and how it maps to the MetaXcan tutorial

| Tutorial section | Step here | What it does |
|---|---|---|
| 3. Download data | `download` | Fetch the AD GWAS from `input.url` (reference bundle: Zenodo 3657902) |
| — | `standardize` | Map columns, derive missing beta/SE/z/p, write GWAS-VCF (`bcftools sort` + tabix) |
| — | `export` | `bcftools query` (+ optional `-i` filter / `-r` regions; autosomes) → MetaXcan input table |
| 4. Harmonize (hg19 → hg38) | `harmonize` | `gwas_parsing.py`: liftover when `build` is not hg38, align alleles to the 1000G/GTEx panel; also writes a harmonized hg38 GWAS-VCF |
| 5. Impute GWAS | `impute` | `gwas_summary_imputation.py`, 22 chromosomes × 10 sub-batches, parallel, failed jobs retried |
| 6. Merge 220 batches | `postprocess` | `gwas_summary_imputation_postprocess.py` |
| 7. Single tissue | `spredixcan` | `SPrediXcan.py` in every GTEx v8 tissue (MASHR models) |
| 8. Multiple tissues | `smultixcan` | `SMulTiXcan.py` across tissues |
| 9. Results | `summarize` | Gene-level tables with Bonferroni significance → `results/<ID>/` |

Every step skips work whose output already exists, so an interrupted run can be restarted.

## Results (tutorial section 9)

`results/<ID>/` contains:

- `spredixcan_eqtl_Whole_Blood.csv`: gene-level S-PrediXcan results for Whole Blood, the same table as
  the tutorial's "Gene level results" example
- `spredixcan_eqtl_significant.tsv`: Bonferroni-significant gene–tissue pairs, all tissues
- `smultixcan_eqtl_significant.tsv` / `smultixcan_eqtl_all.tsv.gz`: S-MultiXcan genes across tissues
- `SUMMARY.md`: counts, thresholds and top genes
- (local only, not in git) `spredixcan_eqtl_all.tsv.gz`: every gene–tissue test

Significance: Bonferroni 0.05 / (N genes × M tissues) for S-PrediXcan and 0.05 / N genes for S-MultiXcan.

Column meanings (S-PrediXcan):

| Column | Meaning |
|---|---|
| `gene` | Ensembl ID |
| `gene_name` | HUGO gene name |
| `zscore` | S-PrediXcan association z-score for the gene |
| `effect_size` | Effect size (when available) |
| `pvalue` | P-value of the z-score |
| `var_g` | Variance of predicted expression, W' G W (W = SNP weights, G = covariance) |
| `pred_perf_r2` | R² of the tissue model's prediction of measured expression |
| `pred_perf_pval` / `pred_perf_qval` | P-/q-value of that prediction performance |
| `n_snps_used` | GWAS SNPs used in the S-PrediXcan analysis |
| `n_snps_in_cov` | SNPs in the covariance matrix |
| `n_snps_in_model` | SNPs in the model |
| `best_gwas_p` / `largest_weight` | Smallest GWAS p among model SNPs / largest model weight |

## Web app

A small web front end searches public GWAS repositories, downloads a dataset and runs this pipeline on it.

```bash
bash webapp/start.sh          # then open http://localhost:8000 on the machine running it
```

1. **Search**: disease name (English) across FinnGen R12, Pan-UK Biobank (per-ancestry analyses) and the
   GWAS Catalog (studies with full summary statistics; harmonised hg38 files preferred). Filter by ancestry.
2. **Configure**: download URL, genome build (from GWAS Catalog metadata), column mapping (auto-detected,
   editable) and case/control counts are filled in; choose eQTL/sQTL, tissues, full or fast mode.
3. **Run**: writes `config/traits/<ID>.yaml` and queues the run (one pipeline at a time).
4. **Results**: S-MultiXcan and S-PrediXcan gene-level tables with Bonferroni flags and TSV download.

Network hosts used: `ftp.ebi.ac.uk` (GWAS Catalog), `storage.googleapis.com` (FinnGen),
`pan-ukb-us-east-1.s3.amazonaws.com` (Pan-UKBB), `zenodo.org` (reference bundle).

## Setup

```bash
git clone https://github.com/hakyimlab/MetaXcan
git clone https://github.com/hakyimlab/summary-gwas-imputation
conda env create -f envs_metaxcan.yaml        # S-PrediXcan / S-MultiXcan, runs twas.py and the web app
conda env create -f envs_gwasimputation.yaml  # harmonization and imputation
conda env create -f envs_vcftools.yaml        # bcftools / htslib
mkdir -p data && (cd data && curl -L "https://zenodo.org/api/records/3657902/files/sample_data.tar/content" | tar -xf -)
```

The Zenodo bundle (16.7 GB) provides the reference data: 1000G EUR genotype panel and LD blocks,
liftover chains, GTEx v8 MASHR eQTL/sQTL models and covariances. Paths to the interpreters, bcftools
and data live in `config/pipeline.yaml`.

## Adding a disease or dataset

1. Copy `config/traits/TEMPLATE.yaml` to `config/traits/<ID>.yaml` (or use the web app).
2. Set `input.path` / `input.url`, `input.format` (`finngen`, `gwas_catalog`, `igap`, `pan_ukbb`, `bbj`,
   `generic`, `custom` with a full `input.columns` map, or `gwasvcf`), `input.build` and sample sizes.
   `input.columns` can rename any column, including per-variant `n_cases` / `n_controls`.
3. `python twas.py run config/traits/<ID>.yaml`

`python twas.py presets` lists the built-in layouts.

### Useful bcftools commands on the intermediate files

```bash
V=work/<ID>/00_gwasvcf/<ID>.hg38.vcf.gz
bcftools query -i 'FORMAT/LP>7.3' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%LP]\n' $V   # genome-wide significant
bcftools view -r 19:44000000-46500000 $V                                                # APOE region
bcftools merge -m none work/A/02_harmonized/A.hg38.vcf.gz work/B/02_harmonized/B.hg38.vcf.gz
```

`compare` writes `results/cross_trait/`: a gene × dataset S-MultiXcan p-value matrix, genes significant
in ≥ 2 datasets, and a `bcftools merge` of the harmonized GWAS-VCFs (same hg38 alleles).

## Optional AD candidate evidence recheck (local extension)

After an upstream drug screen, a frozen candidate list can be checked against AD
genetic gene-set enrichment and CMap/neural LINCS expression reversal. The checks
remain separate, distinguish missing from negative results, and do not claim causal
targets or independent efficacy when reusing the screening reference.

```bash
# Edit a local copy of the example; no OpenAI/DeepSeek key is needed.
cp config/rechecks/ad_template.yaml config/rechecks/local_ad.yaml
python twas.py recheck config/rechecks/local_ad.yaml
# A new output directory is required for a changed/repeated run.
python twas.py recheck config/rechecks/local_ad.yaml --output /path/to/new-output
```

The example is intended for the pipeline hook's `${results}` variable. For the
standalone command, replace that variable with the actual full TWAS result path.
Inputs may instead explicitly choose `genetics.mode: gwas_mapped`; this uses
exploratory locus-mapped gene overrepresentation and is not a completed TWAS.

To opt into the last step of `twas.py run`, add this to a trait YAML:

```yaml
candidate_recheck:
  enabled: true
  config: ${root}/config/rechecks/local_ad.yaml
```

The hook exposes `${results}`, `${work}`, and `${trait_id}` to the recheck config.
It can also be invoked alone with `--steps candidate_recheck` after upstream data
are available. Existing traits without this setting skip the step. See
[the method and local test plan](docs/ad_candidate_recheck_design.md) and
[configuration/output details](docs/ad_candidate_recheck_usage.md).

## Notes

- **APOE region**: the APOE locus (chr19 ~44.9 Mb) dominates AD GWAS and produces many significant
  genes through LD (APOE, APOC1, TOMM40, NECTIN2, ...). Interpret genes in this region jointly. For a
  sensitivity analysis without it, add to the trait config (and use a new `id`):
  `vcf_filter: 'CHROM!="19" || POS<43900000 || POS>45900000'` (coordinates are in the input file's
  build: this window is hg38; for hg19 inputs such as Kunkle 2019 use 44400000–46400000).
- Raw rows whose alleles are not A/C/G/T strings (e.g. indels coded `I`/`D`) are dropped at `standardize`.
- GWAS-VCF convention: ALT = effect allele. The harmonized VCF uses panel REF/ALT.
- Ancestry: the reference is European. A trait whose `ancestry` differs is refused unless it supplies a
  matching `reference:` block or sets `allow_ancestry_mismatch: true`.
