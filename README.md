# TWAS pipeline: GWAS-VCF + S-PrediXcan / S-MultiXcan

A config-driven transcriptome-wide association study (TWAS) pipeline built on
[MetaXcan](https://github.com/hakyimlab/MetaXcan) and
[summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation),
with GTEx v8 MASHR expression (eQTL) and splicing (sQTL) models.

Each disease is one YAML file in `config/traits/`. Any GWAS layout is first converted to
**GWAS-VCF** ([Lyon et al. 2021](https://doi.org/10.1186/s13059-020-02248-0), the MRC-IEU OpenGWAS format)
and then handled with `bcftools`. This gives every dataset the same representation, so new diseases
(neurodermatitis, pulmonary fibrosis, ...) only need a new config file.

Traits included:

| Config | Disease | Source | Build |
|---|---|---|---|
| `CAD_CARDIoGRAMplusC4D.yaml` | Coronary artery disease | CARDIoGRAMplusC4D (MetaXcan tutorial data) | hg19 |
| `LSC_FinnGen_R12.yaml` | Lichen simplex chronicus (neurodermatitis) | FinnGen R12 `L12_LICHENSIMPLEX` | hg38 |
| `IPF_FinnGen_R12.yaml` | Idiopathic pulmonary fibrosis | FinnGen R12 | hg38 |
| `TEMPLATE.yaml` | Template for new traits | | |

## Steps

```
raw GWAS (any layout, any build)
  │ download      fetch input.url if the file is missing
  │ standardize   map columns (preset + overrides), derive missing beta/SE/z/p, drop rows without
  │               ACGT alleles → bcftools sort + tabix → work/<id>/00_gwasvcf/<id>.<build>.vcf.gz
  │ export        bcftools query (+ optional -i filter / -r regions, autosomes by default)
  │               → MetaXcan input table
  │ harmonize     gwas_parsing.py: liftover to hg38 when needed, align alleles to the 1000G/GTEx panel
  │               → harmonized table + harmonized hg38 GWAS-VCF (for cross-trait merging)
  │ impute        summary-statistics imputation, 22 chromosomes × 10 sub-batches, run in parallel
  │ postprocess   merge imputed batches
  │ spredixcan    S-PrediXcan in every tissue (eqtl and/or sqtl)
  │ smultixcan    S-MultiXcan across tissues
  ▼ summarize     Bonferroni tables + results/<id>/SUMMARY.md
```

Every step skips work whose output already exists, so an interrupted run can simply be restarted.

## Usage

```bash
python twas.py run config/traits/CAD_CARDIoGRAMplusC4D.yaml                  # everything
python twas.py run config/traits/LSC_FinnGen_R12.yaml --steps standardize export harmonize
python twas.py run config/traits/IPF_FinnGen_R12.yaml --from spredixcan      # resume from a step
python twas.py run config/traits/*.yaml                                      # several traits
python twas.py compare config/traits/LSC_FinnGen_R12.yaml config/traits/CAD_CARDIoGRAMplusC4D.yaml
python twas.py presets                                                       # known input layouts
bash run_detached.sh config/traits/CAD_CARDIoGRAMplusC4D.yaml                # background run
```

`compare` writes `results/cross_trait/`: a gene × trait S-MultiXcan p-value matrix, genes significant
in ≥ 2 traits, and a `bcftools merge` of the harmonized GWAS-VCFs (one sample column per trait,
same hg38 alleles), which is useful for colocalization, MR and cross-trait look-ups.

### Adding a disease

1. Copy `config/traits/TEMPLATE.yaml` to `config/traits/<ID>.yaml`.
2. Set `input.path` (or `input.url`), `input.format` (`cardiogram`, `finngen`, `gwas_catalog`,
   `pan_ukbb`, `bbj`, `generic`, or `gwasvcf` for files that are already GWAS-VCF), `input.build`,
   and sample sizes. Use `input.columns` to rename any column the preset gets wrong.
3. `python twas.py run config/traits/<ID>.yaml`

### Useful bcftools commands on the intermediate files

```bash
V=work/<ID>/00_gwasvcf/<ID>.hg38.vcf.gz
bcftools query -i 'FORMAT/LP>7.3' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%LP]\n' $V   # genome-wide significant
bcftools view -r 6:28000000-34000000 $V                                                # MHC region
bcftools merge -m none work/A/02_harmonized/A.hg38.vcf.gz work/B/02_harmonized/B.hg38.vcf.gz
```

### Ancestry

The bundled reference (1000G LD, GTEx v8 models) is European. A trait whose `ancestry` differs from
`reference.ancestry` is refused unless you give it a matching `reference:` block (LD blocks,
genotype panel, models) or set `allow_ancestry_mismatch: true`. For East Asian GWAS (e.g. Biobank Japan,
Chinese cohorts), use EAS LD blocks and an EAS genotype panel at minimum.

## Setup

```bash
git clone https://github.com/hakyimlab/MetaXcan
git clone https://github.com/hakyimlab/summary-gwas-imputation
conda env create -f envs_metaxcan.yaml        # S-PrediXcan / S-MultiXcan, runs twas.py
conda env create -f envs_gwasimputation.yaml  # harmonization and imputation
conda env create -f envs_vcftools.yaml        # bcftools / htslib
mkdir -p data && (cd data && curl -L "https://zenodo.org/api/records/3657902/files/sample_data.tar/content" | tar -xf -)
```

Paths to the interpreters, bcftools and data live in `config/pipeline.yaml`.

## Notes

- Raw rows whose alleles are not A/C/G/T strings (e.g. indels coded `I`/`D`) are dropped at
  `standardize`, because they cannot be matched to the reference panel.
- GWAS-VCF convention: ALT = effect allele, `ES`/`EZ` refer to ALT. Before harmonization REF is the
  study's other allele, not a reference-genome check. The harmonized VCF uses panel REF/ALT.
- Bonferroni: S-PrediXcan 0.05 / (gene–tissue tests); S-MultiXcan 0.05 / genes.

## Results

`results/<ID>/SUMMARY.md` and tables:
`spredixcan_<type>_significant.tsv`, `smultixcan_<type>_significant.tsv`, `smultixcan_<type>_all.tsv.gz`
(full per-tissue S-PrediXcan tables are kept locally in `results/<ID>/spredixcan_<type>_all.tsv.gz`).

Column meanings (S-PrediXcan): `gene` Ensembl ID; `gene_name` HUGO symbol; `zscore`; `pvalue`;
`var_g` predicted expression variance (W'GW); `pred_perf_r2` / `pred_perf_pval` model performance;
`n_snps_used` / `n_snps_in_cov` / `n_snps_in_model` SNP counts.
