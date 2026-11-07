# TWAS of coronary artery disease with S-PrediXcan / S-MultiXcan

Transcriptome-wide association study (TWAS) of the CARDIoGRAMplusC4D coronary artery disease (CAD)
GWAS (additive model, N = 184,305, 60,801 cases), run with
[MetaXcan](https://github.com/hakyimlab/MetaXcan) and
[summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation)
using GTEx v8 MASHR expression models (49 tissues).

Input data: the MetaXcan tutorial bundle, [Zenodo 3657902](https://zenodo.org/records/3657902) (`sample_data.tar`, 16.7 GB).

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_harmonize.sh` | Parse the GWAS, lift over hg19 → hg38, align to the 1000G reference panel |
| 2 | `scripts/02_impute_all.sh` | Summary-statistics imputation: 22 chromosomes × 10 sub-batches = 220 jobs (`NJOBS` in parallel, re-runnable) |
| 3 | `scripts/03_postprocess.sh` | Merge the 220 imputed batches with the harmonized GWAS |
| 4 | `scripts/04_spredixcan.sh` | S-PrediXcan in each of the 49 GTEx v8 tissues |
| 5 | `scripts/05_smultixcan.sh` | S-MultiXcan: joint test across tissues |
| 6 | `scripts/06_summarize.py` | Bonferroni significance and summary tables in `results/` |

`scripts/run_all.sh` runs every step. Paths are in `scripts/config.sh`.

Significance thresholds: S-PrediXcan 0.05 / (number of gene–tissue tests); S-MultiXcan 0.05 / (number of genes).

## Setup

```bash
git clone https://github.com/hakyimlab/MetaXcan
git clone https://github.com/hakyimlab/summary-gwas-imputation
conda env create -f envs_metaxcan.yaml
conda env create -f envs_gwasimputation.yaml
mkdir -p data && cd data && curl -L "https://zenodo.org/api/records/3657902/files/sample_data.tar/content" | tar -xf - && cd ..
bash scripts/run_all.sh
```

## Results

See [`results/SUMMARY.md`](results/SUMMARY.md). Main tables:

- `results/spredixcan_significant.tsv`: Bonferroni-significant gene–tissue pairs
- `results/spredixcan_Whole_Blood.tsv`: all Whole Blood results
- `results/smultixcan_significant.tsv`: Bonferroni-significant genes across tissues
- `results/*.tsv.gz`: full result tables

Column meanings (S-PrediXcan): `gene` Ensembl ID; `gene_name` HUGO symbol; `zscore` association z-score;
`pvalue`; `var_g` predicted expression variance (W'GW); `pred_perf_r2` / `pred_perf_pval` model prediction
performance; `n_snps_used` / `n_snps_in_cov` / `n_snps_in_model` SNP counts.
