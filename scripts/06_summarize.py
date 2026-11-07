"""Summarize S-PrediXcan / S-MultiXcan output into small tables under results/.

Significance follows the workflow PDF: Bonferroni 0.05 / (N genes x M tissues)
for S-PrediXcan, and 0.05 / N genes for S-MultiXcan.
"""
import glob
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output")
RES = os.path.join(ROOT, "results")
TRAIT = "CARDIoGRAM_C4D_CAD_ADDITIVE"
os.makedirs(RES, exist_ok=True)

# ---- S-PrediXcan, all tissues ----
frames = []
for f in sorted(glob.glob(os.path.join(OUT, "spredixcan", "eqtl", f"{TRAIT}__PM__*.csv"))):
    tissue = re.search(r"__PM__(.*)\.csv$", f).group(1)
    d = pd.read_csv(f)
    d.insert(0, "tissue", tissue)
    frames.append(d)
sp = pd.concat(frames, ignore_index=True)
sp = sp.dropna(subset=["pvalue"])
n_tests = len(sp)
thr = 0.05 / n_tests
sp = sp.sort_values("pvalue")
sp.to_csv(os.path.join(RES, "spredixcan_all_tissues.tsv.gz"), sep="\t", index=False)
sig = sp[sp.pvalue < thr]
sig.to_csv(os.path.join(RES, "spredixcan_significant.tsv"), sep="\t", index=False)
wb = sp[sp.tissue == "Whole_Blood"]
wb.to_csv(os.path.join(RES, "spredixcan_Whole_Blood.tsv"), sep="\t", index=False)

# ---- S-MultiXcan ----
sm = pd.read_csv(os.path.join(OUT, "smultixcan", "eqtl", f"{TRAIT}_smultixcan.txt"), sep="\t")
sm = sm.dropna(subset=["pvalue"]).sort_values("pvalue")
sm_thr = 0.05 / len(sm)
sm.to_csv(os.path.join(RES, "smultixcan.tsv.gz"), sep="\t", index=False)
sm_sig = sm[sm.pvalue < sm_thr]
sm_sig.to_csv(os.path.join(RES, "smultixcan_significant.tsv"), sep="\t", index=False)

# ---- Markdown summary ----
top_sp = (sig.groupby(["gene", "gene_name"])
             .agg(n_tissues=("tissue", "nunique"), best_p=("pvalue", "min"),
                  best_tissue=("tissue", "first"))
             .reset_index().sort_values("best_p"))
lines = [
    f"# TWAS results: {TRAIT} (coronary artery disease)",
    "",
    "## S-PrediXcan (GTEx v8 MASHR, all tissues)",
    f"- Tissues: {sp.tissue.nunique()}; gene-tissue tests: {n_tests:,}",
    f"- Bonferroni threshold: {thr:.2e}",
    f"- Significant gene-tissue pairs: {len(sig):,} ({top_sp.shape[0]} unique genes)",
    f"- Significant in Whole_Blood: {(wb.pvalue < thr).sum()}",
    "",
    "Top 30 genes (by best p-value across tissues):",
    "",
    "| gene | n tissues sig. | best tissue | best p |",
    "|---|---|---|---|",
]
for r in top_sp.head(30).itertuples():
    lines.append(f"| {r.gene_name} | {r.n_tissues} | {r.best_tissue} | {r.best_p:.2e} |")
lines += [
    "",
    "## S-MultiXcan (joint across tissues)",
    f"- Genes tested: {len(sm):,}; Bonferroni threshold: {sm_thr:.2e}",
    f"- Significant genes: {len(sm_sig)}",
    "",
    "Top 30 genes:",
    "",
    "| gene | p-value | n tissues used | best single-tissue p |",
    "|---|---|---|---|",
]
for r in sm.head(30).itertuples():
    lines.append(f"| {r.gene_name} | {r.pvalue:.2e} | {r.n_indep if hasattr(r, 'n_indep') else ''} | {r.p_i_best:.2e} |")
with open(os.path.join(RES, "SUMMARY.md"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
