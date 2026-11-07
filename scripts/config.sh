# Shared paths — edit to match your machine
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
GWAS_TOOLS=$ROOT/summary-gwas-imputation/src
METAXCAN=$ROOT/MetaXcan/software
DATA=$ROOT/data/data
OUTPUT=$ROOT/output
TRAIT=CARDIoGRAM_C4D_CAD_ADDITIVE
PY_IMP=${PY_IMP:-/opt/conda/envs/gwasimputation/bin/python}
PY_MX=${PY_MX:-/opt/conda/envs/MetaXcan/bin/python}
NJOBS=${NJOBS:-4}
