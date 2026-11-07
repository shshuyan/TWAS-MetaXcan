#!/bin/bash
set -euo pipefail; D=$(dirname "$0")
bash $D/01_harmonize.sh
bash $D/02_impute_all.sh
bash $D/03_postprocess.sh
bash $D/04_spredixcan.sh
bash $D/05_smultixcan.sh
/opt/conda/envs/MetaXcan/bin/python $D/06_summarize.py
