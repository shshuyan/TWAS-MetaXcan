#!/bin/bash
# Run all 22 x 10 = 220 imputation jobs, NJOBS in parallel (re-runnable: finished jobs are skipped)
set -euo pipefail; source "$(dirname "$0")/config.sh"
for c in $(seq 1 22); do for b in $(seq 0 9); do echo "$c $b"; done; done \
 | xargs -P $NJOBS -n 2 bash $(dirname "$0")/02_impute_one.sh
