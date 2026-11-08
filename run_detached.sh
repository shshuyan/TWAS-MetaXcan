#!/bin/bash
# Run the pipeline for one or more traits in the background, surviving logout.
#   bash run_detached.sh config/traits/AD_Bellenguez2022.yaml [more traits...]
cd "$(dirname "$0")"
PY=${PY:-/opt/conda/envs/MetaXcan/bin/python}
LOG=logs/run_$(date +%Y%m%d_%H%M%S).log; mkdir -p logs
setsid nohup $PY twas.py run "$@" > "$LOG" 2>&1 < /dev/null &
echo "started, log: $LOG"
