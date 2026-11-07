#!/bin/bash
cd /home/user/twas
while [ ! -f data/.download_done ]; do
  if ! pgrep -f "records/3657902" >/dev/null; then echo "download process died"; exit 1; fi
  sleep 60
done
echo "download done $(date)"
bash scripts/run_all.sh && echo "PIPELINE_DONE $(date)"
