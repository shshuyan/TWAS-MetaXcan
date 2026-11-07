#!/bin/bash
# Start the TWAS web app on http://localhost:${PORT:-8000}
cd "$(dirname "$0")/.."
PY=${PY:-/opt/conda/envs/MetaXcan/bin/python}
exec $PY -m uvicorn webapp.app:app --host ${HOST:-0.0.0.0} --port ${PORT:-8000}
