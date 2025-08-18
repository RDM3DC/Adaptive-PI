#!/usr/bin/env bash
set -euo pipefail
# Phase 16 uses sharded logs; see run_phase16_sharded.py
python src/run_phase16_sharded.py \  --tmax 50000 \  --tstep 0.20 \  --method finufft \  --zeros_csv data/zeros_up_to_50000.csv \  --tmin_eval 30000 \  --shard_glob "data/shards/logp_shard*.npy" \  --outdir outputs_phase16
