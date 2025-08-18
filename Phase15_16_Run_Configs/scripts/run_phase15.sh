#!/usr/bin/env bash
set -euo pipefail
python src/run_phase15.py \
  --n_primes 1000000000 \
  --tmax 22000 \  --tstep 0.25 \  --method finufft \  --zeros_csv data/zeros_up_to_22000.csv \  --tmin_eval 20000 \  --outdir outputs_phase15
