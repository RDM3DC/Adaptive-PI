
# Phase 15–16 NUFFT Upgrade

**Goal:** Scale H(t)=Re ∑ w_j e^{i t x_j} with **1B–10B primes** and many t values efficiently.

## Why NUFFT
H(t) is a nonuniform Fourier transform with nodes x_j=log p_j and weights w_j=1/sqrt(p_j).  
Using **FINUFFT** (or cuFINUFFT) drops cost from naive O(N·T) to near O(N + T log T).

## Install (CPU)
```bash
pip install finufft numpy
```
For **GPU**, see cuFINUFFT docs.

## Run (Phase 15 example)
```bash
python src/run_phase15.py --n_primes 1000000000 \
  --tmax 22000 --tstep 0.25 \
  --method finufft \
  --zeros_csv data/zeros_up_to_22000.csv --tmin_eval 20000 \
  --outdir outputs_phase15
```

## SLURM (template)
See `scripts/slurm_phase15.sbatch` and adjust cores/memory/partitions.

## Notes
- For 10B primes (Phase 16), use sharded prime logs and sum partial NUFFT blocks.
- Use double precision; consider Kahan/Neumaier compensated reductions if doing direct blocks.
- Keep JSON logs: counts, timings, seeds, hardware; save checksums for reproducibility.
