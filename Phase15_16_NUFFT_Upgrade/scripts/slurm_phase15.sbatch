
#!/bin/bash
#SBATCH -J phase15_nufft
#SBATCH -o phase15_%j.out
#SBATCH -e phase15_%j.err
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=0
#SBATCH --time=24:00:00
# For GPU nodes, uncomment the next lines and adjust partition:
# SBATCH --gres=gpu:1
# SBATCH -p gpu

module purge
# module load anaconda/3  # or your env module
# source activate your_env

# CPU finufft (pip install finufft) or GPU cuFINUFFT (see their docs)
python src/run_phase15.py --n_primes 1000000000 --tmax 22000 --tstep 0.25 --method finufft --zeros_csv data/zeros_up_to_22000.csv --tmin_eval 20000 --outdir outputs_phase15
