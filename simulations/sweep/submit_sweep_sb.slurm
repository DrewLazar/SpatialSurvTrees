#!/bin/bash
#SBATCH --job-name=sweep_sb
#SBATCH --output=logs/sb_%j.out
#SBATCH --error=logs/sb_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# spBayesSurv arm of the robustness sweep: 12 MCMC fits per node.
#
# Longer walltime and more memory than the tree arm: survregbayes holds
# a 1043 x 1043 distance matrix plus 2000 retained samples per fit.
# Time one fit interactively first --
#   system.time(run_sweep_replicate(0, 1, fit_sb = TRUE))
# -- then set --time to roughly 12 x that, plus 50% headroom.
#
# NODE is passed by launch_sweep.sh via --export=ALL,NODE=$n
#
# Writes rep_cfX_seedY_sb.rds into the same sweep_out directory as the
# tree arm; the suffix keeps the two from colliding.

mkdir -p logs sweep_out

module load R            # EDIT: your cluster's R module

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "host $(hostname)  node index ${NODE}  start $(date)"
Rscript run_sweep_node_sb.R ${NODE}
echo "done $(date)"
