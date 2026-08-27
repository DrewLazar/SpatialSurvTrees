#!/bin/bash
#SBATCH --job-name=sweep_tree
#SBATCH --output=logs/tree_%j.out
#SBATCH --error=logs/tree_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Tree arm of the robustness sweep: 12 replicates per node, each fitting
# a spatial-only tree, a clinical tree and a two-stage spatial tree.
#
# NODE is passed by launch_sweep.sh via --export=ALL,NODE=$n
#
# Single-threaded on purpose: the work is parallelized across jobs, and
# letting BLAS spawn threads inside each job oversubscribes the node.

mkdir -p logs sweep_out

module load R            # EDIT: your cluster's R module

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "host $(hostname)  node index ${NODE}  start $(date)"
Rscript run_sweep_node.R ${NODE}
echo "done $(date)"
