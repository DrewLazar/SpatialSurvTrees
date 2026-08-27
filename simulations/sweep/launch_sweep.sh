#!/bin/bash
# Submit the robustness sweep. Two arms, submitted independently.
#
#   bash launch_sweep.sh tree   # 240 tree fits over 20 nodes
#   bash launch_sweep.sh sb     # 240 spBayes fits over 20 nodes
#   bash launch_sweep.sh both
#
# Node counts must match N_NODES in the corresponding runner.

ARM=${1:-both}

if [ "$ARM" = "tree" ] || [ "$ARM" = "both" ]; then
  for n in $(seq 1 20); do
    sbatch --export=ALL,NODE=$n submit_sweep.slurm
  done
  echo "submitted 20 tree jobs"
fi

if [ "$ARM" = "sb" ] || [ "$ARM" = "both" ]; then
  for n in $(seq 1 20); do
    sbatch --export=ALL,NODE=$n submit_sweep_sb.slurm
  done
  echo "submitted 20 spBayes jobs"
fi
