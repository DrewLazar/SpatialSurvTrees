# =====================================================================
# One node's slice of the robustness sweep.
# Invoked as:  Rscript run_sweep_node.R <node_index>
#
# TREE ARM only: 30 seeds x 8 levels = 240 fits.
# The spBayesSurv arm is submitted separately via run_sweep_node_sb.R,
# which writes ..._sb.rds into the same sweep_out directory.
#
# Re-running a slice is safe: completed replicates are skipped.
# =====================================================================

## ---- paths: EDIT THESE for the cluster ----------------------------
# setwd("/home/dmlazar/.../robustness_sweep")

source("sim_sweep_generate.R")
source("sim_sweep_analysis.R")

## ---- sweep grid ----------------------------------------------------
CONF_GRID <- c(0, 0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.5)   # 8 levels
SEEDS     <- 1:30                                         # replicates
N_NODES   <- 20

## ---- which slice am I? --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
node <- as.integer(args[1])
if (is.na(node) || node < 1 || node > N_NODES)
  stop(sprintf("usage: Rscript run_sweep_node.R <1..%d>", N_NODES))

jobs <- expand.grid(confound = CONF_GRID, seed = SEEDS)
jobs <- jobs[order(jobs$confound, jobs$seed), ]
mine <- jobs[seq(node, nrow(jobs), by = N_NODES), , drop = FALSE]

cat(sprintf("tree arm: node %d/%d : %d of %d replicates\n",
            node, N_NODES, nrow(mine), nrow(jobs)))

t0 <- Sys.time()
for (i in seq_len(nrow(mine))) {
  run_sweep_replicate(confound = mine$confound[i],
                      seed     = mine$seed[i],
                      fit_sb   = FALSE,
                      out_dir  = "sweep_out")
}
cat(sprintf("node %d done in %.1f min\n", node,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
