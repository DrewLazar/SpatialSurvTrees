# =====================================================================
# One node's slice of the spBayesSurv arm of the robustness sweep.
# Invoked as:  Rscript run_sweep_node_sb.R <node_index>
#
# Separate from the tree arm (run_sweep_node.R) so the two can be
# submitted independently: the spBayes fits are far slower, so they
# want their own walltime and their own node count.
#
# Writes rep_cfX_seedY_sb.rds to the SAME sweep_out directory. The _sb
# suffix keeps the two arms from colliding, and the collate step reads
# both. Re-running a slice is safe: completed replicates are skipped.
# =====================================================================

## ---- paths: EDIT THESE for the cluster ----------------------------
# setwd("/home/dmlazar/.../robustness_sweep")

source("sim_sweep_generate.R")
source("sim_sweep_analysis.R")

## ---- sweep grid ----------------------------------------------------
CONF_GRID <- c(0, 0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.5)   # 8 levels
SB_SEEDS  <- 1:30                                         # replicates
N_NODES   <- 20                                           # nodes for THIS arm

## ---- which slice am I? --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
node <- as.integer(args[1])
if (is.na(node) || node < 1 || node > N_NODES)
  stop(sprintf("usage: Rscript run_sweep_node_sb.R <1..%d>", N_NODES))

# 8 levels x 30 seeds = 240 spBayes fits, same replication as the tree arm,
# so the two curves carry comparable standard errors.

jobs <- expand.grid(confound = CONF_GRID, seed = SB_SEEDS)
jobs <- jobs[order(jobs$confound, jobs$seed), ]
mine <- jobs[seq(node, nrow(jobs), by = N_NODES), , drop = FALSE]

cat(sprintf("spBayes arm: node %d/%d : %d of %d replicates\n",
            node, N_NODES, nrow(mine), nrow(jobs)))

t0 <- Sys.time()
for (i in seq_len(nrow(mine))) {
  ti <- Sys.time()
  run_sweep_replicate(confound = mine$confound[i],
                      seed     = mine$seed[i],
                      fit_sb   = TRUE,
                      out_dir  = "sweep_out")
  cat(sprintf("  [%d/%d] cf=%.2f seed=%d  %.1f min\n",
              i, nrow(mine), mine$confound[i], mine$seed[i],
              as.numeric(difftime(Sys.time(), ti, units = "mins"))))
}
cat(sprintf("spBayes node %d done in %.1f min\n", node,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
