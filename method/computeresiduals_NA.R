# Residual computation using the Nelson-Aalen cumulative hazard estimator.
#
# Motivation: the Kaplan-Meier version computes e_i = -log(S_KM(t_i)), which
# diverges when the leaf KM curve reaches 0 (e.g. t_i at or beyond the largest
# event time in a small/heavy-censoring leaf). The previous code clamped
# S_t >= 1e-10, which silently produced e_i ~ 23, a large outlier, exactly in
# those leaves. The Nelson-Aalen estimator estimates the cumulative hazard
# Lambda(t) = sum_{t_j <= t} d_j / n_j directly; it is finite everywhere, is the
# natural nonparametric estimator of the object the theory is stated in terms of
# (the cumulative hazard), and removes the need for any clamp.
#
# Drop-in replacement for compute_residuals(). find_leaf() is unchanged.

# Find which leaf a patient falls into. `node` is a data.tree node (walked via
# node$children and isLeaf); `obj` is the R6 DipolarSurvivalTree object that
# carries the kernel method obj$K. This mirrors the class's own predrec().
find_leaf <- function(node, x, obj) {
  if (isLeaf(node)) {
    return(node)
  } else {
    X  <- as.matrix(node$data)
    w0 <- node$opt_w0_mupmdiff$w0
    mupmdiff <- node$opt_w0_mupmdiff$mupmdiff
    split <- (t(mupmdiff) %*% obj$K(X, x)) + w0
    if (split < 0) {
      find_leaf(node$children[[1]], x, obj)
    } else {
      find_leaf(node$children[[2]], x, obj)
    }
  }
}

# Nelson-Aalen cumulative hazard at time t_i, read from a survfit object.
# survfit stores, at each event/censoring time, n.risk and n.event; the
# Nelson-Aalen estimate is the cumulative sum of n.event/n.risk over all
# tabulated times <= t_i. This is finite for every t_i (it simply stops
# accumulating once past the last event), so no clamp is needed.
nelson_aalen_at <- function(km, t_i) {
  # km is a survfit object (the leaf's KMest slot)
  times  <- km$time
  nrisk  <- km$n.risk
  nevent <- km$n.event
  # increments d_j / n_j at each tabulated time
  incr <- ifelse(nrisk > 0, nevent / nrisk, 0)
  # cumulative hazard up to and including t_i
  idx <- which(times <= t_i)
  if (length(idx) == 0) return(0)        # before the first event time
  sum(incr[idx])
}

# Compute e_i = Lambda_NA_leaf(t_i) for all patients.
# `tree` is the data.tree returned by obj$createtree(...); `obj` is the R6
# DipolarSurvivalTree object (carries the kernel obj$K).
compute_residuals <- function(alldata, tree, obj, covariates, time, censor) {
  n <- nrow(alldata)
  e_resid <- numeric(n)
  testX <- as.matrix(alldata[, covariates])
  
  for (i in 1:n) {
    leaf <- find_leaf(tree, testX[i, ], obj)
    km   <- leaf$KMest                    # survfit object stored at the leaf
    t_i  <- alldata[i, time]
    e_resid[i] <- nelson_aalen_at(km, t_i)
  }
  
  return(e_resid)
}