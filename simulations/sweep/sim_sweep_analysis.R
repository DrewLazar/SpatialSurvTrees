# =====================================================================
# ROBUSTNESS SWEEP ANALYSIS : one replicate per call
# ---------------------------------------------------------------------
# Fits spatial-only and two-stage on a (confound, seed) dataset using
# the SAME machinery as the three-zone recovery analysis (identical Dip
# constructors, kappa, pruning constant), then records:
#
#   GLOBAL   3-zone ARI + accuracy  <- primary y-axis (trusted metric)
#   HOTSPOT  AUC + accuracy on the confounded high zone  <- secondary
#   REALIZED cor(tpi, recovered risk)  <- x-axis, matches LeukSurv 0.05
#
# Their divergence at high confounding is the identifiability story:
# two-stage should keep medium/low while losing high specifically.
#
# spBayes is NOT in the loop (MCMC cost). Run it separately at one
# reference confound level if a benchmark point is wanted.
#
# One small .rds per replicate -> crash-safe, resumable, toppable.
#
# USAGE
#   source("sim_sweep_generate.R"); source("sim_sweep_analysis.R")
#   run_sweep_replicate(confound = 0, seed = 1, debug = TRUE)   # smoke
#   for (cf in CONF_GRID) for (sd in 1:3) run_sweep_replicate(cf, sd)
#   res <- collect_sweep()
# =====================================================================

suppressMessages({ library(data.tree); library(survival) })
# Method files: flat alongside these scripts on the cluster, or under
# ../../method/ in the repository layout.
SRC <- if (file.exists("dipole_tree.R")) "." else "../../method"
source(file.path(SRC, "dipole_tree.R"))
source(file.path(SRC, "computeresiduals_NA.R"))
source(file.path(SRC, "bootstrap_pruning.R"))

CONF_GRID <- c(0, 0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.5)   # dense near 0

# compute_residuals (and possibly other sourced helpers) reference `time`
# and `censor` from the global environment rather than their own arguments.
# In the recovery scripts these were defined at top level, which masked the
# issue; here the fitting happens INSIDE run_sweep_replicate, so we must
# also define them globally or a bare `time` resolves to stats::time (a
# function) and Surv() fails with "invalid subscript type 'closure'".
time <- "stop"; censor <- "status"

# ---- scoring helpers -------------------------------------------------
.adj_rand <- function(a,b){ tab<-table(a,b); N<-sum(tab); sc<-function(x)sum(choose(x,2))
  ai<-rowSums(tab); bj<-colSums(tab)
  (sc(as.vector(tab))-sc(ai)*sc(bj)/choose(N,2))/(0.5*(sc(ai)+sc(bj))-sc(ai)*sc(bj)/choose(N,2)) }

.score3 <- function(risk, truth_k){
  km <- kmeans(risk, centers = 3, nstart = 20)$cluster
  ord <- order(tapply(risk, km, mean))
  remap <- match(km, ord)
  list(acc = mean(remap == truth_k), ari = .adj_rand(km, truth_k))
}

.auc <- function(s, th){                      # orientation-free rank AUC
  r <- rank(s); n1 <- sum(th); n0 <- sum(!th)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  a <- (sum(r[th]) - n1*(n1+1)/2) / (n1*n0); max(a, 1-a) }

.acc_hot <- function(risk, th){
  km <- kmeans(risk, centers = 2, nstart = 10)
  ch <- risk > mean(km$centers)
  max(mean(ch == th), 1 - mean(ch == th)) }

.intvar <- function(data, covs){
  mu <- colMeans(data[covs]); mean(rowSums((data[covs]-matrix(mu,nrow(data),length(mu),byrow=TRUE))^2)) }

# centred log relative risk -- same construction as the LeukSurv script
.clr <- function(pt){ rr <- median(pt)/pt; log(rr) - mean(log(rr)) }

# =====================================================================
## Settings matched to the two-, three- and four-zone reruns.
ADAPTIVE  <- TRUE; EPS_PROBS <- 0.23; EPS_FLOOR <- 0.2
KAPPA     <- exp(-2)

run_sweep_replicate <- function(confound, seed, out_dir = "sweep_out",
                                hotspot_zone = "high",
                                prune_c = 3.7, nsize = 20,
                                fit_sb = FALSE,
                                verbose = TRUE, debug = FALSE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tag <- sprintf("cf%.3f_seed%03d%s", confound, seed, if (fit_sb) "_sb" else "")
  out_file <- file.path(out_dir, paste0("rep_", tag, ".rds"))
  if (file.exists(out_file) && !debug) {
    if (verbose) cat("skip (exists):", out_file, "\n"); return(invisible(out_file)) }

  gen <- make_sweep_data(confound = confound, seed = seed,
                         hotspot_zone = hotspot_zone)
  d <- gen$sim_data; n <- nrow(d)
  truth_k <- as.integer(factor(d$zone, levels = c("low","medium","high")))
  is_hot  <- gen$is_hot
  time <- "stop"; censor <- "status"
  quantiles <- c(.25,.75); tolerance <- 1e-2
  cov_cl <- c("num_age","num_wbc","num_tpi")
  cov_sp <- c("num_xcoord","num_ycoord")

  d_sp <- c(dist(d[,cov_sp]))
  eps_sp  <- as.numeric(quantile(d_sp[d_sp != 0], EPS_PROBS))
  ksig_sp <- .intvar(d, cov_sp)

  ## The two arms are disjoint: the spBayes arm does NOT refit the trees,
  ## so a (confound, seed) pair contributes tree metrics from exactly one
  ## row and spBayes metrics from exactly one row. Aggregating over all
  ## rows with na.rm therefore never double counts.
  if (!fit_sb) {

  # ---- spatial-only ---------------------------------------------------
  Dip.sp <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
    d, time, censor, cov_sp, quantiles, tolerance,
    epsilon = eps_sp, kappa = KAPPA, nsize = nsize,
    pureweight = 1, mixedweight = 1,
    Ksigma = ksig_sp, Kconstant = 1, Kpoly_order = 2,metric="kernel",
    gaussweight = .60, polyweight = .40, linearweight = 0,
    gausscovariates_index = 1:2, polycovariates_index = 1:2,
    linearcovariates_index = 1:2, ncovariatestosearch = 2,
    adaptive = ADAPTIVE, epsilon_floor_alpha = EPS_FLOOR, probs = EPS_PROBS)
  tree.sp <- Dip.sp$createtree(1:n)
  tree.sp <- bootstrapPruning(tree.sp, Dip.sp, prune_c, time = time, censor = censor)[[4]]
  pt.sp <- Dip.sp$predicttime(d, tree.sp); risk.sp <- -pt.sp

  # ---- two-stage ------------------------------------------------------
  d_cl <- c(dist(d[,cov_cl]))
  eps_cl <- as.numeric(quantile(d_cl[d_cl != 0], EPS_PROBS))
  Dip.cl <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
    d, time, censor, cov_cl, quantiles, tolerance,
    epsilon = eps_cl, kappa = KAPPA, nsize = 10,
    pureweight = 1, mixedweight = 1,metric="kernel",
    Kconstant = 0, Kpoly_order = 1,
    gaussweight = 0, polyweight = 0, linearweight = 1,
    linearcovariates_index = 1:3, gausscovariates_index = 1:3,
    polycovariates_index = 1:3, ncovariatestosearch = 3,
    adaptive = ADAPTIVE, epsilon_floor_alpha = EPS_FLOOR, probs = EPS_PROBS)
  tree.cl <- Dip.cl$createtree(1:n)
  d$e_resid <- compute_residuals(d, tree.cl, Dip.cl, cov_cl, time, censor)
  Dip.ts <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
    d, "e_resid", censor, cov_sp, quantiles, tolerance,
    epsilon = eps_sp, kappa = KAPPA, nsize = nsize,
    pureweight = 1, mixedweight = 1,metric="kernel",
    Ksigma = ksig_sp, Kconstant = 1, Kpoly_order = 2,
    gaussweight = .60, polyweight = .40, linearweight = 0,
    gausscovariates_index = 1:2, polycovariates_index = 1:2,
    linearcovariates_index = 1:2, ncovariatestosearch = 2,
    adaptive = ADAPTIVE, epsilon_floor_alpha = EPS_FLOOR, probs = EPS_PROBS)
  tree.ts <- Dip.ts$createtree(1:n)
  tree.ts <- bootstrapPruning(tree.ts, Dip.ts, prune_c, time = "e_resid", censor = censor)[[4]]
  pt.ts <- Dip.ts$predicttime(d, tree.ts); risk.ts <- -pt.ts

  } else {
    risk.sp <- risk.ts <- pt.sp <- pt.ts <- NULL
    tree.sp <- tree.cl <- tree.ts <- list(leafCount = NA_integer_)
    d$e_resid <- NA_real_
  }

  # ---- spBayesSurv PH + GRF (optional arm) -----------------------------
  # Same settings as the base-case benchmark in the zone designs, so the
  # sweep values are comparable with Table 5.
  risk.sb <- NULL
  if (fit_sb) {
    suppressPackageStartupMessages({ library(spBayesSurv); library(fields) })
    dsb <- d; dsb$id <- 1:n
    fit.sb <- tryCatch(
      survregbayes(Surv(stop, status) ~ num_age + num_wbc + num_tpi +
                     frailtyprior("grf", id),
                   data = dsb, survmodel = "PH", dist = "loglogistic",
                   mcmc  = list(nburn = 2000, nsave = 2000, nskip = 1, ndisplay = 1e6),
                   prior = list(maxL = 15, a0 = 1, b0 = 1,
                                nknots = 30, nblock = 30, nu = 1),
                   DIST = function(x1, x2) rdist(x1, x2),
                   Coordinates = as.matrix(dsb[, cov_sp]),
                   InitParamMCMC = FALSE),
      error = function(e) { cat("  spBayes error:", conditionMessage(e), "\n"); NULL })
    if (!is.null(fit.sb)) {
      v <- if (!is.null(fit.sb$v)) fit.sb$v else fit.sb$frail
      if (!is.null(v)) risk.sb <- if (is.matrix(v)) rowMeans(v) else v
    }
  }

  # ---- metrics --------------------------------------------------------
  if (!fit_sb) {
    g.sp <- .score3(risk.sp, truth_k);  g.ts <- .score3(risk.ts, truth_k)
    lrr.sp <- .clr(pt.sp); lrr.ts <- .clr(pt.ts)
  } else {
    g.sp <- g.ts <- list(acc = NA_real_, ari = NA_real_)
    lrr.sp <- lrr.ts <- NA_real_
  }
  .sc <- function(a, b) if (all(is.na(a))) NA_real_ else
    cor(a, b, method = "spearman")
  if (!is.null(risk.sb)) {
    g.sb   <- .score3(risk.sb, truth_k)
    lrr.sb <- risk.sb - mean(risk.sb)
  } else {
    g.sb <- list(acc = NA_real_, ari = NA_real_); lrr.sb <- NA
  }

  res <- data.frame(
    confound = confound, seed = seed, hotspot_zone = hotspot_zone, n = n,
    cens_rate = gen$cens_rate,
    # x-axis candidates
    cor_tpi_hot     = gen$cor_tpi_hot,                       # dialled, generative
    cor_tpi_risk_ts = .sc(lrr.ts, d$num_tpi),   # LeukSurv-matching
    cor_tpi_risk_sp = .sc(lrr.sp, d$num_tpi),
    cor_tpi_resid   = .sc(d$e_resid, d$num_tpi),
    within_zone_cor_north = gen$within_zone_cor_north,
    # GLOBAL three-zone recovery (primary)
    ari_ts = g.ts$ari, ari_sp = g.sp$ari,
    acc_ts = g.ts$acc, acc_sp = g.sp$acc,
    # HOTSPOT-local recovery (secondary)
    auc_hot_ts = .auc(risk.ts, is_hot), auc_hot_sp = .auc(risk.sp, is_hot),
    acc_hot_ts = .acc_hot(risk.ts, is_hot), acc_hot_sp = .acc_hot(risk.sp, is_hot),
    leaves_sp = tree.sp$leafCount, leaves_cl = tree.cl$leafCount,
    leaves_ts = tree.ts$leafCount,
    # spBayesSurv arm (NA when fit_sb = FALSE)
    ari_sb = g.sb$ari, acc_sb = g.sb$acc,
    auc_hot_sb = if (!is.null(risk.sb)) .auc(risk.sb, is_hot) else NA_real_,
    acc_hot_sb = if (!is.null(risk.sb)) .acc_hot(risk.sb, is_hot) else NA_real_,
    cor_tpi_risk_sb = if (!is.null(risk.sb)) .sc(lrr.sb, d$num_tpi) else NA_real_,
    fit_sb = fit_sb,
    stringsAsFactors = FALSE)

  if (verbose) {
    if (fit_sb)
      cat(sprintf("[%s] spBayes ARI=%.3f acc=%.3f hotAUC=%.3f cor(tpi,risk)=%+.3f\n",
                  tag, res$ari_sb, res$acc_sb, res$auc_hot_sb, res$cor_tpi_risk_sb))
    else
      cat(sprintf("[%s] global ARI ts=%.3f sp=%.3f | hot AUC ts=%.3f sp=%.3f | cor(tpi,risk_ts)=%+.3f\n",
                  tag, g.ts$ari, g.sp$ari, res$auc_hot_ts, res$auc_hot_sp, res$cor_tpi_risk_ts))
  }

  if (debug) {
    cat("per-zone mean risk, two-stage:\n");   print(round(tapply(risk.ts, d$zone, mean), 3))
    cat("per-zone mean risk, spatial-only:\n"); print(round(tapply(risk.sp, d$zone, mean), 3))
    cat("per-zone residual mean:\n");           print(round(tapply(d$e_resid, d$zone, mean), 3))
    return(invisible(list(gen = gen, d = d, res = res,
                          risk.sp = risk.sp, risk.ts = risk.ts,
                          pt.sp = pt.sp, pt.ts = pt.ts,
                          tree.sp = tree.sp, tree.cl = tree.cl, tree.ts = tree.ts)))
  }
  saveRDS(res, out_file)
  invisible(out_file)
}

collect_sweep <- function(out_dir = "sweep_out") {
  f <- list.files(out_dir, pattern = "^rep_.*\\.rds$", full.names = TRUE)
  if (!length(f)) { cat("no replicate files in", out_dir, "\n"); return(NULL) }
  res <- do.call(rbind, lapply(f, readRDS))
  cat(sprintf("collected %d rows: %d tree, %d spBayes | levels: %s\n",
              nrow(res), sum(!res$fit_sb), sum(res$fit_sb),
              paste(sort(unique(res$confound)), collapse = ", ")))
  n_tab <- table(res$confound, ifelse(res$fit_sb, "spBayes", "tree"))
  cat("replicates per level:\n"); print(n_tab)
  res
}

## Level means and standard errors, computed column by column over the
## rows where each column is present. Tree columns are NA in spBayes rows
## and vice versa, so na.rm is what keeps the two arms separate.
summarise_sweep <- function(res, cols = NULL) {
  if (is.null(cols))
    cols <- grep("^(ari|acc|auc_hot|acc_hot|cor_tpi)", names(res), value = TRUE)
  lv <- sort(unique(res$confound))
  out <- data.frame(confound = lv)
  for (cc in cols) {
    m <- sapply(lv, function(L) mean(res[[cc]][res$confound == L], na.rm = TRUE))
    s <- sapply(lv, function(L) {
      v <- res[[cc]][res$confound == L]; v <- v[!is.na(v)]
      if (length(v) < 2) NA_real_ else sd(v) / sqrt(length(v)) })
    nn <- sapply(lv, function(L) sum(!is.na(res[[cc]][res$confound == L])))
    out[[cc]]              <- m
    out[[paste0(cc, "_se")]] <- s
    out[[paste0(cc, "_n")]]  <- nn
  }
  out
}

# ---- cluster array decoder (optional) --------------------------------
# SEEDS <- 1:10
# idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# ci <- ((idx-1) %% length(CONF_GRID)) + 1
# si <- ((idx-1) %/% length(CONF_GRID)) + 1
# run_sweep_replicate(CONF_GRID[ci], SEEDS[si])
# array range: 1 .. length(CONF_GRID)*length(SEEDS)
