# =====================================================================
# DEGENERATE DESIGN: analysis, metrics, and maps
# ---------------------------------------------------------------------
# Run after sim_degenerate_generate.R, which leaves in the environment:
#   sim_data, z_s, in_high, boundary, g_of_s, grid, g_patient,
#   b_wbc, b_tpi
#
# Fits three analyses and scores all metrics reported in the paper:
#   accuracy, far-boundary accuracy, ARI, standardized MSE, span ratio, AUC
#
# Scoring comes from ../method/metrics.R so every design is scored
# identically. Far-boundary performance is an ACCURACY, matching Table 2.
# Then builds the 2x2 shared-scale recovery figure.
#
# =====================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.tree); library(survival); library(ggplot2)
library(patchwork); library(scales)
source("../method/dipole_tree.R")
source("../method/computeresiduals_NA.R")
source("../method/bootstrap_pruning.R")
source("../method/metrics.R")

stopifnot(exists("sim_data"), exists("in_high"), exists("grid"),
          exists("boundary"), exists("z_s"), exists("g_patient"))

alldata    <- sim_data
alldata$zone <- ifelse(in_high, "high", "low")
truth_high <- in_high
n          <- nrow(alldata)
time <- "stop"; censor <- "status"
quantiles <- c(.25, .75); tolerance <- 1e-2; nsize <- 20
covariates_clinical <- c("num_age", "num_wbc", "num_tpi")
covariates_spatial  <- c("num_xcoord", "num_ycoord")

# adaptive-margin settings, matching the LeukSurv application
ADAPTIVE <- TRUE; EPS_PROBS <- 0.23; EPS_FLOOR <- 0.2

# =====================================================================
# 1. Design-specific scoring inputs
# =====================================================================
# The scoring functions live in ../method/metrics.R. What is specific to
# this design is the standardized true frailty, the far-boundary set, and
# the true frailty span.

zc <- scale(z_s)[, 1]

# Far-boundary set: the interior half, by margin in g. One boundary here,
# so the margin is just the distance to it.
gap <- abs(g_patient - boundary)
far <- gap > median(gap)

span_true <- span(log(z_s / median(z_s)))

# =====================================================================
# 2. Spatial-only tree
# =====================================================================
cat("\nFitting spatial-only ...\n")
X_sp <- alldata[, covariates_spatial]; d_sp <- c(dist(X_sp))
epsilon_sp <- as.numeric(quantile(d_sp[d_sp != 0], EPS_PROBS))
ksigma_sp  <- intvar(alldata, covariates_spatial)

Dip.sp <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  alldata, time, censor, covariates_spatial, quantiles, tolerance,
  epsilon = epsilon_sp, kappa = exp(-2), nsize = nsize,
  pureweight = 1, mixedweight = 1,metric="kernel",
  Ksigma = ksigma_sp, Kconstant = 1, Kpoly_order = 2,
  gaussweight = .60, polyweight = .40, linearweight = 0,
  gausscovariates_index = 1:2, polycovariates_index = 1:2,
  linearcovariates_index = 1:2, ncovariatestosearch = 2,
  adaptive = ADAPTIVE, epsilon_floor_alpha = EPS_FLOOR, probs = EPS_PROBS)

tree.sp <- Dip.sp$createtree(1:n)
tree.sp <- bootstrapPruning(tree.sp, Dip.sp, 3.7, time = time, censor = censor)[[4]]
pt.sp   <- Dip.sp$predicttime(alldata, tree.sp)
risk.sp <- -pt.sp
cat("   spatial-only leaves =", tree.sp$leafCount, "\n")

# =====================================================================
# 3. Two-stage
# =====================================================================
cat("Fitting two-stage ...\n")
X_cl <- alldata[, covariates_clinical]; d_cl <- c(dist(X_cl))
epsilon_cl <- as.numeric(quantile(d_cl[d_cl != 0], EPS_PROBS))

Dip.cl <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  alldata, time, censor, covariates_clinical, quantiles, tolerance,
  epsilon = epsilon_cl, kappa = exp(-2), nsize = 10,
  pureweight = 1, mixedweight = 1,
  Kconstant = 0, Kpoly_order = 1,
  gaussweight = 0, polyweight = 0, linearweight = 1,
  linearcovariates_index = 1:3, gausscovariates_index = 1:3,
  polycovariates_index = 1:3, ncovariatestosearch = 3, metric="kernel",
  adaptive = ADAPTIVE, epsilon_floor_alpha = EPS_FLOOR, probs = EPS_PROBS)

tree.cl <- Dip.cl$createtree(1:n)
alldata$e_resid <- compute_residuals(alldata, tree.cl, Dip.cl,
                                     covariates_clinical, time, censor)
cat("   clinical leaves =", tree.cl$leafCount, "\n")

Dip.ts <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  alldata, "e_resid", censor, covariates_spatial, quantiles, tolerance,
  metric="kernel",epsilon = epsilon_sp, kappa = exp(-2), nsize = nsize,
  pureweight = 1, mixedweight = 1,
  Ksigma = ksigma_sp, Kconstant = 1, Kpoly_order = 2,
  gaussweight = .60, polyweight = .40, linearweight = 0,
  gausscovariates_index = 1:2, polycovariates_index = 1:2,
  linearcovariates_index = 1:2, ncovariatestosearch = 2,
  adaptive = ADAPTIVE, epsilon_floor_alpha = EPS_FLOOR, probs = EPS_PROBS)

tree.ts       <- Dip.ts$createtree(1:n)
tree.ts.prune <- bootstrapPruning(tree.ts, Dip.ts, 3.7,
                                  time = "e_resid", censor = censor)[[4]]
pt.ts   <- Dip.ts$predicttime(alldata, tree.ts.prune)
risk.ts <- -pt.ts
cat("   spatial leaves =", tree.ts.prune$leafCount, "\n")

# =====================================================================
# 4. spBayesSurv PH + GRF
# =====================================================================
library(spBayesSurv); library(fields)
dat_sb    <- sim_data; dat_sb$id <- 1:n
locs      <- as.matrix(dat_sb[, c("num_xcoord", "num_ycoord")])
m         <- 30
prior_sb  <- list(maxL = 15, a0 = 1, b0 = 1, nknots = m, nblock = m, nu = 1)
mcmc_sb   <- list(nburn = 2000, nsave = 2000, nskip = 1, ndisplay = 1000)
cor.dist  <- function(x1, x2) rdist(x1, x2)

cat("Fitting spBayesSurv PH+GRF (slow) ...\n")
set.seed(1)
fit.sb <- tryCatch(
  survregbayes(Surv(stop, status) ~ num_age + num_wbc + num_tpi +
                 frailtyprior("grf", id),
               data = dat_sb, survmodel = "PH", dist = "loglogistic",
               mcmc = mcmc_sb, prior = prior_sb, DIST = cor.dist,
               Coordinates = locs, InitParamMCMC = FALSE),
  error = function(e) { cat("  spBayesSurv error:", conditionMessage(e), "\n"); NULL })

risk.sb <- NULL
if (!is.null(fit.sb)) {
  vpost <- if (!is.null(fit.sb$v)) {
    if (is.matrix(fit.sb$v)) rowMeans(fit.sb$v) else fit.sb$v
  } else if (!is.null(fit.sb$frail)) {
    if (is.matrix(fit.sb$frail)) rowMeans(fit.sb$frail) else fit.sb$frail
  } else NULL
  if (is.null(vpost)) { cat("  frailty samples not found; names(fit.sb):\n"); print(names(fit.sb)) }
  else risk.sb <- vpost
}

# =====================================================================
# 5. Score everything
# =====================================================================
cat("\n===== Zone recovery (2-means threshold) =====\n")
s.sp <- score_binary(risk.sp, lrr(pt.sp), "spatial-only",
                     truth_high, far, zc, span_true)
s.ts <- score_binary(risk.ts, lrr(pt.ts), "two-stage",
                     truth_high, far, zc, span_true)
if (!is.null(risk.sb))
  s.sb <- score_binary(risk.sb, risk.sb - mean(risk.sb), "spBayes PH+GRF",
                       truth_high, far, zc, span_true)

cat("\n============ TABLE ROWS (paper order) ============\n")
table_header()
table_row(s.sp, "Spatial-only"); table_row(s.ts, "Two-stage")
if (exists("s.sb")) table_row(s.sb, "spBayesSurv PH+GRF")
cat("\ntrue frailty log-risk span (5th-95th) =", round(span_true, 3), "\n")

# =====================================================================
# 6. 2x2 recovery figure, shared relative-risk scale
# =====================================================================
cap      <- function(x, lo = 0.3, hi = 3) pmin(pmax(x, lo), hi)
rr_scale <- scale_color_gradient2(low = "blue", mid = "gray70", high = "red",
                                  midpoint = 1, name = "Relative\nRisk",
                                  limits = c(0.3, 3), oob = squish)
bt <- list(theme_minimal(), coord_fixed(xlim = c(0, 1), ylim = c(0, 1)),
           labs(x = "Easting", y = "Northing"))
bl <- geom_contour(data = grid, aes(num_xcoord, num_ycoord, z = g),
                   breaks = boundary, color = "black", linewidth = 0.8)

md  <- sim_data
eta <- b_wbc * md$num_wbc + b_tpi * md$num_tpi
md$rrA <- cap({ rr <- z_s * exp(eta); rr / median(rr) })
md$rrB <- cap(median(pt.sp) / pt.sp)
md$rrC <- cap(median(pt.ts) / pt.ts)

panel <- function(col, ttl)
  ggplot(md, aes(num_xcoord, num_ycoord, color = .data[[col]])) +
  geom_point(size = 1.5, alpha = .85) + rr_scale + bl + labs(title = ttl) + bt

pA <- panel("rrA", "(A) Simulated total risk")
pB <- panel("rrB", "(B) Spatial-only")
pC <- panel("rrC", "(C) Two-stage")
if (!is.null(risk.sb)) {
  md$rrD <- cap({ r <- exp(risk.sb); r / median(r) })
  pD <- panel("rrD", "(D) spBayesSurv PH+GRF")
} else {
  pD <- ggplot() + theme_void() + labs(title = "(D) spBayesSurv PH+GRF: not fitted")
}

fig <- (pA | pB) / (pC | pD) + plot_layout(guides = "collect")
ggsave("fig_degenerate_2x2_sharedscale.pdf", fig, width = 11, height = 9)
ggsave("fig_degenerate_2x2_sharedscale.png", fig, width = 11, height = 9, dpi = 300)
cat("\nFigure written: fig_degenerate_2x2_sharedscale.pdf\n")

# =====================================================================
# 7. Save
# =====================================================================
objs <- c("sim_data", "z_s", "in_high", "boundary", "g_of_s", "grid",
          "g_patient", "b_wbc", "b_tpi", "risk.sp", "risk.ts",
          "pt.sp", "pt.ts", "s.sp", "s.ts", "span_true")
if (exists("risk.sb") && !is.null(risk.sb)) objs <- c(objs, "risk.sb")
if (exists("s.sb")) objs <- c(objs, "s.sb")
save(list = objs, file = "sim_degenerate_key_objects.RData")
cat("Saved:", paste(objs, collapse = ", "), "\n")
