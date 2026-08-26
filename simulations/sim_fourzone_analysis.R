# =====================================================================
# FOUR-ZONE DESIGN: analysis, metrics, maps, residuals
# ---------------------------------------------------------------------
# Run after sim_fourzone_generate.R, which leaves in the environment:
#   sim_data, z_s, zone, grid, lp1, lp2, lp3, split1_thresh,
#   split2_thresh, split3_thresh, right_side, z_levels,
#   b_wbc, b_tpi
#
# Settings match the two- and three-zone reruns:
#   adaptive = TRUE, probs = 0.23, epsilon_floor_alpha = 0.2,
#   kappa = exp(-2) on all three trees, .60/.40 kernel mixture on both
#   spatial fits so spatial-only and two-stage differ ONLY in the
#   clinical adjustment.
#
# Metrics: accuracy, far-boundary accuracy, ARI, standardized MSE, span
# ratio. AUC is omitted: four classes do not define a single ROC curve.
#
# Hyperparameters are held fixed rather than tuned per design.
# =====================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.tree); library(survival); library(ggplot2)
library(patchwork); library(scales)
source("../method/dipole_tree.R")
source("../method/computeresiduals_NA.R")
source("../method/bootstrap_pruning.R")
source("../method/metrics.R")

stopifnot(exists("sim_data"), exists("zone"), exists("grid"),
          exists("split1_thresh"), exists("split2_thresh"),
          exists("split3_thresh"), exists("z_s"))

ZLEV    <- c("low", "medium_low", "medium_high", "high")
alldata <- sim_data
n       <- nrow(alldata)
truth_k <- as.integer(factor(alldata$zone, levels = ZLEV))
time <- "stop"; censor <- "status"
quantiles <- c(.25, .75); tolerance <- 1e-2; nsize <- 20
covariates_clinical <- c("num_age", "num_wbc", "num_tpi")
covariates_spatial  <- c("num_xcoord", "num_ycoord")

ADAPTIVE <- TRUE; EPS_PROBS <- 0.23; EPS_FLOOR <- 0.2
KAPPA    <- exp(-2)

# =====================================================================
# 1. Design-specific scoring inputs
# =====================================================================
# The scoring functions live in ../method/metrics.R. What is specific to
# this design is the far-boundary set, the standardized true frailty, and
# the true frailty span.

# Far-boundary set. Three boundaries: split 1 everywhere, split 2 on the
# right side, split 3 on the left. Each margin is scaled by the SD of its
# own linear predictor so the three are comparable; a patient's margin is
# the distance to whichever boundary is nearer.
m1 <- abs(lp1 - split1_thresh) / sd(lp1)
m2 <- abs(lp2 - split2_thresh) / sd(lp2); m2[!right_side] <- Inf
m3 <- abs(lp3 - split3_thresh) / sd(lp3); m3[ right_side] <- Inf
margin <- pmin(m1, m2, m3)
far <- margin > median(margin)
cat(sprintf("far set: n=%d | zone shares:", sum(far)))
print(round(prop.table(table(alldata$zone[far])), 3))

zc        <- scale(z_s)[, 1]
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
  epsilon = epsilon_sp, kappa = KAPPA, nsize = nsize,
  pureweight = 1, mixedweight = 1,
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
  epsilon = epsilon_cl, kappa = KAPPA, nsize = 10,
  pureweight = 1, mixedweight = 1,
  Kconstant = 0, Kpoly_order = 1,
  gaussweight = 0, polyweight = 0, linearweight = 1,
  linearcovariates_index = 1:3, gausscovariates_index = 1:3,
  polycovariates_index = 1:3, ncovariatestosearch = 3,
  adaptive = ADAPTIVE, epsilon_floor_alpha = EPS_FLOOR, probs = EPS_PROBS)

tree.cl <- Dip.cl$createtree(1:n)
alldata$e_resid <- compute_residuals(alldata, tree.cl, Dip.cl,
                                     covariates_clinical, time, censor)
cat("   clinical leaves =", tree.cl$leafCount, "\n")

Dip.ts <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  alldata, "e_resid", censor, covariates_spatial, quantiles, tolerance,
  epsilon = epsilon_sp, kappa = KAPPA, nsize = nsize,
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
dat_sb   <- sim_data; dat_sb$id <- 1:n
locs     <- as.matrix(dat_sb[, c("num_xcoord", "num_ycoord")])
m        <- 30
prior_sb <- list(maxL = 15, a0 = 1, b0 = 1, nknots = m, nblock = m, nu = 1)
mcmc_sb  <- list(nburn = 2000, nsave = 2000, nskip = 1, ndisplay = 1000)
cor.dist <- function(x1, x2) rdist(x1, x2)

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
  if (is.null(vpost)) { cat("  frailty samples not found\n"); print(names(fit.sb)) }
  else risk.sb <- vpost
}

# =====================================================================
# 5. Metrics
# =====================================================================
cat("\n===== Zone recovery (4-means, clusters ordered by mean risk) =====\n")
s.sp <- score_multiclass(risk.sp, lrr(pt.sp), "spatial-only",
                         truth_k, far, zc, span_true, k = 4)
s.ts <- score_multiclass(risk.ts, lrr(pt.ts), "two-stage",
                         truth_k, far, zc, span_true, k = 4)
if (!is.null(risk.sb))
  s.sb <- score_multiclass(risk.sb, risk.sb - mean(risk.sb), "spBayes PH+GRF",
                           truth_k, far, zc, span_true, k = 4)

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
bl <- list(
  geom_contour(data = grid, aes(num_xcoord, num_ycoord, z = lp1),
               breaks = split1_thresh, color = "black", linewidth = 0.8),
  geom_contour(data = grid[!is.na(grid$lp2_masked), ],
               aes(num_xcoord, num_ycoord, z = lp2_masked),
               breaks = split2_thresh, color = "black", linewidth = 0.8),
  geom_contour(data = grid[!is.na(grid$lp3_masked), ],
               aes(num_xcoord, num_ycoord, z = lp3_masked),
               breaks = split3_thresh, color = "black", linewidth = 0.8))

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
ggsave("fig_fourzone_2x2_sharedscale.pdf", fig, width = 11, height = 9)
ggsave("fig_fourzone_2x2_sharedscale.png", fig, width = 11, height = 9, dpi = 300)
cat("\nFigure written: fig_fourzone_2x2_sharedscale.pdf\n")

# =====================================================================
# 7. Residual diagnostics
# =====================================================================
d <- alldata
d$zone <- factor(d$zone, levels = ZLEV)
cat("\n===== Residual diagnostics =====\n")
cat(sprintf("%-13s %8s %8s %12s %12s\n",
            "zone", "mean", "median", "cor(north)", "R2(wbc,tpi)"))
for (zz in ZLEV) {
  sub <- d[d$zone == zz, ]
  r2  <- summary(lm(e_resid ~ num_wbc + num_tpi, data = sub))$r.squared
  cat(sprintf("%-13s %8.3f %8.3f %12.3f %12.3f\n",
              zz, mean(sub$e_resid), median(sub$e_resid),
              cor(sub$e_resid, sub$num_ycoord), r2))
}
cat("\nAdjacent-pair Mann-Whitney (higher frailty -> smaller residuals):\n")
for (k in 1:3) {
  pr  <- c(ZLEV[k+1], ZLEV[k])
  sub <- d[d$zone %in% pr, ]; sub$zone <- factor(sub$zone, levels = pr)
  pv  <- wilcox.test(e_resid ~ zone, data = sub, alternative = "less")$p.value
  cat(sprintf("  %-12s < %-12s : p = %.3g\n", pr[1], pr[2], pv))
}
r2_pool <- summary(lm(e_resid ~ zone + num_wbc + num_tpi, data = d))$r.squared -
           summary(lm(e_resid ~ zone, data = d))$r.squared
cat(sprintf("Clinical covariates explain %.2f%% of residual variance beyond zone\n",
            100 * r2_pool))

# tpi versus wbc share of the within-zone clinical R2, per zone
cat("\nWithin-zone R2 split (tpi alone vs wbc alone):\n")
for (zz in ZLEV) {
  sub <- d[d$zone == zz, ]
  r2t <- summary(lm(e_resid ~ num_tpi, data = sub))$r.squared
  r2w <- summary(lm(e_resid ~ num_wbc, data = sub))$r.squared
  cat(sprintf("  %-13s tpi %.3f   wbc %.3f\n", zz, r2t, r2w))
}

# =====================================================================
# 8. Residuals vs tpi by frailty zone
# =====================================================================
p_rtpi <- ggplot(d, aes(num_tpi, e_resid, color = zone)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(low = "blue", medium_low = "cyan3",
                                medium_high = "orange", high = "red")) +
  labs(x = "tpi", y = "NA residual", title = NULL) +
  theme_minimal()
ggsave("fig_resid_vs_tpi_fourzone.pdf", p_rtpi, width = 5.5, height = 4)
cat("Figure written: fig_resid_vs_tpi_fourzone.pdf\n")

# =====================================================================
# 9. Save
# =====================================================================
objs <- c("sim_data", "alldata", "z_s", "zone", "grid", "lp1", "lp2", "lp3",
          "split1_thresh", "split2_thresh", "split3_thresh", "right_side",
          "b_wbc", "b_tpi", "risk.sp", "risk.ts", "pt.sp", "pt.ts",
          "s.sp", "s.ts", "span_true", "margin", "far", "p_rtpi")
if (exists("risk.sb") && !is.null(risk.sb)) objs <- c(objs, "risk.sb")
if (exists("s.sb")) objs <- c(objs, "s.sb")
save(list = objs, file = "sim_fourzone_key_objects.RData")
cat("Saved:", paste(objs, collapse = ", "), "\n")
