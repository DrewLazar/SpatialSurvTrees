# ============================================================================
# leuksurv_realdata.R
#
# Complete real-data analysis for the spatial survival tree paper, in one script.
# Sections run top to bottom; the two slow steps (spBayesSurv fit, tree fits)
# save their results so you can comment them out and reload on a rerun.
#
# SECTIONS
#   0. Setup: paths, libraries, tree class + helpers, shared plotting scale
#   1. Data prep: load + clean LeukSurv, descriptive stats, overall KM curve
#   2. spBayesSurv PH+GRF benchmark fit            (SLOW ~2h; skip if saved)
#   3. spBayesSurv summary: Table 1 numbers + in-sample C-index
#   4. Two-stage kernel dipole tree: clinical -> NA residuals -> spatial -> map
#   5. Spatial-only baseline: coordinates-only tree on raw times
#   6. Figures: three-panel comparison (smooth | two-stage over spatial-only)
# ============================================================================


# ============================================================================
# 0. SETUP
# ============================================================================
# Paths are relative to application/. Set the working directory there
# before sourcing, or open the file in RStudio and use Session > Set
# Working Directory > To Source File Location.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

CODE_DIR    <- "../method"
DATA_FILE   <- "../data/LeukSurv.csv"
OUT_DIR     <- "outputs"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(scales); library(RColorBrewer)
  library(spBayesSurv); library(coda); library(patchwork)
})

source(file.path(CODE_DIR, "dipole_tree.R"))
source(file.path(CODE_DIR, "computeresiduals_NA.R"))   # Nelson-Aalen residuals
source(file.path(CODE_DIR, "bootstrap_pruning.R"))
source(file.path(CODE_DIR, "kappacrossvalidation.R"))  # supplies intvarfun()

set.seed(22)                                           # reproducibility
`%||%` <- function(a, b) if (is.null(a)) b else a

# shared diverging RdBu scale (identical across all maps)
RISK_LIMITS <- c(-0.30, 0.30)
risk_colors <- rev(brewer.pal(11, "RdBu"))             # blue (low) -> light -> red (high)
diverging_fill <- function(name)
  scale_fill_gradientn(colors = risk_colors, limits = RISK_LIMITS,
                       oob = scales::squish, na.value = "white", name = name)
white_bg <- theme(panel.background = element_rect(fill = "white", color = NA),
                  plot.background  = element_rect(fill = "white", color = NA))
risk_point_map <- function(df, fill_col, legend_title)
  ggplot(df, aes(x = num_xcoord, y = num_ycoord, fill = .data[[fill_col]])) +
  geom_point(shape = 21, color = "grey40", stroke = 0.25, size = 2.0, alpha = 0.95) +
  diverging_fill(legend_title) +
  labs(x = "Easting", y = "Northing") +
  theme_minimal(base_size = 12) + white_bg + coord_fixed()


# ============================================================================
# 1. DATA PREP + DESCRIPTIVE STATS + OVERALL KM
# ============================================================================
LeukSurv <- read.csv(DATA_FILE)
LeukSurv$fac_sex <- as.integer(as.character(LeukSurv$fac_sex) == "M")
names(LeukSurv)[names(LeukSurv) == "time"]  <- "stop"
names(LeukSurv)[names(LeukSurv) == "event"] <- "status"    # confirm censor col name
alldata <- na.omit(LeukSurv)
n <- nrow(alldata)
time <- "stop"; censor <- "status"
quantiles <- c(.25, .75); tolerance <- 1e-2

events <- sum(alldata$status == 1)
cat("=== LeukSurv descriptive statistics ===\n")
cat("n              :", n, "\n")
cat("events (deaths):", events, sprintf("(%.1f%%)\n", 100*events/n))
cat("censored       :", n - events, sprintf("(%.1f%%)\n", 100*(n-events)/n))
cat("time range     :", paste(range(alldata$stop), collapse = " to "), "days\n")
cat("mean time      :", round(mean(alldata$stop), 1), "days\n")
cat("median (raw)   :", median(alldata$stop), "days\n\n")

km <- survfit(Surv(stop, status) ~ 1, data = alldata)
cat("=== Overall Kaplan-Meier ===\n"); print(km)
cat("\nKM median table:\n"); print(summary(km)$table)

p_km <- survminer::ggsurvplot(
  km, data = alldata, conf.int = TRUE, censor = FALSE,
  xlab = "Time (days)", ylab = "Survival probability",
  ggtheme = theme_minimal(base_size = 12), palette = "#2166AC")$plot
ggsave(file.path(OUT_DIR, "leuksurv_km_overall.pdf"), p_km, width = 7, height = 5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "leuksurv_km_overall.png"), p_km, width = 7, height = 5, dpi = 300)


# ============================================================================
# 2. spBayesSurv PH + GRF BENCHMARK FIT   (SLOW ~2h)
#    Comment out this block and load the saved fit on reruns.
# ============================================================================
FIT_FILE <- file.path(OUT_DIR, "fit_ph_leuksurv.rds")
if (!file.exists(FIT_FILE)) {
  data(LeukSurv)                                        # packaged version for the fit
  dbench <- LeukSurv; dbench$id <- 1:nrow(dbench)
  coords_b <- as.matrix(dbench[, c("xcoord", "ycoord")])
  cat("\n=== Fitting spBayesSurv PH + GRF (~2 hours) ===\n")
  t0 <- proc.time()
  fit_ph <- survregbayes(
    formula = Surv(time, cens) ~ age + sex + wbc + tpi + frailtyprior("grf", id),
    data = dbench, survmodel = "PH",
    prior = list(nknots = 150, nblock = 1043),
    mcmc  = list(nburn = 4000, nsave = 4000, nskip = 1, ndisplay = 500),
    Coordinates = coords_b, InitParamMCMC = TRUE)
  cat("Elapsed:", round((proc.time() - t0)["elapsed"]/60, 1), "min\n")
  saveRDS(fit_ph, FIT_FILE)
} else {
  cat("\nLoading saved spBayesSurv fit (skipping the 2h refit).\n")
  fit_ph <- readRDS(FIT_FILE)
}


# ============================================================================
# 3. spBayesSurv PH WITHOUT the GRF frailty  -- comparison fit  (FAST)
#    LPML / DIC / WAIC are interpretable only as DIFFERENCES between models,
#    so the criteria reported for the PH + GRF fit need a referent. This fits
#    the same model minus only the spatial term. It also gives boundary-aware
#    evidence that the GRF earns its place, which a credible interval for a
#    variance cannot: a variance is bounded below by zero, so its interval can
#    never cover zero however weak the spatial signal.
#
#    Identical to block 2 except: frailtyprior("grf", id) dropped from the
#    formula, Coordinates dropped, and nknots / nblock dropped (both configure
#    the full-scale approximation of the spatial process). Data, PH link,
#    Bernstein baseline defaults and MCMC length are held fixed, so the two
#    fits differ only by the GRF.
# ============================================================================
FIT_FILE_NF <- file.path(OUT_DIR, "fit_ph_leuksurv_nofrailty.rds")
if (!file.exists(FIT_FILE_NF)) {
  data(LeukSurv)
  dbench <- LeukSurv; dbench$id <- 1:nrow(dbench)
  cat("\n=== Fitting spBayesSurv PH, no frailty (comparison model) ===\n")
  t0 <- proc.time()
  fit_ph_nf <- survregbayes(
    formula = Surv(time, cens) ~ age + sex + wbc + tpi,
    data = dbench, survmodel = "PH",
    mcmc  = list(nburn = 4000, nsave = 4000, nskip = 1, ndisplay = 500),
    InitParamMCMC = TRUE)
  cat("Elapsed:", round((proc.time() - t0)["elapsed"]/60, 1), "min\n")
  saveRDS(fit_ph_nf, FIT_FILE_NF)
} else {
  cat("\nLoading saved no-frailty fit.\n")
  fit_ph_nf <- readRDS(FIT_FILE_NF)
}

# ---- the three criteria from both fits -------------------------------------
grab <- function(fit) {
  out <- c(LPML = NA_real_, DIC = NA_real_, WAIC = NA_real_)
  for (nm in names(out)) if (!is.null(fit[[nm]])) out[nm] <- as.numeric(fit[[nm]])[1]
  if (all(is.na(out))) {                      # field names vary across versions
    sm <- summary(fit)
    for (nm in names(out)) if (!is.null(sm[[nm]])) out[nm] <- as.numeric(sm[[nm]])[1]
  }
  out
}
crit <- as.data.frame(rbind(`PH + GRF`       = grab(fit_ph),
                            `PH, no frailty` = grab(fit_ph_nf)))
crit$dLPML <- crit$LPML - crit$LPML[2]      # + favours the spatial model
crit$dDIC  <- crit$DIC  - crit$DIC[2]       # - favours the spatial model
crit$dWAIC <- crit$WAIC - crit$WAIC[2]      # - favours the spatial model

cat("\n=== Model comparison: does the GRF earn its place? ===\n")
print(round(crit, 1))
cat("\nLPML higher is better; DIC and WAIC lower is better.\n")
cat("Differences of roughly 5 or more are usually taken as meaningful;\n")
cat("differences of 2 or less are within noise.\n")

# ---- do the clinical effects move when the GRF is dropped? ------------------
# If they shift much, the spatial term was absorbing clinical signal, which is
# directly relevant to the entanglement argument in the paper.
cat("\n=== Clinical coefficients, with and without the GRF ===\n")
cf <- function(fit) { b <- fit$beta; if (is.matrix(b)) rowMeans(b) else as.numeric(b) }
print(round(cbind(`PH + GRF` = cf(fit_ph), `no frailty` = cf(fit_ph_nf)), 4))


# ============================================================================
# 4. spBayesSurv SUMMARY: Table 1 numbers + in-sample C-index
# ============================================================================
data(LeukSurv); dbench <- LeukSurv; nb <- nrow(dbench)
cat("\n=== spBayesSurv PH + GRF posterior summary ===\n")
print(summary(fit_ph))

get_betamean <- function(b) if (is.matrix(b)) rowMeans(b) else as.numeric(b)
beta_scaled <- if (!is.null(fit_ph$beta.scaled)) get_betamean(fit_ph$beta.scaled) else NULL
beta_raw    <- if (!is.null(fit_ph$beta))        get_betamean(fit_ph$beta)        else NULL
X_scaled    <- if (!is.null(fit_ph$X.scaled)) as.matrix(fit_ph$X.scaled) else NULL
X_raw       <- if (!is.null(fit_ph$X))        as.matrix(fit_ph$X)        else NULL
lin_pred <- function(Xmat, beta) {
  if (is.null(Xmat) || is.null(beta)) return(NULL)
  if (ncol(Xmat) == length(beta)) return(as.numeric(Xmat %*% beta))
  if (ncol(Xmat) == length(beta) + 1) return(as.numeric(Xmat[, -1, drop = FALSE] %*% beta))
  if (ncol(Xmat) + 1 == length(beta)) return(as.numeric(Xmat %*% beta[-1]))
  NULL
}
xb <- lin_pred(X_scaled, beta_scaled) %||% lin_pred(X_raw, beta_raw)
vb <- fit_ph$v
vhat_b <- if (is.matrix(vb)) (if (nrow(vb) == nb) rowMeans(vb) else colMeans(vb)) else as.numeric(vb)
surv_b <- Surv(dbench$time, dbench$cens)
c_full <- concordance(surv_b ~ I(xb + vhat_b), reverse = TRUE)$concordance
c_clin <- concordance(surv_b ~ xb,             reverse = TRUE)$concordance
c_spat <- concordance(surv_b ~ vhat_b,         reverse = TRUE)$concordance
cat("\n=== In-sample Harrell C-index (spBayesSurv) ===\n")
cat(sprintf("  clinical + spatial frailty : %.4f\n", c_full))
cat(sprintf("  clinical covariates only   : %.4f\n", c_clin))
cat(sprintf("  spatial frailty only       : %.4f\n", c_spat))
cat(sprintf("  spatial contribution       : %+.4f\n", c_full - c_clin))


# ============================================================================
# 5. TWO-STAGE KERNEL DIPOLE TREE
# ============================================================================
PRUNE_CLINICAL <- FALSE
GW <- 0.45; PW <- 0.55                                  # spatial mixture weights
ALPHA_SPATIAL  <- 1.8

# --- clinical stage (linear kernel) ---
covariates_clinical <- c("num_age", "fac_sex", "num_wbc", "num_tpi")
distX_clin   <- c(dist(alldata[, covariates_clinical]))
epsilon_clin <- quantile(distX_clin[distX_clin != 0], probs = 0.23)
Dipolar.clinical <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  traindata = alldata, time = time, censor = censor,
  covariates = covariates_clinical,
  quantiles = quantiles, tolerance = tolerance,
  epsilon = epsilon_clin, kappa = exp(1), nsize = 10,
  pureweight = 1, mixedweight = 1,
  Ksigma = 1, Kconstant = 0, Kpoly_order = 1,
  ncovariatestosearch = 4,
  gaussweight = 0, polyweight = 0, linearweight = 1,
  gausscovariates_index = c(1,2,3,4), polycovariates_index = c(1,2,3,4),
  linearcovariates_index = c(1,2,3,4),
  probs = 0.23, adaptive = TRUE, epsilon_floor_alpha = 0.2)
clinicaltree <- Dipolar.clinical$createtree(1:n)
if (PRUNE_CLINICAL) clinicaltree <- bootstrapPruning(clinicaltree, Dipolar.clinical, 2.2)[[4]]
cat("\nClinical tree leaves:", clinicaltree$leafCount, "\n")

alldata$e_resid <- compute_residuals(alldata, clinicaltree, Dipolar.clinical,
                                     covariates_clinical, time, censor)
cat("Residual summary:\n"); print(summary(alldata$e_resid))

# --- spatial stage (Gaussian + polynomial mixture) ---
covariates_spatial <- c("num_xcoord", "num_ycoord")
distX_sp   <- c(dist(alldata[, covariates_spatial]))
epsilon_sp <- quantile(distX_sp[distX_sp != 0], probs = 0.23)
ksigma_sp  <- intvarfun(alldata, covariates_spatial)
Dipolar.spatial <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  traindata = alldata, time = "e_resid", censor = "status",
  covariates = covariates_spatial,
  quantiles = quantiles, tolerance = tolerance,
  epsilon = epsilon_sp, kappa = exp(1), nsize = 10, metric = "kernel",
  pureweight = 1, mixedweight = 1,
  Ksigma = ksigma_sp, Kconstant = 1, Kpoly_order = 2,
  ncovariatestosearch = 2,
  gaussweight = GW, polyweight = PW, linearweight = 0,
  gausscovariates_index = c(1,2), polycovariates_index = c(1,2),
  linearcovariates_index = c(1,2),
  probs = 0.23, adaptive = TRUE, epsilon_floor_alpha = 0.2)
spatialtree <- Dipolar.spatial$createtree(1:n)
predictedresid <- Dipolar.spatial$predicttime(alldata, spatialtree)
cat("Spatial tree (unpruned) CI =",
    round(Dipolar.spatial$cindex(alldata$e_resid, predictedresid, alldata$status), 4),
    " leaves =", spatialtree$leafCount, "\n")

spatialtree.prune <- bootstrapPruning(spatialtree, Dipolar.spatial, ALPHA_SPATIAL)[[4]]
predictedresid.prune <- Dipolar.spatial$predicttime(alldata, spatialtree.prune)
CI_sp.prune <- Dipolar.spatial$cindex(alldata$e_resid, predictedresid.prune, alldata$status)
cat("Spatial tree (pruned)   CI =", round(CI_sp.prune, 4),
    " leaves =", spatialtree.prune$leafCount, "\n")

alldata$spatial_pred    <- predictedresid.prune
resid_median            <- median(alldata$e_resid)
alldata$spatial_rr      <- resid_median / alldata$spatial_pred
alldata$centered_log_rr <- log(alldata$spatial_rr) - mean(log(alldata$spatial_rr))

p_tree <- risk_point_map(alldata, "centered_log_rr", "Log relative risk\n(centered)")
ggsave(file.path(OUT_DIR, "leuksurv_twostage_riskmap.pdf"), p_tree, width = 7, height = 6, device = cairo_pdf)


# ============================================================================
# 6. SPATIAL-ONLY BASELINE (coordinates-only tree on RAW survival times)
# ============================================================================
Dipolar.spatialonly <- DipolarSurvivalTree_MixedGaussPolyLinearKernel_recpom$new(
  traindata = alldata, time = "stop", censor = "status",
  covariates = covariates_spatial,
  quantiles = quantiles, tolerance = tolerance,
  epsilon = epsilon_sp, kappa = exp(1), nsize = 10, metric = "kernel",
  pureweight = 1, mixedweight = 1,
  Ksigma = ksigma_sp, Kconstant = 1, Kpoly_order = 2,
  ncovariatestosearch = 2,
  gaussweight = GW, polyweight = PW, linearweight = 0,
  gausscovariates_index = c(1,2), polycovariates_index = c(1,2),
  linearcovariates_index = c(1,2),
  probs = 0.23, adaptive = TRUE, epsilon_floor_alpha = 0.2)
sptree <- Dipolar.spatialonly$createtree(1:n)
sptree.prune <- bootstrapPruning(sptree, Dipolar.spatialonly, ALPHA_SPATIAL)[[4]]
pred_raw.prune <- Dipolar.spatialonly$predicttime(alldata, sptree.prune)
CI_raw.prune <- Dipolar.spatialonly$cindex(alldata$stop, pred_raw.prune, alldata$status)
cat("\nSpatial-only tree (pruned) CI =", round(CI_raw.prune, 4),
    " leaves =", sptree.prune$leafCount, "\n")

alldata$so_pred         <- pred_raw.prune
raw_median              <- median(alldata$stop)
alldata$so_rr           <- raw_median / alldata$so_pred
alldata$so_centered_lrr <- log(alldata$so_rr) - mean(log(alldata$so_rr))

p_so <- risk_point_map(alldata, "so_centered_lrr", "Log relative risk\n(centered)")
ggsave(file.path(OUT_DIR, "leuksurv_spatialonly_riskmap.pdf"), p_so, width = 7, height = 6, device = cairo_pdf)

# --- confounding-mechanism check (numbers reported in the paper) ------------
# Compare clinical covariates between the region the spatial-only method calls
# high risk vs the rest. Higher wbc / deprivation in the "high" region is the
# evidence that spatial-only is attributing clinical burden to location.
alldata$so_high <- alldata$so_centered_lrr > 0
cat("\n=== Clinical covariates by spatial-only risk region ===\n")
print(aggregate(cbind(num_age, num_wbc, num_tpi) ~ so_high, data = alldata, mean))


# ============================================================================
# 7. FIGURES: three-panel comparison (final design)
#    row 1: (a) smooth  (b) two-stage   [these AGREE]
#    row 2:        (c) spatial-only      [this DIVERGES]
#    one shared bottom legend (common scale), descriptive subtitles.
# ============================================================================
data(LeukSurv); df <- LeukSurv; coords <- as.matrix(df[, c("xcoord","ycoord")])
vv <- fit_ph$v
vhat <- if (is.matrix(vv)) (if (nrow(vv) == nrow(df)) rowMeans(vv) else colMeans(vv)) else as.numeric(vv)
spatial_risk <- vhat - mean(vhat)
cat("\nSmooth surface 2/98 pct:", round(quantile(spatial_risk, c(.02,.98)), 3),
    " (limits", paste(RISK_LIMITS, collapse=","), ")\n")
ii <- akima::interp(coords[,1], coords[,2], spatial_risk, duplicate = "mean", nx = 200, ny = 200)
surf <- data.frame(expand.grid(easting = ii$x, northing = ii$y), risk = as.vector(ii$z))
p_smooth <- ggplot(surf, aes(easting, northing, fill = risk)) +
  geom_raster(interpolate = TRUE) +
  diverging_fill("Posterior log-frailty\n(centered)") +
  labs(x = "Easting", y = "Northing") +
  theme_minimal(base_size = 12) + white_bg + coord_fixed()

# shared legend title (all three panels share one common scale)
LEG_TITLE <- "Centered log risk"
relabel <- function(p, tag, subtitle)
  p + labs(tag = tag, subtitle = subtitle) +
  guides(fill = guide_colourbar(
    title = LEG_TITLE, title.position = "top", title.hjust = 0.5,
    barwidth = unit(9, "cm"), barheight = unit(0.4, "cm"))) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.title  = element_text(size = 10),
        plot.margin   = margin(2, 2, 2, 2))

pa <- relabel(p_smooth, "(a)", "spBayesSurv (smooth)")
pb <- relabel(p_tree,   "(b)", "Two-stage tree")
pc <- relabel(p_so,     "(c)", "Spatial-only")

top_row    <- pa + pb + plot_layout(widths = c(1, 1))
bottom_row <- plot_spacer() + pc + plot_spacer() + plot_layout(widths = c(0.35, 1, 0.35))
three_panel <- (top_row / bottom_row) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"),
        panel.spacing = unit(0.4, "lines"))

ggsave(file.path(OUT_DIR, "leuksurv_comparison_3panel.pdf"),
       three_panel, width = 13, height = 12, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "leuksurv_comparison_3panel.png"),
       three_panel, width = 13, height = 12, dpi = 300)

cat("\nDone. Two-stage residual C =", round(CI_sp.prune, 4),
    "; spatial-only C =", round(CI_raw.prune, 4), "\n")


# ============================================================================
# 8. DIAGNOSTICS REPORTED IN THE PAPER
#
# Append to leuksurv_realdata.R after section 7. Everything it needs is in
# scope by then: alldata (with e_resid, spatial_pred, centered_log_rr, so_pred,
# so_centered_lrr), clinicaltree, Dipolar.clinical, spatial_risk, vhat,
# CI_sp.prune, covariates_clinical, n, and find_leaf() from
# computeresiduals_NA.R.
#
# Also REPLACES lines 313-319 of the original script, which used
#   alldata$so_high <- alldata$so_centered_lrr > 0
# selecting 856 of 1,043 patients rather than a region. See 8.6.
#
# No new package dependencies.
# ============================================================================

# ---------------------------------------------------------------------------
# 8.1  Sec 5.1: is tpi spatially structured?
#
# The claim the paper needs is that deprivation clusters geographically, so that
# clinical-spatial confounding is possible. A linear correlation with the
# coordinates only detects a planar trend and misses clustering, so three
# measures are reported: the planar one for reference, Moran's I, and the share
# of tpi variance that location explains.
# ---------------------------------------------------------------------------
cat("\n=== Sec 5.1: spatial structure in tpi ===\n")
cat("  districts:", length(unique(alldata$fac_district)), "\n")

tpi <- alldata$num_tpi
SST <- sum((tpi - mean(tpi))^2)

cat(sprintf("  tpi vs easting   Pearson %+.3f | tpi vs northing Pearson %+.3f\n",
            cor(tpi, alldata$num_xcoord), cor(tpi, alldata$num_ycoord)))
cat(sprintf("  planar multiple R (tpi on both coordinates): %.3f\n",
            sqrt(summary(lm(num_tpi ~ num_xcoord + num_ycoord, data = alldata))$r.squared)))

D <- as.matrix(dist(alldata[, c("num_xcoord", "num_ycoord")]))
diag(D) <- Inf
nn <- function(k) t(apply(D, 1, function(dd) order(dd)[1:k]))

z <- tpi - mean(tpi)
for (k in c(5, 10, 20)) {
  idx <- nn(k)
  lag <- rowMeans(matrix(z[idx], nrow = n))      # row-standardised weights
  cat(sprintf("  Moran's I, k = %-2d nearest neighbours: %.3f\n",
              k, sum(z * lag) / sum(z^2)))
}

for (k in c(10, 25, 50)) {
  idx <- nn(k)
  fit <- rowMeans(matrix(tpi[idx], nrow = n))    # leave-one-out by construction
  cat(sprintf("  variance of tpi explained by location, kNN k = %-2d: R2 = %.3f\n",
              k, 1 - sum((tpi - fit)^2) / SST))
}

fit_d <- ave(tpi, alldata$fac_district)
cat(sprintf("  variance of tpi lying between the %d districts: R2 = %.3f\n",
            length(unique(alldata$fac_district)), 1 - sum((tpi - fit_d)^2) / SST))
# The between-district share is the figure quoted in Sec 5.1. Moran's I and the
# kNN variance shares are reported alongside it as a robustness check: all three
# say location explains roughly a quarter of the variation in tpi, whereas the
# planar correlation sees almost none of it, deprivation being clustered rather
# than trended.

# ---------------------------------------------------------------------------
# 8.2  Sec 5.1: how much tpi signal survives each stage
# ---------------------------------------------------------------------------
cat("\n=== Sec 5.1: tpi correlations ===\n")
for (nm in c("e_resid", "centered_log_rr", "so_centered_lrr")) {
  cat(sprintf("  tpi vs %-16s Pearson %+.4f   Spearman %+.4f\n", nm,
              cor(tpi, alldata[[nm]]), cor(tpi, alldata[[nm]], method = "spearman")))
}

# ---------------------------------------------------------------------------
# 8.3  Sec 5.3: spatial leaf structure
# ---------------------------------------------------------------------------
cat("\n=== Sec 5.3: spatial leaf structure ===\n")
lt <- table(alldata$spatial_pred)
cat("  two-stage leaves:", length(lt),
    " sizes:", paste(sort(as.vector(lt)), collapse = ", "), "\n")
cat("  spatial-only leaves:", length(unique(alldata$so_pred)), "\n")

# ---------------------------------------------------------------------------
# 8.4  Sec 5.3: recovered amplitude
# ---------------------------------------------------------------------------
cat("\n=== Sec 5.3: recovered amplitude (5th-95th percentile span) ===\n")
sp <- function(v) as.numeric(diff(quantile(v, c(0.05, 0.95))))
cat(sprintf("  two-stage    %.3f\n", sp(alldata$centered_log_rr)))
cat(sprintf("  spatial-only %.3f\n", sp(alldata$so_centered_lrr)))
cat(sprintf("  spBayesSurv  %.3f\n", sp(spatial_risk)))

# ---------------------------------------------------------------------------
# 8.5  Sec 5.3: concordance bound implied by the benchmark frailty variance
#
# Under Theorem 1 the residuals are a common increasing transform of Exp(z(s))
# variables, so a pair whose log-frailties differ by D is ordered correctly with
# probability 1/(1 + exp(-|D|)). With v(s) ~ N(0, s2) at two locations,
# D ~ N(0, 2*s2).
# ---------------------------------------------------------------------------
cat("\n=== Sec 5.3: concordance bound from the frailty variance ===\n")
cbound <- function(s2, nsim = 2e6)
  mean(1 / (1 + exp(-abs(rnorm(nsim, 0, sqrt(2 * s2))))))
FRAILTY_VAR <- 0.067    # posterior mean, Table 1 -- confirm against summary(fit_ph)
cat(sprintf("  frailty variance %.3f -> bound %.4f\n",
            FRAILTY_VAR, cbound(FRAILTY_VAR)))
cat(sprintf("  observed two-stage residual concordance: %.4f\n", CI_sp.prune))
cat(sprintf("  (variance of the posterior mean frailties, for reference: %.4f)\n",
            var(vhat)))
# The bound uses the model's frailty VARIANCE PARAMETER, not var(vhat): the
# posterior means are shrunk, so var(vhat) understates the contrast.

# ---------------------------------------------------------------------------
# 8.6  Sec 5.3: the southern band the spatial-only map elevates
#
# Selected by rule rather than hard-coded, so it survives a refit. Every
# elevated leaf is printed first, so a shift is visible rather than silent.
# ---------------------------------------------------------------------------
cat("\n=== Sec 5.3: spatial-only southern band vs rest ===\n")
lf <- aggregate(cbind(so_centered_lrr, num_ycoord) ~ so_pred,
                data = alldata, FUN = mean)
names(lf) <- c("so_pred", "lrr", "mean_north")
lf$n <- as.vector(table(alldata$so_pred)[as.character(lf$so_pred)])
cat("  elevated spatial-only leaves:\n")
print(lf[lf$lrr > 0, ], row.names = FALSE, digits = 3)

BAND_NORTH_MAX <- 0.35
band_leaves <- lf$so_pred[lf$lrr > 0 & lf$mean_north < BAND_NORTH_MAX]
cat("  band leaves (elevated, mean northing <", BAND_NORTH_MAX, "):",
    paste(band_leaves, collapse = ", "), "\n")

alldata$so_high <- alldata$so_pred %in% band_leaves
cat("  band n =", sum(alldata$so_high), " rest n =", sum(!alldata$so_high), "\n")
for (v in c("num_wbc", "num_tpi", "num_age")) {
  a <- alldata[[v]][alldata$so_high]; b <- alldata[[v]][!alldata$so_high]
  cat(sprintf("  %-8s mean %7.2f vs %7.2f | median %7.2f vs %7.2f | MW p = %.2g\n",
              v, mean(a), mean(b), median(a), median(b), wilcox.test(a, b)$p.value))
}
cat(sprintf("  two-stage log RR over band: mean %+.3f (rest %+.3f); %.0f%% below average\n",
            mean(alldata$centered_log_rr[alldata$so_high]),
            mean(alldata$centered_log_rr[!alldata$so_high]),
            100 * mean(alldata$centered_log_rr[alldata$so_high] < 0)))

# ---------------------------------------------------------------------------
# 8.7  Sec 4.1: Nelson-Aalen vs Kaplan-Meier residuals
#
# Mirrors compute_residuals() exactly, reusing find_leaf(), but reads
# -log S_KM(t_i) off the leaf's KMest instead of the Nelson-Aalen increments.
# ---------------------------------------------------------------------------
cat("\n=== Sec 4.1: Nelson-Aalen vs Kaplan-Meier residuals ===\n")
km_at <- function(km, t_i) {
  idx <- which(km$time <= t_i)
  if (length(idx) == 0) 1 else km$surv[max(idx)]
}
testX <- as.matrix(alldata[, covariates_clinical])
e_km  <- numeric(n)
for (i in 1:n) {
  leaf <- find_leaf(clinicaltree, testX[i, ], Dipolar.clinical)
  e_km[i] <- -log(km_at(leaf$KMest, alldata[i, time]))
}
cat("  divergent (non-finite) KM residuals:", sum(!is.finite(e_km)), "of", n, "\n")
ok <- is.finite(e_km)
cat(sprintf("  Spearman vs Nelson-Aalen, finite only : %.4f\n",
            cor(e_km[ok], alldata$e_resid[ok], method = "spearman")))
cat(sprintf("  Spearman vs Nelson-Aalen, all patients: %.4f\n",
            cor(e_km, alldata$e_resid, method = "spearman")))
cat(sprintf("  Nelson-Aalen residual range: %.4f to %.4f\n",
            min(alldata$e_resid), max(alldata$e_resid)))


# ===========================================================================
# SAVE: artifacts for leuksurv_figure_4panel.R
# ---------------------------------------------------------------------------
# alldata_all_models.rds carries LeukSurv plus every fitted column, so the
# figure script runs from that file alone. twostage_session.RData supplies
# surf and spatial_risk, needed only for the smooth benchmark panel; without
# them the figure drops to three panels.
# ===========================================================================
saveRDS(alldata, "alldata_all_models.rds")
save(surf, spatial_risk, file = "twostage_session.RData")
cat("\nSaved: alldata_all_models.rds, twostage_session.RData\n")
