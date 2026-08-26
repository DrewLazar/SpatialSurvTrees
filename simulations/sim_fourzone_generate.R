# =====================================================================
# FOUR-ZONE DESIGN  (exogeneity condition, mean-balanced spatial tpi)
# ---------------------------------------------------------------------
# Built from sim_fourzone_generate.R (the base placement, not the
# variant2 swap): the extreme frailty levels sit on either side of the
# long left-side horizontal boundary, creating a sharp low-adjacent-to-
# high step that a smooth frailty surface cannot represent.
#
# Three splits: a parabola+bump separating left from right, then a
# horizontal split within each side.
#
#   * multiplicative frailty z(s) * lambda_clin(t|x), exponential baseline
#   * crossing-hazard age effect (violates PH, consistent with Theorem 1)
#   * tpi: north-south field, mean-balanced within each of the four zones
#   * closed-form inversion; censoring solved to a target fraction
#
# Frailty levels are in a uniform ratio FRAILTY_RATIO between adjacent
# zones, set below. For
# exponential residuals sd = mean, so the separability of two zones
# depends only on their ratio, not their levels:
#   d = (1 - 1/r) / sqrt((1 + 1/r^2)/2)
# so every adjacent pair carries the same separability and the zone count
# is a complexity increase rather than a change in contrast.
#
# STRENGTH. b_tpi is set so that
#     b_tpi * SD(tpi within zone) = CONFOUND_RATIO * log(FRAILTY_RATIO)
# matching the convention used in the two- and three-zone designs. The
# realized ratio is printed below.
#
# tpi is loaded on northing only. The zone boundaries are largely
# horizontal, so an easting component would be stripped by the
# within-zone centering and contribute nothing.
# =====================================================================

library(survival)
library(ggplot2)

set.seed(20260628)

DATA_FILE      <- "../data/LeukSurv.csv"
CONFOUND_RATIO <- 1.0
TARGET_CENS    <- 0.10
FRAILTY_RATIO  <- 2.5    # uniform ratio between adjacent zones

# ---- 1. Real LeukSurv coordinates ------------------------------------
LeukSurv <- read.csv(DATA_FILE)
n <- nrow(LeukSurv)
easting  <- LeukSurv$num_xcoord
northing <- LeukSurv$num_ycoord

# ---- 2. Four-zone geometry -------------------------------------------
hotspot_x <- 0.25; hotspot_y <- 0.60
bump <- exp(-((easting - hotspot_x)^2 + (northing - hotspot_y)^2) / (2 * 0.15^2))

lp1 <- -1.2*easting - 0.8*northing + 1.2*northing^2 + 1.2*bump
split1_thresh <- median(lp1)
right_side <- (lp1 < split1_thresh)

lp2 <- -1.2*northing + 0.24*easting^2
split2_thresh <- median(lp2[right_side])
above_line_right <- (lp2 < split2_thresh)

lp3 <- -northing + 0.15*easting^2
split3_thresh <- median(lp3[!right_side])
above_line_left <- (lp3 < split3_thresh)

zone <- rep(NA_character_, n)
zone[!right_side & above_line_left]  <- "low"           # upper-left
zone[!right_side & !above_line_left] <- "high"          # lower-left
zone[right_side & above_line_right]  <- "medium_low"    # upper-right
zone[right_side & !above_line_right] <- "medium_high"   # lower-right
zone <- factor(zone, levels = c("low","medium_low","medium_high","high"))
cat("Zone counts:\n"); print(table(zone))

# ---- 3. Frailty levels: uniform ratio between adjacent zones ---------
# Levels are FRAILTY_RATIO^(-1.5, -0.5, 0.5, 1.5), so every adjacent pair
# has the same ratio and the geometric mean of the four is one. Raising
# the ratio widens the contrast without disturbing that property.
z_levels <- FRAILTY_RATIO^c(-1.5, -0.5, 0.5, 1.5)
names(z_levels) <- levels(zone)
z_s <- z_levels[as.character(zone)]
names(z_s) <- NULL
frailty_contrast <- log(FRAILTY_RATIO)

cohen_d <- (1 - 1/FRAILTY_RATIO) / sqrt((1 + 1/FRAILTY_RATIO^2)/2)
cat(sprintf("\nFrailty levels (ratio %.2f): %s\n",
            FRAILTY_RATIO, paste(sprintf("%.3f", z_levels), collapse = ", ")))
cat(sprintf("  adjacent-pair separability (Cohen's d) = %.3f\n", cohen_d))
cat(sprintf("  full log-risk span = %.3f\n", log(max(z_levels)/min(z_levels))))

# ---- 4. Covariates ---------------------------------------------------
age <- rbinom(n, 1, 0.5)          # binary crossing-hazard carrier
wbc <- rnorm(n, 0, 1)             # ordinary PH

tpi_raw <- 3.0 * northing + rnorm(n, 0, 0.4)
tpi <- tpi_raw
for (zz in levels(zone)) {
  idx <- zone == zz
  tpi[idx] <- tpi_raw[idx] - mean(tpi_raw[idx])
}

# ---- 5. Clinical hazard ----------------------------------------------
lambda0 <- 0.02
b_wbc   <- 0.60
gamma0  <- 2.00
t_cross <- 18

tpi_sd_within <- mean(sapply(levels(zone), function(zz) sd(tpi[zone == zz])))
b_tpi <- CONFOUND_RATIO * frailty_contrast / tpi_sd_within

eta_fixed <- b_wbc * wbc + b_tpi * tpi

# ---- CHECKS -----------------------------------------------------------
cat("\n--- tpi mean-balance per zone (want all ~0) ---\n")
for (zz in levels(zone)) {
  idx <- zone == zz
  cat(sprintf("  %-12s mean tpi = %+.4f  (n=%d)\n", zz, mean(tpi[idx]), sum(idx)))
}
cat(sprintf("\n  cor(z_s, tpi)      = %+.3f   [want ~0: mean-balance]\n", cor(z_s, tpi)))
cat(sprintf("  cor(northing, tpi) = %+.3f   [want substantial]\n", cor(northing, tpi)))
cat(sprintf("  cor(easting,  tpi) = %+.3f   [want weak]\n", cor(easting, tpi)))
cat("  per-zone cor(northing, tpi):\n")
for (zz in levels(zone)) {
  idx <- zone == zz
  cat(sprintf("    %-12s %+.3f\n", zz, cor(northing[idx], tpi[idx])))
}

cat("\n--- Confounding strength ---\n")
cat(sprintf("      b_tpi                         = %.3f  (set, not tuned)\n", b_tpi))
cat(sprintf("      within-zone SD of tpi         = %.3f\n", tpi_sd_within))
cat(sprintf("      within-zone clinical gradient = %.3f\n", b_tpi * tpi_sd_within))
cat(sprintf("      adjacent-zone frailty contrast= %.3f\n", frailty_contrast))
cat(sprintf("      realized confound ratio       = %.3f  (target %.2f)\n",
            b_tpi * tpi_sd_within / frailty_contrast, CONFOUND_RATIO))

# ---- 6. Closed-form inversion (piecewise crossing-age) ---------------
k_co <- 2 * gamma0 / t_cross
H1   <- (exp(gamma0) / k_co) * (1 - exp(-2 * gamma0))
U <- runif(n); target <- -log(U)
scale_i <- z_s * lambda0 * exp(eta_fixed)
rhs <- target / scale_i
Tevent <- numeric(n)
for (i in 1:n) {
  if (age[i] == 0) {
    Tevent[i] <- rhs[i]
  } else if (rhs[i] <= H1) {
    Tevent[i] <- -(1 / k_co) * log(1 - (k_co / exp(gamma0)) * rhs[i])
  } else {
    Tevent[i] <- t_cross + (rhs[i] - H1) / exp(-gamma0)
  }
}

# ---- 7. Censoring ----------------------------------------------------
foltime <- as.numeric(quantile(Tevent, 0.95))

cens_frac <- function(rate) {
  set.seed(99)
  Cc <- rexp(n, rate = rate)
  mean(!(Tevent <= pmin(Cc, foltime)))
}
rate_cens <- uniroot(function(r) cens_frac(r) - TARGET_CENS,
                     interval = c(1e-5, 1e-1))$root
cat(sprintf("\ncensoring rate %.6f (mean censoring time %.0f) -> %.3f censored\n",
            rate_cens, 1 / rate_cens, cens_frac(rate_cens)))

set.seed(99)
Ccens  <- rexp(n, rate = rate_cens)
status <- as.integer(Tevent <= pmin(Ccens, foltime))
stop_t <- pmin(Tevent, Ccens, foltime)

# ---- 8. Assemble -----------------------------------------------------
sim_data <- data.frame(
  nid = 1:n, status = status, start = 0, stop = stop_t,
  num_xcoord = easting, num_ycoord = northing,
  num_age = age, num_wbc = wbc, num_tpi = tpi,
  zone = as.character(zone), z_true = z_s
)
cat("\n--- Dataset summary ---\n")
cat("n =", n, " events =", sum(status),
    " censoring rate =", round(mean(1 - status), 3), "\n")
cat(sprintf("  end of follow-up (95th pct of event times) = %.1f\n", foltime))
cat(sprintf("  administrative censoring: %d of %d censored (%.0f%%)\n",
            sum(status == 0 & stop_t >= foltime - 1e-8), sum(status == 0),
            100 * sum(status == 0 & stop_t >= foltime - 1e-8) / max(sum(status == 0), 1)))
for (zz in levels(zone)) {
  idx <- sim_data$zone == zz
  cat(sprintf("  zone %-12s: n=%d  median time=%.2f  events=%d  cens=%.3f\n",
              zz, sum(idx), median(sim_data$stop[idx]), sum(sim_data$status[idx]),
              mean(1 - sim_data$status[idx])))
}

# ---- 9. Grid + boundaries for plotting -------------------------------
grid <- expand.grid(num_xcoord = seq(0, 1, length.out = 400),
                    num_ycoord = seq(0, 1, length.out = 400))
gb <- exp(-((grid$num_xcoord - hotspot_x)^2 + (grid$num_ycoord - hotspot_y)^2) /
            (2 * 0.15^2))
grid$lp1 <- -1.2*grid$num_xcoord - 0.8*grid$num_ycoord +
             1.2*grid$num_ycoord^2 + 1.2*gb
grid$lp2 <- -1.2*grid$num_ycoord + 0.24*grid$num_xcoord^2
grid$lp3 <- -grid$num_ycoord + 0.15*grid$num_xcoord^2
grid$lp2_masked <- ifelse(grid$lp1 <  split1_thresh, grid$lp2, NA)  # right side
grid$lp3_masked <- ifelse(grid$lp1 >= split1_thresh, grid$lp3, NA)  # left side

cat("\nObjects: sim_data, z_s, zone, grid, lp1, lp2, lp3,\n")
cat("         split1_thresh, split2_thresh, split3_thresh, right_side,\n")
cat("         z_levels, b_wbc, b_tpi, frailty_contrast\n")
