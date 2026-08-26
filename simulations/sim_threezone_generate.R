# =====================================================================
# THREE-ZONE DESIGN (Prop 1 mean-balance, non-degenerate)
# ---------------------------------------------------------------------
# Sibling of the two-zone (sim_zperpx_generate_v2.R) and four-zone
# designs. Same machinery: real LeukSurv coordinates, crossing-hazard
# age, north-south tpi gradient mean-balanced within zones, exponential
# censoring.
#
#   * Frailty z(S): two-split three-zone construction. Split 1 (parabola
#     + bump) divides medium (left) from the right; split 2 (gentle
#     horizontal, right side only) divides the right into high (below)
#     and low (above).
#   * tpi: north-south gradient, mean-balanced WITHIN each zone.
#     Centering removes the between-zone means; the within-zone N-S
#     slope survives, so tpi stays spatially structured while remaining
#     mean-independent of the frailty.
#   * age: independent binary crossing-hazard carrier (non-PH).
#   * wbc: ordinary independent PH covariate.
#
# STRENGTH. b_tpi is set so that
#     b_tpi * SD(tpi within zone) = CONFOUND_RATIO * log(2)
# where log(2) is the frailty contrast between adjacent zones (levels
# are in a uniform ratio of two). CONFOUND_RATIO = 1 puts the spurious
# within-zone gradient on equal footing with the true frailty step. This
# is a stated design parameter, not a value tuned to a result, and the
# realized ratio is printed below. The same convention is used in the
# two- and four-zone designs.
# =====================================================================

library(survival)
library(ggplot2)

set.seed(20260701)

DATA_FILE      <- "../data/LeukSurv.csv"
CONFOUND_RATIO <- 1.0

# ---- 1. Spatial support: real LeukSurv coordinates -------------------
LeukSurv <- read.csv(DATA_FILE)
n <- nrow(LeukSurv)
easting  <- LeukSurv$num_xcoord
northing <- LeukSurv$num_ycoord

# ---- 2. Frailty z(s): two-split three-zone construction --------------
hotspot_x <- 0.25; hotspot_y <- 0.60
lp1_of_s <- function(x, y) {
  -1.2*x - 0.8*y + 1.2*y^2 +
    1.2*exp(-((x - hotspot_x)^2 + (y - hotspot_y)^2) / (2 * 0.15^2))
}
lp2_of_s <- function(x, y) { -1.2*y + 0.24*x^2 }

lp1 <- lp1_of_s(easting, northing)
split1_thresh <- median(lp1)
right_side <- (lp1 < split1_thresh)
lp2 <- lp2_of_s(easting, northing)
split2_thresh <- median(lp2[right_side])
above_line <- (lp2 < split2_thresh)

zone <- rep("medium", n)
zone[right_side & !above_line] <- "high"
zone[right_side &  above_line] <- "low"
zone <- factor(zone, levels = c("low", "medium", "high"))

z_levels <- c(low = 0.6, medium = 1.2, high = 2.4)
z_s <- z_levels[as.character(zone)]
frailty_contrast <- log(2)          # adjacent-zone ratio
cat("Zone counts:\n"); print(table(zone))

# ---- 3. Clinical covariates ------------------------------------------
age <- rbinom(n, 1, 0.5)
wbc <- rnorm(n, 0, 1)

# tpi: north-south field, centered within each zone
tpi_raw <- 2.0 * (northing - 0.5) + rnorm(n, 0, 0.35)
tpi <- tpi_raw
for (zz in levels(zone)) {
  idx <- zone == zz
  tpi[idx] <- tpi_raw[idx] - mean(tpi_raw[idx])
}

# ---- 4. Clinical hazard ----------------------------------------------
lambda0 <- 0.02
b_wbc   <- 0.60
gamma0  <- 2.00
t_cross <- 18
foltime <- NA          # set after Tevent is drawn: 95th percentile

tpi_sd_within <- mean(sapply(levels(zone), function(zz) sd(tpi[zone == zz])))
b_tpi <- CONFOUND_RATIO * frailty_contrast / tpi_sd_within

eta_fixed <- b_wbc * wbc + b_tpi * tpi

# ---- CHECKS -----------------------------------------------------------
cat("\n--- Independence / structure checks ---\n")
cat(sprintf("(i)   mean within-zone cor(tpi, northing) = %+.3f  [want strong]\n",
            mean(sapply(levels(zone), function(zz){
              idx <- zone == zz; cor(northing[idx], tpi[idx]) }))))
cat("(ii)  tpi mean by zone (want ~equal -> balanced):\n")
print(round(tapply(tpi, zone, mean), 4))
cat(sprintf("(iii) cor(as.integer(zone), tpi) = %+.3f  [want ~0]\n",
            cor(as.integer(zone), tpi)))

cat("\n--- Confounding strength ---\n")
cat(sprintf("      b_tpi                         = %.3f  (set, not tuned)\n", b_tpi))
cat(sprintf("      within-zone SD of tpi         = %.3f\n", tpi_sd_within))
cat(sprintf("      within-zone clinical gradient = %.3f\n", b_tpi * tpi_sd_within))
cat(sprintf("      adjacent-zone frailty contrast= %.3f\n", frailty_contrast))
cat(sprintf("      realized confound ratio       = %.3f  (target %.2f)\n",
            b_tpi * tpi_sd_within / frailty_contrast, CONFOUND_RATIO))

# ---- 5. Closed-form inversion (piecewise crossing-age) ---------------
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

# ---- 6. Censoring ----------------------------------------------------
# End of follow-up at the 95th percentile of the realized event times, so
# administrative censoring binds for about 5% of patients rather than
# never. The remaining censoring is the random Exp(1/2000) draw.
foltime <- as.numeric(quantile(Tevent, 0.95))

# Censoring rate solved to hit a target censoring fraction, so that the
# designs are comparable despite their different frailty levels and
# follow-up windows. The rate is a consequence of the target, not a
# separately chosen parameter.
TARGET_CENS <- 0.16
cens_frac <- function(rate) {
  set.seed(99)
  Cc <- rexp(n, rate = rate)
  mean(!(Tevent <= pmin(Cc, foltime)))
}
rate_cens <- uniroot(function(r) cens_frac(r) - TARGET_CENS,
                     interval = c(1e-5, 1e-1))$root
cat(sprintf("censoring rate %.6f (mean censoring time %.0f) -> %.3f censored\n",
            rate_cens, 1/rate_cens, cens_frac(rate_cens)))

set.seed(99)
Ccens <- rexp(n, rate = rate_cens)
status <- as.integer(Tevent <= pmin(Ccens, foltime))
stop_t <- pmin(Tevent, Ccens, foltime)

# ---- 7. Assemble -----------------------------------------------------
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
  cat(sprintf("  zone %-6s: n=%d  median time=%.2f  events=%d  cens=%.3f\n",
              zz, sum(idx), median(sim_data$stop[idx]), sum(sim_data$status[idx]),
              mean(1 - sim_data$status[idx])))
}

# ---- 8. Grid + boundaries for plotting -------------------------------
grid <- expand.grid(num_xcoord = seq(0, 1, length.out = 400),
                    num_ycoord = seq(0, 1, length.out = 400))
grid$lp1 <- lp1_of_s(grid$num_xcoord, grid$num_ycoord)
grid$lp2 <- lp2_of_s(grid$num_xcoord, grid$num_ycoord)
grid$lp2_masked <- ifelse(grid$lp1 < split1_thresh, grid$lp2, NA)

in_high <- zone == "high"

cat("\nObjects: sim_data, z_s, zone, in_high, grid, lp1, lp2,\n")
cat("         split1_thresh, split2_thresh, lp1_of_s, lp2_of_s, z_levels,\n")
cat("         right_side, b_wbc, b_tpi, frailty_contrast\n")43