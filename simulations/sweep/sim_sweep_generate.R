# =====================================================================
# ROBUSTNESS SWEEP GENERATOR : three-zone base + confounded high zone
# ---------------------------------------------------------------------
# confound = 0 reproduces sim_threezone_generate_v2.R exactly: same
# geometry, same frailty levels, same tpi field, same b_tpi rule
# (confound_ratio * log(2) / within-zone SD), same censoring target and
# 95th-percentile follow-up, same crossing-hazard engine, same order of
# random draws. Only the hotspot shift below moves across the sweep.

#
# As `confound` rises, tpi is shifted upward for HIGH-zone patients, so
# tpi progressively predicts high-zone membership and the exogeneity
# condition z(S) _|_ X degrades in a controlled way. Nothing else moves.
#
# Hotspot = high zone: all three methods recover it at confound 0, so
# any separation in the hotspot-local curves is attributable to the
# confounding alone.
#
# USAGE: source, then  gen <- make_sweep_data(confound = 0.5, seed = 3)
# =====================================================================

library(survival)

# Flat alongside this script on the cluster, or ../../data/ in the repo.
LEUK_PATH <- if (file.exists("LeukSurv.csv")) "LeukSurv.csv" else "../../data/LeukSurv.csv"
HOTSPOT_X <- 0.25; HOTSPOT_Y <- 0.60
Z_LEVELS  <- c(low = 0.6, medium = 1.2, high = 2.4)

lp1_of_s <- function(x, y) {
  -1.2*x - 0.8*y + 1.2*y^2 +
    1.2*exp(-((x - HOTSPOT_X)^2 + (y - HOTSPOT_Y)^2) / (2 * 0.15^2))
}
lp2_of_s <- function(x, y) { -1.2*y + 0.24*x^2 }

make_sweep_data <- function(confound = 0.0, seed = 1,
                            hotspot_zone = "high",
                            b_wbc = 0.60,
                            confound_ratio = 1.0, target_cens = 0.16,
                            lambda0 = 0.02, gamma0 = 2.00, t_cross = 18,
                            leuk_path = LEUK_PATH, verbose = FALSE) {
  set.seed(seed)

  LeukSurv <- read.csv(leuk_path)
  n <- nrow(LeukSurv)
  easting  <- LeukSurv$num_xcoord
  northing <- LeukSurv$num_ycoord

  # ---- three-zone frailty (identical to the validated design) --------
  lp1 <- lp1_of_s(easting, northing)
  split1_thresh <- median(lp1)
  right_side <- (lp1 < split1_thresh)
  lp2 <- lp2_of_s(easting, northing)
  split2_thresh <- median(lp2[right_side])
  above_line <- (lp2 < split2_thresh)
  zone <- rep("medium", n)
  zone[right_side & !above_line] <- "high"
  zone[right_side &  above_line] <- "low"
  zone <- factor(zone, levels = c("low","medium","high"))
  z_s  <- Z_LEVELS[as.character(zone)]
  is_hot <- zone == hotspot_zone

  # ---- clinical covariates -------------------------------------------
  age <- rbinom(n, 1, 0.5)
  wbc <- rnorm(n, 0, 1)

  # tpi: N-S gradient, mean-balanced WITHIN zones (base, confound = 0),
  # then a graded shift added to the hotspot zone.
  tpi_raw <- 2.0 * (northing - 0.5) + rnorm(n, 0, 0.35)   # as sim_threezone_generate_v2.R
  tpi <- tpi_raw
  for (zz in levels(zone)) {
    idx <- zone == zz
    tpi[idx] <- tpi_raw[idx] - mean(tpi_raw[idx])
  }
  # b_tpi is fixed from the UNCONFOUNDED field, so the covariate effect is
  # identical at every confounding level and only the hotspot shift moves.
  tpi_sd_within <- mean(sapply(levels(zone), function(zz) sd(tpi[zone == zz])))
  b_tpi <- confound_ratio * log(2) / tpi_sd_within

  tpi[is_hot] <- tpi[is_hot] + confound * sd(tpi)

  # ---- multiplicative frailty + crossing-hazard age -------------------
  eta_fixed <- b_wbc * wbc + b_tpi * tpi
  k_co <- 2 * gamma0 / t_cross
  H1   <- (exp(gamma0) / k_co) * (1 - exp(-2 * gamma0))
  U <- runif(n); target <- -log(U)
  rhs <- target / (z_s * lambda0 * exp(eta_fixed))
  Tevent <- numeric(n)
  for (i in 1:n) {
    if (age[i] == 0)          Tevent[i] <- rhs[i]
    else if (rhs[i] <= H1)    Tevent[i] <- -(1/k_co)*log(1 - (k_co/exp(gamma0))*rhs[i])
    else                      Tevent[i] <- t_cross + (rhs[i] - H1)/exp(-gamma0)
  }

  # End of follow-up at the 95th percentile of event times, and a censoring
  # rate solved once (at this replicate's confound = 0 event times) so the
  # censoring MECHANISM is identical across the sweep. The realized
  # fraction drifts slightly as confounding shortens high-zone survival,
  # which is a consequence of the confounding rather than a moving part.
  foltime <- as.numeric(quantile(Tevent, 0.95))
  cens_frac <- function(rate) {
    set.seed(99)
    mean(!(Tevent <= pmin(rexp(n, rate = rate), foltime)))
  }
  rate_cens <- tryCatch(
    uniroot(function(r) cens_frac(r) - target_cens, interval = c(1e-6, 1))$root,
    error = function(e) 1/2000)

  set.seed(99)
  Ccens  <- rexp(n, rate = rate_cens)
  status <- as.integer(Tevent <= pmin(Ccens, foltime))
  stop_t <- pmin(Tevent, Ccens, foltime)

  sim_data <- data.frame(
    nid = 1:n, status = status, start = 0, stop = stop_t,
    num_xcoord = easting, num_ycoord = northing,
    num_age = age, num_wbc = wbc, num_tpi = tpi,
    zone = as.character(zone), z_true = z_s, is_hot = is_hot
  )

  # ---- diagnostics ----------------------------------------------------
  cor_tpi_hot <- cor(as.numeric(is_hot), tpi)
  within_zone_cor_north <- mean(sapply(levels(zone), function(zz){
    idx <- zone == zz; cor(northing[idx], tpi[idx]) }))
  zone_tpi_means <- tapply(tpi, zone, mean)

  grid <- expand.grid(num_xcoord = seq(0, 1, length.out = 400),
                      num_ycoord = seq(0, 1, length.out = 400))
  grid$lp1 <- lp1_of_s(grid$num_xcoord, grid$num_ycoord)
  grid$lp2 <- lp2_of_s(grid$num_xcoord, grid$num_ycoord)
  grid$lp2_masked <- ifelse(grid$lp1 < split1_thresh, grid$lp2, NA)

  if (verbose) {
    cat(sprintf("confound=%.2f seed=%d  cor(tpi,hot)=%+.3f  within-zone cor(tpi,N)=%+.3f  cens=%.3f\n",
                confound, seed, cor_tpi_hot, within_zone_cor_north, mean(1-status)))
    cat("  zone counts: "); print(table(zone))
    cat("  zone tpi means: "); print(round(zone_tpi_means, 3))
  }

  list(sim_data = sim_data, zone = zone, z_s = z_s, is_hot = is_hot,
       hotspot_zone = hotspot_zone, z_levels = Z_LEVELS, grid = grid,
       split1_thresh = split1_thresh, split2_thresh = split2_thresh,
       confound = confound, seed = seed,
       b_tpi = b_tpi, foltime = foltime, rate_cens = rate_cens,
       cor_tpi_hot = cor_tpi_hot,
       within_zone_cor_north = within_zone_cor_north,
       zone_tpi_means = zone_tpi_means,
       cens_rate = mean(1 - status))
}
