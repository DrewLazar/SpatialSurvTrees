# =====================================================================
# DEGENERATE DESIGN  (S _|_ X : full independence, two-zone frailty)
# ---------------------------------------------------------------------
# The baseline / degenerate case for the exogeneity ladder. Here ALL
# clinical covariates (age, wbc, tpi) are generated INDEPENDENTLY of
# location, so S _|_ X holds in full (not merely mean-balanced). Under
# this strong condition the clinical covariates carry NO spatial
# information, so a spatial-only tree is not biased by them -- it only
# pays a variance penalty. The two-stage method's stage-1 residualization
# therefore has no confound to remove and is expected to COINCIDE with
# the spatial-only tree (possibly slightly worse, from stage-1 variance).
#
# This is the LOWER RUNG of the ladder: it establishes that the clinical
# stage is redundant under full exogeneity, which is a correctness check
# and motivates the non-degenerate designs (z(S) _|_ X with X spatially
# structured) where the two-stage method earns its keep.
#
# Kept consistent with the two-zone z(S)-mean-balanced design (same
# frailty boundary, same crossing-hazard age, same closed-form inversion,
# same censoring). The ONLY change is that tpi is drawn independently of
# location rather than as a spatial field.
# =====================================================================

# Paths are relative to simulations/. Set the working directory there
# before sourcing, or open the file in RStudio and use Session > Set
# Working Directory > To Source File Location.

library(survival)
library(ggplot2)

set.seed(20260628)

# ---- 1. Spatial support: real LeukSurv coordinates -------------------
LeukSurv <- read.csv("../data/LeukSurv.csv")
n <- nrow(LeukSurv)
easting  <- LeukSurv$num_xcoord
northing <- LeukSurv$num_ycoord

# ---- 2. Frailty z(s): parabola + bump (E-W split), same as z-perp-X --
hotspot_x <- 0.25; hotspot_y <- 0.60
g_of_s <- function(x, y) {
  -3*x - 2*y + 3*y^2 +
    3*exp(-((x - hotspot_x)^2 + (y - hotspot_y)^2) / (2 * 0.15^2))
}
g_patient <- g_of_s(easting, northing)
boundary  <- median(g_patient)
in_high   <- g_patient > boundary
z_high <- 1.6; z_low <- 0.7
z_s <- ifelse(in_high, z_high, z_low)
cat("Zone counts:  high =", sum(in_high), "  low =", sum(!in_high), "\n")

# ---- 3. Clinical covariates: ALL INDEPENDENT OF LOCATION -------------
# This is the ONLY structural difference from the z-perp-X design:
# tpi is NOT a spatial field here. Every covariate is drawn independently
# of the coordinates, so S _|_ X holds fully.
age <- rbinom(n, 1, 0.5)          # binary crossing-hazard carrier
wbc <- rnorm(n, 0, 1)             # ordinary PH
tpi <- rnorm(n, 0, 1)             # ordinary PH, INDEPENDENT of location

# ---- CHECKS: confirm full independence S _|_ X -----------------------
cat("\n--- Independence checks (all should be ~0) ---\n")
cat(sprintf("  cor(in_high, age)      = %+.3f\n", cor(in_high, age)))
cat(sprintf("  cor(in_high, wbc)      = %+.3f\n", cor(in_high, wbc)))
cat(sprintf("  cor(in_high, tpi)      = %+.3f\n", cor(in_high, tpi)))
cat(sprintf("  cor(northing, tpi)     = %+.3f   [~0: tpi NOT spatial here]\n",
            cor(northing, tpi)))
cat(sprintf("  cor(easting,  tpi)     = %+.3f\n", cor(easting, tpi)))

# ---- 4. Clinical hazard: same as z-perp-X ----------------------------
lambda0 <- 0.02
b_wbc   <- 0.60
b_tpi   <- 1.10       # kept identical to z-perp-X for comparability, even
                      # though tpi is non-spatial here (so it cannot confound
                      # the spatial tree -- it is just independent noise).
gamma0  <- 2.00
t_cross <- 18
foltime <- 500
eta_fixed <- b_wbc * wbc + b_tpi * tpi

# ---- 5. Closed-form inversion (same piecewise age structure) ---------
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
Ccens <- rexp(n, rate = 1 / 2000)
status <- as.integer(Tevent <= pmin(Ccens, foltime))
stop_t <- pmin(Tevent, Ccens, foltime)

# ---- 7. Assemble -----------------------------------------------------
sim_data <- data.frame(
  nid = 1:n, status = status, start = 0, stop = stop_t,
  num_xcoord = easting, num_ycoord = northing,
  num_age = age, num_wbc = wbc, num_tpi = tpi,
  zone = ifelse(in_high, "high", "low"), z_true = z_s
)
cat("\n--- Dataset summary ---\n")
cat("n =", n, " events =", sum(status),
    " censoring rate =", round(mean(1 - status), 3), "\n")
for (zz in c("high", "low")) {
  idx <- sim_data$zone == zz
  cat(sprintf("  zone %-4s: n=%d  median time=%.2f  events=%d\n",
              zz, sum(idx), median(sim_data$stop[idx]), sum(sim_data$status[idx])))
}

# ---- 8. Visuals: true frailty and (non-spatial) tpi ------------------
grid <- expand.grid(num_xcoord = seq(0,1,length.out=300),
                    num_ycoord = seq(0,1,length.out=300))
grid$g <- g_of_s(grid$num_xcoord, grid$num_ycoord)

p_frailty <- ggplot() +
  geom_point(data = sim_data, aes(num_xcoord, num_ycoord, color = zone),
             size = 1.4, alpha = 0.85) +
  scale_color_manual(values = c(high="red", low="blue"), name="True zone") +
  geom_contour(data = grid, aes(num_xcoord, num_ycoord, z = g),
               breaks = boundary, color = "black", linewidth = 1) +
  labs(title = "True frailty z(s): E-W boundary (degenerate S _|_ X)",
       x="Easting", y="Northing") +
  coord_fixed(xlim=c(0,1), ylim=c(0,1)) + theme_minimal()

# tpi map -- should look like SPATIAL NOISE (no gradient), unlike z-perp-X
p_tpi <- ggplot() +
  geom_point(data = sim_data, aes(num_xcoord, num_ycoord, color = num_tpi),
             size = 1.4, alpha = 0.85) +
  scale_color_gradient2(low="blue", mid="gray80", high="red",
                        midpoint = median(sim_data$num_tpi), name="tpi") +
  geom_contour(data = grid, aes(num_xcoord, num_ycoord, z = g),
               breaks = boundary, color = "black", linewidth = 1) +
  labs(title = "tpi: independent of location (no spatial gradient)",
       x="Easting", y="Northing") +
  coord_fixed(xlim=c(0,1), ylim=c(0,1)) + theme_minimal()

cat("\nPlots: p_frailty (true zones), p_tpi (should be spatial NOISE, no gradient).\n")
cat("Objects: sim_data, z_s, in_high, boundary, g_of_s, grid, g_patient\n")
