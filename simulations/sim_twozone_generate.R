# =====================================================================
# SPATIALLY-STRUCTURED-X DESIGN  (Prop 1 mean-balance, non-degenerate)
# ---------------------------------------------------------------------
# GEOMETRY, stated precisely, because it determines what the design can
# demonstrate.
#
# The frailty boundary is the level set of g: a curve running roughly
# NORTH-SOUTH that separates a high-frailty WEST from a low-frailty EAST.
# So "crossing the boundary" means moving east-west, and "travelling
# along the boundary" means moving north-south.
#
# tpi is centered within each frailty zone, which is what secures the
# Proposition 1 mean-balance condition E[z(S)|X] = const. Centering
# removes exactly the ACROSS-boundary (east-west) component of tpi: that
# component is what would otherwise make tpi differ between zones. What
# survives centering is the ALONG-boundary (north-south) component.
#
# The confounding available to mislead a spatial-only analysis is
# therefore the WITHIN-ZONE north-south gradient, and only that. A tpi
# field loaded on easting contributes nothing after centering.
#
#   * Frailty z(S): parabola + bump boundary, high west / low east.
#   * tpi: a NORTH-SOUTH field, so its surviving within-zone gradient
#     competes with the frailty for the spatial tree's attention.
#   * age: independent binary crossing-hazard carrier (non-PH).
#   * wbc: ordinary independent PH covariate.
#
# STRENGTH. Whether the spatial-only tree is actually misled depends on
# how the within-zone clinical gradient compares with the between-zone
# frailty contrast. b_tpi is therefore set so that
#
#     b_tpi * SD(tpi within zone)  =  CONFOUND_RATIO * log(z_high/z_low)
#
# CONFOUND_RATIO = 1 puts the spurious within-zone gradient on equal
# footing with the true frailty step. This is a stated design parameter,
# not a value tuned to a result, and the realized ratio is printed below.
#
# CHECKS (measured, not assumed):
#   (i)   cor(in_high, tpi) ~ 0     -> mean-balance across zones
#   (ii)  cor(northing, tpi) strong -> tpi genuinely spatial, N-S
#   (iii) realized confound ratio   -> the strength of the confounding
# =====================================================================

library(survival)
library(ggplot2)

set.seed(202606262)

DATA_FILE      <- "../data/LeukSurv.csv"
CONFOUND_RATIO <- 1.0    # within-zone clinical gradient / frailty contrast

# ---- 1. Spatial support: real LeukSurv coordinates -------------------
LeukSurv <- read.csv(DATA_FILE)
n <- nrow(LeukSurv)
easting  <- LeukSurv$num_xcoord
northing <- LeukSurv$num_ycoord

# ---- 2. Frailty z(s): parabola + bump (high west / low east) ---------
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
frailty_contrast <- log(z_high / z_low)
cat("Zone counts:  high =", sum(in_high), "  low =", sum(!in_high), "\n")

# ---- 3. Clinical covariates ------------------------------------------
age <- rbinom(n, 1, 0.5)          # independent, crossing-hazard carrier
wbc <- rnorm(n, 0, 1)             # independent, ordinary PH

# tpi: NORTH-SOUTH field. No easting term: an easting component would be
# stripped by the within-zone centering below and contribute nothing.
tpi_raw <- 3.0 * northing + rnorm(n, 0, 0.4)
tpi <- tpi_raw
tpi[in_high]  <- tpi_raw[in_high]  - mean(tpi_raw[in_high])
tpi[!in_high] <- tpi_raw[!in_high] - mean(tpi_raw[!in_high])

# ---- 4. Clinical hazard ----------------------------------------------
lambda0 <- 0.02
b_wbc   <- 0.60
gamma0  <- 2.00
t_cross <- 18
foltime <- 500

# b_tpi set from the design parameter, not chosen by hand
tpi_sd_within <- mean(c(sd(tpi[in_high]), sd(tpi[!in_high])))
b_tpi <- CONFOUND_RATIO * frailty_contrast / tpi_sd_within

eta_fixed <- b_wbc * wbc + b_tpi * tpi

# ---- CHECKS -----------------------------------------------------------
cat("\n--- Independence / structure checks ---\n")
cat(sprintf("(i)   cor(in_high, tpi)   = %+.3f   [want ~0: mean-balanced]\n",
            cor(in_high, tpi)))
cat(sprintf("      cor(in_high, age)   = %+.3f   cor(in_high, wbc) = %+.3f\n",
            cor(in_high, age), cor(in_high, wbc)))
cat(sprintf("(ii)  cor(northing, tpi)  = %+.3f   [want strong: N-S spatial]\n",
            cor(northing, tpi)))
cat(sprintf("      cor(easting,  tpi)  = %+.3f   [want weak: no E-W component]\n",
            cor(easting, tpi)))
cat(sprintf("(iii) mean tpi  high-zone = %+.3f   low-zone = %+.3f  [want ~equal]\n",
            mean(tpi[in_high]), mean(tpi[!in_high])))
cat(sprintf("      tpi ~ zone Mann-Whitney p = %.3f  [want LARGE: balanced]\n",
            wilcox.test(tpi ~ in_high)$p.value))

cat("\n--- Confounding strength ---\n")
cat(sprintf("      b_tpi                        = %.3f  (set, not tuned)\n", b_tpi))
cat(sprintf("      within-zone SD of tpi        = %.3f\n", tpi_sd_within))
cat(sprintf("      within-zone clinical gradient= %.3f\n", b_tpi * tpi_sd_within))
cat(sprintf("      between-zone frailty contrast= %.3f\n", frailty_contrast))
cat(sprintf("      realized confound ratio      = %.3f  (target %.2f)\n",
            b_tpi * tpi_sd_within / frailty_contrast, CONFOUND_RATIO))

# ---- 5. Closed-form inversion (piecewise age structure) ---------------
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

# ---- 8. Visuals ------------------------------------------------------
grid <- expand.grid(num_xcoord = seq(0,1,length.out=300),
                    num_ycoord = seq(0,1,length.out=300))
grid$g <- g_of_s(grid$num_xcoord, grid$num_ycoord)

p_frailty <- ggplot() +
  geom_point(data = sim_data, aes(num_xcoord, num_ycoord, color = zone),
             size = 1.4, alpha = 0.85) +
  scale_color_manual(values = c(high="red", low="blue"), name="True zone") +
  geom_contour(data = grid, aes(num_xcoord, num_ycoord, z = g),
               breaks = boundary, color = "black", linewidth = 1) +
  labs(title = "True frailty z(s): high west, low east",
       x="Easting", y="Northing") +
  coord_fixed(xlim=c(0,1), ylim=c(0,1)) + theme_minimal()

p_tpi <- ggplot() +
  geom_point(data = sim_data, aes(num_xcoord, num_ycoord, color = num_tpi),
             size = 1.4, alpha = 0.85) +
  scale_color_gradient2(low="blue", mid="gray80", high="red",
                        midpoint = median(sim_data$num_tpi), name="tpi") +
  geom_contour(data = grid, aes(num_xcoord, num_ycoord, z = g),
               breaks = boundary, color = "black", linewidth = 1) +
  labs(title = "tpi: north-south gradient, within-zone balanced",
       x="Easting", y="Northing") +
  coord_fixed(xlim=c(0,1), ylim=c(0,1)) + theme_minimal()

cat("\nPlots: p_frailty, p_tpi.\n")
cat("Objects: sim_data, z_s, in_high, boundary, g_of_s, grid, g_patient,",
    "b_wbc, b_tpi\n")
