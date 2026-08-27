# =====================================================================
# FIGURE: the confounding mechanism
# ---------------------------------------------------------------------
# Two panels showing tpi across the study region at the low and high
# ends of the sweep. At low confounding tpi is mean-balanced across the
# three frailty zones; at high confounding it is elevated in the hot
# zone specifically, which it therefore predicts, while the other two
# zones are untouched.
#
# The realized cor(tpi, hot-zone membership) is printed in each panel
# title, so the figure reports the confounding it actually achieved
# rather than the dialled parameter.
#
# Both panels use the SAME seed, so the only difference between them is
# the hotspot shift.
# =====================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2); library(patchwork); library(scales)
source("sim_sweep_generate.R")
SEED     <- 1
CONF_LOW <- 0.0     # first level of CONF_GRID
CONF_HI  <- 2.5     # last level of CONF_GRID

A <- make_sweep_data(confound = CONF_LOW, seed = SEED, verbose = TRUE)
B <- make_sweep_data(confound = CONF_HI,  seed = SEED, verbose = TRUE)

cat(sprintf("\nrealized cor(tpi, hot): low %.3f | high %.3f\n",
            A$cor_tpi_hot, B$cor_tpi_hot))

# ---- boundary contours ----------------------------------------------
grid <- expand.grid(num_xcoord = seq(0, 1, length.out = 400),
                    num_ycoord = seq(0, 1, length.out = 400))
grid$lp1 <- lp1_of_s(grid$num_xcoord, grid$num_ycoord)
grid$lp2 <- lp2_of_s(grid$num_xcoord, grid$num_ycoord)
grid$lp2_masked <- ifelse(grid$lp1 < A$split1_thresh, grid$lp2, NA)

bl <- list(
  geom_contour(data = grid, aes(num_xcoord, num_ycoord, z = lp1),
               breaks = A$split1_thresh, color = "black", linewidth = 0.9),
  geom_contour(data = grid[!is.na(grid$lp2_masked), ],
               aes(num_xcoord, num_ycoord, z = lp2_masked),
               breaks = A$split2_thresh, color = "black", linewidth = 0.9))

# ---- shared colour scale across both panels -------------------------
tpi_all <- c(A$sim_data$num_tpi, B$sim_data$num_tpi)
LIM <- as.numeric(quantile(abs(tpi_all), 0.98))
cat(sprintf("colour limits: +/- %.2f\n", LIM))

panel <- function(gen, ttl) {
  d <- gen$sim_data
  d$hot <- gen$is_hot
  ggplot(d, aes(num_xcoord, num_ycoord, color = num_tpi)) +
    # hot-zone points are outlined so the confounded region is visible
    geom_point(data = subset(d, hot), aes(num_xcoord, num_ycoord),
               color = "black", size = 2.6, inherit.aes = FALSE) +
    geom_point(size = 1.6, alpha = 0.9) +
    bl +
    scale_color_gradient2(low = "blue", mid = "gray90", high = "red",
                          midpoint = 0, limits = c(-LIM, LIM),
                          oob = squish, name = "tpi") +
    labs(title = ttl, x = "Easting", y = "Northing") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 11)
}

pA <- panel(A, sprintf("(A) Low confounding  [cor(tpi,high) = %.2f]",  A$cor_tpi_hot))
pB <- panel(B, sprintf("(B) High confounding  [cor(tpi,high) = %.2f]", B$cor_tpi_hot))

fig <- (pA | pB) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("fig_contamination_sidebyside.pdf", fig, width = 11, height = 4.6)
ggsave("fig_contamination_sidebyside.png", fig, width = 11, height = 4.6, dpi = 300)
cat("Figure written: fig_contamination_sidebyside.pdf\n")

# ---- numbers for the caption ----------------------------------------
cat("\n=== mean tpi by zone ===\n")
for (nm in c("low confounding", "high confounding")) {
  gen <- if (nm == "low confounding") A else B
  cat(sprintf("  %s:\n", nm))
  print(round(tapply(gen$sim_data$num_tpi, gen$sim_data$zone, mean), 3))
}
cat(sprintf("\nhot zone: %s (n = %d)\n",
            unique(as.character(A$sim_data$zone[A$is_hot])), sum(A$is_hot)))
cat(sprintf("hot-zone tpi shift, high panel: %+.3f\n",
            mean(B$sim_data$num_tpi[B$is_hot]) - mean(A$sim_data$num_tpi[A$is_hot])))
