# =====================================================================
# FIGURE: the covariate tpi across the two designs
# ---------------------------------------------------------------------
# Panel (A) degenerate design: tpi independent of location.
# Panel (B) exogeneity design: tpi carries a north-south gradient that
#           survives the within-zone centering.
#
# Both generators are sourced into their own environment, so their
# identically-named objects do not collide. Nothing else is needed:
# each generator defines sim_data, grid, boundary and g_of_s.
#
# tpi is standardized within each design before plotting. The two
# designs use different b_tpi, so the raw scales are not comparable;
# standardizing puts both panels in SD units on a shared colour scale.
# =====================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2); library(patchwork); library(scales)

GEN_DEGENERATE <- "sim_degenerate_generate.R"
GEN_EXOGENEITY <- "sim_twozone_generate.R"

run_gen <- function(path) {
  e <- new.env()
  cat("sourcing", path, "...\n")
  sys.source(path, envir = e)
  stopifnot(all(c("sim_data", "grid", "boundary") %in% ls(e)))
  list(dat = e$sim_data, grid = e$grid, boundary = e$boundary)
}

A <- run_gen(GEN_DEGENERATE)
B <- run_gen(GEN_EXOGENEITY)

# standardize tpi within each design
A$dat$tpi_sd <- as.numeric(scale(A$dat$num_tpi))
B$dat$tpi_sd <- as.numeric(scale(B$dat$num_tpi))

# shared symmetric colour limits, clipped at the 2nd/98th percentile so a
# few outliers do not flatten the middle of the scale
lim <- quantile(abs(c(A$dat$tpi_sd, B$dat$tpi_sd)), 0.98)
lim <- round(lim, 1)
cat(sprintf("colour limits: +/- %.1f SD\n", lim))

tpi_panel <- function(P, ttl) {
  ggplot() +
    geom_point(data = P$dat,
               aes(num_xcoord, num_ycoord, color = tpi_sd),
               size = 1.4, alpha = 0.85) +
    geom_contour(data = P$grid, aes(num_xcoord, num_ycoord, z = g),
                 breaks = P$boundary, color = "black", linewidth = 0.9) +
    scale_color_gradient2(low = "blue", mid = "gray90", high = "red",
                          midpoint = 0, limits = c(-lim, lim),
                          oob = squish, name = "tpi (SD)") +
    labs(title = ttl, x = "Easting", y = "Northing") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 11)
}

pA <- tpi_panel(A, "(A) Degenerate design")
pB <- tpi_panel(B, "(B) Exogeneity design")

fig <- (pA | pB) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("fig_tpi_contrast.pdf", fig, width = 11, height = 4.5)
ggsave("fig_tpi_contrast.png", fig, width = 11, height = 4.5, dpi = 300)
cat("Figure written: fig_tpi_contrast.pdf\n")

# ---- numbers for the caption ----------------------------------------
cat("\n=== tpi structure by design ===\n")
for (nm in c("A", "B")) {
  P <- get(nm)
  cat(sprintf("  %s  cor(northing, tpi) = %+.3f | cor(easting, tpi) = %+.3f\n",
              ifelse(nm == "A", "degenerate", "exogeneity"),
              cor(P$dat$num_ycoord, P$dat$num_tpi),
              cor(P$dat$num_xcoord, P$dat$num_tpi)))
  cat(sprintf("     mean tpi by zone: high %+.3f  low %+.3f\n",
              mean(P$dat$num_tpi[P$dat$zone == "high"]),
              mean(P$dat$num_tpi[P$dat$zone == "low"])))
}
