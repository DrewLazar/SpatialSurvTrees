# =====================================================================
# LeukSurv comparison figure, with a study-region context panel
# ---------------------------------------------------------------------
# Panel (a) shows the 24 administrative districts with the 1,043 patient
# locations on top, coloured by censoring. Panels (b)-(d) are the three
# risk maps as before, on their common colour scale.
#
# The district boundaries ship with spBayesSurv as nwengland.bnd, in the
# SAME scaled coordinate system as LeukSurv's xcoord/ycoord -- that is
# how survregbayes matches districts to ICAR frailties -- so the
# polygons drop straight onto the existing axes with no transform.
#
# Needs: BayesX (for read.bnd). Install with
#   install.packages("BayesX")
#
# Run after loading twostage_session.RData, as with the earlier figure
# script: needs alldata, surf, spatial_risk, OUT_DIR.
# =====================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2); library(scales); library(RColorBrewer); library(patchwork)
library(BayesX)

ALL_MODELS   <- "alldata_all_models.rds"
SESSION_FILE <- "twostage_session.RData"   # optional: supplies surf, spatial_risk
RISK_LIMITS <- c(-0.30, 0.30)

# ---------------------------------------------------------------------
# 0. Data
# ---------------------------------------------------------------------
# alldata_all_models.rds carries LeukSurv plus every fitted column, so
# the script runs from that file alone. twostage_session.RData is loaded
# if present: it supplies `surf` and `spatial_risk`, needed ONLY for the
# smooth benchmark panel. Without them the figure is built with three
# panels instead of four.
if (!exists("surf") && file.exists(SESSION_FILE)) {
  load(SESSION_FILE)
  cat("loaded", SESSION_FILE, "\n")
}
if (!exists("alldata")) {
  alldata <- readRDS(ALL_MODELS)
  cat("loaded alldata from", ALL_MODELS, "-", nrow(alldata), "rows\n")
}
need <- c("spatial_pred", "centered_log_rr", "so_pred", "so_centered_lrr",
          "num_xcoord", "num_ycoord", "status")
stopifnot(all(need %in% names(alldata)))

HAVE_SMOOTH <- exists("surf") && exists("spatial_risk")
if (!HAVE_SMOOTH)
  cat("note: surf / spatial_risk not in the environment;",
      "the smooth benchmark panel will be omitted.\n",
      "     load twostage_session.RData first to include it.\n")
if (!exists("OUT_DIR")) OUT_DIR <- "outputs"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)


# ---------------------------------------------------------------------
# 1. District boundaries -> a data frame ggplot can draw
# ---------------------------------------------------------------------
bnd <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))
bnd_df <- do.call(rbind, lapply(seq_along(bnd), function(k) {
  m <- bnd[[k]]
  data.frame(x = m[, 1], y = m[, 2], piece = k)
}))
cat(sprintf("boundary: %d polygons, %d vertices\n", length(bnd), nrow(bnd_df)))
cat(sprintf("  boundary range  E %.3f-%.3f  N %.3f-%.3f\n",
            min(bnd_df$x), max(bnd_df$x), min(bnd_df$y), max(bnd_df$y)))
cat(sprintf("  patient  range  E %.3f-%.3f  N %.3f-%.3f\n",
            min(alldata$num_xcoord), max(alldata$num_xcoord),
            min(alldata$num_ycoord), max(alldata$num_ycoord)))

# Shared extent covering BOTH the district outline and the patients: the
# outline reaches past northing 1, so clipping at 1 would cut the
# northern tip off. All four panels use these limits so the axes match.
XLIM <- range(c(bnd_df$x, alldata$num_xcoord)) + c(-0.02, 0.02)
YLIM <- range(c(bnd_df$y, alldata$num_ycoord)) + c(-0.02, 0.02)
cat(sprintf("  shared extent   E %.3f-%.3f  N %.3f-%.3f\n",
            XLIM[1], XLIM[2], YLIM[1], YLIM[2]))
# If those two ranges disagree the boundary is in a different coordinate
# system and the panel would be misleading -- stop rather than plot it.
stopifnot(min(bnd_df$x) < min(alldata$num_xcoord) + 0.1,
          max(bnd_df$x) > max(alldata$num_xcoord) - 0.1)

# ---------------------------------------------------------------------
# 2. Panel (a): study region, districts, patients by censoring
# ---------------------------------------------------------------------
alldata$obs <- factor(ifelse(alldata$status == 1, "Death", "Censored"),
                      levels = c("Death", "Censored"))
cat(sprintf("censored: %d of %d (%.1f%%)\n",
            sum(alldata$status == 0), nrow(alldata),
            100 * mean(alldata$status == 0)))

white_bg <- theme(panel.background = element_rect(fill = "white", color = NA),
                  plot.background  = element_rect(fill = "white", color = NA))

p_region <- ggplot() +
  geom_polygon(data = bnd_df, aes(x, y, group = piece),
               fill = "grey97", color = "grey55", linewidth = 0.3) +
  geom_point(data = alldata,
             aes(num_xcoord, num_ycoord, color = obs, shape = obs),
             size = 1.1, alpha = 0.75) +
  scale_color_manual(values = c(Death = "grey25", Censored = "#0072B2"),
                     name = NULL) +
  scale_shape_manual(values = c(Death = 16, Censored = 17), name = NULL) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Easting", y = "Northing") +
  coord_fixed(xlim = XLIM, ylim = YLIM) +
  theme_minimal(base_size = 12) + white_bg +
  # legend inside the panel, top right, where the region is empty; kept
  # out of the collected guides so it stays with its own map
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "grey75"),
        legend.margin = margin(3, 5, 3, 5),
        legend.key.size = unit(0.9, "lines"),
        legend.text = element_text(size = 9))

# ---------------------------------------------------------------------
# 3. Panels (b)-(d): the three risk maps, unchanged
# ---------------------------------------------------------------------
risk_colors <- rev(brewer.pal(11, "RdBu"))
shared_fill <- scale_fill_gradientn(colors = risk_colors, limits = RISK_LIMITS,
                                    oob = scales::squish, na.value = "white",
                                    name = "Centered log risk")

base_map <- function(p)
  p + shared_fill + labs(x = "Easting", y = "Northing") +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  theme_minimal(base_size = 12) + white_bg +
  coord_fixed(xlim = XLIM, ylim = YLIM)

point_map <- function(df, col)
  base_map(
    ggplot(df, aes(num_xcoord, num_ycoord, fill = .data[[col]])) +
      geom_point(shape = 21, color = "grey40", stroke = 0.25,
                 size = 2.0, alpha = 0.95))

p_tree   <- point_map(alldata, "centered_log_rr")
p_so     <- point_map(alldata, "so_centered_lrr")
if (HAVE_SMOOTH)
  p_smooth <- base_map(ggplot(surf, aes(easting, northing, fill = risk)) +
                         geom_raster(interpolate = TRUE))

# ---------------------------------------------------------------------
# 4. Assemble
# ---------------------------------------------------------------------
tag <- function(p, t, sub)
  p + labs(tag = t, subtitle = sub) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 11),
        plot.margin = margin(2, 2, 2, 2))

# 2x2 grid. guides = "collect" would sweep up the censoring legend too,
# so the region panel is marked guides = "keep": its legend stays inside
# its own map while the three risk panels share one colour bar.
pa <- tag(p_region, "(a)", "Study region") + plot_layout(guides = "keep")

if (HAVE_SMOOTH) {
  fig <- ((pa | tag(p_smooth, "(b)", "spBayesSurv (smooth)")) /
            (tag(p_tree, "(c)", "Two-stage tree") |
               tag(p_so,   "(d)", "Spatial-only"))) +
    plot_layout(guides = "collect")
  fname <- "leuksurv_comparison_4panel"; H <- 11
} else {
  fig <- (pa | tag(p_tree, "(b)", "Two-stage tree") |
            tag(p_so, "(c)", "Spatial-only")) +
    plot_layout(guides = "collect")
  fname <- "leuksurv_comparison_3panel_region"; H <- 5
}

ggsave(file.path(OUT_DIR, paste0(fname, ".pdf")), fig,
       width = 11, height = H, device = cairo_pdf)
ggsave(file.path(OUT_DIR, paste0(fname, ".png")), fig,
       width = 11, height = H, dpi = 300)
cat("Figure written:", paste0(fname, ".pdf"), "\n")

# ---------------------------------------------------------------------
# 5. Do the recovered boundaries fall near district edges?
#    The Introduction argues that care is delivered through discrete
#    administrative units; this is the cheapest check of that claim.
# ---------------------------------------------------------------------
bpts <- as.matrix(bnd_df[, c("x", "y")])
d_to_border <- apply(as.matrix(alldata[, c("num_xcoord", "num_ycoord")]), 1,
                     function(p) min(sqrt(colSums((t(bpts) - p)^2))))
alldata$d_border <- d_to_border

# patients whose two-stage leaf differs from that of their nearest neighbour
xy  <- as.matrix(alldata[, c("num_xcoord", "num_ycoord")])
D   <- as.matrix(dist(xy)); diag(D) <- Inf
nn  <- apply(D, 1, which.min)
edge <- alldata$spatial_pred != alldata$spatial_pred[nn]

cat("\n=== tree boundaries vs district boundaries ===\n")
cat(sprintf("  patients at a recovered leaf edge: %d\n", sum(edge)))
cat(sprintf("  median distance to a district border: edge %.4f vs interior %.4f\n",
            median(alldata$d_border[edge]), median(alldata$d_border[!edge])))
cat(sprintf("  Mann-Whitney p = %.3g\n",
            wilcox.test(alldata$d_border[edge], alldata$d_border[!edge])$p.value))
cat("  (smaller distance at leaf edges would support the administrative-boundary\n")
cat("   argument of Section 1; a null result is also worth knowing)\n")