# =====================================================================
# COMBINED RECOVERY FIGURE: four analyses x three designs
# ---------------------------------------------------------------------
# Replaces fig_twozone_2x2_sharedscale, fig_threezone_2x2 and
# fig_fourzone_2x2_sharedscale with one 4x3 grid: analyses down the
# rows, designs across the columns, everything on a common
# relative-risk scale. Three columns rather than four leaves each panel
# a third larger, and the four analyses of a given design stack
# directly above one another.
#
# The usual furniture is stripped: no axis labels, no tick labels, one
# shared legend. Rows and columns are labelled once. The maps are read
# for pattern, not for coordinates.
#
# Reads the saved objects from the three analysis runs. Point SAVE_DIR
# at wherever those .RData files live.
# =====================================================================
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
SAVE_DIR <- "."
library(ggplot2); library(patchwork); library(scales)

FILES <- c(
  "Two zones"   = "sim_twozone_key_objects.RData",
  "Three zones" = "sim_threezone_key_objects.RData",
  "Four zones"  = "sim_fourzone_key_objects.RData")

# ---------------------------------------------------------------------
# 1. Load each design into its own environment
# ---------------------------------------------------------------------
load_design <- function(f) {
  e <- new.env()
  load(file.path(SAVE_DIR, f), envir = e)
  need <- c("sim_data", "z_s", "grid", "pt.sp", "pt.ts", "b_wbc", "b_tpi")
  miss <- setdiff(need, ls(e))
  if (length(miss)) stop(f, " is missing: ", paste(miss, collapse = ", "))
  e
}
D <- lapply(FILES, load_design)
names(D) <- names(FILES)

# ---------------------------------------------------------------------
# 2. Shared scale and boundary contours
# ---------------------------------------------------------------------
cap <- function(x, lo = 0.3, hi = 3) pmin(pmax(x, lo), hi)

# Same mapping as the other map figures in the paper: linear, midpoint
# at a relative risk of 1, clipped to [0.3, 3].
rr_scale <- scale_color_gradient2(
  low = "blue", mid = "gray70", high = "red", midpoint = 1,
  name = "Relative\nRisk", limits = c(0.3, 3), oob = squish)

# Each design stores its boundaries differently: the two-zone design has
# a single level set of g, the three- and four-zone designs have two and
# three masked splits.
boundary_layers <- function(e) {
  g <- e$grid
  L <- list()
  if (!is.null(g$g))
    L <- c(L, list(geom_contour(data = g, aes(num_xcoord, num_ycoord, z = g),
                                breaks = e$boundary, color = "black",
                                linewidth = 0.35)))
  if (!is.null(g$lp1))
    L <- c(L, list(geom_contour(data = g, aes(num_xcoord, num_ycoord, z = lp1),
                                breaks = e$split1_thresh, color = "black",
                                linewidth = 0.35)))
  if (!is.null(g$lp2_masked)) {
    gg <- g[!is.na(g$lp2_masked), ]
    L <- c(L, list(geom_contour(data = gg, aes(num_xcoord, num_ycoord, z = lp2_masked),
                                breaks = e$split2_thresh, color = "black",
                                linewidth = 0.35)))
  }
  if (!is.null(g$lp3_masked)) {
    gg <- g[!is.na(g$lp3_masked), ]
    L <- c(L, list(geom_contour(data = gg, aes(num_xcoord, num_ycoord, z = lp3_masked),
                                breaks = e$split3_thresh, color = "black",
                                linewidth = 0.35)))
  }
  L
}

# ---------------------------------------------------------------------
# 3. One small panel
# ---------------------------------------------------------------------
bare <- theme_minimal(base_size = 8) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA,
                                    linewidth = 0.3),
        plot.title = element_text(size = 8, hjust = 0.5,
                                  margin = margin(b = 2)),
        plot.margin = margin(1, 1, 1, 1))

panel <- function(e, rr, title = NULL, ylab = NULL) {
  d <- e$sim_data; d$rr <- rr
  p <- ggplot(d, aes(num_xcoord, num_ycoord, color = rr)) +
    geom_point(size = 0.55, alpha = 0.85) +
    boundary_layers(e) + rr_scale +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + bare
  if (!is.null(title)) p <- p + labs(title = title)
  if (!is.null(ylab))
    p <- p + labs(y = ylab) +
    theme(axis.title.y = element_text(size = 8.5, angle = 90,
                                      margin = margin(r = 2)))
  p
}

# ---------------------------------------------------------------------
# 4. Build the grid
# ---------------------------------------------------------------------
ROWS <- c("Simulated risk", "Spatial-only", "Two-stage", "spBayesSurv")

# risk on the common relative-risk scale, one list per design
RR <- lapply(D, function(e) {
  eta <- e$b_wbc * e$sim_data$num_wbc + e$b_tpi * e$sim_data$num_tpi
  list(
    cap({ v <- e$z_s * exp(eta); v / median(v) }),
    cap(median(e$pt.sp) / e$pt.sp),
    cap(median(e$pt.ts) / e$pt.ts),
    if (!is.null(e$risk.sb)) cap({ v <- exp(e$risk.sb); v / median(v) })
    else rep(NA_real_, nrow(e$sim_data)))
})

# rows are analyses, columns are designs
rows <- list()
for (r in 1:4) {
  ps <- lapply(seq_along(D), function(k)
    panel(D[[k]], RR[[k]][[r]],
          title = if (r == 1) names(D)[k] else NULL,
          ylab  = if (k == 1) ROWS[r] else NULL))
  rows[[r]] <- Reduce(`|`, ps)
}

fig <- (rows[[1]] / rows[[2]] / rows[[3]] / rows[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.width = unit(1.4, "cm"),
        legend.key.height = unit(0.30, "cm"),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 7.5))

ggsave("fig_recovery_grid.pdf", fig, width = 6.6, height = 7.4)
ggsave("fig_recovery_grid.png", fig, width = 6.6, height = 7.4, dpi = 400)
cat("Figure written: fig_recovery_grid.pdf\n")
cat("Check legibility at 0.95\\textwidth before committing:\n")
cat("  each panel lands near 1.7 x 1.4 in on the page.\n")