# =====================================================================
# FIGURE: Nelson-Aalen residuals against tpi, by frailty zone
# ---------------------------------------------------------------------
# Three panels: two-, three- and four-zone designs, on a COMMON residual
# scale so the panels can be read against one another.
#
# Rebuilt from the saved analysis objects rather than from the p_rtpi
# plots stored alongside them, because those were built per design with
# independent y-limits.
#
# Each RData supplies `alldata` with e_resid, num_tpi and zone.
# =====================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2); library(patchwork)

FILES <- c(two   = "sim_twozone_key_objects.RData",
           three = "sim_threezone_key_objects.RData",
           four  = "sim_fourzone_key_objects.RData")

ZLEV <- list(
  two   = c("low", "high"),
  three = c("low", "medium", "high"),
  four  = c("low", "medium_low", "medium_high", "high"))

# Blue (lowest frailty) through red (highest), consistent across panels.
ZCOL <- list(
  two   = c(low = "blue", high = "red"),
  three = c(low = "blue", medium = "orange", high = "red"),
  four  = c(low = "blue", medium_low = "cyan3",
            medium_high = "orange", high = "red"))

get_design <- function(f, lev) {
  e <- new.env(); load(f, envir = e)
  stopifnot("alldata" %in% ls(e))
  d <- e$alldata
  stopifnot(all(c("e_resid", "num_tpi", "zone") %in% names(d)))
  d$zone <- factor(d$zone, levels = lev)
  if (any(is.na(d$zone)))
    stop("zone labels in ", f, " do not match: ",
         paste(unique(e$alldata$zone), collapse = ", "))
  d[, c("e_resid", "num_tpi", "zone")]
}

D <- Map(get_design, FILES, ZLEV[names(FILES)])
names(D) <- names(FILES)

# ---- common residual scale ------------------------------------------
# Clipped at the 99th percentile across all three designs: a handful of
# very large residuals would otherwise flatten every fitted line.
YMAX <- as.numeric(quantile(unlist(lapply(D, `[[`, "e_resid")), 0.99))
YMAX <- ceiling(YMAX)
cat(sprintf("common y-limit: 0 to %.0f (99th pct across designs)\n", YMAX))
for (nm in names(D))
  cat(sprintf("  %-5s zone: n=%d, %.1f%% of residuals above the limit\n",
              nm, nrow(D[[nm]]), 100 * mean(D[[nm]]$e_resid > YMAX)))

panel <- function(nm, ttl) {
  d <- D[[nm]]
  ggplot(d, aes(num_tpi, e_resid, color = zone)) +
    geom_point(size = 0.55, alpha = 0.40) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
    scale_color_manual(values = ZCOL[[nm]], name = "Frailty zone") +
    coord_cartesian(ylim = c(0, YMAX)) +
    labs(title = ttl, x = "tpi", y = "NA residual") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text  = element_text(size = 8))
}

pA <- panel("two",   "(A) Two zones")
pB <- panel("three", "(B) Three zones")
pC <- panel("four",  "(C) Four zones")

fig <- pA | pB | pC
ggsave("fig_resid_vs_tpi_threepanel.pdf", fig, width = 12, height = 4.6)
ggsave("fig_resid_vs_tpi_threepanel.png", fig, width = 12, height = 4.6, dpi = 300)
cat("Figure written: fig_resid_vs_tpi_threepanel.pdf\n")

# ---- fitted slopes, for the caption ---------------------------------
cat("\n=== residual ~ tpi slope by zone (caption numbers) ===\n")
for (nm in names(D)) {
  d <- D[[nm]]
  cat(sprintf("  %s zones:\n", nm))
  for (zz in levels(d$zone)) {
    sub <- d[d$zone == zz, ]
    fit <- lm(e_resid ~ num_tpi, data = sub)
    cat(sprintf("    %-12s slope %+.3f   R2 %.3f   mean resid %.3f\n",
                zz, coef(fit)[2], summary(fit)$r.squared, mean(sub$e_resid)))
  }
}
cat("\ntpi range by design (panels narrow as zone count rises,\n")
cat("because tpi is centered within each zone):\n")
for (nm in names(D))
  cat(sprintf("  %-5s  %.2f to %.2f\n", nm,
              min(D[[nm]]$num_tpi), max(D[[nm]]$num_tpi)))
