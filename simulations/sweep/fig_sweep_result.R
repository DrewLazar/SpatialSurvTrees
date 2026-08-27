# =====================================================================
# ROBUSTNESS SWEEP: collate and plot
# ---------------------------------------------------------------------
# Standalone. Reads the per-replicate .rds files directly, so it needs
# nothing from the fitting scripts -- point SWEEP_DIR at the downloaded
# folder and run.
#
# Produces
#   sweep_all.rds      every replicate, one row each
#   sweep_summary.csv  level means, standard errors and counts
#   fig_sweep_result.pdf
#   a block of numbers for the Section 6.3 text
#
# The two arms are disjoint: tree columns are NA in spBayes rows and
# vice versa, so every aggregate is taken with na.rm and the counts are
# reported alongside so a thinned level is visible.
# =====================================================================

SWEEP_DIR <- "sweep_out"
library(ggplot2); library(patchwork)

# ---- 1. collate ------------------------------------------------------
# The collated frame is committed as sweep_all.rds, so the figure can be
# rebuilt without the per-replicate files. If sweep_out/ is present the
# replicates are recollated from it instead.
f <- list.files(SWEEP_DIR, pattern = "^rep_.*\\.rds$", full.names = TRUE)
if (length(f)) {
  res2 <- do.call(rbind, lapply(f, readRDS))
  saveRDS(res2, "sweep_all.rds")
} else if (file.exists("sweep_all.rds")) {
  res2 <- readRDS("sweep_all.rds")
  cat("using committed sweep_all.rds\n")
} else {
  stop("no replicate files in ", SWEEP_DIR, " and no sweep_all.rds")
}

cat(sprintf("collected %d rows: %d tree, %d spBayes\n",
            nrow(res2), sum(!res2$fit_sb), sum(res2$fit_sb)))
cat("replicates per level:\n")
print(table(round(res2$confound, 3), ifelse(res2$fit_sb, "spBayes", "tree")))

# ---- 2. summarise ----------------------------------------------------
METRICS <- c("ari_ts","ari_sp","ari_sb", "acc_ts","acc_sp","acc_sb",
             "auc_hot_ts","auc_hot_sp","auc_hot_sb",
             "cor_tpi_risk_ts","cor_tpi_risk_sp","cor_tpi_risk_sb",
             "cor_tpi_resid","cor_tpi_hot","cens_rate")
METRICS <- intersect(METRICS, names(res2))

lv <- sort(unique(res2$confound))
S <- data.frame(confound = lv)
for (m in METRICS) {
  S[[m]]              <- sapply(lv, function(L) mean(res2[[m]][res2$confound == L], na.rm = TRUE))
  S[[paste0(m,"_se")]] <- sapply(lv, function(L) {
    v <- res2[[m]][res2$confound == L]; v <- v[!is.na(v)]
    if (length(v) < 2) NA_real_ else sd(v)/sqrt(length(v)) })
  S[[paste0(m,"_n")]]  <- sapply(lv, function(L) sum(!is.na(res2[[m]][res2$confound == L])))
}
write.csv(S, "sweep_summary.csv", row.names = FALSE)

# x axis is the REALIZED coupling, averaged over replicates at each level
S$x <- S$cor_tpi_hot

# ---- 3. plotting helpers --------------------------------------------
COL <- c("Two-stage" = "#1b6ca8", "Spatial-only" = "#c0392b",
         "spBayesSurv" = "#7d3c98", "Two-stage residuals" = "#148f77")

# majority-zone baseline: zones are 25/50/25, the majority being the
# middle one, so a map that resolves neither extreme scores 0.50
MAJ <- 0.50

band <- function(nm, m, se, lab) {
  data.frame(x = S$x, y = S[[m]], lo = S[[m]] - S[[se]],
             hi = S[[m]] + S[[se]], method = lab)
}

panel <- function(dd, ylab, ttl, hline = NULL, hlab = NULL) {
  p <- ggplot(dd, aes(x, y, color = method, fill = method))
  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline, linetype = 2, color = "grey45")
  p <- p +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.4) +
    scale_color_manual(values = COL, name = NULL) +
    scale_fill_manual(values = COL, guide = "none") +
    labs(x = "Confounding, cor(tpi, hot zone)", y = ylab, title = ttl) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(size = 11))
  if (!is.null(hline) && !is.null(hlab))
    p <- p + annotate("text", x = min(dd$x), y = hline, label = hlab,
                      hjust = 0, vjust = -0.6, size = 2.9, color = "grey35")
  p
}

# ---- 4. panels -------------------------------------------------------
dA <- rbind(band("ari_ts","ari_ts","ari_ts_se","Two-stage"),
            band("ari_sp","ari_sp","ari_sp_se","Spatial-only"),
            band("ari_sb","ari_sb","ari_sb_se","spBayesSurv"))
pA <- panel(dA, "Adjusted Rand index", "(A) Global partition")

dB <- rbind(band("acc_ts","acc_ts","acc_ts_se","Two-stage"),
            band("acc_sp","acc_sp","acc_sp_se","Spatial-only"),
            band("acc_sb","acc_sb","acc_sb_se","spBayesSurv"))
pB <- panel(dB, "Zone-recovery accuracy", "(B) Accuracy",
            hline = MAJ, hlab = "majority-zone baseline")

dC <- rbind(band("auc_hot_ts","auc_hot_ts","auc_hot_ts_se","Two-stage"),
            band("auc_hot_sp","auc_hot_sp","auc_hot_sp_se","Spatial-only"),
            band("auc_hot_sb","auc_hot_sb","auc_hot_sb_se","spBayesSurv"))
pC <- panel(dC, "Hotspot AUC", "(C) Confounded zone only")

dD <- rbind(band("cor_tpi_risk_ts","cor_tpi_risk_ts","cor_tpi_risk_ts_se","Two-stage"),
            band("cor_tpi_risk_sp","cor_tpi_risk_sp","cor_tpi_risk_sp_se","Spatial-only"),
            band("cor_tpi_risk_sb","cor_tpi_risk_sb","cor_tpi_risk_sb_se","spBayesSurv"),
            band("cor_tpi_resid","cor_tpi_resid","cor_tpi_resid_se","Two-stage residuals"))
pD <- panel(dD, "cor(tpi, recovered risk)", "(D) tpi contamination")

fig <- (pA | pB) / (pC | pD) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("fig_sweep_result.pdf", fig, width = 11, height = 8.5)
ggsave("fig_sweep_result.png", fig, width = 11, height = 8.5, dpi = 300)
cat("Figure written: fig_sweep_result.pdf\n")

# ---- 5. numbers for the text ----------------------------------------
cat("\n=========== NUMBERS FOR SECTION 6.3 ===========\n")
cat("\nrealized coupling by level:\n")
print(round(S$x, 3))

show <- function(m, lab) {
  cat(sprintf("\n%s:\n", lab))
  cat(sprintf("  first %.3f (se %.3f)   last %.3f (se %.3f)   n per level %d\n",
              S[[m]][1], S[[paste0(m,"_se")]][1],
              S[[m]][nrow(S)], S[[paste0(m,"_se")]][nrow(S)],
              S[[paste0(m,"_n")]][1]))
  print(round(S[[m]], 3))
}
for (m in c("ari_ts","ari_sp","ari_sb","acc_ts","acc_sp","acc_sb",
            "auc_hot_ts","auc_hot_sp","auc_hot_sb",
            "cor_tpi_risk_ts","cor_tpi_risk_sp","cor_tpi_risk_sb","cor_tpi_resid"))
  show(m, m)

# two-stage advantage on ARI, in standard errors, at each level
cat("\ntwo-stage minus spatial-only ARI, in SEs of the difference:\n")
d  <- S$ari_ts - S$ari_sp
sd_ <- sqrt(S$ari_ts_se^2 + S$ari_sp_se^2)
print(round(d / sd_, 2))

# accuracy against the majority baseline
cat("\naccuracy vs majority baseline (0.50), one-sample t on replicates:\n")
for (L in lv) {
  a <- res2$acc_sp[res2$confound == L]; a <- a[!is.na(a)]
  b <- res2$acc_ts[res2$confound == L]; b <- b[!is.na(b)]
  cat(sprintf("  cf=%.2f  spatial-only p=%.3f   two-stage p=%.3f\n", L,
              t.test(a, mu = MAJ)$p.value, t.test(b, mu = MAJ)$p.value))
}

# hotspot AUC crossing between two-stage and spatial-only
cat("\nhotspot AUC, two-stage minus spatial-only:\n")
dd <- S$auc_hot_ts - S$auc_hot_sp
print(round(setNames(dd, round(S$x, 3)), 3))
k <- which(dd[-1] * dd[-length(dd)] < 0)[1]
if (!is.na(k)) {
  x0 <- S$x[k] + (S$x[k+1] - S$x[k]) * dd[k] / (dd[k] - dd[k+1])
  cat(sprintf("  curves cross at realized coupling %.3f\n", x0))
} else cat("  no crossing in the swept range\n")
