# =====================================================================
# metrics.R
#
# Scoring functions shared by the simulation analyses. Extracted verbatim
# from the per-design scripts so that every design is scored identically.
#
# CONVENTION: far-boundary performance is reported as ACCURACY among the
# half of observations farthest from a true boundary, matching Table 2,
# where higher is better in every column except standardized MSE. The
# field is named far_acc rather than fbe, since fbe previously held an
# error in one script and an accuracy in the others.
# =====================================================================


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

# Mean squared distance from the covariate centroid. Used to set the
# Gaussian kernel bandwidth Ksigma.
intvar <- function(data, covs) {
  mu <- colMeans(data[covs])
  mean(rowSums((data[covs] - matrix(mu, nrow(data), length(mu),
                                    byrow = TRUE))^2))
}

# Adjusted Rand index between two partitions.
adj_rand <- function(a, b) {
  tab <- table(a, b); N <- sum(tab); sc <- function(x) sum(choose(x, 2))
  ai <- rowSums(tab); bj <- colSums(tab)
  (sc(as.vector(tab)) - sc(ai)*sc(bj)/choose(N,2)) /
    (0.5*(sc(ai) + sc(bj)) - sc(ai)*sc(bj)/choose(N,2))
}

# AUC via the rank-sum identity, orientation-free: a score that ranks the
# zones perfectly backwards scores the same as one that ranks them
# perfectly forwards, since sign is not identified by the tree.
auc_score <- function(s, th) {
  r <- rank(s); n1 <- sum(th); n0 <- sum(!th)
  a <- (sum(r[th]) - n1*(n1+1)/2) / (n1*n0); max(a, 1-a)
}

# 5th-to-95th percentile span, the recovered-contrast measure.
span <- function(v) as.numeric(diff(quantile(v, c(0.05, 0.95))))

# Predicted times to centred log relative risk, the common scale on which
# span ratios compare like with like against the true frailty span.
lrr <- function(pt) { v <- log(median(pt) / pt); v - mean(v) }

# Standardize a risk score and orient it to agree in sign with the
# standardized true frailty zc.
std_to <- function(x, zc) {
  x <- scale(x)[, 1]
  if (cor(x, zc) < 0) x <- -x
  x
}

# Binary call: 2-means on the risk score, split at the midpoint of the
# two centres, then flipped if it disagrees with truth more than half the
# time (tree risk has no identified sign).
thr_kmeans <- function(risk) {
  mean(kmeans(risk, centers = 2, nstart = 10)$centers)
}

call_from_thr <- function(risk, thr, th) {
  ch <- risk > thr
  if (mean(ch == th) < 0.5) ch <- !ch
  ch
}


# ---------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------

# Two-zone designs. truth_high is logical; far is a logical vector marking
# the half of observations farthest from the true boundary; zc is the
# standardized true frailty; span_true is span(log(z_s / median(z_s))).
score_binary <- function(risk, risk_lrr, name,
                         truth_high, far, zc, span_true,
                         verbose = TRUE) {
  ck  <- call_from_thr(risk, thr_kmeans(risk), truth_high)
  acc <- max(mean(ck == truth_high), 1 - mean(ck == truth_high))
  out <- list(call    = ck,
              acc     = acc,
              far_acc = mean((ck == truth_high)[far]),
              ari     = adj_rand(ck, truth_high),
              mse     = mean((std_to(risk, zc) - zc)^2),
              span    = span(risk_lrr) / span_true,
              auc     = auc_score(risk, truth_high))
  if (verbose)
    cat(sprintf("  %-16s acc=%.3f far=%.3f ARI=%.3f MSE=%.3f span=%.3f AUC=%.3f\n",
                name, out$acc, out$far_acc, out$ari, out$mse,
                out$span, out$auc))
  out
}

# Three- and four-zone designs. k-means on the risk score with clusters
# relabelled in order of mean risk, so they align with low < ... < high.
# AUC is not defined for a multiclass partition and is returned as NA.
score_multiclass <- function(risk, risk_lrr, name,
                             truth_k, far, zc, span_true, k,
                             verbose = TRUE) {
  km    <- kmeans(risk, centers = k, nstart = 20)$cluster
  ord   <- order(tapply(risk, km, mean))
  remap <- match(km, ord)
  out <- list(call    = remap,
              acc     = mean(remap == truth_k),
              far_acc = mean((remap == truth_k)[far]),
              ari     = adj_rand(remap, truth_k),
              mse     = mean((std_to(risk, zc) - zc)^2),
              span    = span(risk_lrr) / span_true,
              auc     = NA_real_)
  if (verbose)
    cat(sprintf("  %-16s acc=%.3f far=%.3f ARI=%.3f MSE=%.3f span=%.3f\n",
                name, out$acc, out$far_acc, out$ari, out$mse, out$span))
  out
}


# ---------------------------------------------------------------------
# Table output, in the column order of Table 2
# ---------------------------------------------------------------------

table_header <- function() {
  cat(sprintf("%-22s %8s %9s %8s %9s %8s %8s\n",
              "Method", "Acc", "Far-bdry", "ARI", "Std.MSE", "Span", "AUC"))
}

table_row <- function(s, nm) {
  auc <- if (is.na(s$auc)) "     ---" else sprintf("%8.3f", s$auc)
  cat(sprintf("%-22s %8.3f %9.3f %8.3f %9.3f %8.2f %s\n",
              nm, s$acc, s$far_acc, s$ari, s$mse, s$span, auc))
}
