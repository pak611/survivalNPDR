# ── sim_utils.R ───────────────────────────────────────────────────────────────
# Shared simulation helpers used by:
#   main_effect_size_sweep_rsurv.r
#   interaction_effect_size_sweep_rsurv.r
#   test_k/test_k_lognormal.R
# Source this file at the top of each script:
#   source(file.path(root, "scripts", "simulation", "sim_utils.R"))
# ─────────────────────────────────────────────────────────────────────────────

library(rsurv)

# ── Data generation ───────────────────────────────────────────────────────────

#' Simulate lognormal AFT survival data with main effects.
#'
#' The first n_main features are causal with alternating +/- beta.
#' All remaining features are pure noise.
simulate_rsurv_main <- function(n, p, n_main, beta, meanlog, sdlog, max_time, seed) {
  set.seed(seed)
  feat_names   <- paste0("simvar", seq_len(p))
  signal_names <- paste0("simvar", seq_len(n_main))
  X <- as.data.frame(matrix(rnorm(n * p), nrow = n, ncol = p,
                             dimnames = list(NULL, feat_names)))
  betas <- c(rep( beta, ceiling(n_main / 2)),
             rep(-beta, floor(n_main / 2)))
  fmla <- as.formula(paste("~", paste(signal_names, collapse = " + ")))
  set.seed(seed + 1)
  t_event <- rsurv::raftreg(runif(n), formula = fmla, beta = betas,
                             dist = "lnorm", meanlog = meanlog, sdlog = sdlog, data = X)
  time   <- pmin(t_event, max_time)
  status <- as.integer(t_event <= max_time)
  list(data             = cbind(X, time = time, status = status),
       functional_vars  = signal_names)
}

# ── Metric helpers ────────────────────────────────────────────────────────────

#' auPRC from signed scores for positive (functional) and negative features.
compute_auc <- function(pos, neg) {
  if (length(pos) > 0 && length(neg) > 0)
    PRROC::pr.curve(scores.class0 = abs(pos), scores.class1 = abs(neg))$auc.integral
  else NA_real_
}

#' Mean Reciprocal Rank: 1 / rank of the first functional feature.
compute_mrr <- function(ranked, functional_set) {
  r <- which(ranked[!is.na(ranked) & nzchar(ranked)] %in% functional_set)[1]
  if (is.na(r)) NA_real_ else 1 / r
}

#' Precision, Recall, Enrichment Factor, and MCC at top-k% cutoffs.
compute_topk_metrics <- function(ranked, functional_set, total, pcts = c(.01, .05, .10)) {
  ranked <- ranked[!is.na(ranked) & nzchar(ranked)]
  n_pos  <- length(functional_set)
  base_r <- n_pos / total
  out    <- list()
  for (pct in pcts) {
    lbl  <- sprintf("%dpct", round(pct * 100))
    k    <- min(max(1, ceiling(total * pct)), length(ranked))
    hits <- sum(ranked[1:k] %in% functional_set)
    prec <- hits / k
    rec  <- if (n_pos > 0) hits / n_pos else NA_real_
    ef   <- if (base_r > 0) prec / base_r else NA_real_
    TP <- hits; FP <- k - hits; FN <- n_pos - hits; TN <- total - n_pos - FP
    denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc   <- if (denom > 0) (TP * TN - FP * FN) / denom else NA_real_
    out[[paste0("precision_at_", lbl)]] <- prec
    out[[paste0("recall_at_",    lbl)]] <- rec
    out[[paste0("ef_at_",        lbl)]] <- ef
    out[[paste0("mcc_at_",       lbl)]] <- mcc
  }
  out
}

# ── Recording helpers ─────────────────────────────────────────────────────────

#' Append a full feature ranking (all features) to a records data frame.
record_full_ranking <- function(recs, m, f, feats, func_set) {
  feats <- feats[!is.na(feats) & nzchar(feats)]
  if (!length(feats)) return(recs)
  rbind(recs, data.frame(method = m, fold = f, rank = seq_along(feats),
                          feature = feats, is_functional = feats %in% func_set,
                          stringsAsFactors = FALSE))
}

#' Append the top-50 features to a records data frame.
record_top_features <- function(recs, m, f, feats) {
  feats <- feats[!is.na(feats) & nzchar(feats)]
  if (!length(feats)) return(recs)
  rbind(recs, data.frame(method = m, fold = f, rank = seq_along(feats),
                          feature = feats, stringsAsFactors = FALSE))
}
