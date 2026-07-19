# Statistical analysis of feature selection benchmarks.
#
# Reproduces all method comparisons reported in the paper:
#   - Main effects sweep  (N=200, p=1000, 10 functional features, sigma=0.25)
#   - Interaction effects sweep (N=500, p=1000, 10 interaction pairs, sigma=0.25)
#
# For each sweep the script:
#   1. Loads fold-level metrics and aggregates to simulation-replicate level
#      (the independent unit of inference).
#   2. Runs pre-specified pairwise Wilcoxon signed-rank tests with Holm correction.
#   3. Saves results to CSV alongside the input data.
#
# Usage:
#   Rscript analysis_statistical.R
#
# Outputs (written to paper_tables/):
#   main_effect_sweep_wilcoxon_simlevel.csv
#   interaction_sweep_wilcoxon_simlevel.csv

library(dplyr)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# Works whether sourced interactively or run via Rscript.
script_path <- tryCatch(
  normalizePath(sys.frame(0)$ofile),
  error = function(e) NULL
)
if (is.null(script_path)) {
  # Rscript passes the file via commandArgs
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", args)
  if (length(file_flag) > 0) {
    script_path <- normalizePath(sub("^--file=", "", args[file_flag[1]]))
  }
}
if (!is.null(script_path)) {
  BASE <- file.path(dirname(script_path), "../results")
} else {
  BASE <- file.path(getwd(), "results")
}
BASE <- normalizePath(BASE, mustWork = FALSE)

MAIN_FOLD_METRICS <- file.path(BASE, "main_effect_size_sweep_rsurv_fold_metrics.csv")
INT_FOLD_METRICS  <- file.path(BASE, "interaction_effect_size_sweep_rsurv_fold_metrics.csv")

# ---------------------------------------------------------------------------
# Metrics evaluated
# ---------------------------------------------------------------------------
METRICS <- c("ef_at_1pct", "auprc", "mrr")

# ---------------------------------------------------------------------------
# Pre-specified pairwise contrasts
# These are the scientifically motivated comparisons defined before analysis.
# Correcting only over these (rather than all C(n,2) pairs) preserves power.
#
# Main effects:  Is sNPDR competitive with standard methods?
#                Are Ranger and sNPDR-LM-KM clearly worse?
MAIN_PAIRS <- list(
  c("sNPDR",       "Cox"),
  c("sNPDR",       "SVM"),
  c("Ranger",      "sNPDR"),
  c("Ranger",      "Cox"),
  c("sNPDR-LM-KM", "sNPDR")
)

# Interaction effects:  Does sNPDR-LM-KM outperform each competitor?
INT_PAIRS <- list(
  c("sNPDR-LM-KM", "sNPDR"),
  c("sNPDR-LM-KM", "Cox"),
  c("sNPDR-LM-KM", "SVM"),
  c("sNPDR-LM-KM", "mboost"),
  c("sNPDR-LM-KM", "LASSO-Cox"),
  c("sNPDR-LM-KM", "Elastic Net-Cox")
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

#' Aggregate fold-level metrics to simulation-replicate level.
#'
#' This makes the simulation replicate the unit of inference, removing the
#' pseudo-replication that arises from treating CV folds as independent.
aggregate_to_sim_level <- function(fold_df, beta_col, metrics) {
  fold_df %>%
    group_by(method, sim_id, .data[[beta_col]]) %>%
    summarise(across(all_of(metrics), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
}


#' Apply Holm step-down correction to a vector of p-values.
holm_correction <- function(p_values) {
  p.adjust(p_values, method = "holm")
}


#' Run paired Wilcoxon signed-rank tests for pre-specified method pairs.
#'
#' Holm correction is applied globally across all (beta x metric x pair)
#' comparisons, which is the most conservative defensible approach.
#'
#' @param sim_agg  Simulation-level aggregated data frame.
#' @param beta_col Name of the effect-size column (character).
#' @param metrics  Character vector of metrics to compare.
#' @param pairs    List of length-2 character vectors (method_A, method_B).
#' @return Data frame with one row per (beta, metric, pair).
pairwise_wilcoxon <- function(sim_agg, beta_col, metrics, pairs) {
  betas <- sort(unique(sim_agg[[beta_col]]))
  rows  <- list()

  for (beta in betas) {
    b <- sim_agg[sim_agg[[beta_col]] == beta, ]
    for (metric in metrics) {
      for (pair in pairs) {
        mA <- pair[1]; mB <- pair[2]
        a_vals <- b[b$method == mA, c("sim_id", metric)]
        b_vals <- b[b$method == mB, c("sim_id", metric)]
        common <- intersect(a_vals$sim_id, b_vals$sim_id)
        if (length(common) < 5) next

        a_common <- a_vals[a_vals$sim_id %in% common, ][[metric]]
        b_common <- b_vals[b_vals$sim_id %in% common, ][[metric]]
        diff     <- a_common - b_common

        p_raw <- if (all(diff == 0)) 1.0 else
          wilcox.test(diff, alternative = "two.sided", exact = FALSE)$p.value

        rows[[length(rows) + 1]] <- data.frame(
          beta_col_val = beta,
          metric       = metric,
          method_A     = mA,
          method_B     = mB,
          median_A     = round(median(a_common, na.rm = TRUE), 4),
          median_B     = round(median(b_common, na.rm = TRUE), 4),
          mean_A       = round(mean(a_common,   na.rm = TRUE), 4),
          mean_B       = round(mean(b_common,   na.rm = TRUE), 4),
          p_value      = p_raw,
          n_sims       = length(common),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  result <- do.call(rbind, rows)
  names(result)[names(result) == "beta_col_val"] <- beta_col
  result$p_adj_holm <- holm_correction(result$p_value)
  result
}


print_summary <- function(df, beta_col, label) {
  cat(sprintf("\n%s\n %s\n%s\n", strrep("=", 60), label, strrep("=", 60)))
  sig <- df[df$p_adj_holm < 0.05, ]
  sig <- sig[order(sig$p_adj_holm), ]
  if (nrow(sig) > 0) {
    print(sig[, c(beta_col, "metric", "method_A", "method_B",
                  "median_A", "median_B", "p_value", "p_adj_holm")],
          row.names = FALSE)
  } else {
    cat("  No pairwise comparisons significant after Holm correction.\n")
    best_idx <- which.min(df$p_value)
    cat(sprintf("  Smallest raw p: %.4f (%s vs %s, beta=%.2f, metric=%s)\n",
                df$p_value[best_idx], df$method_A[best_idx], df$method_B[best_idx],
                df[[beta_col]][best_idx], df$metric[best_idx]))
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# ── Main effects ─────────────────────────────────────────────────────────────
cat("\nLoading main effect sweep ...\n")
main_fold <- read.csv(MAIN_FOLD_METRICS, stringsAsFactors = FALSE)
cat(sprintf("  %d fold-level rows, %d sims, %d methods, %d effect sizes\n",
            nrow(main_fold), length(unique(main_fold$sim_id)),
            length(unique(main_fold$method)), length(unique(main_fold$beta_main))))

main_sim <- aggregate_to_sim_level(main_fold, "beta_main", METRICS)

# Pairwise Wilcoxon
main_wilcox <- pairwise_wilcoxon(main_sim, "beta_main", METRICS, MAIN_PAIRS)
write.csv(main_wilcox, file.path(BASE, "main_effect_sweep_wilcoxon_simlevel.csv"),
          row.names = FALSE)
print_summary(main_wilcox, "beta_main", "Main effects — pairwise Wilcoxon (sim-level, Holm)")

# ── Interaction effects ───────────────────────────────────────────────────────
cat("\nLoading interaction sweep ...\n")
int_fold <- read.csv(INT_FOLD_METRICS, stringsAsFactors = FALSE)
cat(sprintf("  %d fold-level rows, %d sims, %d methods, %d effect sizes\n",
            nrow(int_fold), length(unique(int_fold$sim_id)),
            length(unique(int_fold$method)), length(unique(int_fold$beta_int))))

int_sim <- aggregate_to_sim_level(int_fold, "beta_int", METRICS)

# Pairwise Wilcoxon
int_wilcox <- pairwise_wilcoxon(int_sim, "beta_int", METRICS, INT_PAIRS)
write.csv(int_wilcox, file.path(BASE, "interaction_sweep_wilcoxon_simlevel.csv"),
          row.names = FALSE)
print_summary(int_wilcox, "beta_int", "Interaction effects — pairwise Wilcoxon (sim-level, Holm)")

cat("\nDone. Output written to results/\n")
