
# test_k_lognormal.R
#
# Sweeps neighborhood size k for BOTH main-effect and interaction detection
# using lognormal AFT simulations, then saves two result CSVs consumed by
# plot_k_combined_lognormal.R.
#
# Usage: Rscript test_k_lognormal.R [sdlog] [p]
#   sdlog  lognormal sigma for interaction sim (default 0.5)
#   p      number of features                  (default 1000)

library(PRROC)
library(rsurv)
library(parallel)

# ── Root detection ────────────────────────────────────────────────────────────
root <- normalizePath(file.path(dirname(normalizePath(
  commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))] |>
    sub("--file=", "", x = _)
)), "../.."), mustWork = FALSE)
if (!file.exists(file.path(root, "sNPDR/DESCRIPTION")))
  root <- getwd()
devtools::load_all(file.path(root, "sNPDR"), quiet = TRUE)
source(file.path(root, "scripts", "simulation", "sim_utils.R"))

# ── Command-line arguments ────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
base_sdlog <- if (length(args) >= 1) as.numeric(args[1]) else 0.5
p          <- if (length(args) >= 2) as.integer(args[2]) else 1000L
p_int      <- if (length(args) >= 3) as.integer(args[3]) else 100L
cat(sprintf("Running with base_sdlog = %.2f, p_main = %d, p_int = %d\n", base_sdlog, p, p_int))

# ── Shared parameters ─────────────────────────────────────────────────────────
set.seed(12345)
num_replicates <- 20
sim_seeds      <- sample(1:10000, num_replicates)

n            <- 200
base_meanlog <- log(100)
max_time     <- 500
k_values     <- c(5, 20, 100, n - 1, n)

n_workers <- min(length(k_values) * num_replicates,
                 max(1L, parallel::detectCores() - 2L))
cat(sprintf("k values: %s | reps: %d | workers: %d\n",
            paste(k_values, collapse = ", "), num_replicates, n_workers))

out_dir <- file.path(root, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helper: run one mclapply sweep ────────────────────────────────────────────
run_sweep <- function(label, sim_fn, snpdr_args, col_name) {
  cat(sprintf("\n── %s k-sweep ──────────────────────────────────\n", label))
  param_grid <- expand.grid(k_idx = seq_along(k_values),
                             rep_i = seq_len(num_replicates))

  res_list <- parallel::mclapply(
    seq_len(nrow(param_grid)), mc.cores = n_workers, mc.set.seed = FALSE,
    FUN = function(task_idx) {
      k <- k_values[param_grid$k_idx[task_idx]]
      i <- param_grid$rep_i[task_idx]

      nbd_method <- if (k == n) "multisurf" else "relieff"
      knn        <- if (k == n) 0L else as.integer(k)

      sim <- sim_fn(i)
      dat <- sim$data

      args_full <- c(list(
        outcome       = c("time_var" = "time", "status_var" = "status"),
        dataset       = dat,
        nbd.method    = nbd_method,
        nbd.metric    = "manhattan",
        knn           = knn,
        msurf.sd.frac = 0.5
      ), snpdr_args)

      mdl    <- do.call(sNPDR::npdr_surv, args_full)
      mdl    <- mdl[!is.na(mdl$beta), ]
      scores <- abs(mdl$beta)
      labels <- as.integer(mdl$feature %in% sim$functional_vars)
      auc    <- PRROC::pr.curve(scores.class0 = scores,
                                weights.class0 = labels,
                                curve = FALSE)$auc.integral

      cat(sprintf("%s  k=%d  rep=%d/%d  auPRC=%.4f\n",
                  label, k, i, num_replicates, auc))
      setNames(data.frame(k, i, auc), c("K", "Replicate", col_name))
    }
  )
  do.call(rbind, res_list)
}

# ── Interaction simulation function ───────────────────────────────────────────
simulate_rsurv_interactions <- function(n, p, n_main, n_int,
                                        beta_main, beta_int,
                                        meanlog, sdlog, max_time, seed) {
  set.seed(seed)
  n_signal_int <- 2 * n_int
  main_names   <- paste0("simvar",   seq_len(n_main))
  int_names    <- paste0("intervar", seq_len(n_signal_int))
  noise_names  <- paste0("noisevar", seq_len(p - n_main - n_signal_int))

  X <- as.data.frame(matrix(rnorm(n * p), nrow = n, ncol = p,
                             dimnames = list(NULL, c(main_names, int_names, noise_names))))

  betas_main <- c(rep( beta_main, ceiling(n_main / 2)),
                  rep(-beta_main, floor(n_main / 2)))
  betas_int  <- rep(c(beta_int, -beta_int), length.out = n_int)
  int_terms  <- paste0("intervar", seq(1, n_signal_int, 2),
                       ":intervar", seq(2, n_signal_int, 2))

  fmla  <- as.formula(paste("~", paste(main_names, collapse = " + "),
                             "+", paste(int_terms,  collapse = " + ")))
  betas <- c(betas_main, betas_int)

  set.seed(seed + 1)
  u       <- if (sdlog == 0) rep(0.5, n) else runif(n)
  t_event <- raftreg(u, formula = fmla, beta = betas, dist = "lnorm",
                     meanlog = meanlog, sdlog = max(sdlog, 1e-9), data = X)

  time   <- pmin(t_event, max_time)
  status <- as.integer(t_event <= max_time)
  list(data            = cbind(X, time = time, status = status),
       functional_vars = c(main_names, int_names))
}

# ── Run main-effect sweep ─────────────────────────────────────────────────────
# Parameters match main_effect_size_sweep_rsurv.r (n_main=10, beta=1.0, sdlog=0.5)
main_results <- run_sweep(
  label      = "Main",
  sim_fn     = function(i) simulate_rsurv_main(
    n, p, n_main = 10, beta = 1.0,
    meanlog = base_meanlog, sdlog = 0.5,
    max_time = max_time, seed = sim_seeds[i]),
  snpdr_args = list(diff.type = "signed"),
  col_name   = "auPRC"
)
main_csv <- file.path(out_dir, "auPRC_K_main_lognormal.csv")
write.csv(main_results, main_csv, row.names = FALSE)
cat("Saved:", main_csv, "\n")

# ── Run interaction sweep ─────────────────────────────────────────────────────
int_results <- run_sweep(
  label      = "Interaction",
  sim_fn     = function(i) simulate_rsurv_interactions(
    n, p_int, n_main = 2, n_int = 10,
    beta_main = 0.5, beta_int = 0.25,
    meanlog = base_meanlog, sdlog = base_sdlog,
    max_time = max_time, seed = sim_seeds[i]),
  snpdr_args = list(diff.type = "absolute", km.impute = TRUE),
  col_name   = "AUC"
)
int_suffix <- sprintf("PR_K_AUC_interactions_lognormal_lmkm_sigma%.2f_p%d.csv", base_sdlog, p_int)
int_csv <- file.path(out_dir, int_suffix)
write.csv(int_results, int_csv, row.names = FALSE)
cat("Saved:", int_csv, "\n")

cat("\nDone. Run plot_k_combined_lognormal.R to generate Figure 7.\n")
