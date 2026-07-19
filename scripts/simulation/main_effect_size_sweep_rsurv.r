# ── main_effect_size_sweep_rsurv.r ───────────────────────────────────────────
# Sweeps beta_main across multiple values using rsurv (true lognormal AFT).
# Outputs CSVs in the same format as main_effect_size_sweep_cox_epistasis.r so
# that plot_ef_vs_effect_size.R can consume them directly.
library(survival)
library(dplyr)
library(tibble)
library(glmnet)
library(caret)
library(survivalsvm)
library(PRROC)
library(R.utils)
library(ggplot2)
library(ranger)
# library(mboost)
library(parallel)
library(rsurv)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("--file=", "", script_args[grep("--file=", script_args)])
if (length(script_file) > 0) {
  root <- normalizePath(file.path(dirname(script_file), "../..")
) } else {
  root <- getwd()
}
devtools::load_all(file.path(root, "sNPDR"), quiet = TRUE)

# ── Parameters ────────────────────────────────────────────────────────────────
n_samples    <- as.integer(Sys.getenv("N_SAMPLES",    "200"))
n_attributes <- as.integer(Sys.getenv("N_ATTR",       "1000"))
n_main       <- as.integer(Sys.getenv("N_MAIN",       "10"))
base_meanlog <- as.numeric(Sys.getenv("BASE_MEANLOG",  as.character(log(100))))
base_sdlog   <- as.numeric(Sys.getenv("BASE_SDLOG",    "0.25"))
max_time     <- as.numeric(Sys.getenv("MAX_TIME",      "125"))
base_seed    <- as.integer(Sys.getenv("BASE_SEED",     "2467"))
snpdr_timeout_sec <- as.numeric(Sys.getenv("SNPDR_TIMEOUT", "7200"))

beta_sweep_str <- Sys.getenv("BETA_MAIN_SWEEP", "")
if (nchar(beta_sweep_str) > 0) {
  beta_sweep <- as.numeric(strsplit(beta_sweep_str, ",")[[1]])
} else {
  beta_sweep <- c(0.00, 0.07, 0.13, 0.20, 0.27, 0.33, 0.40)
}
N_sim <- as.integer(Sys.getenv("N_SIM", "20"))
cat(sprintf("Main effect size sweep (rsurv AFT): beta_main = %s\n",
            paste(beta_sweep, collapse = ", ")))

param_grid <- expand.grid(sim_id = seq_len(N_sim), beta_val = beta_sweep)
n_tasks    <- nrow(param_grid)
n_workers  <- min(n_tasks, max(1L, parallel::detectCores() - 2L))
cat(sprintf("Tasks: %d (%d beta x %d sims), workers: %d\n",
            n_tasks, length(beta_sweep), N_sim, n_workers))

# ── Shared simulation + metric helpers ───────────────────────────────────────
source(file.path(root, "scripts", "simulation", "sim_utils.R"))
drop_cols <- c("time", "status")

# ── Parallel sweep ────────────────────────────────────────────────────────────
results_list <- parallel::mclapply(
  seq_len(n_tasks), mc.cores = n_workers, mc.set.seed = FALSE,
  FUN = function(task_idx) tryCatch({
    beta_val <- param_grid$beta_val[task_idx]
    sim_id   <- param_grid$sim_id[task_idx]
    seed     <- base_seed + (sim_id - 1L) * 997L

    cat(sprintf("\n=== beta_main=%.3g | sim %d/%d | seed=%d ===\n",
                beta_val, sim_id, N_sim, seed))

    sim <- simulate_rsurv_main(n_samples, n_attributes, n_main,
                               beta_val, base_meanlog, base_sdlog, max_time, seed)
    dat             <- sim$data
    functional.vars <- sim$functional_vars

    set.seed(123L + sim_id)
    ev <- sum(dat$status); ce <- sum(dat$status == 0)
    k  <- if (min(ev, ce) >= 10) 10L else max(2L, min(ev, ce))
    k_folds <- if (min(ev, ce) >= k)
      caret::createFolds(dat$status, k = k, list = TRUE, returnTrain = TRUE)
    else
      caret::createFolds(seq_len(nrow(dat)), k = k, list = TRUE, returnTrain = TRUE)

    n <- nrow(dat); p <- ncol(dat) - 2L

    errors <- data.frame(
      method=character(), sim_id=integer(), fold=integer(),
      c_index=numeric(), auprc=numeric(), mrr=numeric(),
      ef_at_1pct=numeric(), ef_at_5pct=numeric(), ef_at_10pct=numeric(),
      precision_at_1pct=numeric(), precision_at_5pct=numeric(), precision_at_10pct=numeric(),
      recall_at_1pct=numeric(), recall_at_5pct=numeric(), recall_at_10pct=numeric(),
      mcc_at_1pct=numeric(), mcc_at_5pct=numeric(), mcc_at_10pct=numeric(),
      runtime=numeric(), n=integer(), p=integer(),
      n_main=integer(), n_int=integer(),
      beta_main=numeric(), beta_int=numeric(),
      censparam=numeric(), lambda=numeric(),
      stringsAsFactors=FALSE)
    top_recs  <- data.frame(method=character(),fold=integer(),rank=integer(),
                             feature=character(),stringsAsFactors=FALSE)
    full_recs <- data.frame(method=character(),fold=integer(),rank=integer(),
                             feature=character(),is_functional=logical(),
                             stringsAsFactors=FALSE)

    append_row <- function(method, fold_idx, mrr_v, auc_v, topk, rt) {
      errors <<- rbind(errors, data.frame(
        method=method, sim_id=sim_id, fold=fold_idx,
        c_index=NA_real_, auprc=auc_v, mrr=mrr_v,
        ef_at_1pct=topk$ef_at_1pct, ef_at_5pct=topk$ef_at_5pct,
        ef_at_10pct=topk$ef_at_10pct,
        precision_at_1pct=topk$precision_at_1pct,
        precision_at_5pct=topk$precision_at_5pct,
        precision_at_10pct=topk$precision_at_10pct,
        recall_at_1pct=topk$recall_at_1pct,
        recall_at_5pct=topk$recall_at_5pct,
        recall_at_10pct=topk$recall_at_10pct,
        mcc_at_1pct=topk$mcc_at_1pct,
        mcc_at_5pct=topk$mcc_at_5pct,
        mcc_at_10pct=topk$mcc_at_10pct,
        runtime=rt, n=n, p=p, n_main=n_main, n_int=0L,
        beta_main=beta_val, beta_int=0.0,
        censparam=NA_real_, lambda=NA_real_,
        stringsAsFactors=FALSE))
    }

    for (fold_idx in seq_along(k_folds)) {
      cat("  Fold", fold_idx, "\n")
      tr_idx     <- k_folds[[fold_idx]]
      train_data <- dat[tr_idx, ]
      test_data  <- dat[-tr_idx, ]
      attr_train <- scale(train_data[, setdiff(names(train_data), drop_cols)])
      attr_test  <- scale(test_data[,  setdiff(names(test_data),  drop_cols)],
                          center = attr(attr_train,"scaled:center"),
                          scale  = attr(attr_train,"scaled:scale"))
      perm       <- sample(ncol(attr_train))
      attr_train <- attr_train[, perm, drop=FALSE]
      attr_test  <- attr_test[,  perm, drop=FALSE]

      ### sNPDR-signed (signed ΔT ~ signed ΔX, best for main effects)
      rt <- system.time({
        mdl_sgn <- tryCatch(withTimeout({
          sNPDR::npdr_surv(
            outcome=c(time_var="time",status_var="status"),
            dataset=train_data, nbd.method="multisurf",
            nbd.metric="manhattan", knn=10, msurf.sd.frac=0.5,
            diff.type="signed", km.impute=TRUE)
        }, timeout=snpdr_timeout_sec, onTimeout="silent"),
        error=function(e){cat("sNPDR-signed err:",conditionMessage(e),"\n");NULL})
      })[3]
      if (!is.null(mdl_sgn) && nrow(mdl_sgn)>0) {
        mdl_sgn$beta <- as.numeric(mdl_sgn$beta)
        mdl_sgn <- mdl_sgn[!is.na(mdl_sgn$beta), ]
        mdl    <- mdl_sgn[order(abs(mdl_sgn$beta), decreasing=TRUE), ]
        mname  <- "sNPDR-signed"
        ranked <- mdl$feature
        idx    <- which(ranked %in% functional.vars)
        full_recs <<- record_full_ranking(full_recs,mname,fold_idx,ranked,functional.vars)
        top_recs  <<- record_top_features(top_recs,mname,fold_idx,head(ranked,50))
        append_row(mname, fold_idx,
                   compute_mrr(ranked,functional.vars),
                   compute_auc(abs(mdl$beta[idx]),abs(mdl$beta[-idx])),
                   compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
        cat(sprintf("    %s MRR=%.3f EF@1%%=%.1f\n", mname,
                    compute_mrr(ranked,functional.vars),
                    compute_topk_metrics(ranked,functional.vars,ncol(attr_train))$ef_at_1pct))
      }
      ### sNPDR-absolute (|ΔT| ~ |ΔX|, best for interaction effects)
      rt <- system.time({
        mdl_abs <- tryCatch(withTimeout({
          sNPDR::npdr_surv(
            outcome=c(time_var="time",status_var="status"),
            dataset=train_data, nbd.method="multisurf",
            nbd.metric="manhattan", knn=10, msurf.sd.frac=0.5,
            diff.type="absolute", km.impute=TRUE)
        }, timeout=snpdr_timeout_sec, onTimeout="silent"),
        error=function(e){cat("sNPDR-absolute err:",conditionMessage(e),"\n");NULL})
      })[3]
      if (!is.null(mdl_abs) && nrow(mdl_abs)>0) {
        mdl_abs$beta <- as.numeric(mdl_abs$beta)
        mdl_abs <- mdl_abs[!is.na(mdl_abs$beta), ]
        mdl    <- mdl_abs[order(mdl_abs$beta, decreasing=TRUE), ]
        mname  <- "sNPDR-absolute"
        ranked <- mdl$feature
        idx    <- which(ranked %in% functional.vars)
        full_recs <<- record_full_ranking(full_recs,mname,fold_idx,ranked,functional.vars)
        top_recs  <<- record_top_features(top_recs,mname,fold_idx,head(ranked,50))
        append_row(mname, fold_idx,
                   compute_mrr(ranked,functional.vars),
                   compute_auc(abs(mdl$beta[idx]),abs(mdl$beta[-idx])),
                   compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
        cat(sprintf("    %s MRR=%.3f EF@1%%=%.1f\n", mname,
                    compute_mrr(ranked,functional.vars),
                    compute_topk_metrics(ranked,functional.vars,ncol(attr_train))$ef_at_1pct))
      }

      ### sNPDR-LASSO (binomial glmnet, joint elastic-net)
      rt <- system.time({
        mdl_lasso <- tryCatch(withTimeout({
          sNPDR::npdr_surv(
            outcome=c(time_var="time",status_var="status"),
            dataset=train_data, nbd.method="relieff",
            nbd.metric="manhattan", knn=10, msurf.sd.frac=0.5,
            regression="binomial", regularize=TRUE, alpha=1)
        }, timeout=snpdr_timeout_sec, onTimeout="silent"),
        error=function(e){cat("sNPDR-LASSO err:",conditionMessage(e),"\n");NULL})
      })[3]
      if (!is.null(mdl_lasso) && nrow(mdl_lasso)>0) {
        mdl_lasso$beta <- as.numeric(mdl_lasso$beta)
        mdl_lasso <- mdl_lasso[!is.na(mdl_lasso$beta), ]
        ranked <- mdl_lasso$feature
        idx    <- which(ranked %in% functional.vars)
        full_recs <<- record_full_ranking(full_recs,"LASSO-sNPDR",fold_idx,ranked,functional.vars)
        top_recs  <<- record_top_features(top_recs,"LASSO-sNPDR",fold_idx,head(ranked,50))
        append_row("LASSO-sNPDR", fold_idx,
                   compute_mrr(ranked,functional.vars),
                   compute_auc(abs(mdl_lasso$beta[idx]),abs(mdl_lasso$beta[-idx])),
                   compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
        cat(sprintf("    LASSO-sNPDR MRR=%.3f EF@1%%=%.1f\n",
                    compute_mrr(ranked,functional.vars),
                    compute_topk_metrics(ranked,functional.vars,ncol(attr_train))$ef_at_1pct))
      }

      ### sReliefF
      # rt <- system.time({
      #   srf <- tryCatch(withTimeout({
      #     sReliefF(
      #       x       = as.data.frame(attr_train),
      #       time    = train_data$time,
      #       status  = train_data$status,
      #       m       = nrow(train_data),
      #       k_pct   = 0.20,
      #       feature_types = rep("continuous", ncol(attr_train)))
      #   }, timeout=snpdr_timeout_sec, onTimeout="silent"),
      #   error=function(e){cat("sReliefF err:",conditionMessage(e),"\n");NULL})
      # })[3]
      # if (!is.null(srf)) {
      #   ranked <- srf$ranking$feature
      #   scores <- srf$scores[ranked]
      #   idx    <- which(ranked %in% functional.vars)
      #   full_recs <<- record_full_ranking(full_recs,"sReliefF",fold_idx,ranked,functional.vars)
      #   top_recs  <<- record_top_features(top_recs,"sReliefF",fold_idx,head(ranked,50))
      #   append_row("sReliefF", fold_idx,
      #              compute_mrr(ranked,functional.vars),
      #              compute_auc(scores[idx],scores[-idx]),
      #              compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
      #   cat(sprintf("    sReliefF MRR=%.3f EF@1%%=%.1f\n",
      #               compute_mrr(ranked,functional.vars),
      #               compute_topk_metrics(ranked,functional.vars,ncol(attr_train))$ef_at_1pct))
      # }

      ### Cox
      rt <- system.time({
        cox_list <- tryCatch(withTimeout({
          lapply(colnames(attr_train), function(x) {
            tryCatch({
              m <- survival::coxph(as.formula(paste0("Surv(time,status)~`",x,"`")),
                                   data=train_data)
              df <- as.data.frame(summary(m)$coefficients)
              df$Feature <- x; df
            }, error=function(e) NULL)
          })
        }, timeout=60, onTimeout="silent"),
        error=function(e){cat("Cox err:",conditionMessage(e),"\n");NULL})
      })[3]
      if (!is.null(cox_list)) {
        cox_list <- Filter(Negate(is.null), cox_list)
        cox.df <- do.call(rbind, cox_list) %>%
          rename(beta="coef", p.value="Pr(>|z|)") %>%
          mutate(p.adj=p.adjust(p.value,"bonferroni",n=ncol(attr_train))) %>%
          arrange(p.value)
        ranked  <- cox.df$Feature
        idx     <- which(cox.df$Feature %in% functional.vars)
        full_recs <<- record_full_ranking(full_recs,"Cox",fold_idx,ranked,functional.vars)
        top_recs  <<- record_top_features(top_recs,"Cox",fold_idx,head(ranked,50))
        append_row("Cox", fold_idx,
                   compute_mrr(ranked,functional.vars),
                   compute_auc(cox.df$beta[idx],cox.df$beta[-idx]),
                   compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
      }

      ### SVM
      rt <- system.time({
        svm_fit <- tryCatch(withTimeout({
          survivalsvm(Surv(train_data$time,train_data$status)~.,
                      data=as.data.frame(attr_train), type="regression",
                      gamma.mu=0.2, opt.meth="quadprog", kernel="lin_kernel")
        }, timeout=60, onTimeout="silent"),
        error=function(e){cat("SVM err:",conditionMessage(e),"\n");NULL})
      })[3]
      if (!is.null(svm_fit)) {
        w       <- as.vector(t(svm_fit$model.fit$SV) %*% svm_fit$model.fit$Beta)
        names(w) <- colnames(attr_train)
        svm.df  <- data.frame(Feature=names(w),Weight=w,stringsAsFactors=FALSE) %>%
          arrange(desc(abs(Weight)))
        ranked  <- svm.df$Feature
        idx     <- which(svm.df$Feature %in% functional.vars)
        full_recs <<- record_full_ranking(full_recs,"SVM",fold_idx,ranked,functional.vars)
        top_recs  <<- record_top_features(top_recs,"SVM",fold_idx,head(ranked,50))
        append_row("SVM", fold_idx,
                   compute_mrr(ranked,functional.vars),
                   compute_auc(svm.df$Weight[idx],svm.df$Weight[-idx]),
                   compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
      }

      ### Ranger — disabled: mclapply+withTimeout+ranger fork-signal conflicts
      # rt <- system.time({
      #   rng <- tryCatch(withTimeout({
      #     ranger::ranger(Surv(time,status)~., data=train_data,
      #                    num.trees=500, mtry=floor(sqrt(p)),
      #                    splitrule="logrank", importance="permutation",
      #                    num.threads=1L)
      #   }, timeout=120, onTimeout="silent"),
      #   error=function(e){cat("Ranger err:",conditionMessage(e),"\n");NULL})
      # })[3]
      # if (!is.null(rng)) {
      #   imp    <- sort(rng$variable.importance, decreasing=TRUE)
      #   ranked <- names(imp)
      #   idx    <- which(ranked %in% functional.vars)
      #   full_recs <<- record_full_ranking(full_recs,"Ranger",fold_idx,ranked,functional.vars)
      #   top_recs  <<- record_top_features(top_recs,"Ranger",fold_idx,head(ranked,50))
      #   append_row("Ranger", fold_idx,
      #              compute_mrr(ranked,functional.vars),
      #              compute_auc(imp[idx],imp[-idx]),
      #              compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
      # }

      ### mboost
      # rt <- system.time({
      #   mb <- tryCatch(withTimeout({
      #     mboost::glmboost(Surv(time,status)~.,
      #                      data=as.data.frame(cbind(as.data.frame(attr_train),
      #                                               time=train_data$time,
      #                                               status=train_data$status)),
      #                      family=mboost::CoxPH(),
      #                      control=mboost::boost_control(mstop=1000,nu=0.1))
      #   }, timeout=120, onTimeout="silent"),
      #   error=function(e){cat("mboost err:",conditionMessage(e),"\n");NULL})
      # })[3]
      # if (!is.null(mb)) {
      #   vi  <- mboost::varimp(mb); vi <- vi[names(vi)!="(Intercept)"]
      #   sc  <- setNames(rep(0,p), colnames(attr_train))
      #   sc[names(vi)] <- as.vector(vi)
      #   mb.df  <- data.frame(Feature=names(sc),Weight=sc,stringsAsFactors=FALSE) %>%
      #     arrange(desc(Weight))
      #   ranked <- mb.df$Feature
      #   idx    <- which(ranked %in% functional.vars)
      #   full_recs <<- record_full_ranking(full_recs,"mboost",fold_idx,ranked,functional.vars)
      #   top_recs  <<- record_top_features(top_recs,"mboost",fold_idx,head(ranked,50))
      #   append_row("mboost", fold_idx,
      #              compute_mrr(ranked,functional.vars),
      #              compute_auc(mb.df$Weight[idx],mb.df$Weight[-idx]),
      #              compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
      # }

      ### LASSO-Cox
      rt <- system.time({
        lasso <- tryCatch(withTimeout({
          X_tr <- as.matrix(attr_train)
          y_tr <- Surv(train_data$time, train_data$status)
          cv   <- glmnet::cv.glmnet(X_tr,y_tr,family="cox",alpha=1,nfolds=5)
          glmnet::glmnet(X_tr,y_tr,family="cox",alpha=1,lambda=cv$lambda.min)
        }, timeout=120, onTimeout="silent"),
        error=function(e){cat("LASSO err:",conditionMessage(e),"\n");NULL})
      })[3]
      if (!is.null(lasso)) {
        coefs <- setNames(as.vector(coef(lasso)), colnames(attr_train))
        coefs[is.na(coefs)] <- 0
        l.df   <- data.frame(Feature=names(coefs),Weight=abs(coefs),
                              stringsAsFactors=FALSE) %>% arrange(desc(Weight))
        ranked <- l.df$Feature; idx <- which(ranked %in% functional.vars)
        full_recs <<- record_full_ranking(full_recs,"LASSO-Cox",fold_idx,ranked,functional.vars)
        top_recs  <<- record_top_features(top_recs,"LASSO-Cox",fold_idx,head(ranked,50))
        append_row("LASSO-Cox", fold_idx,
                   compute_mrr(ranked,functional.vars),
                   compute_auc(l.df$Weight[idx],l.df$Weight[-idx]),
                   compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
      }

      ### Elastic Net-Cox
      # rt <- system.time({
      #   enet <- tryCatch(withTimeout({
      #     X_tr <- as.matrix(attr_train)
      #     y_tr <- Surv(train_data$time, train_data$status)
      #     cv   <- glmnet::cv.glmnet(X_tr,y_tr,family="cox",alpha=0.5,nfolds=5)
      #     glmnet::glmnet(X_tr,y_tr,family="cox",alpha=0.5,lambda=cv$lambda.min)
      #   }, timeout=120, onTimeout="silent"),
      #   error=function(e){cat("Enet err:",conditionMessage(e),"\n");NULL})
      # })[3]
      # if (!is.null(enet)) {
      #   coefs <- setNames(as.vector(coef(enet)), colnames(attr_train))
      #   coefs[is.na(coefs)] <- 0
      #   e.df   <- data.frame(Feature=names(coefs),Weight=abs(coefs),
      #                         stringsAsFactors=FALSE) %>% arrange(desc(Weight))
      #   ranked <- e.df$Feature; idx <- which(ranked %in% functional.vars)
      #   full_recs <<- record_full_ranking(full_recs,"Elastic Net-Cox",fold_idx,ranked,functional.vars)
      #   top_recs  <<- record_top_features(top_recs,"Elastic Net-Cox",fold_idx,head(ranked,50))
      #   append_row("Elastic Net-Cox", fold_idx,
      #              compute_mrr(ranked,functional.vars),
      #              compute_auc(e.df$Weight[idx],e.df$Weight[-idx]),
      #              compute_topk_metrics(ranked,functional.vars,ncol(attr_train)), rt)
      # }

    } # end fold loop

    if (nrow(top_recs)  > 0) { top_recs$beta_val  <- beta_val; top_recs$sim_id  <- sim_id }
    if (nrow(full_recs) > 0) { full_recs$beta_val <- beta_val; full_recs$sim_id <- sim_id }
    list(errors=errors, top_recs=top_recs, full_recs=full_recs)
  }, error = function(e) {
    cat("Task error (task", task_idx, "):", conditionMessage(e), "\n")
    NULL
  })
)

# ── Combine ────────────────────────────────────────────────────────────────────
results_list <- Filter(Negate(is.null), results_list)
errors       <- do.call(rbind, lapply(results_list, `[[`, "errors"))

output_dir   <- file.path(root, "results")
graphics_dir <- file.path(root, "figures")
dir.create(output_dir,   recursive=TRUE, showWarnings=FALSE)
dir.create(graphics_dir, recursive=TRUE, showWarnings=FALSE)

write.csv(errors,
          file.path(output_dir, "main_effect_size_sweep_rsurv_fold_metrics.csv"),
          row.names=FALSE)

# Stage 1: average across folds within each simulation
sim_level <- errors %>%
  group_by(method, n, p, n_main, n_int, beta_main, beta_int, censparam, lambda, sim_id) %>%
  summarise(
    mrr         = mean(mrr,         na.rm=TRUE),
    ef_at_1pct  = mean(ef_at_1pct,  na.rm=TRUE),
    ef_at_5pct  = mean(ef_at_5pct,  na.rm=TRUE),
    ef_at_10pct = mean(ef_at_10pct, na.rm=TRUE),
    auprc       = mean(auprc,       na.rm=TRUE),
    runtime     = mean(runtime,     na.rm=TRUE),
    .groups="drop")

# Stage 2: median/IQR across simulations
results_summary <- sim_level %>%
  group_by(method, n, p, n_main, n_int, beta_main, beta_int, censparam, lambda) %>%
  summarise(
    mean_mrr          = mean(mrr,            na.rm=TRUE),
    sd_mrr            = sd(mrr,              na.rm=TRUE),
    mean_ef_at_1pct   = mean(ef_at_1pct,     na.rm=TRUE),
    sd_ef_at_1pct     = sd(ef_at_1pct,       na.rm=TRUE),
    median_ef_at_1pct = median(ef_at_1pct,   na.rm=TRUE),
    q1_ef_at_1pct     = quantile(ef_at_1pct, 0.25, na.rm=TRUE),
    q3_ef_at_1pct     = quantile(ef_at_1pct, 0.75, na.rm=TRUE),
    mean_ef_at_5pct   = mean(ef_at_5pct,     na.rm=TRUE),
    sd_ef_at_5pct     = sd(ef_at_5pct,       na.rm=TRUE),
    mean_ef_at_10pct  = mean(ef_at_10pct,    na.rm=TRUE),
    sd_ef_at_10pct    = sd(ef_at_10pct,      na.rm=TRUE),
    mean_auprc        = mean(auprc,          na.rm=TRUE),
    sd_auprc          = sd(auprc,            na.rm=TRUE),
    mean_runtime      = mean(runtime,        na.rm=TRUE),
    n_sims            = n(),
    .groups="drop")

print(results_summary)
write.csv(results_summary,
          file.path(output_dir, "main_effect_size_sweep_rsurv_results_summary.csv"),
          row.names=FALSE)

# ── EF@1% line plot ────────────────────────────────────────────────────────────
pal <- c(Cox="#E41A1C", SVM="#FF7F00", Ranger="#6A3D9A",
         "sNPDR-signed"="#1F78B4", "sNPDR-absolute"="#33A02C",
         "LASSO-sNPDR"="#17BECF", "LASSO-Cox"="#F781BF")
shp <- c(Cox=16, SVM=17, Ranger=15,
         "sNPDR-signed"=18, "sNPDR-absolute"=19,
         "LASSO-sNPDR"=6, "LASSO-Cox"=4)
ord <- names(pal)

df <- results_summary
df <- df[df$method %in% ord, ]
df$method <- factor(df$method, levels=ord)
n_sims    <- round(median(results_summary$n_sims, na.rm=TRUE))

p_ef1 <- ggplot(df, aes(x=beta_main, y=median_ef_at_1pct,
                         colour=method, fill=method, group=method, shape=method)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey40", linewidth=0.6) +
  geom_ribbon(aes(ymin=pmax(q1_ef_at_1pct, 0),
                  ymax=q3_ef_at_1pct),
              alpha=0.15, colour=NA) +
  geom_line(linewidth=1.1) + geom_point(size=3) +
  scale_colour_manual(values=pal, name="Method", breaks=ord) +
  scale_fill_manual(values=pal,   name="Method", breaks=ord) +
  scale_shape_manual(values=shp,  name="Method", breaks=ord) +
  scale_x_continuous(breaks=sort(unique(df$beta_main))) +
  scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(.02,.12))) +
  labs(x="Main Effect Size",
       y=sprintf("EF@1%% (median \u00b1 IQR, %d simulations)", n_sims)) +
  theme_bw(base_size=13) +
  theme(panel.grid.minor=element_blank(),
        axis.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13,face="bold"),
        legend.position="right",
        plot.title=element_blank(),
        plot.subtitle=element_blank(),
        panel.border=element_rect(color="black",fill=NA,linewidth=1))

ggsave(file.path(graphics_dir, "main_ef1pct_sweep_rsurv_lineplot.png"),
       plot=p_ef1, width=12, height=7, dpi=300)
cat("Main EF@1% lineplot (rsurv) saved.\n")
cat("All outputs in:", output_dir, "and", graphics_dir, "\n")
