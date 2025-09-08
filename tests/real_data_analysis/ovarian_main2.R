
library(gridExtra)
library(ggplot2)
library(patchwork)
library(survminer)
library(survivalsvm)

datasets <- c("GSE9891", "GSE32062", "GSE13876")
plot_list <- list()

devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")

for (ds in datasets) {
  dataset <- ds
  load(paste0("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/", dataset, "/", dataset, ".rda"))
  gset <- get(dataset)
  expr_data <- exprs(gset)
  feature_data <- fData(gset)
  probe_to_gene <- feature_data[, c("probeset", "gene")]
  colnames(probe_to_gene) <- c("ID", "gene")
  expr_data_df <- as.data.frame(expr_data)
  expr_data_df$ID <- rownames(expr_data_df)
  expr_data_with_genes <- merge(expr_data_df, probe_to_gene, by = "ID", all.x = TRUE)
  expr_data_with_genes$gene <- make.unique(as.character(expr_data_with_genes$gene))
  rownames(expr_data_with_genes) <- expr_data_with_genes$gene
  expr_data_with_genes <- expr_data_with_genes[, !(names(expr_data_with_genes) %in% c("ID", "gene"))]
  pheno_data <- pData(gset)
  pheno_data <- pheno_data[, c("days_to_death", "vital_status")]
  pheno_data$vital_status <- ifelse(pheno_data$vital_status == "deceased", 1, 0)
  sample_names <- colnames(expr_data_with_genes)
  pheno_data <- pheno_data[match(sample_names, rownames(pheno_data)), ]
  combined_data <- cbind(pheno_data, t(expr_data_with_genes))
  dat <- as.data.frame(combined_data)
  dat$time <- dat$days_to_death
  dat$status <- dat$vital_status
  dat <- dat[, setdiff(names(dat), c("days_to_death", "vital_status"))]
  dat <- dat[complete.cases(dat[, c("time", "status")]), ]
  attr_mat <- dat[, setdiff(names(dat), c("time", "status"))]
  attr_mat <- attr_mat[, sapply(attr_mat, is.numeric)]
  variance_values <- apply(attr_mat, 2, var, na.rm = TRUE)
  top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:15000])
  #top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:100])
  attr_mat <- attr_mat[, top_variance_genes]
  dat <- cbind(attr_mat, time = dat$time, status = dat$status)

  model_names <- c("survNPDR", "Cox", "SVM")

  for (model_name in model_names) {
    if (model_name == "survNPDR") {
      survNPDR <- sNPDR::npdr_surv_binomial(outcome = c("time_var" = "time", "status_var" = "status"),
                                            dataset = dat, attr.diff.type = "standard",
                                            nbd.method = "multisurf", nbd.metric = "manhattan",
                                            knn = 20, msurf.sd.frac = 0.5, glmnet.alpha = 1,
                                            glmnet.lower = -Inf, model.type = "binomial")
      threshold <- quantile(survNPDR$p.adj, 0.90)
      survNPDR <- survNPDR %>% filter(p.adj < threshold) %>% arrange(desc(abs(beta)))
      top_features <- survNPDR$Feature[1:30]
      top_coefficients <- survNPDR$beta[1:30]
    } else if (model_name == "Cox") {
  cox.model <- lapply(colnames(dat |> select(-c(time, status))), function(x) {
    tryCatch({
      formula <- as.formula(paste("Surv(time, status) ~ `", x, "`", sep = ""))
      summary(coxph(formula, data = dat))
    }, error = function(e) NULL)
  })
  
  # Remove NULL elements
  cox.model <- cox.model[!sapply(cox.model, is.null)]
  
  # Check if cox.model is empty
  if (length(cox.model) > 0) {
    results_df <- do.call(rbind, lapply(cox.model, function(s) {
      coef_table <- s$coefficients
      data.frame(Feature = rownames(coef_table)[1], Beta = coef_table[1, "coef"], P.Value = coef_table[1, "Pr(>|z|)"])
    }))
    
    # Ensure results_df is not empty before arranging
    if (!is.null(results_df) && nrow(results_df) > 0) {
      results_df <- results_df %>% arrange(desc(abs(Beta)))
      top_features <- results_df$Feature[1:30]
      top_coefficients <- results_df$Beta[1:30]
    } else {
      message("No valid coefficients found for Cox model.")
      next  # Skip to the next model
    }
  } else {
    message("Cox model returned no valid results.")
    next  # Skip to the next model
  }
} else if (model_name == "SVM") {
      svm.model <- survivalsvm(Surv(time, status) ~ ., data = dat, type = "regression", gamma.mu = 1,
                               opt.meth = "quadprog", kernel = "add_kernel")
      weights <- t(svm.model$model.fit$SV) %*% svm.model$model.fit$Beta
      svm.results <- data.frame("Feature" = svm.model$var.names, beta = weights) %>%
        arrange(desc(abs(beta)))
      top_features <- svm.results$Feature[1:30]
      top_coefficients <- svm.results$beta[1:30]
    }

    existing_features <- intersect(top_features, colnames(dat))
    X <- dat %>%
    select(all_of(existing_features)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

    beta_vals <- top_coefficients[match(existing_features, top_features)]
    risk_scores <- X %*% as.numeric(beta_vals)

    dat$group <- cut(risk_scores, breaks = quantile(risk_scores, probs = c(0, 0.4, 0.6, 1)),
                     labels = c("Low", "Intermediate", "High"), include.lowest = TRUE)
    dat_km <- dat %>% filter(group != "Intermediate")

    km_fit <- survfit(Surv(time, status) ~ group, data = dat_km)
    plot_title <- paste0("KM - ", model_name, " (", dataset, ")")
    km_plot <- ggsurvplot(km_fit, data = dat_km, pval = TRUE,
                          title = plot_title,
                          xlab = "Time (days)", ylab = "Survival Probability",
                          legend.title = "Risk Group", legend.labs = c("Low", "High"),
                          font.main = c(14, "bold"), font.x = c(12), font.y = c(12),
                          font.tickslab = c(10), font.legend = c(10),
                          linetype = "solid", size = 1.0, risk.table = FALSE)

    plot_list[[paste0(dataset, "_", model_name)]] <- km_plot$plot
  }
}

# Arrange and save the 3x3 grid
final_grid <- gridExtra::grid.arrange(grobs = plot_list, nrow = 3, ncol = 3)
ggsave("paper_graphics/ovarian_main2_KMPlot.png",
       plot = final_grid, width = 16, height = 12, dpi = 300)
