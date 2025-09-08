

# Testing k for main effects

# Load necessary libraries and your custom functions
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/simSurvData/cox_epistasis")
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")



# Assuming the true positives are known
true_functional_features <- paste0("simvar", 1:50)  # True functional features for this simulation

# Number of replicates
num_replicates <- 20

# Generate a list of random seeds based on a single seed
set.seed(2467)
sim_seeds <- sample(1:10000, num_replicates)
n <- 200
p <- 500
n_main <- 50
n_int <- 2
beta_main <- 0.7
beta_int <- 0.01
censparam <- 1/4
lambda <- 1/1000

# initialize the errors dataframe
errors <- data.frame(K = integer(), Replicate = integer(), AUC = numeric())


# Define multiple values of k (neighbors)
k_values <- c(5, 20, 100, n - 1, n)


for (k in k_values) {
  for (i in 1:num_replicates) {

    browser()

    # Simulate data
    simdata <- simul.int(sim_seeds[i], n = n, p = p,
                         n.main = n_main,
                         n.int = n_int,
                         beta.main = beta_main, 
                         beta.int = beta_int, 
                         censparam = censparam, 
                         lambda = lambda)

    dat <- simdata$data

    # Choose method based on whether k is numeric or "adaptive"
    if (k == n) {
      nbd_method <- "multisurf"
      knn <- NULL  # multisurf does not use fixed k
    } else {
      nbd_method <- "relieff"
      knn <- as.numeric(k)
    }

    # Call survival NPDR
    survNPDR.model <- sNPDR::npdr_surv_binomial(
      outcome = c("time_var" = "time", "status_var" = "status"),
      dataset = dat, 
      attr.diff.type = "standard",
      nbd.method = nbd_method, 
      nbd.metric = "manhattan",
      knn = k,  # NULL for multisurf
      msurf.sd.frac = 0.5, 
      glmnet.alpha = 1, 
      model.type = "binomial", 
      KM.weight = FALSE,
      KM.kernel.type = "gaussian",
      KM.kernel.sigma = 1.5)

    # Generate AUC
    functional.vars <- grep("simvar", colnames(dat), value = TRUE)
    idx_func <- which(survNPDR.model$Feature %in% functional.vars)
    func_betas <- survNPDR.model$beta[idx_func]
    neg_betas <- survNPDR.model$beta[-idx_func]
    pr_curve_survNPDR <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    auc_survNPDR <- pr_curve_survNPDR$auc.integral

    # Save results
    errors <- rbind(errors, data.frame(
      K = as.character(k),  # store "adaptive" as a string
      Replicate = i,
      AUC = auc_survNPDR
    ))

    print(paste("Completed", (i / num_replicates) * 100, "% of the simulations."))
    print(paste("Completed simulations for k =", k))
    print(paste("AUC for k =", k, "is", auc_survNPDR))
  }
}
# View the results
print(errors)

# Optionally save to a file
write.csv(errors, "paper_tables/PR_K_AUC.csv", row.names = FALSE)


# Read in the PR_K_AUC.csv file
errors <- read.csv("paper_tables/PR_K_AUC.csv")


# Replace k = 200 with "multisurf" in the errors data frame
errors$K <- as.factor(errors$K)  # Ensure K is a factor
levels(errors$K)[levels(errors$K) == "200"] <- "multisurf"  # Replace 200 with "multisurf"

# Create the AUC boxplot with K indicated in the legend
ggplot(errors, aes(x = K, y = AUC, fill = K)) +
  geom_boxplot(outlier.shape = NA, alpha = 1.0, color = "black") +  # Added black border to boxes
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -2.5, size = 10, color = "black") +  # Adjusted vjust to bring means even higher
  labs(
    title = NULL,
    x = NULL,  # Remove x-axis label
    y = "AUC",
    fill = "K"  # Legend title
  ) +
  scale_fill_brewer(palette = "Set2") +  # Updated color scheme
  coord_cartesian(ylim = c(0.0, 1.0)) +  # Adjust y-axis range to 0.15 to 0.25
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 36, face = "bold"),  # Increased axis title size
    axis.text = element_text(size = 34),  # Increased axis text size
    axis.text.x = element_blank(),  # Remove X-axis text
    legend.position = "top",  # Overlay legend at 80% width and height of the plot
    legend.direction = "horizontal",
    legend.title = element_text(size = 35, face = "bold"),  # Larger legend title
    legend.text = element_text(size = 32),  # Larger legend text
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

# Save the plot
ggsave("paper_tables/AUC_K_main_overlay_legend.png", width = 12, height = 8, dpi = 300)

# Save the plot
ggsave("paper_tables/AUC_K_main.png", width = 12, height = 8, dpi = 300)
