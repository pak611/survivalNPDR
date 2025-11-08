

# Testing k for main effects

# Load necessary libraries and your custom functions
root <- getwd()
devtools::load_all(file.path(root, "survival_NPDR/cox_epistasis"))
devtools::load_all(file.path(root, "survival_NPDR/sNPDR"))



# Assuming the true positives are known
true_functional_features <- paste0("intervar", seq(3, 42))  # True functional features for this simulation


# Set a single seed for reproducibility
set.seed(12345)

# Number of replicates
num_replicates <- 20

# Generate a list of seeds for simulations
sim_seed <- sample(1:10000, num_replicates)
n <- 200
p <- 500
n_main <- 2
n_int <- 50
beta_main <- 0.001
beta_int <- 0.002
censparam <- 1/4
lambda <- 1/1000
# initialize the errors dataframe
errors <- data.frame(K = integer(), Replicate = integer(), AUC = numeric())


# Define multiple values of k (neighbors)
k_values <- c(5, 20, 100, n - 1, n)

# Loop over k values
for (k in k_values) {
  
  # Loop for running the model 100 times for each k
  for (i in 1:num_replicates) {

    #browser()

    if (k == n) {
      nbd.method <- "multisurf"
      knn <- NULL  # multisurf does not use fixed k
    } else {
      nbd.method <- "relieff"
      knn <- as.numeric(k)
    }
    

    # Simulate data
    simdata <- simul.int(sim_seed[i], n = n, p = p,
                        n.main = n_main,
                        n.int = n_int,
                        beta.main = beta_main, 
                        beta.int = beta_int, 
                        censparam = censparam, 
                        lambda = lambda)

    dat <- simdata$data


    # Call survival NPDR
    survNPDR.model <- sNPDR::npdr_surv_binomial(
            outcome = c("time_var" = "time", "status_var" = "status"),
            dataset = dat, 
            attr.diff.type = "standard",
            nbd.method = nbd.method, 
            nbd.metric = "manhattan",
            knn = k, 
            msurf.sd.frac = 0.5, 
            glmnet.alpha = 1, 
            model.type = "binomial", 
            KM.weight = FALSE,
            KM.kernel.type = "gaussian",
            KM.kernel.sigma = 1.5)

    #browser()

    
    # Generate AUC
    functional.vars <- grep("intervar", colnames(dat), value = TRUE)
    idx_func <- which(survNPDR.model$Feature %in% functional.vars)
    func_betas <- survNPDR.model$beta[idx_func]
    neg_betas <- survNPDR.model$beta[-idx_func]
    pr_curve_survNPDR <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    auc_survNPDR <- pr_curve_survNPDR$auc.integral

    # Save results
    errors <- rbind(errors, data.frame(
      K = k,
      Replicate = i,
      AUC = auc_survNPDR
    ))

    
    # Print percentage completion
    print(paste("Completed", (i / num_replicates) * 100, "% of the simulations."))

    # Print the value of k
    print(paste("Completed simulations for k =", k))
  }
}

# View the results
print(errors)

# Optionally save to a file
write.csv(errors, file.path(root, "paper_tables", "PR_K_AUC_interactions.csv"), row.names = FALSE)

# Read in the errors data frame from the CSV file
errors <- read.csv(file.path(root, "paper_tables", "PR_K_AUC_interactions.csv"))


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
  coord_cartesian(ylim = c(0.15, 0.50)) +  # Adjust y-axis range to 0.15 to 0.25
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
ggsave(file.path(root, "paper_tables", "AUC_K_interactions_overlay_legend.png"), width = 12, height = 8, dpi = 300)



# Read in interaction errors and main errors and plot them side by side
errors_interaction <- read.csv(file.path(root, "paper_tables", "PR_K_AUC_interactions.csv"))
errors_main <- read.csv(file.path(root, "paper_tables", "PR_K_AUC.csv"))

errors_interaction$Type <- "Interaction"
errors_main$Type <- "Main"

# Combine the two data frames
errors_combined <- rbind(errors_interaction, errors_main)
errors_combined$K <- as.factor(errors_combined$K)  # Ensure K is a factor


# Replace k = 200 with "multisurf" in the errors data frame
errors_combined$K <- as.factor(errors_combined$K)  # Ensure K is a factor

levels(errors_combined$K)[levels(errors_combined$K) == "200"] <- "multisurf"  # Replace 200 with "multisurf"

# Create the AUC boxplot with K indicated in the legend
# Create the AUC boxplot with K indicated in the legend and separated by Type
ggplot(errors_combined, aes(x = K, y = AUC, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black", position = position_dodge(width = 0.8)) +  # Boxplot with dodge for separation
  geom_jitter(aes(color = Type), size = 2, alpha = 0.6, position = position_dodge(width = 0.8)) +  # Jitter points for individual data
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               position = position_dodge(width = 0.8), vjust = -1.5, size = 8, color = "black") +  # Display mean values
  labs(
    x = "K (Neighbors)",
    y = "AUC",
    fill = "Effect Type",  # Legend title for fill
    color = "Effect Type"  # Legend title for jitter points
  ) +
  scale_fill_manual(values = c("Interaction" = "lightgreen", "Main" = "lightblue")) +  # Set custom colors for fill
  scale_color_manual(values = c("Interaction" = "lightgreen", "Main" = "lightblue")) +  # Set custom colors for jitter points
  coord_cartesian(ylim = c(0.0, 0.50)) +  # Adjust y-axis range
  theme_minimal() +
  theme(
    axis.title = element_text(size = 30, face = "bold"),  # Bold and larger axis titles
    axis.text = element_text(size = 30),  # Larger axis text
    legend.position = "top",  # Place legend at the top
    legend.direction = "horizontal",  # Make legend horizontal
    legend.title = element_text(size = 30, face = "bold"),  # Larger legend title
    legend.text = element_text(size = 30),  # Larger legend text
    panel.border = element_rect(color = "black", fill = NA, size = 1.2)  # Add a border around the plot
  )

# Save the plot
ggsave(file.path(root, "paper_tables", "AUC_K_main_vs_interaction.png"), width = 12, height = 15, dpi = 300)

