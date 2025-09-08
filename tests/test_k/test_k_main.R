

# Testing k for main effects

# Load necessary libraries and your custom functions
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/simSurvData/cox_epistasis")
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")



# Assuming the true positives are known
true_functional_features <- paste0("simvar", 1:50)  # True functional features for this simulation

# Initialize a data frame to store precision, recall, and k for each replicate
results_df <- data.frame()

# Number of replicates
num_replicates <- 20

# Parameters for the simulation
ran_seed <- 8675309
num_features <- 500
num_samples <- 200
num_functional <- 50
censor <- 0.7

# Define multiple values of k (neighbors)
k_values <- c(5, 20, 100, num_samples - 1) 

# Loop over k values
for (k in k_values) {
  
  # Loop for running the model 100 times for each k
  for (i in 1:num_replicates) {

    #browser()
    
    # Simulate a new dataset for each replicate
    betas <- withr::with_seed(seed = ran_seed, 
                              c(runif(num_functional, min = 0.25, max = 0.55),
                                runif(num_features - num_functional, min = 0.0001, max = 0.001)))

    sim_seed <- 2468 + i  # Ensure a different seed for each replicate
    simdata <- withr::with_seed(seed = sim_seed,
                                sim.survdata(N = num_samples, T = 20, num.data.frames = 1, 
                                             covariate = 1:num_functional, low = 0, high = 1, 
                                             beta = betas, interactions = FALSE, xvars = num_features, mu = 0,
                                             sd = c(rep(0.1, num_functional), rep(1, num_features - num_functional)),
                                             censor = censor))
    
    # Rename functional variables for readability
    new_names <- paste0("simvar", 1:num_functional)
    old_names <- paste0("X", 1:num_functional)
    names(old_names) <- new_names

    # Clean up the dataset by renaming variables
    dat <- simdata$data |> 
      mutate(status = ifelse(failed == TRUE, 1, 0)) |> 
      rename(time = "y") |> 
      select(all_of(c("time", "status", paste0("X", 1:num_features)))) |> 
      rename(all_of(old_names))


    survNPDR.glm.model <- sNPDR::npdr_surv_binomial_glm(
            outcome = c("time_var" = "time", "status_var" = "status"),
            dataset = dat, 
            attr.diff.type = "standard",
            nbd.method = "relieff", 
            nbd.metric = "manhattan",
            knn = k, 
            msurf.sd.frac = 0.5, 
            glmnet.alpha = 0.9, 
            glmnet.lower = -Inf,
            glmnet.lam = "lambda.1se",
            use.glmnet = TRUE,
            model.type = "binomial", 
            KM.weight = FALSE,
            KM.kernel.type = "gaussian",
            KM.kernel.sigma = 1.0)
    
    # Extract features based on the top 10 beta coefficients
    threshold <- quantile(survNPDR.glm.model$beta, 0.85)

    selected_features <- survNPDR.glm.model$Feature[survNPDR.glm.model$beta > threshold]

    #browser()

    
    # Calculate true positives, false positives, false negatives
    true_positives <- sum(selected_features %in% true_functional_features)
    false_positives <- sum(!selected_features %in% true_functional_features)
    false_negatives <- sum(!true_functional_features %in% selected_features)
    
    # Calculate precision and recall
    precision <- true_positives / (true_positives + false_positives)
    recall <- true_positives / (true_positives + false_negatives)
    
    # Store the results in the dataframe
    results_df <- rbind(
      results_df, 
      data.frame(Replicate = i, K = k, Precision = precision, Recall = recall))

    
    # Print percentage completion
    print(paste("Completed", (i / num_replicates) * 100, "% of the simulations."))

    # Print the value of k
    print(paste("Completed simulations for k =", k))
  }
}

# View the results
print(results_df)

# Optionally save to a file
write.csv(results_df, "paper_tables/PR_K.csv", row.names = FALSE)


# Load necessary libraries
library(ggplot2)
library(tidyr)

# Reshape the results_df for plotting
results_long <- results_df %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value")

# Create the plot with larger fonts
ggplot(results_long, aes(x = Metric, y = Value, fill = as.factor(K))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +  # Add jitter for individual points
  scale_fill_brewer(palette = "Set2", name = "K") +  # Use a color palette for K
  labs(
    title = "Precision and Recall Across Different K Values",
    x = "Metric (Precision/Recall)",
    y = "Value"
  ) +
  theme_minimal(base_size = 25) +  # Set a larger base font size
  theme(
    legend.position = "top",  # Place legend at the top
    plot.title = element_text(hjust = 0.5, face = "bold", size = 30),  # Larger and bold title
    axis.title.x = element_text(face = "bold", size = 28),  # Larger x-axis title
    axis.title.y = element_text(face = "bold", size = 28),  # Larger y-axis title
    axis.text.x = element_text(size = 22),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 22),  # Larger y-axis tick labels
    legend.text = element_text(size = 22),  # Larger legend text
    legend.title = element_text(face = "bold", size = 24)  # Larger and bold legend title
  )

# Save the plot
ggsave("paper_tables/PR_K_main.png", width = 12, height = 8, dpi = 300)
