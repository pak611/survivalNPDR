setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("censor_cv.R")

# Load necessary libraries
library(ggplot2)

# Read in the dataframe with AUC results
auc_results <- readRDS("auc_results.rds")

# Reshape the data for easier plotting
library(tidyr)
auc_long <- auc_results %>%
  pivot_longer(
    cols = c(npdr, standard_cox, survival_svm),
    names_to = "algorithm",
    values_to = "AUC"
  )

# Create the plot with smoothed lines
auc_plot <- ggplot(auc_long, aes(x = censor_rate, y = AUC, color = algorithm)) +
  geom_smooth(size = 2, se = FALSE) +  # Smooth the lines and remove confidence interval shading
  geom_point(size = 4) +  # Increase point size
  labs(
    title = "AUC as a Function of Censoring Rate",
    x = "Censoring Rate",
    y = "AUC",
    color = "Algorithm"
  ) +
  theme_classic(base_size = 20) +  # Use a classic theme for a cleaner look
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),  # Increase plot title size
    axis.title = element_text(size = 25),  # Increase axis titles size
    axis.text = element_text(size = 23),  # Increase axis text size
    legend.title = element_text(size = 25),  # Increase legend title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.key.size = unit(1.5, "lines")  # Increase legend key size
  )

# Display the plot
print(auc_plot)


setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("effect_size_cv.R")

# Load necessary libraries
library(ggplot2)

# Read in the dataframe with AUC results
auc_results <- readRDS("auc_results_effect_size.rds")

# Reshape the data for easier plotting
library(tidyr)
auc_long <- auc_results %>%
  pivot_longer(
    cols = c(npdr, standard_cox, survival_svm),
    names_to = "algorithm",
    values_to = "AUC"
  )

# Create the plot
auc_plot <- ggplot(auc_long, aes(x = effect_size, y = AUC, color = algorithm)) +
  geom_line(size = 2) +  # Increase line thickness
  geom_point(size = 4) +  # Increase point size
  labs(
    title = "AUC as a Function of Background Noise",
    x = "Effect Size",
    y = "AUC",
    color = "Algorithm"
  ) +
  theme_minimal(base_size = 20) +  # Increase base font size
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),  # Increase plot title size
    axis.title = element_text(size = 20),  # Increase axis titles size
    axis.text = element_text(size = 18)  # Increase axis text size
  )

# Display the plot
print(auc_plot)


setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("noise_cv.R")

library(tidyr)
library(ggplot2)

# Assuming auc_results is already loaded or created
# Pivot the data to long format
auc_long <- auc_results %>%
  pivot_longer(
    cols = c(npdr, standard_cox, survival_svm),
    names_to = "algorithm",
    values_to = "AUC"
  )

# Create the plot with smoothed lines
auc_plot <- ggplot(auc_long, aes(x = noise_level, y = AUC, color = algorithm)) +
  geom_smooth(size = 2, se = FALSE) +  # Smooth the lines and remove confidence interval shading
  geom_point(size = 4) +  # Increase point size
  labs(
    title = "AUC as a Function of Noise-to-Signal Ratio",
    x = "Noise-to-Signal Ratio",
    y = "AUC",
    color = "Algorithm"
  ) +
  theme_classic(base_size = 20) +  # Use a classic theme for a cleaner look
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),  # Increase plot title size
    axis.title = element_text(size = 25),  # Increase axis titles size
    axis.text = element_text(size = 23),  # Increase axis text size
    legend.title = element_text(size = 25),  # Increase legend title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.key.size = unit(1.5, "lines")  # Increase legend key size
  )

# Display the plot
print(auc_plot)

# setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
# source("test_k_main.R")



# # Load necessary libraries
# library(ggplot2)
# library(dplyr)

# # Read in the precision-recall data
# data <- read.csv("PR_K.csv")


# # Load necessary library
# library(tidyr)

# # Reshape the data into long format
# long_data <- data %>%
#   pivot_longer(cols = c(Precision, Recall), names_to = "Metric", values_to = "Value")

# # Check the structure of the reshaped data
# head(long_data)


# # Create the grouped boxplot
# ggplot(long_data, aes(x = Metric, y = Value, fill = as.factor(K))) +
#   geom_boxplot() +
#   labs(
#        x = "Metric (Precision/Recall)",
#        y = "Value",
#        fill = "K") +
#   theme_minimal(base_size = 20) +  # Increase base font size
#   theme(
#     legend.position = c(0.1, 0.9),  # Position legend at the top left
#     legend.background = element_rect(fill = "white", color = "black"),
#     legend.title = element_text(size = 18),  # Increase legend title size
#     legend.text = element_text(size = 16),   # Increase legend text size
#     plot.title = element_text(size = 24),    # Increase plot title size
#     axis.title = element_text(size = 20),    # Increase axis titles size
#     axis.text = element_text(size = 18),     # Increase axis text size
#     panel.background = element_rect(fill = "gray95", color = NA),  # Light gray panel background
#     panel.border = element_rect(color = "black", fill = NA, size = 2)  # Black border around panel
#   ) +
#   guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

# ggsave("precision_recall_boxplot2.png", width = 10, height = 6)



setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("test_k_interaction.R")



# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read in the precision-recall data
data <- read.csv("precision_recall_k_values_results_with_simulation.csv")


# Load necessary library
library(tidyr)

# Reshape the data into long format
long_data <- data %>%
  pivot_longer(cols = c(Precision, Recall), names_to = "Metric", values_to = "Value")

# Check the structure of the reshaped data
head(long_data)


# Create the grouped boxplot
ggplot(long_data, aes(x = Metric, y = Value, fill = as.factor(K))) +
  geom_boxplot() +
  labs(
       x = "Metric (Precision/Recall)",
       y = "Value",
       fill = "K") +
  theme_minimal(base_size = 20) +  # Increase base font size
  theme(
    legend.position = c(0.1, 0.9),  # Position legend at the top left
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 16),   # Increase legend text size
    plot.title = element_text(size = 24),    # Increase plot title size
    axis.title = element_text(size = 20),    # Increase axis titles size
    axis.text = element_text(size = 18)      # Increase axis text size
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

ggsave("precision_recall_boxplot2.png", width = 10, height = 6)


setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("main_effect_concordance.R")


# read in cross_validation_results.csv
cross_validation_results <- read.csv("cross_validation_results.csv")

print(cross_validation_results)


# Summarize the results
summary_results <- cross_validation_results %>%
  group_by(method) %>%
  summarise(
    mean_c_index = mean(c_index, na.rm = TRUE),
    median_c_index = median(c_index, na.rm = TRUE),
    sd_c_index = sd(c_index, na.rm = TRUE)
  )

  # Summarize the simulation parameters
simulation_parameters <- cross_validation_results %>%
  select(n, p, n_main, n_int, beta_main, beta_int, censparam, lambda) %>%
  distinct()

# Print the simulation parameters
print(simulation_parameters)

# Print the summary results
print(summary_results)


# Calculate mean concordance index for each method
mean_c_index <- cross_validation_results %>%
  group_by(method) %>%
  summarise(mean_c_index = mean(c_index, na.rm = TRUE))

# Create a boxplot comparing the various methods with mean annotations
ggplot(cross_validation_results, aes(x = method, y = c_index, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Comparison of Concordance Index Across Methods",
    x = "Method",
    y = "Concordance Index"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23)  # Increase axis text size
  ) +
  # Add mean concordance index values as text above each boxplot
  geom_text(data = mean_c_index, aes(x = method, y = mean_c_index + 0.15, label = sprintf("Mean: %.3f", mean_c_index)), 
            size = 10, color = "black", fontface = "bold")


setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("interaction_effect_concordance.R")

cross_validation_results <- read.csv("cross_validation_results_interactions.csv")
print(cross_validation_results)


# Summarize the results
summary_results <- cross_validation_results %>%
  group_by(method) %>%
  summarise(
    mean_c_index = mean(c_index, na.rm = TRUE),
    median_c_index = median(c_index, na.rm = TRUE),
    sd_c_index = sd(c_index, na.rm = TRUE)
  )

  # Summarize the simulation parameters
simulation_parameters <- cross_validation_results %>%
  select(n, p, n_main, n_int, beta_main, beta_int, censparam, lambda) %>%
  distinct()

# Print the simulation parameters
print(simulation_parameters)

# Print the summary results
print(summary_results)


# Calculate mean concordance index for each method
mean_c_index <- cross_validation_results %>%
  group_by(method) %>%
  summarise(mean_c_index = mean(c_index, na.rm = TRUE))

# Create a boxplot comparing the various methods with mean annotations
ggplot(cross_validation_results, aes(x = method, y = c_index, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Comparison of Concordance Index Across Methods",
    x = "Method",
    y = "Concordance Index"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23)  # Increase axis text size
  ) +
  # Add mean concordance index values as text above each boxplot
  geom_text(data = mean_c_index, aes(x = method, y = mean_c_index + 0.15, label = sprintf("Mean: %.3f", mean_c_index)), 
            size = 10, color = "black", fontface = "bold")


setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("main_effects_cv.R")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the AUC results
auc_results_loaded <- readRDS("auc_results.rds")

# Rename the columns
colnames(auc_results_loaded) <- c("Iteration", "survNPDR", "survNPDR+GLM", "Cox", "SVM")

# Convert auc_results to long format for ggplot2
auc_results_long <- auc_results_loaded %>%
  pivot_longer(cols = c("survNPDR", "survNPDR+GLM", "Cox", "SVM"), names_to = "Algorithm", values_to = "AUC")

# Calculate average AUC for each algorithm
mean_auc <- auc_results_long %>%
  group_by(Algorithm) %>%
  summarise(mean_AUC = mean(AUC))

# Create the boxplot with mean AUC annotations
ggplot(auc_results_long, aes(x = Algorithm, y = AUC, fill = Algorithm)) +
  geom_boxplot() +
  labs(title = "MAIN EFFECTS",
       x = "Algorithm",
       y = "AUC") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(size = 30),    # Increase plot title size
    axis.title = element_text(size = 30),    # Increase axis titles size
    axis.text = element_text(size = 25),      # Increase axis text size
    panel.background = element_rect(fill = "gray95", color = NA),  # Gray panel background
    panel.border = element_rect(color = "black", fill = NA, size = 2)
  ) +
  # Add mean AUC values as text above each boxplot
  geom_text(data = mean_auc, aes(x = Algorithm, y = mean_AUC + 0.02, label = sprintf("Mean: %.3f", mean_AUC)), 
            size = 10, color = "black", fontface = "bold")

# Save the boxplot
ggsave("auc_results_boxplot.png", width = 10, height = 6)



setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("interactions_cv.R")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the AUC results
auc_results_loaded <- readRDS("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests/auc_results_interactions.rds")
#auc_results_loaded <- readRDS("auc_results_interactions.rds")

# Rename the columns
colnames(auc_results_loaded) <- c("Iteration", "survNPDR", "survNPDR-GLM", "Cox", "SVM")

# Convert auc_results to long format for ggplot2
auc_results_long <- auc_results_loaded %>%
  pivot_longer(cols = c("survNPDR", "survNPDR-GLM", "Cox", "SVM"), names_to = "Algorithm", values_to = "AUC")

# Calculate average AUC for each algorithm
mean_auc <- auc_results_long %>%
  group_by(Algorithm) %>%
  summarise(mean_AUC = mean(AUC))

# Create the boxplot with mean AUC annotations
ggplot(auc_results_long, aes(x = Algorithm, y = AUC, fill = Algorithm)) +
  geom_boxplot() +
  labs(title = "INTERACTIONS",
       x = "Algorithm",
       y = "AUC") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(size = 30),    # Increase plot title size
    axis.title = element_text(size = 30),    # Increase axis titles size
    axis.text = element_text(size = 25),     # Increase axis text size
    panel.background = element_rect(fill = "gray95", color = NA),  # Gray panel background
    panel.border = element_rect(color = "black", fill = NA, size = 2)  # Black border around panel
  ) +
  # Add mean AUC values as text above each boxplot
  geom_text(data = mean_auc, aes(x = Algorithm, y = mean_AUC + 0.02, label = sprintf("Mean: %.3f", mean_AUC)), 
            size = 10, color = "black", fontface = "bold")

# Save the boxplot
ggsave("auc_results_boxplot.png", width = 10, height = 6)




setwd("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests")
source("interaction_effect_size_cv.R")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the AUC results
auc_results_loaded <- readRDS("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests/auc_results_interactions.rds")
#auc_results_loaded <- readRDS("auc_results_interactions.rds")

# Rename the columns
colnames(auc_results_loaded) <- c("Iteration", "survNPDR", "Cox", "SVM")

# Convert auc_results to long format for ggplot2
auc_results_long <- auc_results_loaded %>%
  pivot_longer(cols = c("survNPDR", "Cox", "SVM"), names_to = "Algorithm", values_to = "AUC")

# Calculate average AUC for each algorithm
mean_auc <- auc_results_long %>%
  group_by(Algorithm) %>%
  summarise(mean_AUC = mean(AUC))

# Create the boxplot with mean AUC annotations
ggplot(auc_results_long, aes(x = Algorithm, y = AUC, fill = Algorithm)) +
  geom_boxplot() +
  labs(title = "INTERACTIONS",
       x = "Algorithm",
       y = "AUC") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(size = 30),    # Increase plot title size
    axis.title = element_text(size = 30),    # Increase axis titles size
    axis.text = element_text(size = 30),     # Increase axis text size
    panel.background = element_rect(fill = "gray95", color = NA),  # Gray panel background
    panel.border = element_rect(color = "black", fill = NA, size = 2)  # Black border around panel
  ) +
  # Add mean AUC values as text above each boxplot
  geom_text(data = mean_auc, aes(x = Algorithm, y = mean_AUC + 0.02, label = sprintf("Mean: %.3f", mean_AUC)), 
            size = 10, color = "black", fontface = "bold")

# Save the boxplot
ggsave("auc_results_boxplot.png", width = 10, height = 6)