


# Load necessary libraries and your custom functions
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/simSurvData/cox_epistasis")
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")




# Assuming the true positives are known
true_functional_features <- paste0("intervar", seq(3, 42))  # True functional features for this simulation


# Number of replicates
num_replicates <- 5


# Simulation parameters
sim_seed <- 2467
n <- 200
p <- 500
n_main <- 2
n_int <- 50
beta_main <- 0.001
beta_int <- 0.002
censparam <- 1/4
lambda <- 1/1000


simdata <- simul.int(seed = sim_seed, n = 200, p = 500,
                n.main = n_main,
                n.int = n_int,
                beta.main= beta_main, 
                beta.int = beta_int, 
                censparam = censparam, 
                lambda = lambda)


# Extract the time variable from simdata
time_data <- simdata$data$time

# Load ggplot2
library(ggplot2)
# Create a histogram of survival times with more increments
ggplot(data = data.frame(time = time_data), aes(x = time)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(
        title = "Distribution of Survival Times",
        x = "Time",
        y = "Frequency"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14)
    )
