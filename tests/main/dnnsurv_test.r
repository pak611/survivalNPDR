# Simple DNNSurv test script

# Install required packages if not already installed
if (!requireNamespace("survivalmodels", quietly = TRUE)) install.packages("survivalmodels")
if (!requireNamespace("survcomp", quietly = TRUE)) install.packages("survcomp")

library(survivalmodels)
library(survival)
library(dplyr)

library(reticulate)
use_python("C:/Users/patri/AppData/Local/Programs/Python/Python310/python.exe", required = TRUE)

# Diagnostic: Print which Python reticulate is using and check for pycox
library(reticulate)
cat('Reticulate is using Python at:', py_config()$python, '\n')
if (py_module_available('pycox')) {
  cat('pycox is available in this Python environment.\n')
} else {
  cat('pycox is NOT available in this Python environment.\n')
}

install_pycox(pip = TRUE, install_torch = TRUE)
install_keras(pip = TRUE, install_tensorflow = TRUE)

# Ensure required Python modules are installed
if (reticulate::py_module_available('pip')) {
  if (!reticulate::py_module_available('pycox')) reticulate::py_install('pycox', pip = TRUE)
  if (!reticulate::py_module_available('torch')) reticulate::py_install('torch', pip = TRUE)
  if (!reticulate::py_module_available('torchtuples')) reticulate::py_install('torchtuples', pip = TRUE)
} else {
  warning('pip is not available in the selected Python environment. Please install pip and required modules manually.')
}

# Simulate data (same as in main_effect4.r)
sim_seed <- 2467
n <- 200
p <- 1000
n_main <- 50
n_int <- 2
beta_main <- 0.7
beta_int <- 0.1
censparam <- 1/4
lambda <- 1/100

devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/simSurvData/cox_epistasis")
simdata <- simul.int(sim_seed, n = n, p = p,
                     n.main = n_main, n.int = n_int,
                     beta.main = beta_main, beta.int = beta_int,
                     censparam = censparam, lambda = lambda)
dat <- simdata$data

# Split into train/test (80/20 split)
set.seed(123)
train_idx <- sample(seq_len(nrow(dat)), size = 0.8 * nrow(dat))
train_data <- dat[train_idx, ]
test_data <- dat[-train_idx, ]

# Prepare data for DNNSurv
attr_train <- train_data %>% select(-time, -status) %>% scale()
attr_test <- test_data %>% select(-time, -status) %>% scale(center = attr_train %>% attr("scaled:center"), scale = attr_train %>% attr("scaled:scale"))

# Fit DNNSurv
set.seed(123)
dnn_fit <- deepsurv(
  x = as.matrix(attr_train),
  y = Surv(train_data$time, train_data$status),
  x_test = as.matrix(attr_test),
  y_test = Surv(test_data$time, test_data$status),
  epochs = 100,
  batch_size = 32,
  learning_rate = 1e-3,
  patience = 5,
  verbose = TRUE
)

# Evaluate C-index on test set
risk_scores <- predict(dnn_fit, newdata = as.matrix(attr_test))
library(survcomp)
c.ind.dnn <- concordance.index(x = risk_scores, surv.time = test_data$time, surv.event = test_data$status)
cat("DNNSurv test set C-index:", c.ind.dnn$c.index, "\n")
