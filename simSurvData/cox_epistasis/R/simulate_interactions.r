#' Simulate Interactions for Survival Analysis
#'
#' This function simulates interaction effects for survival analysis using a parametric Cox model with an exponential distribution.
#'
#' @param seed An integer seed for random number generation.
#' @param n An integer specifying the number of observations. Default is 100.
#' @param p An integer specifying the number of predictors. Default is 1000.
#' @param n.main An integer specifying the number of main effect variables. Default is 2.
#' @param n.int An integer specifying the number of interaction effect variables. Default is 2.
#' @param beta.main A numeric value for the effect size of main effects. Default is 2.
#' @param beta.int A numeric value for the effect size of interaction effects. Default is 4.
#' @param censparam A numeric value for the censoring parameter. Default is 1/5.
#' @param lambda A numeric value for the rate parameter of the exponential distribution. Default is 1/20.
#'
#' @return A list containing:
#' \item{data}{A data frame with simulated predictor variables, observed survival times, and censoring status.}
#' \item{info}{A data frame with information on the predictor variables and their effect sizes.}
#'
#' @examples
#' \dontrun{
#'   result <- simul.int(seed = 123, n = 200, p = 500, n.main = 3, n.int = 3)
#'   head(result$data)
#'   head(result$info)
#' }
#'
#' @export
simul.int <- function(seed, n = 100, p = 1000,
                      n.main = 2,
                      n.int = 2,
                      beta.main = 2, 
                      beta.int = 4, 
                      censparam = 1/5, 
                      lambda = 1/20) {

  # Create seeds:
  set.seed(seed + 70)
  seeds <- sample(1:100000, 4)
    
  # Simulate X:
  set.seed(seeds[1])
  data <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Create interactions of variables without effect:
  i <- 1:n.int * 2 - 1 + n.main 
  j <- i + 1
  
  data.int <- data[, i] * data[, j]
  print(paste0("interacting variables are ", i, " and ", j))
  
  # Name the interaction columns descriptively
  colnames(data.int) <- paste0("inter", 1:n.int)
  
  data.pred <- cbind(data, data.int)
  
  # Create Survival times (parametrisches Cox-Modell mit Exponentialverteilung):
  beta.main.vec <- c(rep(beta.main, floor(n.main / 2)), rep(-beta.main, ceiling(n.main / 2)))
  beta.int.vec <- c(rep(beta.int, floor(n.int / 2)), rep(-beta.int, ceiling(n.int / 2)))
  beta <- c(beta.main.vec,          # main effects
            rep(0, p - n.main),     # noise
            beta.int.vec)           # interactions
  linpred <- exp(data.pred %*% beta)
  set.seed(seeds[3])
  runifdata <- runif(n, 0, 1)
  obs.time <- array()
  obs.time <- (-log(runif(n, 0, 1)) / (lambda * linpred))
  
  # Censoring:
  set.seed(seeds[1])
  cens.time <- (-log(runif(n, 0, 1)) / (censparam))
  obs.status <- ifelse(obs.time <= cens.time, 1, 0)
  obs.time <- ifelse(obs.time <= cens.time, obs.time, cens.time)
  
  obs.status[obs.time > 1 / censparam * 2] <- 0
  obs.time[obs.time > 1 / censparam * 2] <- 1 / censparam * 2
  obs.time[which(obs.time == 9)] <- 10 ^ -5
  
  # Summarize:
  simdata <- as.data.frame(cbind(data, obs.time, obs.status))
  colnames(simdata) <- c(paste('X', 1:p, sep = ''), 'time', 'status')
  
  # Rename main effect variables
  main_effect_indices <- 1:n.main
  colnames(simdata)[main_effect_indices] <- paste0("simvar", main_effect_indices)
  
  # Rename interaction variables
  colnames(simdata)[i] <- paste0("intervar", i)
  colnames(simdata)[j] <- paste0("intervar", j)
  
  # Data information:
  info <- as.data.frame(
    rbind(cbind(colnames(simdata)[1:p], beta[1:p])[which(beta[1:p] != 0), ],
          cbind(paste(colnames(simdata)[i], ':', colnames(simdata)[j], sep = ""), beta.int.vec))
  )
  colnames(info) <- c('X', 'Effect size')
  
  out <- list()
  out$data <- simdata
  out$info <- info
  return(out)
}