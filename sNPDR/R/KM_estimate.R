#' Compute Kaplan-Meier Survival Estimates and Merge With Data
#'
#' This function computes Kaplan-Meier survival estimates for a given dataset
#' and adds a column to the dataset with these survival probabilities.
#' @param NN.df A dataframe containing at least two columns: `time` and `outcome`,
#' where `time` is the time to event or censoring and `outcome` is a binary
#' indicator of the event occurrence (1 if the event occurred, 0 if censored).
#' @return A dataframe similar to `NN.df` but with an additional column `surv_prob`
#' indicating the survival probability at each time point.
#' @importFrom dplyr left_join mutate
#' @importFrom survival Surv survfit
#' @importFrom broom tidy
#' @examples
#' # Assuming NN.df has columns `time` and `outcome`
#' # new_df <- computeKM(NN.df)
#' @export
computeKM <- function(NN.df) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("broom", quietly = TRUE)) {
    stop("Package 'broom' is required.")
  }

  #browser()
  
  # Creating the survival object correctly
  surv_obj <- survival::Surv(time = NN.df$time, event = as.numeric(NN.df$status))
  
  # Compute the Kaplan-Meier estimate
  km_estimate <- survival::survfit(surv_obj ~ 1)
  
  # Tidy the km_estimate for merging
  km_data <- broom::tidy(km_estimate) %>%
    dplyr::select(time, estimate) %>%
    dplyr::rename(surv_prob = estimate)
  
  # Merge the Kaplan-Meier estimates with NN.df based on time
  new_df <- NN.df %>%
    dplyr::left_join(km_data, by = "time")
  
  return(new_df)
}

