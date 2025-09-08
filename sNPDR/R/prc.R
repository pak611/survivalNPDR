#' Compute Precision-Recall Curve with Smoothing and Outlier Removal
#'
#' This function computes the precision-recall curve for given scores of two classes,
#' removes zero precision points, eliminates large jumps, and optionally interpolates for smoothing.
#'
#' @param score.class0 A numeric vector of scores for class 0 (negative class).
#' @param score.class1 A numeric vector of scores for class 1 (positive class).
#' @param interpolate A logical value indicating whether to interpolate for smoothing (default: TRUE).
#' @return A data frame with thresholds, precision, and recall values.
#' @export
compute_prc <- function(score.class0, score.class1, interpolate = TRUE) {
  # Combine scores and create labels
  scores <- c(score.class0, score.class1)
  labels <- c(rep(0, length(score.class0)), rep(1, length(score.class1)))

  # Create thresholds
  thresholds <- sort(unique(scores), decreasing = TRUE)

  # Initialize precision and recall lists
  precision <- c()
  recall <- c()

  # Total positives
  total_positives <- sum(labels)

  # Iterate through thresholds
  for (threshold in thresholds) {
    # Classify based on the threshold
    predictions <- ifelse(scores >= threshold, 1, 0)

    # True Positives (TP), False Positives (FP), False Negatives (FN)
    TP <- sum(predictions == 1 & labels == 1)
    FP <- sum(predictions == 1 & labels == 0)
    FN <- total_positives - TP

    # Calculate precision and recall
    prec <- ifelse((TP + FP) > 0, TP / (TP + FP), 1) # Avoid division by zero
    rec <- TP / total_positives

    # Append to lists
    precision <- c(precision, prec)
    recall <- c(recall, rec)
  }

  # Remove zero precision points
  non_zero_indices <- which(precision > 0)
  thresholds <- thresholds[non_zero_indices]
  precision <- precision[non_zero_indices]
  recall <- recall[non_zero_indices]

  # Remove outliers based on changes in precision and recall
  precision_diff <- abs(diff(precision))
  recall_diff <- abs(diff(recall))

  # Use Median Absolute Deviation (MAD) to filter out large jumps
  precision_mad <- mad(precision_diff)
  recall_mad <- mad(recall_diff)

  valid_indices <- which(
    c(TRUE, precision_diff <= 3 * precision_mad) & # Allow 3 MAD tolerance for precision
    c(TRUE, recall_diff <= 3 * recall_mad)        # Allow 3 MAD tolerance for recall
  )

  thresholds <- thresholds[valid_indices]
  precision <- precision[valid_indices]
  recall <- recall[valid_indices]

  # Interpolation for smoothing, if enabled
  if (interpolate) {
    interpolated_recall <- seq(min(recall), max(recall), length.out = 100)
    interpolated_precision <- approx(recall, precision, xout = interpolated_recall, method = "linear")$y

    return(data.frame(
      threshold = approx(recall, thresholds, xout = interpolated_recall, method = "linear")$y,
      recall = round(interpolated_recall, digits = 6), # Ensure consistent precision
      precision = round(interpolated_precision, digits = 6) # Ensure consistent precision
    ))
  }

  # Return without interpolation
  return(data.frame(
    threshold = thresholds,
    recall = round(recall, digits = 6), # Ensure consistent precision
    precision = round(precision, digits = 6) # Ensure consistent precision
  ))
}
