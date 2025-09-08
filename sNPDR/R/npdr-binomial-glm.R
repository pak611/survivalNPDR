npdr_surv_binomial_glm <- function(outcome = c("time_var" = "time", 
                                              "status_var" = "status"),
                                  dataset, 
                                  attr.diff.type = "numeric-abs",
                                  nbd.method     = "relieff", 
                                  nbd.metric     = "manhattan",
                                  knn            = 20, 
                                  msurf.sd.frac  = 0.5, 
                                  use.glmnet     = FALSE,
                                  glmnet.alpha   = 1, 
                                  glmnet.lower   = -Inf, 
                                  lambda.grid    = seq(0.01, 1, by = 0.01),
                                  glmnet.lam     = "lambda.min",
                                  model.type     = "multinomial", 
                                  KM.weight      = TRUE, 
                                  KM.kernel.type = "gaussian", 
                                  KM.kernel.sigma = 1.0) {
    #browser()
    # Record the start time
    start_time <- Sys.time()
    
    # Extract time and status variable names
    time_var <- as.character(outcome["time_var"])
    status_var <- as.character(outcome["status_var"])

    # Ensure dataset is a data frame
    if(!is.data.frame(dataset)) {
        dataset <- as.data.frame(dataset)
    }

    # Extract the attribute matrix by removing outcome columns
    attr_mat <- dataset |> select(-all_of(as.character(outcome)))
    num_samp <- nrow(dataset)

    # Check and adjust knn value if it's too large
    if(knn > num_samp - 1) {
        warning("The number of nearest neighbors `knn` is too large! Setting `knn = nrow(dataset) - 1`.")
        knn <- num_samp - 1
    }

    # Find nearest neighbors using sNPDR package
    knn_edge_list <- sNPDR::nearestNeighbors(attr.mat   = attr_mat, 
                                             sd.frac    = msurf.sd.frac, 
                                             k          = knn,
                                             nbd.method = nbd.method,
                                             nbd.metric = nbd.metric)

    
    Ri.pheno.vals <- dataset[knn_edge_list$Ri_idx, c(time_var, status_var)]
    NN.pheno.vals <- dataset[knn_edge_list$NN_idx, c(time_var, status_var)]

    Ri.KM.est <- computeKM(Ri.pheno.vals) # Compute the Kaplan-Meier Estimate for Ri
    NN.KM.est <- computeKM(NN.pheno.vals) # Compute the Kaplan-Meier Estimate for NN

    # Initialize the npdr design matrix
    num_attr  <- ncol(attr_mat)
    num_pairs <- nrow(knn_edge_list)
    npdr_diff_mat <- matrix(0, nrow = num_pairs, ncol = num_attr)

    # Compute attribute differences and apply KM weights if specified
    for (i in seq(1, num_attr)) {
        attr_vec <- attr_mat[, i]
        attr_Ris <- attr_vec[knn_edge_list$Ri_idx]
        attr_NNs <- attr_vec[knn_edge_list$NN_idx]
        attr_diff_vec <- sNPDR::npdrDiff(attr_Ris, attr_NNs, diff.type = attr.diff.type)
        
        if (KM.weight) {
            if (KM.kernel.type == 'linear') {
                KM.kernel <- Ri.KM.est$surv_prob * NN.KM.est$surv_prob
            } else if (KM.kernel.type == 'gaussian') {
                KM.kernel <- exp((Ri.KM.est$surv_prob - NN.KM.est$surv_prob)^2 / (2 * KM.kernel.sigma^2))
            }
            attr_diff_vec <- attr_diff_vec * KM.kernel
        }
        npdr_diff_mat[, i] <- attr_diff_vec
    }
  
    colnames(npdr_diff_mat) <- colnames(attr_mat)
  
    # Add time and status information to the knn_edge_list
    knn_edge_list <- knn_edge_list |> 
            mutate(Ri_time = dataset[Ri_idx, time_var],
                   NN_time = dataset[NN_idx, time_var],
                   Ri_status = dataset[Ri_idx, status_var],
                   NN_status = dataset[NN_idx, status_var])

    # # Per-Neighborhood Variance Filter
    # # Calculate variance for each feature within each neighborhood
    # for (i in seq_len(num_attr)) {
    #     # Compute variance for the feature within each neighborhood
    #     variances <- tapply(npdr_diff_mat[, i], knn_edge_list$Ri_idx, var, na.rm = TRUE)
    #     threshold <- quantile(variances, 0.5, na.rm = TRUE)

    #     # Create a mask to set values to NA where the variance is below the threshold
    #     for (j in seq_len(num_pairs)) {
    #         if (variances[knn_edge_list$Ri_idx[j]] <= threshold) {
    #             npdr_diff_mat[j, i] <- NA
    #         }
    #     }
    # }


     # Filter out pairs where NN_time < Ri_time and NN_status == 0
        filter_indices <- which(knn_edge_list$NN_time <= knn_edge_list$Ri_time & knn_edge_list$NN_status == 0)
        npdr_diff_mat_full <- knn_edge_list[-filter_indices, ]
        npdr_diff_mat <- npdr_diff_mat[-filter_indices, ]

       

        # Create a class column for the multinomial model
        npdr_diff_mat_full <- npdr_diff_mat_full |> 
            mutate(class = case_when(
                NN_time > Ri_time ~ 1,
                NN_time == Ri_time ~ 1,
                NN_time < Ri_time ~ 0
            ))

        # Merge the npdr_diff_mat_full with the npdr_diff_mat
        npdr_diff_mat_full <- cbind(npdr_diff_mat_full, npdr_diff_mat) |> 
            select(colnames(attr_mat), class)

    
        # Initialize an empty list to store the results
        results_list <- list()

        # Get the feature names
        features <- colnames(npdr_diff_mat_full)[colnames(npdr_diff_mat_full) != "class"]

        if(!use.glmnet) {
        # Iterate over each feature
        for (feature in features) {
            # Subset the data to remove rows with NA in the current feature
            subset_data <- npdr_diff_mat_full |>
                select(class, all_of(feature)) |>
                na.omit()  # Remove rows with NA in the current feature

            # Fit the binomial logistic regression model for the current feature
            binom_fit <- tryCatch(
            {
                glm(class ~ ., data = subset_data, family = binomial(link = "logit"))
            },
            error = function(e) {
                return(NULL)  # Handle cases where the model cannot be fit
            }
                )

                # Check if the model was successfully fitted
                if (!is.null(binom_fit)) {
                    # Extract and format the coefficients
                    coef_summary <- as.data.frame(summary(binom_fit)$coefficients)
                    coef_summary <- coef_summary[2, ]  # Extract the coefficient for the feature (excluding intercept)
                    coef_summary$Feature <- feature

                    # Store the result
                    results_list[[feature]] <- coef_summary
                }
            }

        }

        if(use.glmnet) {

    
                if (glmnet.lam=="lambda.min"){
                # Run glmnet on the diff attribute columns
                npdrNET.model <- glmnet::cv.glmnet(as.matrix(select(npdr_diff_mat_full, -class)) , npdr_diff_mat_full$class,
                                                    alpha = glmnet.alpha, family = "binomial",
                                                    lower.limits = glmnet.lower, type.measure = "class"
                )
                npdrNET.coeffs <- as.matrix(predict(npdrNET.model, 
                                                    type = "coefficients",
                                                    s=npdrNET.model$lambda.min))
                }else if(glmnet.lam=="lambda.1se"){
                # Run glmnet on the diff attribute columns
                npdrNET.model <- glmnet::cv.glmnet(as.matrix(select(npdr_diff_mat_full, -class)), npdr_diff_mat_full$class,
                                        alpha = glmnet.alpha, family = "binomial",
                                        lower.limits = glmnet.lower, type.measure = "class"
                )
                npdrNET.coeffs <- as.matrix(predict(npdrNET.model, 
                                                    type = "coefficients",
                                                    s=npdrNET.model$lambda.1se))

                  # add npdrNET.coeffs to results_list
                coeffs_df <- as.data.frame(as.matrix(npdrNET.coeffs))
                coeffs_df$Feature <- rownames(coeffs_df)
             
               # Rename the columns appropriately
                names(coeffs_df) <- c("beta", "Feature")

                # Swap the order of the two columns
                coeffs_df <- coeffs_df[, c("Feature", "beta")]
                coeffs_df <- coeffs_df[-1,]


                } else{ # numeric lambda value
                npdrNET.model <- glmnet::glmnet(as.matrix(select(npdr_diff_mat_full, -class)), npdr_diff_mat_full$class,
                                        alpha = glmnet.alpha, family = "binomial",
                                        lambda=glmnet.lam, thresh = 1e-14,
                                        lower.limits = glmnet.lower, type.measure = "class"
                )
                #npdrNET.coeffs <- as.matrix(coef(npdrNET.model))
                npdrNET.coeffs <- as.matrix(predict(npdrNET.model, 
                                                    type = "coefficients",
                                                    s=glmnet.lam))
              
                # add npdrNET.coeffs to results_list
                coeffs_df <- as.data.frame(as.matrix(npdrNET.coeffs))
                coeffs_df$Feature <- rownames(coeffs_df)
                # Rename the columns appropriately
                names(coeffs_df) <- c("beta", "Feature")

                # Swap the order of the two columns
                coeffs_df <- coeffs_df[, c("Feature", "beta")]
                coeffs_df <- coeffs_df[-1,]

            }

        }

            
            # Sort the results by the absolute value of beta (importance)
            results_df <- coeffs_df[order(abs(coeffs_df$beta), decreasing = TRUE), ]


    
    # Return the results
    return(results_df)
}
