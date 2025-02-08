#' Calculate optimal tuning parameters using cross-validation
#' @param exposure_data Matrix of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @param Kmax_range Vector of Kmax values to try
#' @param nfolds Number of cross-validation folds
#' @return List with optimal parameters and cv results
calculate_optimal_parameters <- function(exposure_data, outcome_beta, outcome_se,
                                         Kmax_range = 3:8, nfolds = 5) {
  n <- nrow(exposure_data)
  folds <- cut(seq(1, n), breaks = nfolds, labels = FALSE)
  
  cv_results <- matrix(NA, nrow = length(Kmax_range), ncol = nfolds)
  rownames(cv_results) <- paste0("Kmax_", Kmax_range)
  
  for(i in seq_along(Kmax_range)) {
    Kmax <- Kmax_range[i]
    
    for(fold in 1:nfolds) {
      # Split data into training and validation sets
      test_idx <- which(folds == fold)
      train_exposure <- exposure_data[-test_idx,]
      train_outcome_beta <- outcome_beta[-test_idx]
      train_outcome_se <- outcome_se[-test_idx]
      
      # Scale training data
      train_exposure_scaled <- scale(train_exposure)
      
      # Fit model on training data
      fit <- try({
        MVMR_BFA(exposure_data = train_exposure,
                 outcome_beta = train_outcome_beta,
                 outcome_se = train_outcome_se,
                 Kmax = Kmax,
                 plot_diagnostics = FALSE,
                 generate_report = FALSE)
      })
      
      if(!inherits(fit, "try-error")) {
        # Get test data
        test_exposure <- exposure_data[test_idx,]
        test_outcome <- outcome_beta[test_idx]
        
        # Scale test data using training mean and sd
        test_exposure_scaled <- scale(test_exposure, 
                                      center = attr(train_exposure_scaled, "scaled:center"),
                                      scale = attr(train_exposure_scaled, "scaled:scale"))
        
        # Extract parameters needed for factor score calculation
        n_factors <- ncol(fit$factor_scores)
        p <- ncol(exposure_data)
        
        # Get final factor allocation and loadings
        final_dedic <- fit$bfa_fit$dedic[nrow(fit$bfa_fit$dedic), ]
        alpha_post <- fit$bfa_fit$alpha
        
        # Calculate alpha_means
        alpha_means <- matrix(0, nrow = p, ncol = n_factors)
        for(j in 1:p) {
          k <- final_dedic[j]
          if(k > 0) {
            alpha_means[j, k] <- mean(alpha_post[, j])
          }
        }
        
        # Get unique variances
        sigma_post <- colMeans(fit$bfa_fit$sigma)
        sigma_mat <- diag(sigma_post)
        
        # Calculate factor correlation matrix
        R_means <- matrix(1, n_factors, n_factors)
        if(n_factors > 1) {
          R_post <- fit$bfa_fit$R
          k <- 1
          for(j1 in 1:(n_factors-1)) {
            for(j2 in (j1+1):n_factors) {
              R_means[j1,j2] <- R_means[j2,j1] <- mean(R_post[,k])
              k <- k + 1
            }
          }
        }
        
        # Calculate regression weights for factor scores
        B <- R_means %*% t(alpha_means) %*% solve(alpha_means %*% R_means %*% t(alpha_means) + sigma_mat)
        
        # Calculate factor scores for test data
        test_scores <- test_exposure_scaled %*% t(B)
        
        # Calculate prediction error
        pred <- test_scores %*% fit$mvmr_results$model@Estimate
        cv_results[i, fold] <- mean((test_outcome - pred)^2)
      }
    }
  }
  
  # Calculate mean CV error for each Kmax
  mean_cv_error <- rowMeans(cv_results, na.rm = TRUE)
  
  # Find optimal Kmax
  optimal_Kmax <- Kmax_range[which.min(mean_cv_error)]
  
  return(list(
    optimal_Kmax = optimal_Kmax,
    cv_results = cv_results,
    mean_cv_error = mean_cv_error
  ))
}

# Create visualization of cross-validation results
library(ggplot2)

plot_cv_results <- function(cv_results) {
  cv_data <- data.frame(
    Kmax = as.numeric(gsub("Kmax_", "", rownames(cv_results$cv_results))),
    CV_Error = cv_results$mean_cv_error
  )
  
  ggplot(cv_data, aes(x = Kmax, y = CV_Error)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = cv_results$optimal_Kmax, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Cross-validation Error by Maximum Number of Factors",
         x = "Maximum Number of Factors (Kmax)",
         y = "Mean Cross-validation Error")
}

# Plot cross-validation results
cv_plot <- plot_cv_results(cv_results)
print(cv_plot)
