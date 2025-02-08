#' Plot diagnostics for MVMR_BFA results
#' @param mvmr_bfa_results Results from MVMR_BFA function
#' @param diagnostics Diagnostic statistics from calculate_diagnostics
#' @return List of ggplot objects
plot_mvmr_bfa_diagnostics <- function(mvmr_bfa_results, diagnostics) {
  library(ggplot2)
  library(tidyr)
  
  plots <- list()
  
  # Factor loadings heatmap from alpha_means matrix
  p <- nrow(mvmr_bfa_results$bfa_summary$alpha)
  n_factors <- max(mvmr_bfa_results$bfa_fit$dedic[nrow(mvmr_bfa_results$bfa_fit$dedic), ])
  
  # Create alpha_means matrix
  alpha_means <- matrix(0, nrow = p, ncol = n_factors)
  final_dedic <- mvmr_bfa_results$bfa_fit$dedic[nrow(mvmr_bfa_results$bfa_fit$dedic), ]
  alpha_post <- mvmr_bfa_results$bfa_fit$alpha
  
  # Extract factor loadings
  for(i in 1:p) {
    k <- final_dedic[i]
    if(k > 0) {  # if variable loads on a factor
      alpha_means[i, k] <- mean(alpha_post[, i])
    }
  }
  
  # Convert to data frame for plotting
  clean_names <- gsub("alpha:", "", rownames(mvmr_bfa_results$bfa_summary$alpha))
  
  rownames(alpha_means) <- clean_names
  colnames(alpha_means) <- paste0("Factor", 1:n_factors)
  
  loadings_data <- as.data.frame(alpha_means)
  loadings_data$Variable <- rownames(alpha_means)
  
  loadings_long <- tidyr::pivot_longer(
    loadings_data,
    cols = -Variable,
    names_to = "Factor",
    values_to = "Loading"
  )
  
  plots$loadings <- ggplot(loadings_long, aes(x = Factor, y = Variable, fill = Loading)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    labs(title = "Factor Loadings Heatmap") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Trace plots for number of factors
  nfac_data <- data.frame(
    Iteration = 1:length(mvmr_bfa_results$bfa_fit$nfac),
    NumberOfFactors = mvmr_bfa_results$bfa_fit$nfac
  )
  
  plots$nfac_trace <- ggplot(nfac_data, aes(x = Iteration, y = NumberOfFactors)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Trace Plot of Number of Factors",
         x = "Iteration",
         y = "Number of Factors")
  
  # Factor score distributions
  factor_scores_df <- as.data.frame(mvmr_bfa_results$factor_scores)
  factor_scores_long <- tidyr::pivot_longer(
    factor_scores_df,
    cols = everything(),
    names_to = "Factor",
    values_to = "Score"
  )
  
  plots$factor_dist <- ggplot(factor_scores_long, aes(x = Score, fill = Factor)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Factor) +
    theme_minimal() +
    labs(title = "Factor Score Distributions",
         x = "Score",
         y = "Density")
  
  # MCMC trace plots for representative loadings
  alpha_mcmc <- as.data.frame(mvmr_bfa_results$bfa_fit$alpha[,1:min(5, ncol(mvmr_bfa_results$bfa_fit$alpha))])
  colnames(alpha_mcmc) <- gsub("alpha:", "", colnames(alpha_mcmc))
  alpha_long <- tidyr::pivot_longer(
    alpha_mcmc,
    cols = everything(),
    names_to = "Loading",
    values_to = "Value"
  )
  alpha_long$Iteration <- rep(1:nrow(alpha_mcmc), ncol(alpha_mcmc))
  
  plots$alpha_trace <- ggplot(alpha_long, aes(x = Iteration, y = Value, color = Loading)) +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    labs(title = "MCMC Traces for Factor Loadings",
         x = "Iteration",
         y = "Loading Value")
  
  return(plots)
}
