#' Compare sensitivity analysis results
#' @param sensitivity_results Results from perform_sensitivity_analysis
#' @return List of comparison metrics and plots
compare_sensitivity_results <- function(sensitivity_results) {
  library(ggplot2)
  
  comparisons <- list()
  
  # Extract key metrics for each model
  metrics <- lapply(names(sensitivity_results), function(name) {
    if(inherits(sensitivity_results[[name]], "try-error")) {
      return(NULL)
    }
    
    result <- sensitivity_results[[name]]
    
    # Extract parameters from name
    params <- strsplit(name, "_")[[1]]
    Kmax <- as.numeric(params[2])
    model <- params[4]
    Nid <- as.numeric(params[6])
    
    # Get key metrics
    estimates <- result$mvmr_results$model@Estimate
    se <- result$mvmr_results$model@StdError
    n_factors <- ncol(result$factor_scores)
    
    data.frame(
      model_name = name,
      Kmax = Kmax,
      model_type = model,
      Nid = Nid,
      n_factors = n_factors,
      mean_estimate = mean(abs(estimates)),
      mean_se = mean(se)
    )
  })
  
  # Combine all metrics
  metrics_df <- do.call(rbind, metrics[!sapply(metrics, is.null)])
  
  # Create comparison plots
  plots <- list()
  
  # Number of factors vs Kmax
  plots$factors <- ggplot(metrics_df, 
                          aes(x = Kmax, y = n_factors, color = model_type)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Nid) +
    theme_minimal() +
    labs(title = "Number of Factors vs Kmax")
  
  # Mean estimates comparison
  plots$estimates <- ggplot(metrics_df, 
                            aes(x = Kmax, y = mean_estimate, color = model_type)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Nid) +
    theme_minimal() +
    labs(title = "Mean Absolute Estimates vs Kmax")
  
  return(list(
    metrics = metrics_df,
    plots = plots
  ))
}
