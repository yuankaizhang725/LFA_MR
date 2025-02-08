#' Generate summary report for MVMR_BFA analysis
#' @param mvmr_bfa_results Results from MVMR_BFA function
#' @param diagnostics Diagnostic statistics from calculate_diagnostics
#' @return Character string containing formatted report
generate_mvmr_bfa_report <- function(mvmr_bfa_results, diagnostics) {
  report <- ""
  
  # Basic information
  report <- paste0(report, "MVMR-BFA Analysis Report\n",
                   "=======================\n\n")
  
  # Factor Analysis Summary
  report <- paste0(report, "Factor Analysis Results:\n",
                   "-------------------------\n",
                   "Number of factors identified: ", ncol(mvmr_bfa_results$factor_scores), "\n",
                   "Number of variables: ", nrow(mvmr_bfa_results$bfa_summary$alpha), "\n\n")
  
  # Convergence Diagnostics
  report <- paste0(report, "Convergence Diagnostics:\n",
                   "------------------------\n",
                   "MH acceptance rate: ", round(mean(mvmr_bfa_results$convergence$mh_acceptance), 3), "\n",
                   "Factor stability: ", 
                   round(mean(diagnostics$dedication_stability), 3), "\n\n")
  
  # Factor Loading Summary
  report <- paste0(report, "Factor Loading Summary:\n",
                   "---------------------\n")
  loading_summary <- mvmr_bfa_results$bfa_summary$alpha
  report <- paste0(report, "Mean absolute loadings by factor:\n")
  for(i in 1:ncol(loading_summary)) {
    report <- paste0(report, "Factor ", i, ": ", 
                     round(mean(abs(loading_summary[,i])), 3), "\n")
  }
  report <- paste0(report, "\n")
  
  # MVMR Results
  report <- paste0(report, "MVMR Results:\n",
                   "-------------\n")
  mvmr_estimates <- mvmr_bfa_results$mvmr_results$model@Estimate
  mvmr_se <- mvmr_bfa_results$mvmr_results$model@StdError
  
  for(i in seq_along(mvmr_estimates)) {
    report <- paste0(report, 
                     "Factor ", i, ": ",
                     "Estimate = ", round(mvmr_estimates[i], 4),
                     " (SE = ", round(mvmr_se[i], 4), ")\n")
  }
  
  return(report)
}
