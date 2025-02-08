#' Calculate diagnostic statistics for MVMR_BFA
#' @param bfa_fit BFA fit object
#' @param factor_scores Matrix of factor scores
#' @return List of diagnostic statistics
calculate_diagnostics <- function(bfa_fit, factor_scores) {
  # Factor score metrics
  factor_correlations <- cor(factor_scores)
  factor_variances <- apply(factor_scores, 2, var)
  
  # BFA convergence diagnostics
  factor_frequencies <- table(bfa_fit$nfac) / length(bfa_fit$nfac)
  dedication_stability <- apply(bfa_fit$dedic, 2, function(x) {
    length(unique(x)) / length(x)
  })
  
  # Effective sample size for key parameters
  ess_alpha <- effectiveSize(mcmc(bfa_fit$alpha))
  ess_sigma <- effectiveSize(mcmc(bfa_fit$sigma))
  
  # Geweke diagnostics
  geweke_alpha <- geweke.diag(mcmc(bfa_fit$alpha))
  geweke_sigma <- geweke.diag(mcmc(bfa_fit$sigma))
  
  return(list(
    factor_correlations = factor_correlations,
    factor_variances = factor_variances,
    factor_frequencies = factor_frequencies,
    dedication_stability = dedication_stability,
    effective_sample_size = list(alpha = ess_alpha, sigma = ess_sigma),
    geweke_diagnostics = list(alpha = geweke_alpha, sigma = geweke_sigma)
  ))
}
