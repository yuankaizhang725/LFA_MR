library(BayesFM)
library(MendelianRandomization)
library(ggplot2)
library(coda)
library(tidyr)

# Main MVMR_BFA function
#' @title MVMR with Bayesian Factor Analysis
#' @description Implements Multivariate Mendelian Randomization using Bayesian Factor Analysis 
#' to handle multiple correlated exposures
#' 
#' @param exposure_beta Matrix or data frame of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @param Kmax Maximum number of factors to consider
#' @param model Character specifying either "fixed" or "random" effects model
#' @param Nid Integer specifying minimum number of variables per factor
#' @param burnin Integer specifying number of burnin iterations for MCMC
#' @param iter Integer specifying number of MCMC iterations
#' @param plot_diagnostics Logical, whether to generate diagnostic plots
#' @param generate_report Logical, whether to generate summary report
#' @return Object of class mvmr_bfa containing analysis results
MVMR_BFA <- function(exposure_beta, 
                     outcome_beta,
                     outcome_se,
                     Kmax,
                     model = "fixed",
                     Nid = 2,
                     burnin = 500,
                     iter = 1000,
                     nu0 = Kmax+1,
                     kappa = 1/Kmax,
                     kappa0 = 2,
                     xi0 = 1,
                     seed = NULL,
                     plot_diagnostics = TRUE,
                     generate_report = TRUE) {
  
  # Input validation
  validation <- check_mvmr_bfa_inputs(exposure_beta, outcome_beta, outcome_se)
  if (!validation$valid) {
    stop(paste(validation$messages, collapse = "\n"))
  }
  if (length(validation$messages) > 0) {
    warning(paste(validation$messages, collapse = "\n"))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # Scale exposure data
  exposure_scaled <- scale(exposure_beta)
  
  # Perform Bayesian Factor Analysis
  bfa_fit <- befa(exposure_scaled, 
                  Nid = Nid,
                  Kmax = Kmax,
                  nu0 = nu0,
                  kappa = kappa,
                  kappa0 = kappa0,
                  xi0 = xi0,
                  burnin = burnin,
                  iter = iter)
  
  # Post-process BFA results
  bfa_fit <- post.column.switch(bfa_fit)
  bfa_fit <- post.sign.switch(bfa_fit)
  bfa_summary <- summary(bfa_fit)
  
  # Get number of active factors from final iteration
  n_factors <- max(bfa_fit$dedic[nrow(bfa_fit$dedic), ])
  
  # Calculate factor scores using regression method
  p <- ncol(exposure_scaled)
  alpha_post <- bfa_fit$alpha
  alpha_means <- matrix(0, nrow = p, ncol = n_factors)
  
  # Extract final factor allocation
  final_dedic <- bfa_fit$dedic[nrow(bfa_fit$dedic), ]
  
  # Extract factor loadings
  for(i in 1:p) {
    k <- final_dedic[i]
    if(k > 0) {  # if variable loads on a factor
      alpha_means[i, k] <- mean(alpha_post[, i])
    }
  }
  
  # Get unique variances
  sigma_post <- colMeans(bfa_fit$sigma)
  sigma_mat <- diag(sigma_post)
  
  # Get factor correlation matrix from posterior means
  R_means <- matrix(1, n_factors, n_factors)
  if(n_factors > 1) {
    R_post <- bfa_fit$R
    k <- 1
    for(i in 1:(n_factors-1)) {
      for(j in (i+1):n_factors) {
        R_means[i,j] <- R_means[j,i] <- mean(R_post[,k])
        k <- k + 1
      }
    }
  }
  
  # Calculate regression weights for factor scores
  Lambda <- alpha_means  # factor loadings matrix
  Psi <- sigma_mat      # unique variances matrix
  R <- R_means          # factor correlation matrix
  
  B <- R %*% t(Lambda) %*% solve(Lambda %*% R %*% t(Lambda) + Psi)
  
  # Calculate factor scores using regression method
  factor_scores <- exposure_scaled %*% t(B)
  colnames(factor_scores) <- paste0("Factor", 1:n_factors)
  
  # Prepare matrices for MVMR
  bx_matrix <- as.matrix(factor_scores)
  bxse_matrix <- matrix(0, nrow = nrow(factor_scores), ncol = ncol(factor_scores))
  
  # Create MVMR input object
  mvmr_input <- mr_mvinput(bx = bx_matrix,
                           bxse = bxse_matrix, 
                           by = outcome_beta,
                           byse = outcome_se)
  
  # Perform MVMR analysis
  mvmr_results <- mr_mvivw(model = model, mvmr_input)
  
  # Calculate diagnostics
  diagnostics <- calculate_diagnostics(bfa_fit, factor_scores)
  
  # Generate plots if requested
  plots <- NULL
  if (plot_diagnostics) {
    plots <- plot_mvmr_bfa_diagnostics(list(
      bfa_fit = bfa_fit,
      bfa_summary = bfa_summary,
      factor_scores = factor_scores
    ), diagnostics)
  }
  
  # Generate report if requested
  report <- NULL
  if (generate_report) {
    report <- generate_mvmr_bfa_report(list(
      factor_scores = factor_scores,
      bfa_summary = bfa_summary,
      mvmr_results = list(model = mvmr_results),
      convergence = list(
        mh_acceptance = mean(bfa_fit$MHacc)
      )
    ), diagnostics)
  }
  
  # Return comprehensive results
  results <- list(
    factor_scores = factor_scores,
    bfa_summary = bfa_summary,
    bfa_fit = bfa_fit,
    mvmr_results = list(
      model = mvmr_results,
      input = mvmr_input
    ),
    diagnostics = diagnostics,
    plots = plots,
    report = report,
    parameters = list(
      Kmax = Kmax,
      model = model,
      Nid = Nid,
      burnin = burnin,
      iter = iter,
      nu0 = nu0,
      kappa = kappa,
      kappa0 = kappa0,
      xi0 = xi0
    )
  )
  
  class(results) <- "mvmr_bfa"
  return(results)
}

# print method for mvmr_bfa objects
print.mvmr_bfa <- function(x, ...) {
  cat("\nMVMR-BFA Analysis Results\n")
  cat("=======================\n\n")
  
  cat("Factor Analysis Summary:\n")
  cat("- Number of factors:", ncol(x$factor_scores), "\n")
  cat("- Number of variables:", nrow(x$bfa_summary$alpha), "\n")
  cat("- MCMC iterations:", x$parameters$iter, "\n")
  cat("- Burn-in:", x$parameters$burnin, "\n\n")
  
  cat("MVMR Estimates:\n")
  estimates <- x$mvmr_results$model@Estimate
  se <- x$mvmr_results$model@StdError
  
  # Calculate p-values (two-sided test)
  z_scores <- estimates / se
  p_values <- 2 * (1 - pnorm(abs(z_scores)))
  
  # Format output with p-values
  for(i in seq_along(estimates)) {
    cat(sprintf("Factor %d: %.4f (SE: %.4f, p-value: %.4e)\n", 
                i, estimates[i], se[i], p_values[i]))
  }
  
  cat("\nDiagnostics Summary:\n")
  cat("- MH acceptance rate:", round(mean(x$bfa_fit$MHacc), 3), "\n")
  cat("- Mean factor dedication stability:", 
      round(mean(x$diagnostics$dedication_stability), 3), "\n")
}
