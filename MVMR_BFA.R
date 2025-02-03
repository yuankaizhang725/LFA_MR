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
#' @param exposure_data Matrix or data frame of exposure variables
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
MVMR_BFA <- function(exposure_data, 
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
  validation <- check_mvmr_bfa_inputs(exposure_data, outcome_beta, outcome_se)
  if (!validation$valid) {
    stop(paste(validation$messages, collapse = "\n"))
  }
  if (length(validation$messages) > 0) {
    warning(paste(validation$messages, collapse = "\n"))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # Scale exposure data
  exposure_scaled <- scale(exposure_data)
  
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

# Utility functions and methods for MVMR with Bayesian Factor Analysis

#' Check if inputs are valid for MVMR_BFA
#' @param exposure_data Matrix or dataframe of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @return List with validation results and messages
check_mvmr_bfa_inputs <- function(exposure_data, outcome_beta, outcome_se) {
  messages <- list()
  is_valid <- TRUE
  
  # Check data types
  if (!is.matrix(exposure_data) && !is.data.frame(exposure_data)) {
    messages <- c(messages, "exposure_data must be a matrix or data frame")
    is_valid <- FALSE
  }
  
  # Check dimensions
  if (length(outcome_beta) != length(outcome_se)) {
    messages <- c(messages, "outcome_beta and outcome_se must have same length")
    is_valid <- FALSE
  }
  
  if (nrow(exposure_data) != length(outcome_beta)) {
    messages <- c(messages, "Number of rows in exposure_data must match length of outcome vectors")
    is_valid <- FALSE
  }
  
  # Check for missing values
  if (any(is.na(exposure_data))) {
    messages <- c(messages, "Warning: Missing values found in exposure_data")
  }
  
  if (any(is.na(outcome_beta)) || any(is.na(outcome_se))) {
    messages <- c(messages, "Warning: Missing values found in outcome data")
  }
  
  # Check for zero or negative standard errors
  if (any(outcome_se <= 0)) {
    messages <- c(messages, "Error: Standard errors must be positive")
    is_valid <- FALSE
  }
  
  return(list(valid = is_valid, messages = messages))
}

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

#' Plot diagnostics for MVMR_BFA results
#' @param mvmr_bfa_results Results from MVMR_BFA function
#' @param diagnostics Diagnostic statistics from calculate_diagnostics
#' @return List of ggplot objects
plot_mvmr_bfa_diagnostics <- function(mvmr_bfa_results, diagnostics) {
  library(ggplot2)
  library(tidyr)
  
  plots <- list()
  
  # Factor loadings heatmap from alpha_means matrix
  # Get number of exposures and factors
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
  rownames(alpha_means) <- rownames(mvmr_bfa_results$bfa_summary$alpha)
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

#' Perform sensitivity analysis for MVMR_BFA
#' @param exposure_data Matrix of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @param Kmax_range Vector of Kmax values to try
#' @param other parameters as in MVMR_BFA
#' @return List of results for different parameter settings
perform_sensitivity_analysis <- function(exposure_data, outcome_beta, outcome_se,
                                         Kmax_range = 3:8, 
                                         models = c("fixed", "random"),
                                         Nid_range = 2:3) {
  results <- list()
  
  for(Kmax in Kmax_range) {
    for(model in models) {
      for(Nid in Nid_range) {
        key <- paste("Kmax", Kmax, "model", model, "Nid", Nid, sep="_")
        
        results[[key]] <- try({
          MVMR_BFA(exposure_data = exposure_data,
                   outcome_beta = outcome_beta,
                   outcome_se = outcome_se,
                   Kmax = Kmax,
                   model = model,
                   Nid = Nid)
        })
      }
    }
  }
  
  return(results)
}

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


# Print method for mvmr_bfa objects with p-values
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
