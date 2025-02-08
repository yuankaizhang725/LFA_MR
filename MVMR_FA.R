#' Perform parallel analysis to determine number of factors
#' @param data Matrix or data frame of variables
#' @param n.iter Number of iterations for parallel analysis
#' @param centile Percentile to use for comparison
#' @return Number of factors to retain
parallel_analysis <- function(data, n.iter = 1000, centile = 95) {
  # Get dimensions
  n <- nrow(data)
  p <- ncol(data)
  
  # Calculate eigenvalues of original correlation matrix
  orig.r <- cor(data)
  orig.values <- eigen(orig.r)$values
  
  # Generate random data and get eigenvalues
  random.values <- matrix(0, nrow = n.iter, ncol = p)
  for(i in 1:n.iter) {
    random.data <- matrix(rnorm(n * p), nrow = n)
    random.r <- cor(random.data)
    random.values[i,] <- eigen(random.r)$values
  }
  
  # Get percentile values
  pa.values <- apply(random.values, 2, function(x) quantile(x, centile/100))
  
  # Determine number of factors
  n.factors <- sum(orig.values > pa.values)
  
  return(list(
    n.factors = n.factors,
    orig.values = orig.values,
    pa.values = pa.values
  ))
}

#' Factor Analysis for MVMR using ML estimation
#' @param exposure_beta Matrix of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @param n_factor Optional, specific number of factors to extract
#' @param rotation Whether to apply varimax rotation
#' @param parallel_iter Number of iterations for parallel analysis
#' @param parallel_centile Percentile for parallel analysis
#' @param tol Convergence tolerance for ML iterations
#' @return MVMR results with factor analysis
MVMR_FA <- function(exposure_beta, 
                    outcome_beta,
                    outcome_se,
                    n_factor = NULL,
                    rotation = TRUE,
                    parallel_iter = 1000,
                    parallel_centile = 95,
                    min.eigval = 0.7,
                    standardize = TRUE,
                    maxit = 100,
                    tol = 1e-4) {  # tol: stop iteration when change in uniquenesses < tol
  
  # Input validation
  if(!is.matrix(exposure_beta) && !is.data.frame(exposure_beta)) {
    stop("exposure_beta must be a matrix or data frame")
  }
  
  if(length(outcome_beta) != length(outcome_se)) {
    stop("outcome_beta and outcome_se must have same length")
  }
  
  if(nrow(exposure_beta) != length(outcome_beta)) {
    stop("Number of rows in exposure_beta must match length of outcome vectors")
  }
  
  # Standardize data if requested
  if(standardize) {
    exposure_scaled <- scale(exposure_beta)
  } else {
    exposure_scaled <- as.matrix(exposure_beta)
  }
  
  # Calculate correlation matrix
  R <- cor(exposure_scaled)
  p <- ncol(R)
  
  # Determine number of factors
  if(is.null(n_factor)) {
    # Perform parallel analysis only if n_factor not specified
    pa_results <- parallel_analysis(
      exposure_scaled, 
      n.iter = parallel_iter,
      centile = parallel_centile
    )
    
    n_factors <- pa_results$n.factors
    
    # Adjust number of factors if any eigenvalues are too small
    eigen_values <- pa_results$orig.values[1:n_factors]
    if(any(eigen_values < min.eigval)) {
      n_factors <- max(which(eigen_values >= min.eigval))
    }
  } else {
    n_factors <- n_factor
    pa_results <- list(
      orig.values = eigen(R)$values,
      pa.values = rep(NA, ncol(R))
    )
  }
  
  # ML Factor Analysis
  # Initial communality estimates using squared multiple correlations
  diag_inv_R <- solve(R)
  h2 <- 1 - 1/diag(diag_inv_R)
  U2 <- 1 - h2
  Udiag <- diag(U2)
  
  # ML iterations
  for(iter in 1:maxit) {
    # Store old uniquenesses
    old_U2 <- U2
    
    # Adjusted correlation matrix
    Radj <- R - Udiag
    
    # Eigendecomposition of adjusted correlation matrix
    eig <- eigen(Radj)
    
    # Get loadings from first nfactors eigenvectors
    loadings <- eig$vectors[, 1:n_factors] %*% 
      diag(sqrt(pmax(eig$values[1:n_factors], 0)))
    
    # Update communalities and uniquenesses
    h2 <- rowSums(loadings^2)
    U2 <- 1 - h2
    Udiag <- diag(U2)
    
    # Check convergence
    if(max(abs(U2 - old_U2)) < tol) break  # Stop if change in uniquenesses is small
  }
  
  # Apply varimax rotation if requested
  if(rotation) {
    # Normalize
    h2 <- rowSums(loadings^2)
    loadings <- loadings / sqrt(h2)
    
    # Initialize rotated loadings
    rotated <- loadings
    
    # Varimax criterion
    for (i in 1:50) { # Maximum 50 iterations
      d <- 0
      for (j in 1:(n_factors-1)) {
        for (k in (j+1):n_factors) {
          u <- rotated[,j]
          v <- rotated[,k]
          
          a <- sum(u^2 - v^2)
          b <- 2 * sum(u * v)
          
          if (abs(b) > 1e-10) {
            theta <- atan2(b, a) / 4
            c <- cos(theta)
            s <- sin(theta)
            
            # Update columns j and k
            temp <- rotated[,j]
            rotated[,j] <- c * u - s * v
            rotated[,k] <- s * temp + c * v
            
            d <- d + abs(theta)
          }
        }
      }
      if (d < 1e-6) break  # Convergence criterion
    }
    
    # Denormalize
    rotated <- rotated * sqrt(h2)
    loadings <- rotated
  }
  
  # Calculate residual correlation matrix
  residual <- R - (loadings %*% t(loadings) + Udiag)
  
  # Calculate fit measures
  chi_square <- (nrow(R) - 1 - (2*p + 5)/6 - (2*n_factors)/3) * sum(residual^2)
  df <- ((p - n_factors)^2 - p - n_factors)/2
  p_value <- 1 - pchisq(chi_square, df)
  
  # Store ML results
  ml_results <- list(
    uniquenesses = U2,
    communalities = h2,
    residual = residual,
    converged = iter < maxit,
    iterations = iter,
    chi_square = chi_square,
    df = df,
    p_value = p_value
  )
  
  # Calculate factor scores using regression method
  weights <- solve(R) %*% loadings
  scores <- exposure_scaled %*% weights
  
  # Prepare matrices for MVMR
  bx_matrix <- as.matrix(scores)
  bxse_matrix <- matrix(0, nrow = nrow(scores), ncol = ncol(scores))
  
  # Create MVMR input object
  mvmr_input <- mr_mvinput(bx = bx_matrix,
                           bxse = bxse_matrix, 
                           by = outcome_beta,
                           byse = outcome_se)
  
  # Perform MVMR analysis
  mvmr_results <- mr_mvivw(mvmr_input)
  
  # Prepare results object
  results <- list(
    factor_analysis = list(
      loadings = loadings,
      factor_scores = scores,
      n_factors = n_factors,
      parallel_analysis = pa_results,
      ml_results = ml_results,
      rotated = rotation,
      variable_names = names(exposure_beta)
    ),
    mvmr_results = list(
      model = mvmr_results,
      input = mvmr_input
    ),
    parameters = list(
      rotation = rotation,
      parallel_iter = parallel_iter,
      parallel_centile = parallel_centile,
      min.eigval = min.eigval,
      standardize = standardize,
      maxit = maxit,
      tol = tol,
      specified_n_factor = !is.null(n_factor)
    )
  )
  
  class(results) <- "mvmr_fa"
  return(results)
}

# Print method for mvmr_fa objects
print.mvmr_fa <- function(x, ..., digits = 3) {
  cat("\nMVMR-FA Analysis Results (Maximum Likelihood)\n")
  cat("=======================================\n\n")
  
  cat("Factor Analysis Summary:\n")
  cat("- Number of factors:", x$factor_analysis$n_factors, "\n")
  cat("- Number of variables:", nrow(x$factor_analysis$loadings), "\n")
  cat("- Rotation applied:", x$parameters$rotation, "\n")
  cat("- ML iterations:", x$factor_analysis$ml_results$iterations, "\n\n")
  
  cat("Model Fit:\n")
  cat("Chi-square test of model fit:\n")
  cat("- H0: The factor model explains the correlations\n")
  cat("- H1: The factor model does not explain the correlations\n")
  cat("- Chi-square:", round(x$factor_analysis$ml_results$chi_square, 3), "\n")
  cat("- Degrees of freedom:", x$factor_analysis$ml_results$df, "\n")
  cat("- P-value:", format.pval(x$factor_analysis$ml_results$p_value), "\n")
  cat("Note: P-value > 0.05 suggests good model fit\n\n")
  
  cat("Factor Loadings Matrix:\n")
  loadings <- x$factor_analysis$loadings
  if(!is.null(x$factor_analysis$variable_names)) {
    rownames(loadings) <- x$factor_analysis$variable_names
  }
  colnames(loadings) <- paste0("Factor", 1:ncol(loadings))
  print(round(loadings, digits))
  cat("\n")
  
  cat("Communalities:\n")
  communalities <- x$factor_analysis$ml_results$communalities
  names(communalities) <- x$factor_analysis$variable_names
  print(round(communalities, digits))
  cat("\n")
  
  cat("MVMR Estimates:\n")
  estimates <- x$mvmr_results$model@Estimate
  se <- x$mvmr_results$model@StdError
  
  # Calculate p-values
  z_scores <- estimates / se
  p_values <- 2 * (1 - pnorm(abs(z_scores)))
  
  for(i in seq_along(estimates)) {
    cat(sprintf("Factor %d: %.4f (SE: %.4f, p-value: %.4e)\n", 
                i, estimates[i], se[i], p_values[i]))
  }
}

# Plot method for mvmr_fa objects
plot.mvmr_fa <- function(x, exposure_names = NULL, ...) {
  # Create list to store plots
  plots <- list()
  
  # Scree plot with parallel analysis
  scree_data <- data.frame(
    Component = 1:length(x$factor_analysis$parallel_analysis$orig.values),
    Eigenvalue = x$factor_analysis$parallel_analysis$orig.values,
    PA_Threshold = x$factor_analysis$parallel_analysis$pa.values
  )
  
  plots$scree <- ggplot(scree_data, aes(x = Component)) +
    geom_line(aes(y = Eigenvalue, color = "Observed"), size = 1) +
    geom_line(aes(y = PA_Threshold, color = "Parallel Analysis"), size = 1, linetype = "dashed") +
    geom_vline(xintercept = x$factor_analysis$n_factors, linetype = "dotted") +
    theme_minimal() +
    labs(title = "Scree Plot with Parallel Analysis",
         y = "Eigenvalue") +
    scale_color_manual(values = c("Observed" = "blue", "Parallel Analysis" = "red")) +
    theme(legend.title = element_blank())
  
  # Factor loadings heatmap
  loadings_data <- as.data.frame(x$factor_analysis$loadings)
  
  loadings_data$Variable <- x$factor_analysis$variable_names
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)) 
  
  # Factor scores distribution
  scores_data <- as.data.frame(x$factor_analysis$factor_scores)
  scores_long <- tidyr::pivot_longer(
    scores_data,
    cols = everything(),
    names_to = "Factor",
    values_to = "Score"
  )
  
  plots$scores <- ggplot(scores_long, aes(x = Score, fill = Factor)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Factor) +
    theme_minimal() +
    labs(title = "Factor Scores Distribution")
  
  return(plots)
}
