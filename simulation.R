#' @title Generate exposure data
#' 
#' @description Creates simulated exposure data influenced by SNPs, with block structure 
#' and confounding effects.
#' 
#' @param n_snps Integer. Number of SNPs to generate
#' @param n_exposures Integer. Number of exposure variables to generate
#' @param n_blocks Integer. Number of blocks to partition SNPs and exposures into
#' @param n_samples Integer. Number of samples/individuals to simulate
#' @param heritX Numeric. Proportion of variance explained by SNPs for each exposure (default: 0.15)
#' @param compute_stats Logical. Whether to compute summary statistics (default: TRUE)
#' 
#' @return A list containing:
#'   \item{exposures}{Matrix of standardized exposure values}
#'   \item{snps}{Matrix of SNP genotypes}
#'   \item{confounder}{Vector of confounding values}
#'   If compute_stats=TRUE, also includes:
#'   \item{beta_gx}{Matrix of SNP-exposure effect sizes}
#'   \item{se_gx}{Matrix of standard errors for SNP-exposure effects}
#'   \item{snp_blocks}{List of SNP block indices}
#'   \item{exposure_blocks}{List of exposure block indices}
#'   \item{snp_effects}{Matrix of true SNP effects on exposures}
#'   \item{exposure_error}{Matrix of random error terms}
generate_exposure_data <- function(n_snps, n_exposures, n_blocks, n_samples, 
                                   heritX = 0.15, # Variance explained by SNPs for each exposure
                                   compute_stats = TRUE) {
  # Generate SNP data
  maf <- runif(n_snps, 0.01, 0.5)
  snp_matrix <- matrix(nrow = n_samples, ncol = n_snps)
  for(i in 1:n_snps) {
    snp_matrix[,i] <- rbinom(n_samples, 2, maf[i])
  }
  
  # Partition variance
  total_var <- 1
  confounder_prop <- 0.1  # Confounder explains 10% of variance
  genetic_prop <- heritX   # Genetic component explains heritX of remaining variance
  error_prop <- 1 - confounder_prop - genetic_prop
  
  if(error_prop < 0) {
    stop("heritX too large: heritX + confounder proportion must be less than 1")
  }
  
  # Generate scaled confounder
  confounder <- matrix(rnorm(n_samples, 0, sqrt(total_var * confounder_prop)))
  
  # Create blocks of SNPs
  block_boundaries <- seq(1, n_snps, length.out = n_blocks + 1)
  block_boundaries <- round(block_boundaries)
  snp_blocks <- list()
  for(i in 1:n_blocks) {
    snp_blocks[[i]] <- block_boundaries[i]:(block_boundaries[i+1] - 1)
  }
  
  # Create corresponding exposure blocks
  exposure_per_block <- ceiling(n_exposures / n_blocks)
  exposure_blocks <- list()
  for(i in 1:n_blocks) {
    start_idx <- (i-1) * exposure_per_block + 1
    end_idx <- min(i * exposure_per_block, n_exposures)
    exposure_blocks[[i]] <- start_idx:end_idx
  }
  
  # Generate SNP effects on exposures
  snp_effects <- matrix(0, nrow = n_snps, ncol = n_exposures)
  mean_effect_size <- sqrt((total_var * genetic_prop) / (n_snps/n_blocks))  # per block scaling
  
  for(i in 1:n_blocks) {
    for(j in exposure_blocks[[i]]) {
      snp_effects[snp_blocks[[i]], j] <- rnorm(
        length(snp_blocks[[i]]),
        mean = mean_effect_size,
        sd = mean_effect_size/5
      )
    }
  }
  
  # Generate exposures with scaled error
  exposure_error <- matrix(rnorm(n_samples * n_exposures, 0, sqrt(total_var * error_prop)), 
                           nrow = n_samples)
  exposures <- snp_matrix %*% snp_effects + exposure_error
  
  # Add confounder effect
  exposures <- apply(exposures, 2, function(x) x + confounder)
  
  # Standardize exposures to unit variance
  exposures <- scale(exposures)
  
  if(compute_stats) {
    # Compute summary statistics
    beta_gx <- se_gx <- matrix(nrow = n_snps, ncol = n_exposures)
    
    for(i in 1:n_exposures) {
      model <- summary(lm(exposures[,i] ~ snp_matrix))$coefficients[-1,]
      beta_gx[,i] <- model[,1]
      se_gx[,i] <- model[,2]
    }
    
    return(list(
      exposures = exposures,
      snps = snp_matrix,
      confounder = confounder,
      beta_gx = beta_gx,
      se_gx = se_gx,
      snp_blocks = snp_blocks,
      exposure_blocks = exposure_blocks,
      snp_effects = snp_effects,
      exposure_error = exposure_error
    ))
  }
  
  return(list(
    exposures = exposures,
    snps = snp_matrix,
    confounder = confounder
  ))
}


#' @title Generate outcome data with exposure effects and optional pleiotropy
#' 
#' @description Creates simulated outcome data influenced by exposures, with optional 
#' pleiotropic effects from SNPs and confounding.
#' 
#' @param exposures Matrix. Standardized exposure variables
#' @param snps Matrix. SNP genotypes
#' @param confounder Vector. Confounding values
#' @param r2_exposures Numeric. Proportion of variance explained by exposures (default: 0.20)
#' @param n_causal_exposures Integer. Number of causal exposures. If NULL, all exposures are causal (default: NULL)
#' @param pleiotropy_type Character. Type of pleiotropy: "none", "balanced", or "unbalanced" (default: "none")
#' @param pleiotropy_proportion Numeric. Proportion of SNPs with pleiotropic effects (default: 0)
#' @param r2_pleiotropy Numeric. Proportion of variance explained by pleiotropic effects (default: 0)
#' 
#' @return A list containing:
#'   \item{outcome}{Vector of standardized outcome values}
#'   \item{beta_gy}{Vector of SNP-outcome effect sizes}
#'   \item{se_gy}{Vector of standard errors for SNP-outcome effects}
#'   \item{causal_exposures}{Vector of indices for causal exposures}
#'   \item{exposure_effects}{Vector of true exposure effects}
#'   \item{pleiotropic_snps}{Vector of indices for pleiotropic SNPs (NULL if no pleiotropy)}
#'   \item{pleiotropic_effects}{Vector of true pleiotropic effects}
generate_outcome <- function(exposures, snps, confounder, 
                             r2_exposures = 0.20,  # Variance explained by exposures
                             n_causal_exposures = NULL,
                             pleiotropy_type = "none",  # "none", "balanced", or "unbalanced"
                             pleiotropy_proportion = 0,  # Proportion of SNPs with pleiotropic effects
                             r2_pleiotropy = 0  # Variance explained by pleiotropic effects
) {
  
  n_samples <- nrow(exposures)
  n_exposures <- ncol(exposures)
  n_snps <- ncol(snps)
  
  if(is.null(n_causal_exposures)) {
    n_causal_exposures <- n_exposures
  }
  
  # Validate inputs
  if(!pleiotropy_type %in% c("none", "balanced", "unbalanced")) {
    stop("pleiotropy_type must be one of: 'none', 'balanced', or 'unbalanced'")
  }
  if(pleiotropy_proportion < 0 || pleiotropy_proportion > 1) {
    stop("pleiotropy_proportion must be between 0 and 1")
  }
  
  # Calculate variance components
  total_var <- 1
  confounder_prop <- 0.1  # Confounder explains 10% of variance
  exposure_var <- r2_exposures * total_var  # Variance explained by exposures
  pleiotropic_var <- r2_pleiotropy * total_var  # Variance explained by pleiotropic effects
  error_var <- total_var * (1 - confounder_prop - r2_exposures - r2_pleiotropy)  # Remaining variance
  
  # Validate total variance
  if(error_var < 0) {
    stop("Total variance explained (confounder + exposures + pleiotropy) cannot exceed 1")
  }
  
  # Generate exposure effects
  exposure_effects <- rep(0, n_exposures)
  causal_indices <- sample(1:n_exposures, n_causal_exposures)
  mean_effect_size <- sqrt(exposure_var / n_causal_exposures)
  exposure_effects[causal_indices] <- rnorm(
    n_causal_exposures,
    mean = mean_effect_size,
    sd = mean_effect_size/5
  )
  
  # Initialize pleiotropic effects
  pleiotropic_effects <- rep(0, n_snps)
  
  # Generate pleiotropic effects if needed
  if(pleiotropy_type != "none" && pleiotropy_proportion > 0) {
    n_pleiotropic_snps <- round(n_snps * pleiotropy_proportion)
    pleiotropic_snps <- sample(1:n_snps, n_pleiotropic_snps)
    
    # Calculate effect size to achieve desired variance explained
    base_effect_size <- sqrt(pleiotropic_var / n_pleiotropic_snps)
    
    if(pleiotropy_type == "balanced") {
      # Balanced pleiotropy: mean zero effect
      pleiotropic_effects[pleiotropic_snps] <- rnorm(
        n_pleiotropic_snps,
        mean = 0,
        sd = base_effect_size
      )
    } else {  # unbalanced
      # Unbalanced pleiotropy: all positive effects
      pleiotropic_effects[pleiotropic_snps] <- rnorm(
        n_pleiotropic_snps,
        mean = base_effect_size,
        sd = base_effect_size/5
      )
    }
  }
  
  # Generate outcome
  outcome_error <- rnorm(n_samples, 0, sqrt(error_var))
  outcome <- exposures %*% exposure_effects + 
    snps %*% pleiotropic_effects + 
    confounder + 
    outcome_error
  
  # Standardize outcome to unit variance
  outcome <- scale(outcome)
  
  # Compute SNP-outcome associations
  model <- summary(lm(outcome ~ snps))$coefficients[-1,]
  beta_gy <- model[,1]
  se_gy <- model[,2]
  
  return(list(
    outcome = outcome,
    beta_gy = beta_gy,
    se_gy = se_gy,
    causal_exposures = causal_indices,
    exposure_effects = exposure_effects,
    pleiotropic_snps = if(pleiotropy_type != "none") pleiotropic_snps else NULL,
    pleiotropic_effects = pleiotropic_effects
  ))
}
