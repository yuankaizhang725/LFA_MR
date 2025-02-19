# LFA_MR

## Overview

This repository contains the R implementation of integration of latent factor analysis into multivariable Mendelian randomization (MVMR), including MVMR with factor analysis (MVMR-FA) and multivariable MR with Bayesian factor analysis (MVMR-BFA). These methods are designed to handle correlated exposures in Mendelian Randomization studies by leveraging factor analysis to identify latent factors.

## Installation

```{r}
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install LFA_MR from GitHub
devtools::install_github("yuankaizhang725/LFA_MR")
```

## Methods Overview

### MVMR-FA (Multivariable MR with Factor Analysis)
MVMR-FA uses standard factor analysis to extract latent factors from correlated exposures and performs MR on these factors to estimate causal effects.

### MVMR-BFA (Multivariable MR with Bayesian Factor Analysis)
MVMR-BFA employs Bayesian factor analysis, providing sparse loading matrix, especially useful for complex correlation structures.

## Usage Examples

## Basic Usage of MVMR_FA

The `MVMR_FA` function implements multivariable Mendelian randomization with standard factor analysis.

### Simple Usage with Automatic Factor Determination

```{r mvmr_fa_auto}
# MVMR_FA with automatic factor determination through parallel analysis
mvmr_fa_results_auto <- MVMR_FA(
  exposure_beta = betaX,
  outcome_beta = betaY,
  outcome_se = seY,
  # Let parallel analysis determine the number of factors
  n_factor = NULL,
  rotation = TRUE,
  parallel_iter = 1000,
  parallel_centile = 95,
  min.eigval = 0.7,
  standardize = TRUE
)

# View summary of results
print.mvmr_fa(mvmr_fa_results_auto)
```

### Specifying Number of Factors

```{r mvmr_fa_specific}
# MVMR_FA with specified number of factors
mvmr_fa_results <- MVMR_FA(
  exposure_beta = betaX,
  outcome_beta = betaY,
  outcome_se = seY,
  n_factor = 2,            # Specify exactly 2 factors
  rotation = TRUE,         # Apply varimax rotation
  standardize = TRUE,      # Standardize variables
  maxit = 100,             # Maximum iterations for ML estimation
  tol = 1e-4               # Convergence tolerance
)
```

### Visualize factor analysis results

```{r}
plots <- plot.mvmr_fa(mvmr_fa_results_advanced)

# Scree plot
print(plots$scree)

# Factor loadings heatmap
print(plots$loadings)

# Factor scores distribution
print(plots$scores)
```

## Basic Usage of MVMR_BFA

The `MVMR_BFA` function implements multivariable Mendelian randomization with Bayesian factor analysis, providing sparse loading matrix.

### Simple Usage

```{r mvmr_bfa_simple}
mvmr_bfa_results <- MVMR_BFA(
  exposure_beta = betaX,
  outcome_beta = betaY,
  outcome_se = seY,
  Kmax = 5,                  # Maximum factors to consider
  model = "fixed",
  Nid = 2,                   # Minimum variables per factor
  burnin = 2000,
  iter = 10000,
  nu0 = 7,                   # Custom prior degrees of freedom (Kmax+2)
  kappa = 0.15,              # Custom shrinkage parameter
  kappa0 = 3,                # Custom prior correlation parameter
  xi0 = 0.8,                 # Custom prior scale parameter
  seed = 725,
  plot_diagnostics = TRUE,
  generate_report = TRUE
)

# Show results
print.mvmr_bfa(mvmr_bfa_results)

# Plot posterior distributions of causal effects
plot(mvmr_bfa_results, type = "causal")

# Plot factor loadings with credible intervals
plot(mvmr_bfa_results, type = "loadings")
```


## Sensitivity Analysis

Sensitivity analysis helps assess the robustness of results to different modeling choices.

### 1. Set up parameter ranges to test

```{r setup_sensitivity}
Kmax_range <- 3:5            # Test different maximum number of factors
models <- c("fixed", "random") # Test both fixed and random effects models
Nid_range <- 2:3             # Test different minimum variables per factor
```

### 2. Run sensitivity analysis

```{r run_sensitivity}
sensitivity_results <- perform_sensitivity_analysis(
  exposure_beta = betaX,     # Your exposure data matrix
  outcome_beta = betaY,      # Your outcome beta vector
  outcome_se = seY,          # Your outcome standard error vector
  Kmax_range = Kmax_range,
  models = models,
  Nid_range = Nid_range
)
```

### 3. Compare results

```{r compare_sensitivity}
comparison <- compare_sensitivity_results(sensitivity_results)

# View metrics dataframe
print(comparison$metrics)
```

### 4. Plot comparisons

```{r plot_comparisons, fig.width=10, fig.height=6}
# Number of factors vs Kmax plot
print(comparison$plots$factors)

# Mean estimates comparison plot
print(comparison$plots$estimates)
```

### 5. Run cross-validation

```{r cross_validation}
cv_results <- calculate_optimal_parameters(
  exposure_beta = betaX,
  outcome_beta = betaY,
  outcome_se = seY,
  Kmax_range = 3:8,
  nfolds = 5
)

# Check results
print(cv_results$optimal_Kmax)   # Optimal number of maximum factors
print(cv_results$mean_cv_error)  # Mean cross-validation error for each Kmax
```


## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or support, please open an issue on GitHub or contact the repository maintainer at yuankaizhang725@github.com.
