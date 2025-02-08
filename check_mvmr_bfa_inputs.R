#' Check if inputs are valid for MVMR_BFA
#' @param exposure_beta Matrix or dataframe of exposure variables
#' @param outcome_beta Vector of outcome betas
#' @param outcome_se Vector of outcome standard errors
#' @return List with validation results and messages
check_mvmr_bfa_inputs <- function(exposure_beta, outcome_beta, outcome_se) {
  messages <- list()
  is_valid <- TRUE
  
  # Check data types
  if (!is.matrix(exposure_beta) && !is.data.frame(exposure_beta)) {
    messages <- c(messages, "exposure_beta must be a matrix or data frame")
    is_valid <- FALSE
  }
  
  # Check dimensions
  if (length(outcome_beta) != length(outcome_se)) {
    messages <- c(messages, "outcome_beta and outcome_se must have same length")
    is_valid <- FALSE
  }
  
  if (nrow(exposure_beta) != length(outcome_beta)) {
    messages <- c(messages, "Number of rows in exposure_beta must match length of outcome vectors")
    is_valid <- FALSE
  }
  
  # Check for missing values
  if (any(is.na(exposure_beta))) {
    messages <- c(messages, "Warning: Missing values found in exposure_beta")
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
