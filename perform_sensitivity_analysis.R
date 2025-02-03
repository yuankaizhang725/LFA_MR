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