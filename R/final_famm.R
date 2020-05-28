
#------------------------------------------------------------------------------#
# Update Model Formula to Include Weighted Eigenfunctions
#------------------------------------------------------------------------------#
create_formula <- function(data, MFPC, bs, bf_mean, m_mean){

  # Arguments
  # data        : Data that contains all the necessary information provided by
  #                 prepare_gam()
  # MFPC        : List containing MFPC objects for each coponent of the variance
  # bs          : Spline basis function, only tested for "ps" (as in sparseFLMM)
  # bf_mean     : Basis dimension for functional intercept (as in sparseFLMM)
  # m_mean      : Order of penalty for basis function (as in sparseFLMM)


  # Prune down the list according to the included variance components
  fRI_list = list(E = "n_long",
                  B = "subject_long",
                  C = "word_long")
  fRI_list <- fRI_list[intersect(names(MFPC), names(fRI_list))]

  # Extract info of how many weighted Eigenfunctions are to be included
  n <- lapply(MFPC, function(x) length(x$values))

  # Paste the additional formula terms
  # Combine the names of the component with the variable name for the level
  new_terms <- lapply(seq_along(n), function(i){
    term <- paste0("w", names(n)[[i]], "_", 1:n[[i]])
  })
  new_terms <- mapply(function(x, y){
    paste0("s(", x, ", ", paste0(y, collapse = ", "), ", bs = \"pcre\")")
  }, fRI_list, new_terms, SIMPLIFY = FALSE)

  # Create the formula components
  # No scalar intercepts necessary with metric by-variables
  covariates <- paste(attr(data, "gam_terms"), collapse = " + ")
  random_effects <- do.call(paste, list(new_terms, collapse = " + "))

  # Components plus terms that are always the same
  formula <- as.formula(paste0("y_vec ~ 0 + dim + ",
                               "s(t, by = dim, k = ", bf_mean, ", bs = \"", bs,
                               "\", m = ",
                               paste0("c(", paste(m_mean, collapse = ","), ")"),
                               ") + ", covariates, " + ", random_effects))
  formula

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Compute Final Model on All the Dimensions
#------------------------------------------------------------------------------#
final_model <- function(formula, data, final_method, model_list, weight_refit){

  # Arguments
  # formula     : Formula for model estimation on all dimensnions
  # data        : Data that contains all the necessary information provided by
  #                 prepare_gam()
  # final_method: Function used for estimation of final model to allow for
  #                 potential heterscedasticity ("bam", "w_bam", "gamm",
  #                 "gaulss")
  # model_list  : List containing sparseFLMM objects for each dimension


  switch(final_method,

         # Assumption: Homoscedasticity
         "bam" = mgcv::bam(formula = formula, data = data, discrete = TRUE),

         # Assumption: Heteroscedasticity depending on dimension
         "w_bam" = {
           if (!weight_refit) {
             weights <- sapply(model_list, function(x) {
               sig <- grep("^cov_hat_", names(x))
               1 / x[[sig]]$sigmasq
             })
           } else {
             # refit the model to get an update of the sigma^2 estimate
             weights <- refit_for_weights(formula = formula, data = data,
                                          model_list = model_list)
           }
           weights <- rep(weights, times = table(data$dim))
           data$norm_weights <- weights/mean(weights)
           mgcv::bam(formula = formula, data = data, weights = norm_weights,
                     discrete = TRUE)
         },

         # Assumption: Heteroscedasticity depending on dimension
         # Could be expanded to incorporate a general residual variance structure
         # representable with nlme
         "gamm" = mgcv::gamm(formula = formula, data = data,
                             weights = varIdent(form = ~ 1|dim)),

         # Assumption: Gaussian location scale model family
         # Could be expanded to incorporate more covariates
         "gaulss" = mgcv::gam(formula = list (formula, ~ 0 + dim), data = data,
                              family = gaulss())
  )

}
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# Compute Final Model on All the Dimensions
#------------------------------------------------------------------------------#
#' Refit the model under an independence assumption
#'
#' This is an internal function. Refit the model under an independence
#' assumption now with the basis functions from the MFPCA. Goal is to extract an
#' estimate for the error variances.
#'
#' @param formula Formula to fit the final model.
#' @param data Data that contains all the variables specified in formula.
#' @param model_list List containing sparseFLMM objects for each dimension
refit_for_weights <- function(formula, data, model_list){

  # remove the by argument from the smooth terms
  formula_indep <- paste(formula)
  formula_indep[3] <- gsub(", by = dim[.a-z]*[.1-9]*", "", formula_indep[3])
  formula_indep[3] <- sub("0 \\+ dim \\+", "", formula_indep[3])
  formula_indep <- as.formula(paste(formula_indep[2], formula_indep[1],
                                    formula_indep[3]))

  # fit the model on a subdataset and extract the sigma estimate
  weights <- sapply(names(model_list), function (x) {
    dat <- subset(data, dim == x)
    1 / mgcv::bam(formula = formula_indep, data = dat, discrete = TRUE)$sig2
  })
  weights

}
