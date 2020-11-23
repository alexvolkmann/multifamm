#------------------------------------------------------------------------------#
# Implementation of sparseFLMM on Each Dimension
#------------------------------------------------------------------------------#
apply_sparseFLMM <- function(fRI_B, fRI_C, nested, bs, bf_mean, bf_covariates,
                             m_mean, covariate, num_covariates, covariate_form,
                             interaction, which_interaction, bf_covs, m_covs,
                             var_level, use_famm, save_model_famm, one_dim,
                             data, ...){

  # Arguments
  # fRI_B             : Boolean for including functional random intercept
  #                       for individual (B in Cederbaum)
  # fRI_C             : Boolean for including functional random intercept
  #                       for word (C in Cederbaum)
  # bs                : Spline basis function, only tested for "ps" (as in
  #                       sparseFLMM)
  # bf_mean           : Basis dimension for functional intercept (as in
  #                       sparseFLMM)
  # bf_covariates     : Basis dimension for all covariates (as in sparseFLMM)
  # m_mean            : Order of penalty for basis function (as in sparseFLMM)
  # covariate         : Covariate effects (as in sparseFLMM)
  # num_covariates    : Number of covariates included in the model (as in
  #                       sparseFLMM)
  # covariate_form    : Vector of strings for type of covariate (as in
  #                       sparseFLMM)
  # interaction       : TRUE if there are interactions between covariates (as in
  #                       sparseFLMM)
  # which_interaction : Symmetric matrix specifying the interaction terms (as in
  #                       sparseFLMM)
  # bf_covs           : Vector of marginal basis dimensions for fRI covariance
  #                       estimation (as in sparseFLMM)
  # m_covs            : list of marginal orders for the penaltyin fRI covariance
  #                       estimation (as in sparseFLMM)
  # var_level         : Pre-specified level of explained variance on each
  #                       dimension (as in sparseFLMM)
  # use_famm          : Re-estimate the mean in FAMM context (as in sparseFLMM)
  # save_model_famm   : Give out the FAMM model object (as in sparseFLMM)
  # one_dim           : Specify the name of the dimension if sparseFLMM is to be
  #                       computed only on one dimension
  # data              : Data that contains the information with some fixed
  #                       varnames

  # Data.tables changes even outside of function so work with copy
  # Not memory efficient but other variable names outside of function
  data <- data.table::copy(data)

  # Only change variable names if they are specified
  # Also transform information of fRI to arguments for sparseFLMM
  if(fRI_B){
    if(fRI_C){
      # Also nested fRI
      use_RI <- FALSE
      use_simple <- FALSE
    }else{
      # Only fRI for individual
      use_RI <- TRUE
      use_simple <- FALSE
    }
  }else{
    # Independent curves
    # Rather unintuitively, sparseFLMM estimates E as B
    data.table::set(data, j = "subject_long", value = data$n_long)
    use_RI <- TRUE
    use_simple <- TRUE
  }

  # If only one dimension is to be computed, subset the data to the dimension
  if (! is.null(one_dim)) {
    data <- subset(data, dim == one_dim)
    data.table::set(data, j = "dim", value = droplevels(data$dim))
    use_famm <- TRUE
    save_model_famm <- TRUE
  }

  # Apply sparseFLMM to every dimension of the data
  out <- lapply(levels(data$dim), sparseFLMM_one_d, data = data,
                use_RI = use_RI, use_simple = use_simple,
                bs = bs, bf_mean = bf_mean, bf_covariates = bf_covariates,
                m_mean = m_mean, covariate = covariate,
                num_covariates = num_covariates,
                covariate_form = covariate_form, interaction = interaction,
                which_interaction = which_interaction, bf_covs = bf_covs,
                m_covs = m_covs, var_level = var_level, use_famm = use_famm,
                save_model_famm = save_model_famm, nested = nested, ... = ...)
  names(out) <- levels(data$dim)

  out

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Execution of sparseFLMM on One Dimension
#------------------------------------------------------------------------------#
sparseFLMM_one_d <- function(dimension, data, use_RI, use_simple, bs, bf_mean,
                             bf_covariates, m_mean, covariate, num_covariates,
                             covariate_form, interaction, which_interaction,
                             bf_covs, m_covs, var_level, use_famm,
                             save_model_famm, nested, ...){

  # Arguments
  # dimension         : String describing the dimension to analyse
  # data              : Data that contains the information with some fixed
  #                       varnames
  # use_RI            : FALSE - crossed fRI, TRUE - simpler structure (as in
  #                       sparseFLMM)
  # use_simple        : TRUE - only frI for curve, FALSE - also frI for subject
  #                       (as in sparseFLMM)
  # bs                : Spline basis function, only tested for "ps" (as in
  #                       sparseFLMM)
  # bf_mean           : Basis dimension for functional intercept (as in
  #                       sparseFLMM)
  # bf_covariates     : Basis dimension for all covariates (as in sparseFLMM)
  # m_mean            : Order of penalty for basis function (as in sparseFLMM)
  # covariate         : Covariate effects (as in sparseFLMM)
  # num_covariates    : Number of covariates included in the model (as in
  #                       sparseFLMM)
  # covariate_form    : Vector of strings for type of covariate (as in
  #                       sparseFLMM)
  # interaction       : TRUE if there are interactions between covariates (as in
  #                       sparseFLMM)
  # which_interaction : Symmetric matrix specifying the interaction terms (as in
  #                       sparseFLMM)
  # bf_covs           : Vector of marginal basis dimensions for fRI covariance
  #                       estimation (as in sparseFLMM)
  # m_covs            : list of marginal orders for the penaltyin fRI covariance
  #                       estimation (as in sparseFLMM)
  # var_level         : Pre-specified level of explained variance on each
  #                       dimension (as in sparseFLMM)
  # use_famm          : Re-estimate the mean in FAMM context (as in sparseFLMM)
  # save_model_famm   : Give out the FAMM model object (as in sparseFLMM)

  # Subset the dimension data
  data <- subset(data, dim == dimension)

  # sparseFLMM
  cat("--------------------------------------\n")
  cat("sparseFLMM for Dimension", dimension, "\n")
  cat("--------------------------------------\n")
  sparseFLMM::sparseFLMM(curve_info = data,
                         use_RI = use_RI,
                         use_simple = use_simple,
                         bs = bs,
                         bf_mean = bf_mean,
                         bf_covariates = bf_covariates,
                         m_mean = m_mean,
                         covariate = covariate,
                         num_covariates = num_covariates,
                         covariate_form = covariate_form,
                         interaction = interaction,
                         which_interaction = which_interaction,
                         bf_covs = bf_covs,
                         m_covs = m_covs,
                         var_level = var_level,
                         use_famm = use_famm,
                         save_model_famm = save_model_famm,
                         nested = nested,
                         ... = ...)

}
#------------------------------------------------------------------------------#
