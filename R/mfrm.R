################################################################################
################################################################################
##                                                                            ##
##       Implementation of the Multivariate Functional Regression Model       ##
##                                                                            ##
################################################################################
################################################################################


#------------------------------------------------------------------------------#
# MultiVariate Functional Regression
#------------------------------------------------------------------------------#
#' Multivariate Functional Additive Mixed Model Regression
#'
#' This is the main function of the package and fits the multivariate functional
#' additive regression model with potentially nested or crossed functional
#' random intercepts.
#'
#' Expand the method proposed by Fabian Scheipl to incorporate the variance
#' decomposition developed by Cederbaum et al. (2016). To account for the
#' correlation between the dimensions, the MFPCA approach by Happ and Greven
#' (2016) is applied.
#'
#' The data set has to be of the following format:
#' \itemize{
#'   \item y_vec (numeric): vector of response values
#'   \item t (numeric): observation point locations
#'   \item n_long (integer): curve identification
#'   \item subject_long (integer): subject identification (NEEDS TO BE
#'     SPECIFIED)
#'   \item word_long (integer): word identification
#'   \item combi_long (integer): repetition
#'   \item dim (factor): level of the dimension
#'   \item covariate.1 (numeric): potential covariate(s) named with trailing
#'     1,2,3,...
#' }
#'
#' @param data Data.table that contains the information with some fixed variable
#'   names, see Details.
#' @param fRI_B Boolean for including functional random intercept for individual
#'   (B in Cederbaum). Defaults to \code{FALSE}.
#' @param fRI_C Boolean for including functional random intercept
#'   for word (C in Cederbaum). Defaults to \code{FALSE}.
#' @param nested \code{TRUE} to specify a model with nested functional random
#'   intercepts for the first and second grouping variable and a smooth error
#'   curve. Defaults to \code{FALSE}.
#' @param bs Spline basis function, only tested for "ps" (as in sparseFLMM).
#' @param bf_mean Basis dimension for functional intercept (as in sparseFLMM).
#' @param bf_covariates Basis dimension for all covariates (as in sparseFLMM).
#' @param m_mean Order of penalty for basis function (as in sparseFLMM).
#' @param covariate Covariate effects (as in sparseFLMM).
#' @param num_covariates Number of covariates included in the model (as in
#'   sparseFLMM).
#' @param covariate_form Vector of strings for type of covariate (as in
#'   sparseFLMM).
#' @param interaction TRUE if there are interactions between covariates (as in
#'   sparseFLMM). Defaults to \code{FALSE}.
#' @param which_interaction Symmetric matrix specifying the interaction terms
#'   (as in sparseFLMM).
#' @param bf_covs Vector of marginal basis dimensions for fRI covariance
#'   estimation (as in sparseFLMM).
#' @param m_covs List of marginal orders for the penaltyin fRI covariance
#'   estimation (as in sparseFLMM).
#' @param var_level Pre-specified level of explained variance on each
#'   dimension (as in sparseFLMM).
#' @param use_famm Re-estimate the mean in FAMM context (as in sparseFLMM) -
#'   overwritten by one_dim.
#' @param save_model_famm Give out the FAMM model object (as in sparseFLMM) -
#'   overwritten by one_dim.
#' @param one_dim Specify the name of the dimension if sparseFLMM is to be
#'   computed only on one dimension.
#' @param mfpc_cutoff Pre-specified level of explained variance of results of
#'   MFPCA (1 for previous versions of mvfr).
#' @param number_mfpc List containing the number of mfPCs needed for each
#'   variance component e.g. list("E" = x, "B" = y).
#' @param mfpc_cut_method Method to determine the level of explained variance
#'   \itemize{
#'     \item total_disp: total dispersion (trace(\eqn{\Sigma})).
#'     \item unidim: separate on each dimension.
#'   }
#' @param final_method Function used for estimation of final model to allow for
#'   potential heteroscedasticity ("w_bam", "bam", "gamm", "gaulss").
#' @param ... Additional arguments to be passed to (mainly) the underlying
#'   sparseFLMM function.
#' @return A list with four elements
#'   \itemize{
#'     \item the final multivariate FAMM
#'     \item the sparseFLMM output for each of the dimensions
#'     \item the MFPC output
#'     \item the data used to fit the model.}
#' @export
#' @import data.table
#' @import funData
mvfr <- function(data, fRI_B = FALSE, fRI_C = FALSE, nested = FALSE, bs = "ps",
                 bf_mean = 8, bf_covariates = 8, m_mean = c(2,3),
                 covariate = FALSE, num_covariates = NULL,
                 covariate_form = NULL, interaction = FALSE,
                 which_interaction = matrix(NA), bf_covs, m_covs,
                 var_level = 0.99, use_famm = FALSE, save_model_famm = FALSE,
                 one_dim = NULL, mfpc_cutoff = 0.95, number_mfpc = NULL,
                 mfpc_cut_method = c("total_disp", "unidim"),
                 final_method = c("w_bam", "bam", "gamm", "gaulss"), ...){

  # Match arguments that are chosen from list of options
  final_method <- match.arg(final_method)
  mfpc_cut_method <- match.arg(mfpc_cut_method)

  # model_list is a list that contains a sparse_flmm object for each dimension
  model_list <- apply_sparseFLMM(fRI_B = fRI_B,
                    fRI_C = fRI_C, nested = nested, bs = bs, bf_mean = bf_mean,
                    bf_covariates = bf_covariates, m_mean = m_mean,
                    covariate = covariate,
                    num_covariates = num_covariates,
                    covariate_form = covariate_form,
                    interaction = interaction,
                    which_interaction = which_interaction,
                    bf_covs = bf_covs, m_covs = m_covs, var_level = var_level,
                    use_famm = use_famm, save_model_famm = save_model_famm,
                    one_dim = one_dim, data = data, ...)

  # Return the sparseFLMM object if model is to be computed on only one
  # dimension
  if (! is.null(one_dim)) {
    return(model_list[[1]])
  }

  # mfpca_info is a list that contains three elements E, B, and C (or less)
  # In every element there are
  #   - multiFunData object
  #   - list containing all the uniExpansions information
  #   - number of possible components to extract
  mfpca_info <- prepare_mfpca(model_list = model_list, fRI_B = fRI_B)

  # So far only unweigted MFPCA
  cat("--------------------------------------\n")
  cat("Compute MFPCA\n")
  cat("--------------------------------------\n")
  MFPC <- lapply(mfpca_info, function(x){
    tryCatch(
      MFPCA::MFPCA(mFData = x$mFData,
            M = x$M_comp,
            uniExpansions = x$uniExpansions),
      error = function(e) return(NULL))
  })

  # Prune the MFPC to prespecified cutoff level of variance explained
  MFPC <- prune_mfpc(MFPC = MFPC, mfpc_cutoff = mfpc_cutoff,
                     model_list = model_list, mfpc_cut_method = mfpc_cut_method,
                     number_mfpc = number_mfpc)

  # Attach weighted MFPC to data
  data <- attach_wfpc(MFPC = MFPC, data = data)

  # Prepare data and model terms for the final model fitting
  data <- prepare_gam(data = data, bs = bs, bf_covariates = bf_covariates,
                      m_mean = m_mean, covariate = covariate,
                      num_covariates = num_covariates,
                      covariate_form = covariate_form,
                      interaction = interaction,
                      which_interaction = which_interaction)

  # Create the model formula
  formula <- create_formula(data = data, MFPC = MFPC, bs = bs,
                            bf_mean = bf_mean, m_mean = m_mean)

  # Final model
  cat("--------------------------------------\n")
  cat(paste("Compute final GAM with function", final_method, "\n"))
  cat("--------------------------------------\n")

  m2 <- final_model(formula = formula, data = data, final_method = final_method,
                    model_list = model_list)

  # List of results as output
  out <- list("model" = m2,
              "model_indep" = model_list,
              "mfpc" = MFPC,
              "data" = data)
  out

}
#------------------------------------------------------------------------------#


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
  data <- copy(data)

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
    # Rather unintuitively, Jona estimates E as B
    data[, subject_long := n_long]
    use_RI <- TRUE
    use_simple <- TRUE
  }

  # If only one dimension is to be computed, subset the data to the dimension
  if (! is.null(one_dim)) {
    data <- subset(data, dim == one_dim)
    data[, dim := droplevels(dim)]
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


#------------------------------------------------------------------------------#
# Prepare Information Necessary for MFPCA
#------------------------------------------------------------------------------#
prepare_mfpca <- function(model_list, fRI_B){

  # Arguments
  # model_list    : List containing sparseFLMM objects for each dimension
  # fRI_B         : Boolean for including functional random intercept for
  #                   individual (B in Cederbaum)

  # Extract the necessary information from the model list
  model_info <- extract_mfpca_info(model_list = model_list, fRI_B = fRI_B)

  # Create list containing the information on the maximum of FPC
  M_comp <- lapply(model_info, function (x) {
    sum(sapply(x, function(y) y$N))
  })

  # Inflate the fPCs when there are different numbers of fPCs with zero curves
  model_info <- inflate(model_info)

  # Determine which model terms are necessary
  # Look at the scores whether there are missing values
  missing_dim <-lapply(model_info, function(x){
    lapply(lapply(x, "[[", 3), function(y) any(is.na(y)))
  })
  # If the scores are missing for all the dimensions, remove the model term
  keep_terms <- sapply(missing_dim, function(x) !all(unlist(x)))
  model_info <- model_info[keep_terms]
  M_comp <- M_comp[keep_terms]

  # Reconstruct the different functional Random Intercepts
  # Scores %*% Eigcts and saved
  model_info <- lapply(model_info, function(x){
    lapply(x, function(y){
      y$RI <- y[[3]] %*% t(y[[2]])
      y
    })
  })

  # Construct multiFunData objects for each variance component
  multiFun_comp <- lapply(model_info, function(x){
    multiFunData(lapply(x, function(y){
      funData(argvals = list(y$grid), X = y$RI)
    }))
  })

  # Create list containing the information of uniExpansion
  uniExpansions <- lapply(model_info, function(x){
    lapply(x, function(y){
      list(type = "given",
           functions = funData(argvals = list(y$grid),
                               X = t(y[[2]])),
           scores = y[[3]],
           ortho = TRUE)
    })
  })

  out <- mapply(function(x, y, z){
    list(mFData = x,
         uniExpansions = y,
         M_comp = z)
  }, multiFun_comp, uniExpansions, M_comp, SIMPLIFY = FALSE)
  out

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Extract Information Necessary for Preparing the MFPCA
#------------------------------------------------------------------------------#
extract_mfpca_info <- function(model_list, fRI_B){

  # Arguments
  # model_list    : List containing sparseFLMM objects for each dimension
  # fRI_B         : Boolean for including functional random intercept for
  #                   individual (B in Cederbaum)

  # Extract the name necessary to find the results of the decomposition
  method_name <- lapply(model_list,
                        function(x) names(x)[grep("^fpc_hat_", names(x))])

  # Create a list containing the elements of interest per dimension
  out <- mapply(function(x,y){

    # Extract information on grid
    grid <- x$my_grid

    # Extract information on estimated eigenfunctions
    # Correct for unintuitively handling of independent curves of Jona
    E_f <- if(fRI_B) x[[y]]$phi_E_hat_grid else x[[y]]$phi_B_hat_grid
    B_f <- if(fRI_B) x[[y]]$phi_B_hat_grid else x[[y]]$phi_E_hat_grid
    C_f <- x[[y]]$phi_C_hat_grid

    # Extract information on estimated scores
    E_s <- if(fRI_B) x[[y]]$xi_E_hat else x[[y]]$xi_B_hat
    B_s <- if(fRI_B) x[[y]]$xi_B_hat else x[[y]]$xi_E_hat
    C_s <- x[[y]]$xi_C_hat

    # Extract information on number of fPCs
    N_E <- if(fRI_B) x[[y]]$N_E else x[[y]]$N_B
    N_B <- if(fRI_B) x[[y]]$N_B else x[[y]]$N_E
    N_C <- x[[y]]$N_C

    # Also extract number of PCs
    list(E = list(grid = grid, E_f = E_f, E_s = E_s, N = N_E),
         B = list(grid = grid, B_f = B_f, B_s = B_s, N = N_B),
         C = list(grid = grid, C_f = C_f, C_s = C_s, N = N_C))

  }, model_list, method_name, SIMPLIFY = FALSE)

  # Invert the list to combine the results for all the dimensions
  out <- apply(do.call(rbind, out), 2, as.list)
  out

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Inflate the List of fPCs with Zero cCrves
#------------------------------------------------------------------------------#
inflate <- function(model_info){

  # Arguments
  # model_info    : List containing extracted information from function
  #                   extract_mfpca_info

  lapply(model_info, function(x) {

    # For each variance component extract info of number of fpcs
    N <- sapply(x, function(y) {
      y$N
    })

    # If there is the same number on all dimensions then do nothing
    if (length(unique(N)) == 1) {
      return(x)
    } else {

      # If there are different numbers inflate the second and third components
      # Position 2: Eigenfunctions
      # Position 3: Eignescores
      m <- max(N)

      # Unequal and maximum is 0 needs to be handled differently
      # cbind not possible but create new matrix
      if(m == 1) {

        # Number of groups needed for correct form of score matrix
        i <- which.max(N)
        n_groups <- nrow(x[[i]][[3]])

        # Inflation step
        x <- lapply(x, function(y) {
          if(dim(y[[2]])[2] != m) {
            y[[2]] <- matrix(0, ncol = 1, nrow = length(y[[1]]))
            y[[3]] <- matrix(0, ncol = 1, nrow = n_groups)
            }
          y
        })

      } else {

        # Inflation step
        x <- lapply(x, function(y) {

          # Make sure to only change the matrix if necessary
          if(dim(y[[2]])[2] != m) {
            y[[2]] <- cbind(y[[2]], matrix(0, ncol = m - ncol(y[[2]]),
                                           nrow = nrow(y[[2]])))
            y[[3]] <- cbind(y[[3]], matrix(0, ncol = m - ncol(y[[3]]),
                                           nrow = nrow(y[[3]])))
            }
          y
        })
      }
      return(x)
    }
  })


}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Prune the MFPC object to include only a prespecified level of explained var
#------------------------------------------------------------------------------#
prune_mfpc <- function(MFPC, mfpc_cutoff, model_list, mfpc_cut_method,
                       number_mfpc){

  # Arguments
  # MFPC            : List containing MFPC objects for each variance component
  # mfpc_cutoff     : Pre-specified level of explained variance of results of
  #                     MFPCA
  # model_list      : List containing sparseFLMM objects for each dimension
  # mfpc_cut_method : Method to determine the level of explained variance
  #                     total_disp: total dispersion (trace(\Sigma))
  #                     unidim: separate on each dimension
  # number_mfpc     : List containing the number of mfPCs needed for each
  #                     variance component

  # No pruning
  if (mfpc_cutoff == 1) return(MFPC)

  # If number of mfPCs is fixed by user prune the MFPC object accordingly
  if (!is.null(number_mfpc)) {
    return(prune_components(MFPC = MFPC, number_mfpc = number_mfpc))
  }

  # Extract estimates of sigma squared integrated
  sigma_sq <- sapply(model_list, function (x) {
    sig <- grep("^cov_hat_", names(x))
    x[[sig]]$sigmasq_int
  })

  # Extract the eigenvalues for each variance component
  # Reweight with norm on single dimension if necessary
  values <- switch(mfpc_cut_method,
    "total_disp" = {
        lapply(MFPC, "[[", "values")
    },
    "unidim" = {

        # Reweight the multivariate Eigenvalues
        tmp <- lapply(MFPC, function (x){
          values <- x$values
          lapply(x$functions@.Data, function (y) {
            funData::norm(y) * values
          })
        })

        # Change the order to have a list of dimensions
        lapply(seq_along(model_list), function (x) {
          lapply(tmp, "[[", x)
        })

    })

  # Compute the number of fPCs depending on the level of explained variance
  # according to the method chosen
  number_mfpc <-  switch(mfpc_cut_method,
    "total_disp" = {
        compute_var(sigma_sq = sigma_sq, values = values,
                    mfpc_cutoff = mfpc_cutoff)
    },
    "unidim" = {

        # Compute the number of fPCs on each dimension
        tmp <- mapply(function (x, y){
          compute_var(sigma_sq = x, values = y, mfpc_cutoff = mfpc_cutoff)
        }, sigma_sq, values, SIMPLIFY = FALSE)

        # Use the maximum of each variance component
        # Even if there is a different variance decomposition on the dimensions,
        # it is assured that the prespecified level is maintained
        tmp <- lapply(seq_along(MFPC), function (x) {
          max(sapply(tmp, "[[", x))
        })
        names(tmp) <- names(MFPC)
        tmp

    })

  # Prune the MFPC object according to the number of mfPCs necessary
  prune_components(MFPC = MFPC, number_mfpc = number_mfpc)

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Attach Weighted Functional Principal Components to the Data
#------------------------------------------------------------------------------#
compute_var <- function(sigma_sq, values, mfpc_cutoff){

  # Arguments
  # sigma_sq        : Vector containing the estimated variances on each
  #                     dimension
  # values          : List containing the multivariate Eigenvalues for each
  #                     variance component
  # mfpc_cutoff     : Pre-specified level of explained variance of results of
  #                     MFPCA


  # Compute the full variance using the total dispersion
  tot_var <- sum(sigma_sq) + sum(unlist(values))

  # Compute the variance that needs to be explained after subtracting
  # the error variance
  var_needed <- tot_var * mfpc_cutoff - sum(sigma_sq)

  # Cumulative sum of largest eigenvalues
  cum_var <- cumsum(sort(unlist(values), decreasing = TRUE))

  # How many fPCs are needed in total
  needed <- sort(unlist(values), decreasing = TRUE)[
   1:min(which(cum_var > var_needed))]

  # Determine how many fPCs are needed on each component
  number_fpc <- mapply(function(x, y) {
    y <- length(needed[grepl(x, names(needed))])
   }, names(values), values, SIMPLIFY = FALSE)

  number_fpc

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Attach Weighted Functional Principal Components to the Data
#------------------------------------------------------------------------------#
prune_components <- function(MFPC, number_mfpc){

  # Arguments
  # MFPC        : List containing MFPC objects for each variance component
  # number_mfpc : List containing the number of mfPCs needed for each variance
  #                 component

  # Drop list element completely
  if (any(sapply(number_mfpc, function (x) x == 0))) {
    pos <- which(sapply(number_mfpc, function (x) x == 0))
    MFPC[[pos]] <- NULL
    number_mfpc[[pos]] <- NULL
  }

  # Prune the components
  mapply(function(x, y){
    x$values <- x$values[1:y]
    x$functions@.Data <- lapply(x$functions@.Data, function (z) {
      z@X <- matrix(z@X[1:y, ], nrow = y)
      z
    })
    x$scores <- matrix(x$scores[, 1:y], ncol = y)
    x$vectors <- if (y == 1) x$vectors else matrix(x$vectors[1:y, 1:y],
                                                   nrow = y)
    x$normFactors <- x$normFactors[1:y]
    x
  }, MFPC, number_mfpc, SIMPLIFY = FALSE)

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Attach Weighted Functional Principal Components to the Data
#------------------------------------------------------------------------------#
attach_wfpc <- function(MFPC, data){

  # Arguments
  # MFPC        : List containing MFPC objects for each coponent of the variance
  # data        : Data where wfPC are attached to


  # Eigenfunctions have to be evaluated on the actual observed observation times
  interpol <- linear_interpol(MFPC, data)

  # Lists of Eigenvalues and argvals per variance component
  values <- lapply(MFPC, function(x) x <- x$values)
  argvals <- lapply(interpol, function(x){
    lapply(x, function(y) y <- unlist(argvals(y)))
    })

  # List of weighted eigenfunctions
  w_interpol <- mapply(function(x, y, z){
    mapply(function(u, v) {
      wE <- t(u@X * sqrt(y))
      rownames(wE) <- v
      wE
      }, x, z, SIMPLIFY = FALSE)
  }, interpol, values, argvals, SIMPLIFY = FALSE)

  # Ensure correct order in data
  # Ordered according to the dimensions
  # Rest of ordering is only for convenience
  data <- data[order(data$dim, data$n_long, data$t), ]

  # Order the weigthed eigenfunctions according to the observations in the data
  # Combine all the dimensions so that the matrices have equal nrow than data
  data_list <- split(data, data$dim)
  w_interpol <- lapply(w_interpol, function(x){
    out <- mapply(function(y, z){
        y[match(z$t, rownames(y)), ]
      }, x, data_list, SIMPLIFY = FALSE)
    if(dim(x[[1]])[2] == 1) out <- lapply(out, matrix, ncol = 1)
    do.call(rbind, out)
    })

  # Rename the elements so they can be distinguished when attached to the data
  w_interpol <- lapply(seq_along(w_interpol), function(i){
    colnames(w_interpol[[i]]) <- paste0("w", names(w_interpol)[[i]], "_",
                                                   1:ncol(w_interpol[[i]]))
    w_interpol[[i]]
  })

  # Combine the original data with the weighted Eigenfunctions
  data <- cbind(data, do.call(cbind, w_interpol))
  data

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Evaluate Eigenfunctions on the Observed Data Times
#------------------------------------------------------------------------------#
linear_interpol <- function(MFPC, data){

  # Arguments
  # MFPC        : List containing MFPC objects for each coponent of the variance
  # data        : Data where wfPC are attached to


  # Which observations have to be interpolated
  # Keep track of the order to order the X values appropriately
  obs_t <- lapply(split(data$t, data$dim), unique)
  grid <- unlist(argvals(MFPC[[1]]$functions), recursive = FALSE)
  arg_vals <- mapply(function(x,y){
    unique(c(y, x))
  }, obs_t, grid, SIMPLIFY = FALSE)
  arg_order <- lapply(arg_vals, order)

  # Attach NAs to the observed values to later interpolate them
  # arg_l helps with mapply function as it has the same structure as obs_x
  obs_x <- lapply(MFPC, function(x) X(x$functions))
  arg_l <- lapply(obs_x, function(x) x <- arg_vals)
  obs_x <- mapply(function(x, y){
    mapply(function(u, v){
      u <- cbind(u, matrix(NA, nrow = nrow(u), ncol = length(v)-ncol(u)))
      u
    }, x, y, SIMPLIFY = FALSE)
  }, obs_x, arg_l, SIMPLIFY = FALSE)

  # Reorder the argvals and xvalues
  arg_vals <- mapply(function(x, y){
    x[y]
  }, arg_vals, arg_order, SIMPLIFY = FALSE)
  arg_l <- lapply(obs_x, function(x) x <- arg_order)
  obs_x <- mapply(function(x, y){
    mapply(function(u, v){
      u <- u[, v]
      u
    }, x, y, SIMPLIFY = FALSE)
  }, obs_x, arg_l, SIMPLIFY = FALSE)

  # Convert the data back to funData object to apply the interpolation
  arg_l <- lapply(obs_x, function(x) x <- arg_vals)
  out <- mapply(function(x, y){
    mapply(function(u, v){
      funData(argvals = list(u), X = matrix(v, ncol = length(u)))
    }, x, y, SIMPLIFY = FALSE)
  }, arg_l, obs_x, SIMPLIFY = FALSE)

  # Apply the interpolation
  out <- lapply(out, function(x) lapply(x, approxNA))
  out

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Adapt funData Function approxNA to Consider x Values
#------------------------------------------------------------------------------#
# Code stems from R package funData Version 1.1
# Only adaptation: argument 'x = unlist(object@argvals)'
approxNA <- function (object) {
  funData(object@argvals, t(zoo::na.approx(t(object@X),
                                           x = unlist(object@argvals),
                                           na.rm = FALSE)))
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Prepare Data Set for Modeling with GAM
#------------------------------------------------------------------------------#
prepare_gam <- function(data, bs, bf_covariates, m_mean, covariate,
                        num_covariates, covariate_form, interaction,
                        which_interaction){

  # Arguments
  # data              : Data that contains all the covariates
  # bs                : Spline basis function, only tested for "ps" (as in
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


  # Transform data in order to have factor variables
  facs <- c("n_long", "subject_long", "word_long")
  facs <- facs[facs %in% colnames(data)]
  setDT(data)[, (facs):= lapply(.SD, as.factor), .SDcols=facs]

  # Rest of function is only necessary to create the interaction variables for
  # the covariates
  if(!covariate) {
    # No covariates
    return(data)
  }

  # Interactions of covariates with dimension
  terms <- vector()
  for (i in 1:num_covariates) {
    if (covariate_form[i] == "by") {

      # Create the necessary interaction variables
      form_temp <- as.formula(paste0("~ 0 + dim:covariate.", i))
      tmp <- model.matrix(form_temp, data = data)
      colnames(tmp) <- sub("\\:covariate", "", colnames(tmp))
      data <- cbind(data, tmp)

      # Create the necessary gam term
      terms <- c(terms,
                 paste0("s(t, k = ", bf_covariates, ", bs = \"", bs, "\", m = ",
                        paste0("c(", paste(m_mean, collapse = ","), ")"),
                        ", by = ", colnames(tmp), ")"))

    } else {
      cat("Experimental usage of covariable of type smooth\n")
      terms <- c(terms,
                 paste0("ti(t, covariate.", i, ", by = dim, k = ",
                        bf_covariates, ", bs = \"", bs, "\", m = ",
                        paste0("c(", paste(m_mean, collapse = ","), ")"),
                        ", mc = c(0, 1))"))
    }
  }

  # Include Interaction between covariates
  if(interaction) {
    for (i in 1:num_covariates) {
      for (k in 1:num_covariates) {
        if (which_interaction[i, k] & (i < k)) {

          # Create the necessary interaction variables
          form_temp <- as.formula(paste0("~ 0 + dim:covariate.", i,
                                         ":covariate.", k))
          tmp <- model.matrix(form_temp, data = data)
          colnames(tmp) <- sub("\\:covariate", "\\.inter", colnames(tmp))
          colnames(tmp) <- sub("\\:covariate", "", colnames(tmp))
          data <- cbind(data, tmp)

          # Create the necessary gam term
          terms <- c(terms,
                     paste0("s(t, k = ", bf_covariates, ", bs = \"", bs,
                            "\", m = ",
                            paste0("c(", paste(m_mean, collapse = ","), ")"),
                            ", by = ", colnames(tmp), ")"))
        }
      }
    }
  }

  # Output is the transformed data and the covariate terms as attributes
  attr(data, "gam_terms") <- terms
  data


}
#------------------------------------------------------------------------------#


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
final_model <- function(formula, data, final_method, model_list){

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
      weights <- sapply(model_list, function(x) {
        sig <- grep("^cov_hat_", names(x))
        1 / x[[sig]]$sigmasq
      })
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
