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
#' @param m_covs List of marginal orders for the penalty in fRI covariance
#'   estimation (as in sparseFLMM).
#' @param var_level Pre-specified level of explained variance on each
#'   dimension (as in sparseFLMM).
#' @param use_famm Re-estimate the mean in FAMM context (as in sparseFLMM) -
#'   overwritten by one_dim.
#' @param save_model_famm Give out the FAMM model object (as in sparseFLMM) -
#'   overwritten by one_dim.
#' @param one_dim Specify the name of the dimension if sparseFLMM is to be
#'   computed only on one dimension.
#' @param mfpc_weight TRUE if the estimated univariate error variance is to be
#'   to be used as weights in the scalar product of the MFPCA.
#' @param mfpc_cutoff Pre-specified level of explained variance of results of
#'   MFPCA (1 for previous versions of mvfr).
#' @param number_mfpc List containing the number of mfPCs needed for each
#'   variance component e.g. list("E" = x, "B" = y).
#' @param mfpc_cut_method Method to determine the level of explained variance
#'   \itemize{
#'     \item total_var: (weighted) sum of variation over the dimensions.
#'     \item unidim: separate on each dimension.
#'   }
#' @param final_method Function used for estimation of final model to allow for
#'   potential heteroscedasticity ("w_bam", "bam", "gamm", "gaulss").
#' @param weight_refit Get the weights for the weighted bam by first refitting
#'   the model under an independence assumption but with mfpc basis functions.
#'   Defaults to FALSE.
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
#'
#' @examples
#' \dontrun{
#' # subset of the phonetic data (very small subset, no meaningful results can
#' # be expected and no random effects other than the random smooth should be
#' # included in the model)
#'
#' data(phonetic_subset)
#'
#' m <- multiFAMM(data = phonetic_subset, covariate = TRUE, num_covariates = 2,
#'                covariate_form = c("by", "by"), interaction = TRUE,
#'                which_interaction = matrix(c(FALSE, TRUE, TRUE, FALSE),
#'                nrow = 2, ncol = 2), bf_covs = c(5), m_covs = list(c(2, 3)),
#'                mfpc_cut_method = "total_var", final_method = "w_bam")
#' }
multiFAMM <- function(data, fRI_B = FALSE, fRI_C = FALSE, nested = FALSE,
                      bs = "ps", bf_mean = 8, bf_covariates = 8,
                      m_mean = c(2,3), covariate = FALSE, num_covariates = NULL,
                      covariate_form = NULL, interaction = FALSE,
                      which_interaction = matrix(NA), bf_covs, m_covs,
                      var_level = 0.99, use_famm = FALSE,
                      save_model_famm = FALSE, one_dim = NULL,
                      mfpc_weight = FALSE, mfpc_cutoff = 0.95,
                      number_mfpc = NULL,
                      mfpc_cut_method = c("total_var", "unidim"),
                      final_method = c("w_bam", "bam", "gamm", "gaulss"),
                      weight_refit = FALSE, ...){

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
  mfpca_info <- prepare_mfpca(model_list = model_list, fRI_B = fRI_B,
                              mfpc_weight = mfpc_weight)

  # So far only unweigted MFPCA or weighted using the estimated error vars
  cat("--------------------------------------\n")
  cat("Compute MFPCA\n")
  cat("--------------------------------------\n")

  MFPC <- conduct_mfpca(mfpca_info, mfpc_weight)

  # Prune the MFPC to prespecified cutoff level of variance explained
  MFPC <- prune_mfpc(MFPC = MFPC, mfpc_cutoff = mfpc_cutoff,
                     model_list = model_list, mfpc_cut_method = mfpc_cut_method,
                     number_mfpc = number_mfpc, mfpca_info = mfpca_info)

  # Attach weighted MFPC to data
  data <- attach_wfpc(MFPC = MFPC, data = data)

  # Prepare data and model terms for the final model fitting
  data <- prepare_gam(data = data, bs = bs, bf_covariates = bf_covariates,
                      m_mean = m_mean, covariate = covariate,
                      num_covariates = num_covariates,
                      covariate_form = covariate_form,
                      interaction = interaction,
                      which_interaction = which_interaction,
                      nested = nested)

  # Create the model formula
  formula <- create_formula(data = data, MFPC = MFPC, bs = bs,
                            bf_mean = bf_mean, m_mean = m_mean)

  # Final model
  cat("--------------------------------------\n")
  cat(paste("Compute final GAM with function", final_method, "\n"))
  cat("--------------------------------------\n")

  m2 <- final_model(formula = formula, data = data, final_method = final_method,
                    model_list = model_list, weight_refit = weight_refit)

  # List of results as output
  out <- list("model" = m2,
              "model_indep" = model_list,
              "mfpc" = MFPC,
              "data" = data)
  out

}
#------------------------------------------------------------------------------#

