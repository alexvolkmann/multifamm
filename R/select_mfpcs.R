
#------------------------------------------------------------------------------#
# Prune the MFPC object to include only a prespecified level of explained var
#------------------------------------------------------------------------------#
#' Prune the MFPC object to include only a prespecified level of explained var
#'
#' This is an internal function contained in the multiFAMM function. This
#' function takes the MFPCA object and decides how many functional principal
#' components are to be included in the model.
#'
#' @param MFPC List containing MFPC objects for each variance component as given
#'   by the function conduct_mfpca()
#' @param model_list List containing sparseFLMM objects for each dimension as
#'   given by the output of apply_sparseFLMM()
#' @param mfpca_info Object containing all the neccessary information for the
#'   MFPCA. List as given by the output of prepare_mfpca().
#' @inheritParams multiFAMM
#'
prune_mfpc <- function(MFPC, mfpc_cutoff, model_list, mfpc_cut_method,
                       number_mfpc, mfpca_info){

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

  # Weight the sigma squared integrated if the total variation is used as a
  # cutoff criterion
  if (mfpc_cut_method == "total_var") {
    sigma_sq <- sigma_sq * mfpca_info[[1]]$weights
  }

  # Extract the eigenvalues for each variance component
  # Reweight with norm on single dimension if necessary
  values <- switch(mfpc_cut_method,
                   "total_var" = {
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
                         "total_var" = {
                           compute_var(sigma_sq = sigma_sq, values = values,
                                       mfpc_cutoff = mfpc_cutoff)
                         },
                         "unidim" = {

                           # Compute the number of fPCs on each dimension
                           tmp <- mapply(function (x, y){
                             compute_var(sigma_sq = x, values = y,
                                         mfpc_cutoff = mfpc_cutoff)
                           }, sigma_sq, values, SIMPLIFY = FALSE)

                           # Use the maximum of each variance component
                           # Even if there is a different variance decomposition
                           # on the dimensions,
                           # it is assured that the prespecified level is
                           # maintained
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


  # Compute the full variance using the total variation
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
