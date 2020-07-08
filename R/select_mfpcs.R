
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
  # Set negative eigenvalues to 0
  values <- lapply(MFPC, "[[", "values")
  if (any(unlist(values) < 0)) {
    warns <- lapply(values, function (x) {
      sum(x < 0)
    })
    for (i in seq_along(warns)) {
      if (warns[[i]] > 0) {
        warning(paste0(warns[[i]], " multivariate negative eigenvalues in ",
                       names(warns)[i], " set to 0."))
      }
    }
    values <- lapply(values, function (x) {
      ifelse(x < 0, 0, x)
    })
  }

  # Reweight with squared norm on single dimension if necessary
  if (mfpc_cut_method == "total_var") {
    norms_sq <- c(rep(1, times = length(unlist(values))))
    names(norms_sq) <- names(unlist(values))
  } else {
    # Use fundatas norm function on each dimension
    norms_sq <- lapply(MFPC, function (x){
      lapply(x$functions@.Data, funData::norm, squared = TRUE)
    })
    # Change the order to have a list of dimensions
    norms_sq <- lapply(seq_along(model_list), function (x) {
        unlist(lapply(norms_sq, "[[", x))
    })
    # Set norms corresponding to negative eigenvalues to 0
    if (any(is.nan(unlist(norms_sq)))) {
      norms_sq <- lapply(norms_sq, function (x) {
        x[which(is.nan(x))] <- 0
        x
        })
    }
  }

  # Compute the number of fPCs depending on the level of explained variance
  # according to the method chosen
  number_mfpc <-  switch(mfpc_cut_method,
                         "total_var" = {
                           compute_var(sigma_sq = sigma_sq, values = values,
                                       norms_sq = norms_sq,
                                       mfpc_cutoff = mfpc_cutoff)
                         },
                         "unidim" = {

                           # Compute the number of fPCs on each dimension
                           tmp <- mapply(function (x, y){
                             compute_var(sigma_sq = x, values = values,
                                         norms_sq = y,
                                         mfpc_cutoff = mfpc_cutoff)
                           }, sigma_sq, norms_sq, SIMPLIFY = FALSE)

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
#' Compute the Number of FPCs needed
#'
#' This is an internal function. The function takes all the information needed
#' to calculate how many FPCs are needed to reach the pre-specified cutoff
#' level.
#'
#' @param sigma_sq Vector containing the estimated variances on each dimension.
#' @param values List containing the multivariate Eigenvalues for each variance
#'   component.
#' @param norms_sq Vector containing the squared norms to be used as weights on
#'   the Eigenvalues.
#' @param mfpc_cutoff Pre-specified level of explained variance of results of
#'   MFPCA.
#'
compute_var <- function(sigma_sq, values, norms_sq, mfpc_cutoff){

  # Compute the full variance using the total variation from MFPCA
  tot_var <- sum(sigma_sq) + as.numeric(unlist(values) %*% unlist(norms_sq))

  # Compute the variance that needs to be explained after subtracting
  # the error variance
  var_needed <- tot_var * mfpc_cutoff - sum(sigma_sq)

  # Cumulative sum of largest eigenvalues, keeping the original order
  val_order <- order(unlist(values), decreasing = TRUE)
  cum_var <- cumsum(unlist(values)[val_order] * norms_sq[val_order])

  # How many fPCs are needed in total
  needed <- unlist(values)[val_order][1:min(which(cum_var > var_needed))]

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
