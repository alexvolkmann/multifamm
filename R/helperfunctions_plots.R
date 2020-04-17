################################################################################
################################################################################
##                                                                            ##
##                          Helperfunctions for Plots                         ##
##                                                                            ##
################################################################################
################################################################################


# Fct Extract Data for FPC Plot -------------------------------------------

#' Extract Data for FPC Plot
#'
#' This is an internal function. It gives as an output the data set necessary to
#' plot the eigenfunctions. It is meant to be general enough to also include
#' more than two dimensions.
#'
#' @param model The multifamm model from which to extract the estimates.
#' @param mcomp The extracted model components as given by the internal function
#' extract_components().
#' @param component Which FPC component to extract. Possible values ("B", "C",
#' "E")
#' @param dimlabels The dimensions of the model.
#' @param two_d Is a two dimensional representation (e.g. x and y coordinate) a
#' better representation for the data. Defaults to FALSE (one variable).
#' @param multi Is the model a multifamm or a sparseFLMM model. Defaults to TRUE
#' (multifamm model).
#' @param m_fac Multiplication factor to represent the difference from the
#' overall mean.
#' @importFrom magrittr %>%
fpc_plot_helper <- function(model, mcomp, component, dimlabels, two_d = FALSE,
                            multi = TRUE, m_fac = 2) {

  # Extract the mean function (all covariates are 0.5 and scalar intercept is
  # already added), eigenfunctions and eigenvalues
  if (missing(model)) {

    # mcomp is given
    meanfun <- multifamm:::funData2DataFrame(fundata = mcomp$meanfun)
    phi <- multifamm:::funData2DataFrame(fundata = mcomp$eigenfcts[[component]])
    lambda <- mcomp$eigenvals[[component]]

  } else {

    # model is given
    meanfun <- predict_mean(model = model, multi = multi, dimlabels = dimlabels)
    meanfun <- multifamm:::funData2DataFrame(fundata = meanfun)
    phi <- multifamm:::funData2DataFrame(fundata =
                                           model$mfpc[[component]]$functions)
    lambda <- model$mfpc[[component]]$values

  }

  # Split the eigenfunctions so that it is a list to apply over
  phi <- split(phi, phi$obs)

  # Construct the data
  # Use a vector to show the difference between mean and FPCs
  add <- do.call(c, mapply(function (phi_oned, lambda_oned) {
    phi_oned$y * sqrt(lambda_oned)
  }, phi, lambda, SIMPLIFY = FALSE))
  label <- paste0("{psi[", 1:length(lambda), "]^", component, "}(t)")
  dat <- data.frame(
    t = rep(meanfun$t, times = length(lambda)),
    val = rep(meanfun$y, times = length(lambda)),
    dim = factor(rep(meanfun$dim, times = length(lambda)),
                 labels = dimlabels),
    effect = factor(rep(1:length(lambda), each = length(meanfun$t)),
                    labels = label))
  dat$val_p <- dat$val + m_fac * add
  dat$val_m <- dat$val - m_fac * add

  # if a two dimensional way of representing the data is more natural, then
  # rearrange the data set
  if (two_d) {
    helpdat <- data.frame(
      do.call(rbind, strsplit(as.character(dat$dim), ".", fixed = TRUE)))
    names(helpdat) <- c("location", "axis")
    dat <- cbind(dat, helpdat)
    dat <- dat %>%
      tidyr::gather(key = type, value = value, val, val_p, val_m) %>%
      tidyr::unite(tmp, axis, type) %>%
      select(-dim) %>%
      tidyr::spread(tmp, value)
  }

  dat

}



# Fct Extract Data for Covariate Plot -------------------------------------

#' Extract Data for Covariate Plot
#'
#' This is an internal function. It gives as an output the data set necessary to
#' plot the covariate effects. It is meant to be general enough to also include
#' more than two dimensions.
#'
#' @param int_include Include the estimate and uncertainty of the scalar
#' intercept in the values of the functional intercept. Defaults to TRUE.
#' @param m_fac Multiplication factor to represent the confidence interval of
#' the estimates.
#' @inheritParams fpc_plot_helper
#' @importFrom magrittr %>%
covariate_plot_helper <- function(model, mcomp, dimlabels, int_include = TRUE,
                                  two_d = FALSE, multi = TRUE, m_fac = 1.96) {

  if (!multi) {
    stop("Not yet implemented for univariate models.")
  }

  # Extract the covariate effects
  if (missing(model)) {

    # mcomp is given

    # Handle intercept differently whether the scalar intercept is included with
    # the functional intercept or not
    if (int_include) {
      # Add scalar and functional intercept
      mcomp$cov_preds$fit[[2]] <- mcomp$cov_preds$fit[[2]] +
        mcomp$cov_preds$fit[[1]]
      mcomp$cov_preds$se.fit[[2]] <- mcomp$cov_preds$se.fit[[2]] +
        mcomp$cov_preds$se.fit[[1]]
    }
    mcomp$cov_preds$fit[[1]] <- NULL
    mcomp$cov_preds$se.fit[[1]] <- NULL

  } else {

    # model is given
    # --- NOT IMPLEMENTED ---
    stop("Not implemented for model. Maybe use extract_components() first.")

  }

  # Collapse the list of multiFunData objects to
  data_list <- lapply(mcomp$cov_preds, function (x) {
    out_outer <- lapply(seq_along(x), function (y) {
      out_inner <- multifamm:::funData2DataFrame(x[[y]])
      out_inner$cov <- y - 1
      out_inner
    })
    out_outer <- do.call(rbind, out_outer)
  })

  # Construct the data
  dat <- data_list$fit %>%
    dplyr::mutate(dim = factor(dim, labels = dimlabels),
                  obs = NULL,
                  effect = factor(cov,
                                  labels = paste0("f[", unique(cov), "](t)")),
                  y_p = y + m_fac*data_list$se.fit$y,
                  y_m = y - m_fac*data_list$se.fit$y)

  # if a two dimensional way of representing the data is more natural, then
  # rearrange the data set
  if (two_d) {

    # --- NOT IMPLEMENTED ---

  }

  dat

}



# Fct Extract Data for Covariance Surface Plot ----------------------------

#' Extract Data for Covariance Surface Plot
#'
#' This is an internal function. It gives as an output the data set necessary to
#' plot the covariance surface. Probably not apt for more than two dimensions.
#'
#' @inheritParams fpc_plot_helper
covariance_surf_plot_helper <- function(mcomp, component, dimlabels) {

  # Extract the matrix of Eigenfunctions
  mat <- do.call(rbind, lapply(mcomp$eigenfcts[[component]],
                               function (x) t(x@X)))

  # Construct a diagonal matric containing the Eigenvalues
  dia <- diag(mcomp$eigenvals[[component]],
              nrow = length(mcomp$eigenvals[[component]]))

  # Compute the Auto- and Cross-covariance
  cov_mat <- mat %*% dia %*% t(mat)

  # Extract the grid on which the fPCs are computed
  grid <- unlist(mcomp$eigenfcts[[component]]@.Data[[1]]@argvals)

  # Name the rows and columns for matrix transformation to data.frame
  rownames(cov_mat) <- colnames(cov_mat) <- rep(grid, times = 2)

  # Create data.frame for plotting
  dat <- data.frame(
    row_dim = rep(rep(dimlabels, each = length(grid)),
                  times = length(dimlabels)*length(grid)),
    col_dim = rep(dimlabels,
                  each = length(dimlabels) * length(grid) * length(grid)),
    row = as.numeric(rownames(cov_mat)[row(cov_mat)]),
    col = as.numeric(colnames(cov_mat)[col(cov_mat)]),
    value = c(cov_mat)
  )

  # Reverse the order of levels for more interpretable plot
  dat$col_dim <- factor(dat$col_dim, levels = rev(levels(dat$col_dim)))
  dat

}
