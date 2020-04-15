################################################################################
################################################################################
##                                                                            ##
##                          Helperfunctions for Plots                         ##
##                                                                            ##
################################################################################
################################################################################

#-------------------------------------------------------------------------------
# Function to Extract the FPCs and Output a Data.Frame
#-------------------------------------------------------------------------------
#' Extract the FPCs and Output a Data.Frame
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
#-------------------------------------------------------------------------------
