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
    meanfun <- predict_mean(model = model, multi = multi, dimnames = dimlabels)
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
#' Note that the final data set might be differently ordered depending on
#' whether the model or the model components are used as a starting point.
#'
#' @param int_include Include the estimate and uncertainty of the scalar
#' intercept in the values of the functional intercept. Defaults to TRUE.
#' @param m_fac Multiplication factor to represent the confidence interval of
#' the estimates. Defaults to 1.96.
#' @inheritParams fpc_plot_helper
#' @importFrom magrittr %>%
covariate_plot_helper <- function(model, mcomp, dimlabels, int_include = TRUE,
                                  multi = TRUE, m_fac = 1.96) {

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

    # Collapse the list of multiFunData objects to
    data_list <- lapply(mcomp$cov_preds, function (x) {
      out_outer <- lapply(seq_along(x), function (y) {
        out_inner <- multifamm:::funData2DataFrame(x[[y]])
        out_inner$cov <- y - 1
        out_inner
      })
      out_outer <- do.call(rbind, out_outer)
    })

  } else {

    # model is given

    cov_preds <- multifamm:::predict_covs(model = model, method = "mul",
                                          type = "terms", unconditional = FALSE)
    data_list <- lapply(cov_preds, function (x) {

      # Adjust the data by removing/adding columns
      x <- x[, ! grepl(paste(c("subject_long", "word_long", "n_long"),
                    collapse = "|"), colnames(x))]
      sc_ints <- which(colnames(x) %in% paste0("dim", dimlabels))
      if (int_include) {
        # Add scalar and functional intercept
        fu_ints <- which(colnames(x) %in% paste0("s(t):dim", dimlabels))
        x[, fu_ints] <- x[, fu_ints] + x[, sc_ints]
      } else {
        # Otherwise remove the scalar intercept from the data set
        print("Note: The scalar intercept is not included in the output.\n")
      }
      x <- x[, -sc_ints]

      # Transform the data
      fu_ints <- which(colnames(x) %in% paste0("s(t):dim", dimlabels))
      colnames(x)[fu_ints] <- paste0(colnames(x[, fu_ints]), ".0")
      x <- x[, order(colnames(x))]
      number_covs <- ncol(x) / length(dimlabels)
      x <- data.frame(t = seq(0, 1, length.out = 100),
                      y = as.vector(x),
                      dim = rep(sort(dimlabels), each = number_covs*100),
                      cov = rep(seq_len(number_covs)-1, each = 100))
      x
    })

  }

  # Construct the data
  dat <- data_list$fit %>%
    dplyr::mutate(dim = factor(dim, labels = dimlabels),
                  obs = NULL,
                  effect = factor(cov,
                                  labels = paste0("f[", unique(cov), "](t)")),
                  y_p = y + m_fac*data_list$se.fit$y,
                  y_m = y - m_fac*data_list$se.fit$y,
                  se = data_list$se.fit$y)

  dat

}



# Fct Transform Covariate Helper Data For 2D Case -------------------------

#' Transform Covaraiate Effect Plot Data for 2D Case
#'
#' This is an internal function. It gives as an output the data set necessary to
#' plot 2D trajectories of the effect plots.
#'
#' If a vector is supplied to covs, then the sum of the first elements up to the
#' second to last element is the given in the variable ... and the whole sum is
#' given in the variable ...
#'
#' @param data Data containing the covariate effect information as provided by
#'   covariate_plot_helper().
#' @param covs Covariate effects to be plotted. Can also be a vector if e.g. an
#'   interaction is to be plotted. Then, the interaction should be the last
#'   element of the vector. Defaults to the intercept.
#' @param indicator Two element vector containing the suffix indicating the 2D
#'   axis (first for x, then for y axis). Defaults to labels ending with ".x"
#'   and ".y".
#' @inheritParams covariate_plot_helper
#' @importFrom magrittr %>%
covariate_plot_helper_2d_transform <- function(data, covs = 0L,
                                               indicator = c("\\.x", "\\.y"),
                                               m_fac = 1.96) {

  # Select the covariate
  if (length(covs) == 1) {
    dat <- data %>%
      dplyr::filter(cov == covs) %>%
      dplyr::select(-cov, -effect)

    # Extract the dimension info
    dat$axis <- as.factor(ifelse(grepl(indicator[1], dat$dim), "x", "y"))
    dat$dim <- as.factor(gsub(paste(indicator, collapse = "|"), "", dat$dim))

    # Reshape the data to wide format
    dat <- dat %>%
      dplyr::rename(val = y, val_m = y_m, val_p = y_p) %>%
      tidyr::pivot_wider(names_from = axis,
                         values_from = c(val, val_m, val_p, se))
    return(dat)

  } else {

    # recursively apply this function for one covariate
    dat_list <- lapply(covs, covariate_plot_helper_2d_transform, data = data,
                       indicator = indicator)

    # define a function that adds multiple covariate effects
    sum_over_effects <- function (dat_list, covs) {
      dat <- dat_list[[length(covs)]]
      dat[, grep("val_[xy]", names(dat))] <- Reduce("+",
        lapply(dat_list[seq_along(covs)], function(x){
          x[, grep("val_[xy]", names(dat))]
        }))
      dat <- dat[, -grep("val_[mp]", names(dat))]
      dat
    }

    # create two data frames with the sum without the last covariate effect and
    # one with all (call it interaction)
    dat <- sum_over_effects(dat_list, covs = covs[-length(covs)])
    dat_int <- sum_over_effects(dat_list, covs = covs)
    dat <- dat %>%
      dplyr::right_join(dat_int, by = c("t", "dim"),
                             suffix = c("", "_int")) %>%
      select(-se_x, -se_y) %>%
      mutate(val_m_x_int = val_x_int - m_fac * se_x_int,
             val_m_y_int = val_y_int - m_fac * se_y_int,
             val_p_x_int = val_x_int + m_fac * se_x_int,
             val_p_y_int = val_y_int + m_fac * se_y_int) %>%
      as.data.frame()
    dat

  }

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




# Fct Create Data to Plot Covariate Effect Comparison Uni-Mul -------------

#' Create Data to Plot Covariate Effect Comparison Uni-Mul
#'
#' This is an internal function. It takes a data set of estimated covariate
#' effects of a multivariate model and compares it to univariately estimated
#' covariate effects.
#'
#' @param dat_m Helper data from multivariate model.
#' @param aco_pr Predictions from first univariate model.
#' @param epg_pr Predictions from second univariate model. Can be NULL.
#' @param mul_level Effect level of the multivariate model to be plotted.
#' @param uni_effect Number of the univariate effect to be plotted.
#' @importFrom magrittr %>%
covariate_comp_plot_helper <- function(dat_m, aco_pr, epg_pr, mul_level,
                                       uni_effect) {

  # Have null object if epg_pr is missing
  if (is.null(epg_pr)) {
    epg_pr <- list("fit" = NULL, "se.fit" = NULL)
  }

  # Remove unneccessary covariate effects from multivariate data set
  dat_m <- dat_m %>%
    dplyr::filter(effect %in% mul_level)

  # Create data set for the univariate covariate effects
  dat_u <- list("fit" = vector("list", length(uni_effect)),
                "se.fit" = vector("list", length(uni_effect)))
  for (i in seq_along(uni_effect)) {
    dat_u$fit[[i]] <- multifamm:::predictedUnivar2DataFrame(
      aco_pr = aco_pr$fit, epg_pr = epg_pr$fit, effect = uni_effect[i])
    dat_u$se.fit[[i]] <- multifamm:::predictedUnivar2DataFrame(
      aco_pr = aco_pr$se.fit, epg_pr = epg_pr$se.fit, effect = uni_effect[i])
  }
  dat_u <- lapply(dat_u, function(x) do.call(rbind, x))

  # Combine the two data sources
  dat <- cbind(dat_m,
               y_uni = dat_u$fit$y,
               y_uni_p = dat_u$fit$y + 1.96*dat_u$se.fit$y,
               y_uni_m = dat_u$fit$y - 1.96*dat_u$se.fit$y,
               se_uni = dat_u$se.fit$y)
  dat$se_dif <- factor((dat$se - dat$se_uni)<0, labels = c("worse", "better"))

  # Create the variables needed for the shaded areas
  grid <- seq(0, 1, length.out = 100)
  jumps <- grid[2] - grid[1]
  box_grid <- seq(min(grid) + jumps, max(grid - jumps),
                  length.out = length(grid) - 1)
  dat$max_val <- max(dat[, c("y_p", "y_uni_p")], na.rm = TRUE) + 0.05
  dat$min_val <- min(dat[, c("y_m", "y_uni_m")], na.rm = TRUE) - 0.05
  dat$x_start <- rep(c(min(grid), box_grid), times = 2)
  dat$x_end <- rep(c(box_grid, max(grid)), times = 2)

  dat

}



# Fct Coverage Area Plot Helper -------------------------------------------

#' Coverage Area Plot Helper Function
#'
#' This function is an internal function. The function takes one single effect
#' function from the simulated effects and gives a data.frame which can be used
#' to plot the coverage area of all simulation runs.
#'
#' @param ylim Two element numeric vector specifying the area for which the
#'   coverage is to be plotted. Defaults to c(-0.5, 0.5).
#' @param by Step size for the grid of y values to check for coverage.
#'   Defaults to 1/99.
#' @inheritParams create_coverage_array
coverage_area_helper <- function(sim_curves, effect_index, ylim = c(-0.5, 0.5),
                                 by = 1/99, m_fac = 1.96){

  niter <- length(sim_curves)

  # Create upper and lower bounds of each simulation run
  # Extract original curve
  # Sum up functional and scalar intercept if necessary
  if (effect_index %in% c(1, 2)) {
    sim_up <- lapply(sim_curves, function (it) {
      it$fit[[1]] + it$fit[[2]] + m_fac * (it$se.fit[[1]] + it$se.fit[[2]])
    })
    sim_do <- lapply(sim_curves, function (it) {
      it$fit[[1]] + it$fit[[2]] - m_fac * (it$se.fit[[1]] + it$se.fit[[2]])
    })
  } else {
    sim_up <- lapply(sim_curves, function (it) {
      it$fit[[effect_index]] + m_fac * (it$se.fit[[effect_index]])
    })
    sim_do <- lapply(sim_curves, function (it) {
      it$fit[[effect_index]] - m_fac * (it$se.fit[[effect_index]])
    })
  }

  # Create the grid points that are to be checked
  xseq <- seq(0, 1, length.out = 100)
  yseq <- seq(from = ylim[1], to = ylim[2], by = by)

  # Evaluate each y position for each x position and each dimension
  dat <- lapply(seq_along(sim_do[[1]]), function (dim) {
    x_list <- lapply(seq_len(100), function (xpos){
      y_list <- lapply(yseq, function (yval) {
        covered <- mapply(function (lim_up, lim_do){
          yval > lim_do[[dim]]@X[, xpos] & yval < lim_up[[dim]]@X[, xpos]
        }, lim_up = sim_up, lim_do = sim_do, SIMPLIFY = FALSE)
        data.frame(dim = paste0("dim", dim),
                   t = xseq[xpos],
                   ypos = yval,
                   cov = sum(unlist(covered))/niter)
      })
      do.call(rbind, y_list)
    })
    do.call(rbind, x_list)
  })
  do.call(rbind, dat)

}
