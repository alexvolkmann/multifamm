################################################################################
################################################################################
##                                                                            ##
##                      Helperfunctions to Extract Components                 ##
##                                                                            ##
################################################################################
################################################################################


# Aim:
# Functions to extract model components and predict covariate or random effects.

# library(mgcv)
# library(funData)
# library(parallel)

#------------------------------------------------------------------------------#
# Extract Model Components to be Compared
#------------------------------------------------------------------------------#
#' Extract Model Components to be Compared
#'
#' This is an internal function that helps to compare different models. The
#' models resulting from a multiFAMM() call are typically very big. This
#' function extracts the main information from a model so that a smaller R
#' object can be saved.
#'
#' So far the grid is fixed to be on [0,1].
#'
#' @param model multiFAMM model object from which to extract the information.
#' @param dimnames Names of the dimensions of the model.
#'
#' @return A list with the following elements
#'   \itemize{
#'     \item \code{error_var}: A list containing the following elements
#'       \itemize{
#'       \item \code{model_weights}: Model weights used in the final multiFAMM.
#'       \item \code{modelsig2}: Estimate of sigma squared in the final model.
#'       \item \code{uni_vars}: Univariate estimates of sigma squared.}
#'     \item \code{eigenvals}: List containing the estimated eigenvalues.
#'     \item \code{fitted_curves}: multiFunData object containing the fitted
#'       curves.
#'     \item \code{eigenfcts}: multiFunData object containing the estimated
#'       eigenfunctions.
#'     \item \code{cov_preds}: multiFunData object containing the estimated
#'       covariate effects.
#'     \item \code{ran_preds}: List containing multiFunData objects of the
#'       predicted random effects.
#'     \item \code{scores}: List containing matrices of the estimated scores.
#'     \item \code{meanfun}: multiFunData object containing the estimated mean
#'       function.
#'     \item \code{var_info}: List containing all eigenvalues and univariate
#'       norms before the MFPC pruning step
#'       \itemize{
#'       \item \code{eigenvals}: Vector of all multivariate eigenvalues.
#'       \item \code{uni_norms}: List of univariate norms of all
#'         eigenfunctions.}}
extract_components <- function (model, dimnames) {

  # Fix a grid
  grid <- seq(0, 1, length.out = 100)


  # Error variance
  modelweights <- unique(model$model$weights)
  modelsig2 <- model$model$sig2
  uni_vars <- sapply(model$model_indep, function (x) {
    x[[grep("^cov_hat_", names(x))]]$sigmasq_int
  })


  # Eigenvalues
  eigenvals <- lapply(model$mfpc, function (x) x$values)


  # Fitted curves
  fitted_curves <- predict_fitted(model = model, grid = grid)
  fitted_curves <- funData::multiFunData(lapply(split(fitted_curves,
                                                      fitted_curves$dim),
                                                function (x) {
    funData::funData(argvals = grid, X = matrix(x$y_fit, ncol = length(grid),
                                                byrow = TRUE))
  }))


  # Eigenfunctions
  eigenfcts <- lapply(model$mfpc, function (x) x$functions)


  # Predictions of covariable effects
  cov_preds <- predict_covs(model = model, method = "mul", type = "terms",
                            unconditional = FALSE, grid = grid)
  # Convert the predictions into a list of multiFunData objects for easier
  # comparison
  # Every covariable is a separate multiFunData object
  cov_preds <- lapply(cov_preds, function (x) {
    use <- ! grepl(paste(c("subject_long", "word_long", "n_long"),
                         collapse = "|"), colnames(x))
    x <- x[, use]
    colnames(x) <- gsub(paste(paste0("dim", dimnames), collapse = "|"), "",
                        colnames(x))
    z <- lapply(unique(colnames(x)), function (y) {
      x[, which(colnames(x) == y)]
    })
    names(z) <- c("int", unique(colnames(x))[-1])
    z <- lapply(z, function (y) {
      funData::multiFunData(lapply(seq_along(dimnames), function (u) {
        funData::funData(argvals = grid, X = matrix(y[, u], nrow = 1))
      }))
    })
  })


  # Functional Random Effects
  ran_preds <- lapply(names(model$mfpc), predict_re, model = model,
                      dimnames = dimnames, grid = grid)
  names(ran_preds) <- names(model$mfpc)
  # Convert the predictions into a list of multiFunData objects for easier
  # comparison
  # Every variance component is a separate multiFunData object
  ran_preds <- lapply(ran_preds, function (x) {
    x <- split(x, f = x$dim)
    funData::multiFunData(lapply(x, function (y) {
      funData::funData(argvals = grid, X = matrix(y$pred, ncol = length(grid),
                                                  byrow = TRUE))
    }))
  })


  # Scores
  # Use linear model for estimation of scores
  # Append all random effects of one variance component into one vector
  ran_ef <- lapply(ran_preds, function (x) {
    do.call(c, lapply(x, function (y) {
      as.vector(t(y@X))
    }))
  })
  # Create list of unit matrices according to levels of variance component
  xid <- lapply(ran_preds, function (x) {
    diag(nrow(x@.Data[[1]]@X))
  })
  # Extract the Eigenfunction as matrices
  xef <- lapply(eigenfcts, function (x) {
    lapply(x, function (y) {
      t(y@X)
    })
  })
  # Construct the design matrix for the linear model
  # Kronecker product for every dimension and stacked
  X <- mapply(function (x, y) {
    do.call(rbind, lapply(x, function (z) {
      y %x% z
    }))
  }, xef, xid, SIMPLIFY = FALSE)
  # Compute score via linear model using QR decomposition
  scores <- mapply(function (x, y, z){
    matrix(qr.coef(qr(x), y), byrow = TRUE, ncol = ncol(z[[1]]))
  }, X, ran_ef, xef, SIMPLIFY = FALSE)


  # Mean function (covs set to 0.5)
  # Helpful for the plots for the FPCs
  meanfun <- predict_mean(model = model, multi = TRUE, dimnames = dimnames)


  # Variance information
  # Can be used to compute the total variation in the data / on one dimension
  var_info <- model$var_info

  # Output
  comps <- list("error_var" = list("modelweights" = modelweights,
                                   "modelsig2" = modelsig2,
                                   "uni_vars" = uni_vars),
                "eigenvals" = eigenvals,
                "fitted_curves" = fitted_curves,
                "eigenfcts" = eigenfcts,
                "cov_preds" = cov_preds,
                "ran_preds" = ran_preds,
                "scores" = scores,
                "meanfun" = meanfun,
                "var_info" = var_info)

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Predict Models for Plotting
#------------------------------------------------------------------------------#
predict_covs <- function (model, method = c("mul", "uni"),
                           type = c("terms", "iterms"), unconditional = FALSE,
                           grid = seq(0, 1, length.out = 100), ...) {

  # Arguments:
  # model         Model to be predicted
  # method        Model was estimated using the "mul"tivariate approach or the
  #                 "uni"variate approach
  # type          Type of prediction (as in predict.gam)
  # unconditional Include smoothing parameter uncertainty (as in predict.gam)
  # grid          Grid of evaluated predicted values


  method <- match.arg(method)
  type <- match.arg(type)

  switch(method,
         "mul" = {

           # Use first row so that there are reasonable values for all the
           # variables
           dat <- model$model$model[1, ]

           # Set all covariable and interaction values to 1
           num <- sapply(dat, is.numeric)
           dim <- grepl("^dim", names(dat))
           dat[, names(num[num & dim])] <- 1

           # Blow up data set
           dat <- dat[rep(1, times = length(grid)), ]
           dat$t <- grid

           # Predict data set
           pred <- predict(model$model, newdata = dat, type = type,
                           se.fit = TRUE, unconditional = unconditional)

           # Predict the intercept for the other dimensions
           pos <- grep("^dim$", colnames(pred$fit))
           colnames(pred$fit)[pos] <- colnames(pred$se.fit)[pos] <- paste0(
             "dim", dat$dim[pos])
           dims <- levels(model$model$model$dim)[
             !levels(model$model$model$dim) %in% levels(droplevels(dat$dim))]
           for (i in dims) {
             dat$dim <- i
             p <- predict(model$model, newdata = dat, type = type,
                          se.fit = TRUE, unconditional = unconditional)

             # Attach intercept and functional intercept to fit
             pred$fit[, paste0("s(t):dim", i)] <- p$fit[, paste0("s(t):dim", i)]
             intercept <- matrix(p$fit[, "dim"], ncol = 1,
                                 dimnames = list(rownames(pred$fit),
                                                 paste0("dim", i)))
             pred$fit <- cbind(pred$fit, intercept)

             # Attach intercept and functional intercept to se.fit
             pred$se.fit[, paste0("s(t):dim", i)] <- p$se.fit[
               ,paste0("s(t):dim", i)]
             intercept <- matrix(p$se.fit[, "dim"], ncol = 1,
                                 dimnames = list(rownames(pred$se.fit),
                                                 paste0("dim", i)))
             pred$se.fit <- cbind(pred$se.fit, intercept)
           }

           pred

         },
         "uni" = {

           # Use first row so that there are reasonable values for all the
           # variables
           dat <- model[[grep("^fpc_famm_hat", names(model))]]$
             famm_estim$model[1, ]

           # Set all covariable and interaction values to 1
           name <- c("^covariate", "^inter_")
           name <- grepl(paste(name, collapse = "|"), names(dat))
           dat[, name] <- 1

           # Blow up data set
           dat <- dat[rep(1, times = length(grid)), ]
           dat$yindex.vec <- grid

           # Predict data set
           mgcv::predict.bam(model[[grep("^fpc_famm_hat",
                                         names(model))]]$famm_estim,
                             newdata = dat, type = type, se.fit = TRUE,
                             unconditional = unconditional)

         })

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Predict Random Effects for Comparison
#------------------------------------------------------------------------------#
predict_re <- function (model, component = c("E", "B", "C"), dimnames,
                        grid = seq(0, 1, length.out = 100)) {

  # Arguments
  # model     : Model from which the functional Random Effects are to be
  #               predicted
  # component : Variance component to be plotted
  # dimnames  : Names of dimensions specified in the model
  # grid      : Grid used for prediction

  component <- match.arg(component)

  # Identify the different levels for which the Random Effects exist
  # Combine to a data.frame
  fac_var <- switch(component, "E" = "n_long", "B" = "subject_long",
                    "C" = "word_long")
  var_levels <- unique(model$model$model[, fac_var])
  dat <- expand.grid(grid, dimnames, var_levels)

  # Extract Eigenfunctions and Eigenvalues from model and compute the weighted
  # Eigenfunctions
  wE <- lapply(model$mfpc[[component]]$functions@.Data, function (x) {
    t(x@X * sqrt(model$mfpc[[component]]$values))
  })

  # Take advantage of the structure of the data.frame
  dat <- cbind(dat, do.call(rbind, wE))
  names(dat) <- c("t", "dim", fac_var, paste0("w", component, "_",
                                              1:ncol(wE[[1]])))

  # For prediction, all the model components are necessary
  newdat <- model$model$model[rep(1, times = nrow(dat)), ]
  newdat[, names(dat)] <- dat

  # Predict Random Effects and attach to data.frame
  pred <- mgcv::predict.bam(model$model, newdata = newdat,
                            type = "terms")
  pred <- pred[, grep(fac_var, colnames(pred))]
  dat <- cbind(dat, pred)

  dat

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Predict Fitted Curves for Comparison
#------------------------------------------------------------------------------#
predict_fitted <- function (model, grid = seq(0, 1, length.out = 100)) {

  # Arguments
  # model     : Model from which the fitted values are to be predicted (on grid)
  # grid      : Grid used for prediction


  # Construct a new data.frame based on the unique functions
  curve_index <- which(!duplicated(model$model$model[, c("n_long", "dim")]))
  newdat <- model$model$model[rep(curve_index, each = length(grid)), ]

  # Number of unique functions for repeating the measurement times
  n <- length(unique(newdat$n_long))
  newdat$t <- rep(grid, times = n)

  # Extract Eigenfunctions and Eigenvalues from model and compute the weighted
  # Eigenfunctions
  wE <- lapply(model$mfpc, function (u) {
    w <- u$values
    lapply(u$functions@.Data, function (v) {
      x <- t(v@X * sqrt(w))
      matrix(x[rep(1:nrow(x), times = n), ], nrow = nrow(x) * n)
    })
  })

  # Combine the weighted Eigenfunctions of the dimensions
  wE <- lapply(wE, do.call, what = rbind)

  # Replace the weighted Eigenfunctions in newdat
  for (i in seq_along(wE)) {
    nam <- paste0("w", names(wE)[i], "_", 1:ncol(wE[[i]]))
    newdat[nam] <- wE[[i]]
  }

  # Predict fitted values and attach to data.frame
  pred <- mgcv::predict.bam(model$model, newdata = newdat,
                            type = "link")
  newdat <- cbind(newdat, y_fit = pred)

  newdat

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Predict The Mean Function For the FPC Plots
#------------------------------------------------------------------------------#
#' Predict The Mean Function For the FPC Plots
#'
#' This is an internal function that helps to interpret the FPCs. Extract the
#' mean function for all covariates set to 0.5. This is useful if combined with
#' the estimated FPCs because one can then add and subtract suitable multiples
#' from this function.
#'
#' @param model multiFAMM model or list of univariate models for which to
#'   predict the mean.
#' @param multi Indicator if it is a multiFAMM model (TRUE) or a list of
#'   univariate models.
#' @param dimnames Vector of strings containing the names of the dimensions.
#'
#' @return A multiFunData object.
predict_mean <- function(model, multi = TRUE, dimnames = c("aco", "epg")) {

  # Handle multivariate models and univariate models differently due to their
  # structure

  if (multi == TRUE) {

    # Use first observation to get the structure of the data.frame
    newdat <- model$model$model[1,]

    # All the covariates have to be set to 0.5
    # BUT only if they are on the same dimension as is computed in the moment
    newdat <- newdat[rep(1, times = 100*length(dimnames)), ]
    newdat$dim <-factor(rep(dimnames, each = 100))
    newdat$t <- rep(seq(0, 1, length.out = 100), times = length(dimnames))
    change_ind <- sapply(newdat$dim, function (x) {
      grepl(paste0("dim", x), colnames(newdat))
    })
    cov_ind <- grepl("^dim.", names(newdat))
    for (i in 1:nrow(newdat)) {
      newdat[i, change_ind[, i]] <- 0.5
      newdat[i, !change_ind[, i] & cov_ind] <- 0
    }

    # Predict the gam terms
    out <- mgcv::predict.bam(model$model, newdata = newdat, type = "terms")
    out <- rowSums(out[, grepl("dim", colnames(out))])

  } else {
    out <- sapply(model, function (x) {

      # Use first observation to get the structure of the data.frame
      newdat <- x$fpc_famm_hat_tri_constr$famm_estim$model[1, ]
      newdat[, grep("^covariate|^inter", names(newdat))] <- 0.5
      newdat <- newdat[rep(1, times = 100), ]
      newdat$yindex.vec <- seq(0, 1, length.out = 100)

      # Predict the gam terms
      out <- mgcv::predict.bam(x$fpc_famm_hat_tri_constr$famm_estim,
                               newdata = newdat, type = "terms")
      rowSums(out[, !grepl("^s\\(id\\_", colnames(out))]) +
        x$fpc_famm_hat_tri_constr$intercept
    })
  }

  # Output as multiFunData
  fundata_list <- lapply(seq_along(dimnames), function (i){
    funData::funData(argvals = seq(0, 1, length.out = 100),
                     X = matrix(out[((i-1)*100+1):(i*100)], ncol = 100,
                                nrow = 1, byrow = TRUE))
  })
  funData::multiFunData(fundata_list)

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Extract Model Components to be Compared from Univariate Model
#------------------------------------------------------------------------------#
extract_components_uni <- function (model) {

  # Arguments
  # model       : Model from which to extract

  grid <- seq(0, 1, length.out = 100)

  # Extract the available variance components
  var_comp <- model[[grep("^fpc_hat_", names(model))]][
    paste0("N_", c("E", "B", "C"))]
  names(var_comp) <- c("E", "B", "C")
  var_comp <- var_comp[sapply(var_comp, function (x) x > 0)]


  # Error variance
  error_var <- model[[grep("^cov_hat_", names(model))]]$sigmasq_int


  # Eigenvalues
  eigenvals <- lapply(names(var_comp), function (x) {
    model[[grep("^fpc_hat_", names(model))]][[paste0("nu_", x, "_hat")]]
  })
  names(eigenvals) <- names(var_comp)


  # Fitted curves
  fitted_curves <- predict_fitted_uni(model = model, grid = grid,
                                      var_comp = var_comp)
  fitted_curves <- funData::funData(argvals = grid,
                                    X = matrix(fitted_curves$y_fit,
                                                      ncol = length(grid),
                                                      byrow = TRUE))


  # Eigenfunctions
  eigenfcts <- lapply(names(var_comp), function (x) {
    y <- model[[grep("^fpc_hat_", names(model))]][[
      paste0("phi_", x, "_hat_grid")]]
    y <- funData::funData(argvals = grid, X = t(y))
    y
  })
  names(eigenfcts) <- names(var_comp)


  # Predictions of covariable effects
  cov_preds <- model[[grep("^fpc_famm_", names(model))]]
  covs <- lapply(cov_preds[grep("^famm_cb_", names(cov_preds))],
                 function (x) {
    funData::funData(argvals = grid, X = t(x$value))
  })
  covs <- c(int = funData::funData(argvals = grid,
                                   X = t(rep(cov_preds[["intercept"]],
                                             times = length(grid)))),
            covs)
  ses <- lapply(cov_preds[grep("^famm_cb_", names(cov_preds))],
                function (x) {
                  funData::funData(argvals = grid, X = t(x$se))
                })
  ses <- c(int = funData::funData(argvals = grid,
                                  X = t(rep(sqrt(vcov(cov_preds$famm_estim)[
                                    "(Intercept)", "(Intercept)"]),
                                    times = length(grid)))),
           ses)
  cov_preds <- list(fit = covs, se.fit = ses)


  # Functional Random Effects
  ran_preds <- lapply(names(var_comp), function (x) {
    y <- model[[grep("^fpc_famm_hat", names(model))]][[
      paste0("famm_predict_", x)]]
    y <- funData::funData(argvals = grid, X = y)
    y
  })
  names(ran_preds) <- names(var_comp)


  # Scores
  scores <- lapply(names(var_comp), function (x) {
    model[[grep("^fpc_famm_hat", names(model))]][[
      paste0("xi_", x, "_hat_famm")]]
  })
  names(scores) <- names(var_comp)

  # Output
  comps <- list("error_var" = error_var,
                "eigenvals" = eigenvals,
                "fitted_curves" = fitted_curves,
                "eigenfcts" = eigenfcts,
                "cov_preds" = cov_preds,
                "ran_preds" = ran_preds,
                "scores" = scores)

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Predict Fitted Curves for Comparison from Univariate Model
#------------------------------------------------------------------------------#
predict_fitted_uni <- function (model, grid = seq(0, 1, length.out = 100),
                                var_comp) {

  # Arguments
  # model     : Model from which the fitted values are to be predicted (on grid)
  # grid      : Grid used for prediction
  # var_comp  : List containing the information of which variance components are
  #               included in the model and how many fPCs there are (e.g.
  #               list("E" = 6, "B" = 4))


  # Construct a new data.frame based on the unique functions
  curve_index <- which(!duplicated(model[[
    grep("^fpc_famm_hat_", names(model))]]$famm_estim$model[, "id_n.vec"]))
  newdat <- model[[
    grep("^fpc_famm_hat_", names(model))]]$famm_estim$model[
      rep(curve_index, each = length(grid)), ]

  # Number of unique functions for repeating the measurement times
  n <- length(unique(newdat$id_n.vec))
  newdat$yindex.vec <- rep(grid, times = n)

  # Extract Eigenfunctions and Eigenvalues from model and compute the weighted
  # Eigenfunctions
  wE <- lapply(names(var_comp), function (x) {
    y <- model[[grep("fpc_hat_", names(model))]]
    w <- y[[grep(paste0("nu_", x, "_hat"), names(y))]]
    X <- t(y[[grep(paste0("phi_", x, "_hat_grid"), names(y))]])
    wX <- t(X * sqrt(w))
    matrix(wX[rep(1:nrow(wX), times = n), ], nrow = nrow(wX) * n)
  })
  names(wE) <- names(var_comp)

  # Replace the weighted Eigenfunctions in newdat
  for (i in seq_along(wE)) {
    nam <- paste0("phi_", names(wE)[i], "_hat_grid.PC", 1:ncol(wE[[i]]))
    newdat[nam] <- wE[[i]]
  }

  # Predict fitted values and attach to data.frame
  pred <- mgcv::predict.bam(model[[grep("^fpc_famm_hat_",
                                        names(model))]]$famm_estim,
                            newdata = newdat, type = "link")
  newdat <- cbind(newdat, y_fit = as.vector(pred))

  newdat

}
#------------------------------------------------------------------------------#
