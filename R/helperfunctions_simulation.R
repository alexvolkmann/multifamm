################################################################################
################################################################################
##                                                                            ##
##                       Helperfunctions for the Simulation                   ##
##                                                                            ##
################################################################################
################################################################################


# Aim:
# Functions to extract model components and predict covariate or random effects.

library(mgcv)
library(funData)
library(parallel)

#------------------------------------------------------------------------------#
# Extract Model Components to be Compared
#------------------------------------------------------------------------------#
extract_components <- function (model, dimnames) {
  
  # Arguments
  # model       : Model from which to extract
  # dimnames    : Names of dimensions in model
  
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
  fitted_curves <- multiFunData(lapply(split(fitted_curves, fitted_curves$dim), 
                                       function (x) {
    funData(argvals = grid, X = matrix(x$y_fit, ncol = length(grid), 
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
      multiFunData(lapply(seq_along(dimnames), function (u) {
        funData(argvals = grid, X = matrix(y[, u], nrow = 1))
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
    multiFunData(lapply(x, function (y) {
      funData(argvals = grid, X = matrix(y$pred, ncol = length(grid), 
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
  
  # Output
  comps <- list("error_var" = list("modelweights" = modelweights,
                                   "modelsig2" = modelsig2,
                                   "uni_vars" = uni_vars),
                "eigenvals" = eigenvals,
                "fitted_curves" = fitted_curves,
                "eigenfcts" = eigenfcts,
                "cov_preds" = cov_preds,
                "ran_preds" = ran_preds,
                "scores" = scores)
  
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
           predict.bam(model[[grep("^fpc_famm_hat", names(model))]]$famm_estim,
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
  pred <- predict.bam(model$model, newdata = newdat, 
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
  pred <- predict.bam(model$model, newdata = newdat, 
                      type = "link")
  newdat <- cbind(newdat, y_fit = pred)
  
  newdat
  
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
  fitted_curves <- funData(argvals = grid, X = matrix(fitted_curves$y_fit,
                                                      ncol = length(grid),
                                                      byrow = TRUE))
  
  
  # Eigenfunctions
  eigenfcts <- lapply(names(var_comp), function (x) {
    y <- model[[grep("^fpc_hat_", names(model))]][[
      paste0("phi_", x, "_hat_grid")]]
    y <- funData(argvals = grid, X = t(y))
    y
  })
  names(eigenfcts) <- names(var_comp)
  
  
  # Predictions of covariable effects
  cov_preds <- model[[grep("^fpc_famm_", names(model))]]
  covs <- lapply(cov_preds[grep("^famm_cb_", names(cov_preds))],
                 function (x) {
    funData(argvals = grid, X = t(x$value))
  })
  covs <- c(int = funData(argvals = grid, X = t(rep(cov_preds[["intercept"]], 
                                                     times = length(grid)))),
            covs)
  ses <- lapply(cov_preds[grep("^famm_cb_", names(cov_preds))],
                function (x) {
                  funData(argvals = grid, X = t(x$se))
                })
  ses <- c(int = funData(argvals = grid, 
                         X = t(rep(sqrt(vcov(cov_preds$famm_estim)[
                           "(Intercept)", "(Intercept)"]), 
                                                    times = length(grid)))),
            ses)
  cov_preds <- list(fit = covs, se.fit = ses)
  
  
  # Functional Random Effects
  ran_preds <- lapply(names(var_comp), function (x) {
    y <- model[[grep("^fpc_famm_hat", names(model))]][[
      paste0("famm_predict_", x)]]
    y <- funData(argvals = grid, X = y)
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
  pred <- predict.bam(model[[grep("^fpc_famm_hat_", names(model))]]$famm_estim,
                      newdata = newdat, type = "link")
  newdat <- cbind(newdat, y_fit = as.vector(pred))
  
  newdat
  
}
#------------------------------------------------------------------------------#
