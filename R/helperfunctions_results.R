################################################################################
################################################################################
##                                                                            ##
##                   Helperfunctions for the Presented Results                ##
##                                                                            ##
################################################################################
################################################################################


# Aim:
# Functions to evaluate or prepare the resulting objects for evaluation. Also
# contains plot functions.


# library(sparseFLMM)
# library(gridExtra)
# library(stargazer)
# library(ggplot2)
# library(viridis)
# library(MFPCA)
# library(grid)


#------------------------------------------------------------------------------#
# Root Relative Mean Squared Error for Scalar Estimates
#------------------------------------------------------------------------------#
rrMSE <- function (theta_true, theta_estim) {

  # Arguments
  # theta_true  : True component
  # theta_estim : Estimated component

  sqrt(mean((theta_true - theta_estim)^2) / mean(theta_true^2))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Multivariate Root Relative Mean Squared Error for Functions
#------------------------------------------------------------------------------#
mrrMSE <- function (fun_true, fun_estim, flip = FALSE) {

  # Arguments
  # fun_true  : True function
  # fun_estim : Estimated function
  # flip      : Are estimated functions to be flipped?

  if (flip == TRUE) {
    fun_estim <- flipFuns(refObject = fun_true, newObject = fun_estim)
  }

  sqrt(mean(norm(fun_true - fun_estim)) / mean(norm(fun_true)))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Univariate Root Relative Mean Squared Error for Functions
#------------------------------------------------------------------------------#
urrMSE <- function (fun_true, fun_estim, flip = FALSE) {

  # Arguments
  # fun_true  : True function
  # fun_estim : Estimated function
  # flip      : Are estimated functions to be flipped?

  if (flip == TRUE) {
    fun_estim <- flipFuns(refObject = fun_true, newObject = fun_estim)
  }

  lapply(seq_along(fun_true@.Data), function (x) {
    sqrt(mean(norm(fun_true@.Data[[x]] - fun_estim@.Data[[x]])) /
           mean(norm(fun_true@.Data[[x]])))
  })

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Flip Estimated Eigenscores According to the Flipping of Eigenfunctions
#------------------------------------------------------------------------------#
flip_scores <- function (fun_true, fun_estim, score_estim) {

  # Arguments
  # fun_true    : True eigenfunctions
  # fun_estim   : Estimated eigenfunctions
  # score_estim : Estimated scores

  # Flip the estimated eigenfunctions
  flipped <- flipFuns(refObject = fun_true, newObject = fun_estim)

  # Information of which eigenfunctions were flipped on each dimension
  flips <- lapply(seq_along(flipped@.Data), function (x) {
    flipped@.Data[[x]]@X[, 1] != fun_estim@.Data[[x]]@X[, 1]
  })

  # If the same eigenfunctions were flipped on each dimension, flip the scores
  if(length(unique(flips)) == 1) {
    flips <- flips[[1]]
    score_estim[, flips] <- score_estim[, flips] * (-1)
  } else {
    warning("No flipping because different flips on the dimensions.")
  }

  # Output
  score_estim

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Flip Estimated Eigenscores According to the Flipping of Univariate
# Eigenfunctions
#------------------------------------------------------------------------------#
flip_scores_uni <- function (fun_true, fun_estim, score_estim) {

  # Arguments
  # fun_true    : True eigenfunctions
  # fun_estim   : Estimated eigenfunctions
  # score_estim : Estimated scores

  # Flip the estimated eigenfunctions
  flipped <- flipFuns(refObject = fun_true, newObject = fun_estim)

  # Information of which eigenfunctions were flipped on each dimension
  flips <- flipped@X[, 1] != fun_estim@X[, 1]

  score_estim[, flips] <- score_estim[, flips] * (-1)

  # Output
  score_estim

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Extract the Estimated Eigenfunctions for One Specific Variance Component
#------------------------------------------------------------------------------#
extract_Eigenfct_sim <- function (number, component, m_true_comp, eigenfcts,
                                  flip = TRUE) {

  # Arguments
  # number      : Number of Eigenfunction to be plotted
  # component   : Name of Variance Component
  # m_true_comp : True model components
  # eigenfcts   : Simulation results of eigenfcts

  # Flip the eigenfunctions if necessary
  eig <- lapply(eigenfcts$mul, function (x) {
    if (flip) {
      flipFuns(refObject = m_true_comp$eigenfcts[[component]],
               newObject = x[[component]])
    } else {
      x[[component]]
    }
  })

  # Extract the observations and create the X-matrix
  estims <- lapply(eig, function (x) {
    lapply(x@.Data, function (y){
      y@X[number, ]
    })
  })
  estims <- do.call(rbind, estims)
  estims <- lapply(1:ncol(estims), function (x) do.call(rbind, estims[, x]))

  # Concatenate the true eigenfunction
  estims <- lapply(seq_along(estims), function (x) {
    rbind(estims[[x]],
          m_true_comp$eigenfcts[[component]]@.Data[[x]]@X[number, ])
  })

  # Construct the multiFunData Object
  multiFunData(lapply(estims, function (x) {
    funData(argvals = getArgvals(eigenfcts$mul[[1]][[component]][[1]]),
            X = x)
  }))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Extract the Estimated Eigenfunctions for One Specific Variance Component for
# a different number of estimated fcts
#------------------------------------------------------------------------------#
extract_Eigenfct_sim_dim <- function (number, component, m_true_comp, eigenfcts,
                                  flip = TRUE) {

  # Arguments
  # number      : Number of Eigenfunction to be plotted
  # component   : Name of Variance Component
  # m_true_comp : True model components
  # eigenfcts   : Simulation results of eigenfcts

  # Flip the eigenfunctions if necessary
  eig <- lapply(eigenfcts$filled, function (x) {
    if (flip) {
      flipFuns(refObject = m_true_comp$eigenfcts[[component]],
               newObject = x[[component]])
    } else {
      x[[component]]
    }
  })

  # Extract the observations and create the X-matrix
  estims <- lapply(eig, function (x) {
    lapply(x@.Data, function (y){
      y@X[number, ]
    })
  })
  estims <- do.call(rbind, estims)
  estims <- lapply(1:ncol(estims), function (x) do.call(rbind, estims[, x]))

  # Concatenate the true eigenfunction
  estims <- lapply(seq_along(estims), function (x) {
    rbind(estims[[x]],
          m_true_comp$eigenfcts[[component]]@.Data[[x]]@X[number, ])
  })

  # Construct the multiFunData Object
  multiFunData(lapply(estims, function (x) {
    funData(argvals = getArgvals(eigenfcts$mul[[1]][[component]][[1]]),
            X = x)
  }))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Extract the Estimated Eigenfunctions for One Specific Variance Component of
# the Univariate Estimation
#------------------------------------------------------------------------------#
extract_Eigenfct_sim_uni <- function (number, component, eigenfcts) {

  # Arguments
  # number      : Number of Eigenfunction to be plotted
  # component   : Name of Variance Component
  # m_true_comp : True model components
  # eigenfcts   : Simulation results of eigenfcts that have been filled up

  # Flip the eigenfunctions if necessary
  eig <- lapply(eigenfcts$filled_uni, function (x) {
    out <- lapply(seq_along(x[[component]]), function (y) {
      flipFuns(refObject = x[[component]][[y]]$tru,
               newObject = x[[component]][[y]]$est)
    })
    names(out) <- names(x[[component]])
    out
  })

  # Extract the observations and create the X-matrix
  estims <- lapply(eig, function (x) {
    lapply(x, function (y){
      y@X[number, ]
    })
  })
  estims <- do.call(rbind, estims)
  estims <- lapply(1:ncol(estims), function (x) do.call(rbind, estims[, x]))

  # Concatenate the true eigenfunction
  estims <- lapply(seq_along(estims), function (x) {
    rbind(estims[[x]],
          eigenfcts$filled_uni[[1]][[component]][[x]]$tru@X[number, ])
  })

  # Construct the multiFunData Object
  multiFunData(lapply(estims, function (x) {
    funData(argvals =
              getArgvals(eigenfcts$filled_uni[[1]][[component]][[1]]$tru),
            X = x)
  }))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Extract the Estimated Covariate Effect for One Specific Model Term
#------------------------------------------------------------------------------#
extract_Covfct_sim <- function (term, m_true_comp, cov_preds) {

  # Arguments
  # term        : Name of covariate effect term
  # m_true_comp : True model components
  # eigenfcts   : Simulation results of covariate effects

  covs <- lapply(cov_preds$mul, function (x) x$fit[[term]])

  # Extract the observations and create the X-matrix
  estims <- lapply(covs, function (x) {
    lapply(x@.Data, function (y){
      y@X
    })
  })
  estims <- do.call(rbind, estims)
  estims <- lapply(1:ncol(estims), function (x) do.call(rbind, estims[, x]))

  # Concatenate the true eigenfunction
  estims <- lapply(seq_along(estims), function (x) {
    rbind(estims[[x]],
          m_true_comp$cov_preds$fit[[term]]@.Data[[x]]@X)
  })

  # Construct the multiFunData Object
  multiFunData(lapply(estims, function (x) {
    funData(argvals = getArgvals(cov_preds$mul[[1]]$fit[[term]][[1]]),
            X = x)
  }))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Extract the Estimated Covariate Effect for One Specific Model Term of the
# Univariate Estimation
#------------------------------------------------------------------------------#
extract_Covfct_sim_uni <- function (term, term_uni, m_true_comp, cov_preds) {

  # Arguments
  # term        : Name of covariate effect term
  # term_uni    : Name of covariate in univariate estimation ("int",
  #                 "famm_cb_mean", "famm_cb_covariate.1", "famm_cb_inter_1_2")
  # m_true_comp : True model components
  # eigenfcts   : Simulation results of covariate effects

  covs <- lapply(cov_preds$uni, function (x) {
    lapply(x, function (y) y$fit[[term_uni]])
  })

  # Extract the observations and create the X-matrix
  estims <- lapply(covs, function (x) {
    lapply(x, function (y){
      y@X
    })
  })
  estims <- do.call(rbind, estims)
  estims <- lapply(1:ncol(estims), function (x) do.call(rbind, estims[, x]))

  # Concatenate the true eigenfunction
  estims <- lapply(seq_along(estims), function (x) {
    rbind(estims[[x]],
          m_true_comp$cov_preds$fit[[term]]@.Data[[x]]@X)
  })

  # Construct the multiFunData Object
  multiFunData(lapply(estims, function (x) {
    funData(argvals = getArgvals(cov_preds$mul[[1]]$fit[[term]][[1]]),
            X = x)
  }))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Compute the True Values to Compare with Fitted Values
#------------------------------------------------------------------------------#
compute_fitted_sim <- function (fitted_cu, I = 10, J = 16, reps = 5) {

  # Arguments
  # fitted_cu   : Object saved from the simulation
  # I           : Number of subjects
  # J           : Number of words
  # reps        : Number of repetitions


  reps_B <- rep(reps*J, times = I)
  reps_C <- rep(reps, times = J)

  # For Random Intercept of Subject
  if ("B" %in% names(fitted_cu$tru[[1]]$re)) {
    re_B_true <- lapply(seq_along(fitted_cu$tru), function (x) {
      multiFunData(lapply(fitted_cu$tru[[x]]$re$B, function (y) {
        funData(argvals = getArgvals(y),
                X = y@X[rep(1:nrow(y@X), times = reps_B), ])
      }))
    })
  } else {
    # Zero object
    re_B_true <- lapply(seq_along(fitted_cu$tru), function (x) {
      argvals <- getArgvals(fitted_cu$tru[[x]]$mu[[1]])
      multiFunData(lapply(seq_along(fitted_cu$tru[[x]]$mu), function (y) {
        funData(argvals = argvals,
                X = matrix(0, nrow = nObs(fitted_cu$tru[[x]]$mu),
                           ncol = length(unlist(argvals))))
      }))
    })
  }

  # Has not been tested for C in names()
  if ("C" %in% names(fitted_cu$tru[[1]]$re)) {
    re_C_true <- lapply(seq_along(fitted_cu$tru), function (x) {
      multiFunData(lapply(fitted_cu$tru[[x]]$re$C, function (y) {
        funData(argvals = getArgvals(y),
                X = y@X[rep(rep(1:nrow(y@X), times = reps_C), times = I), ])
      }))
    })
  } else {
    # Zero object
    re_C_true <- lapply(seq_along(fitted_cu$tru), function (x) {
      argvals <- getArgvals(fitted_cu$tru[[x]]$mu[[1]])
      multiFunData(lapply(seq_along(fitted_cu$tru[[x]]$mu), function (y) {
        funData(argvals = argvals,
                X = matrix(0, nrow = nObs(fitted_cu$tru[[x]]$mu),
                           ncol = length(unlist(argvals))))
      }))
    })
  }

  mapply(function (x, y, z) {x$mu + x$re$E + y + z},
         fitted_cu$tru, re_B_true, re_C_true, SIMPLIFY = FALSE)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Transform a funData object to an Data.Frame
#------------------------------------------------------------------------------#
funData2DataFrame <- function(fundata, multifun = TRUE) {

  # Arguments
  # fundata   : Fundata object to be converted to a data frame
  # multifun  : Is the object a multiFunData object - DEPRECATED -

  # Automatic checking if the funData object belongs to class multiFunData
  multifun <- "multiFunData" %in% class(fundata)

  data_list <-if (multifun == TRUE) {
    lapply(fundata@.Data, function (x) {
      list(argvals = unlist(x@argvals), x = x@X)
    })
  } else {
    list(argvals = unlist(fundata@argvals), x = fundata@X)
  }

  dat <- data.frame(
    t = do.call(c, lapply(data_list, function (x) {
      rep(x$argvals, times = nrow(x$x))
    })),
    y = do.call(c, lapply(data_list, function (x) {
      as.vector(t(x$x))
    })),
    dim = rep(seq_along(data_list), times = sapply(data_list, function (x) {
      length(x$x)
    })),
    obs = do.call(c, lapply(data_list, function (x) {
      rep(seq_len(nrow(x$x)), each = length(x$argvals))
    })))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Evaluate the model components in the simulation
#------------------------------------------------------------------------------#
sim_eval_components <- function (folder, m_true_comp, label_cov,
                                 weighted = TRUE) {

  # Arguments
  # folder      : Folder with saved objects from the simulation
  # m_true_comp : True model components used for the simulation
  # label_cov   : Labels for the covariates
  # weighted    : Was the weighted bam used to compute the models

  if (weighted  == FALSE) {
    w <- "_bam"
  } else {
    w <- NULL
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  load(paste0(folder, "error_var", w,".Rdata"))

  true_sig <- m_true_comp$error_var$modelsig2 /
    m_true_comp$error_var$modelweights

  dat_err <- do.call(rbind, lapply(seq_along(error_var$mul),
                                   function (x, true) {
    sigma_hat <- error_var$mul[[x]]$modelsig2 / error_var$mul[[x]]$modelweights
    if (weighted == FALSE) {sigma_hat <- rep(sigma_hat, times = 2)}
    data.frame(it = rep(x, times = length(sigma_hat)),
               hat = sigma_hat,
               no = factor(1:length(sigma_hat),
                           labels = c("sigma[ACO]^2",
                                      "sigma[EPG]^2")),
               true = true,
               y = mapply(function (x, y) {
                 rrMSE(theta_true = x, theta_estim = y)
               }, true, sigma_hat),
               comp = factor("sigma^2"))
  }, true = true_sig))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  load(paste0(folder, "eigenvals", w,".Rdata"))

  dat_val <- do.call(rbind, lapply(names(m_true_comp$eigenvals), function (x) {
    do.call(rbind, lapply(seq_along(eigenvals$mul), function (y, true) {
      vals <- eigenvals$mul[[y]][[x]]
      data.frame(it = rep(y, times = length(vals)),
                 hat = vals,
                 no = factor(1:length(vals),
                             labels = paste0("upsilon[", seq_along(vals),
                                             "]^", x)),
                 true = true,
                 y = mapply(function (u, v) {
                   rrMSE(theta_true = u, theta_estim = v)
                 }, true, vals),
                 comp = "Eigenvalues")
    }, true = m_true_comp$eigenvals[[x]]))
  }))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  load(paste0(folder, "eigscores", w,".Rdata"))
  load(paste0(folder, "eigenfcts", w,".Rdata"))

  eigscores$mul_flip <- lapply(seq_along(eigscores$mul), function (y) {
    scores <- lapply(names(eigscores$mul[[1]]), function (x) {
      flip_scores(fun_true = m_true_comp$eigenfcts[[x]],
                  fun_estim = eigenfcts$mul[[y]][[x]],
                  score_estim = eigscores$mul[[y]][[x]])
    })
    names(scores) <- names(eigscores$mul[[1]])
    scores
  })

  dat_sco <- do.call(rbind, lapply(names(eigscores$mul_flip[[1]]),
                                   function (x) {
    do.call(rbind, lapply(seq_along(eigscores$mul_flip), function (y) {
      scores <- eigscores$mul_flip[[y]][[x]]
      true <- eigscores$tru[[y]][[x]]
      data.frame(it = rep(y, times = ncol(scores)),
                 no = factor(1:ncol(scores),
                             labels = paste0("rho[",seq_len(ncol(scores)),
                                             "]^", x)),
                 y = sapply(1:ncol(scores), function (z) {
                   rrMSE(theta_true = true[, z], theta_estim = scores[, z])
                 }),
                 comp = "Scores")
    }))
  }))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Eigenfunctions
  dat_fun <- do.call(rbind, lapply(names(eigenfcts$mul[[1]]), function (x) {
    do.call(rbind, lapply(seq_along(eigenfcts$mul), function (y) {
      fcts <- eigenfcts$mul[[y]][[x]]
      data.frame(it = rep(y, times = nObs(fcts)),
                 no = factor(1:nObs(fcts),
                             labels = paste0("psi[",seq_len(nObs(fcts)),
                                             "]^", x)),
                 y = sapply(1:nObs(fcts), function (z) {
                   mrrMSE(fun_true = extractObs(m_true_comp$eigenfcts[[x]],
                                                obs = z),
                          fun_estim = extractObs(fcts, obs = z), flip = TRUE)
                 }),
                 comp = "Eigenfunctions")
    }))
  }))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Covariance operator
  auto_cross <- list(c(1,1), c(1,2), c(2,2))
  true_covs <- lapply(names(m_true_comp$eigenvals), function (x) {
    lapply(auto_cross, function (y) {
      funData(argvals = list(seq(0,1,length.out = 100),
                             seq(0,1,length.out = 100)),
              X = array(covSurv(vals = m_true_comp$eigenvals[[x]],
                          fcts = m_true_comp$eigenfcts[[x]],
                          dim1 = y[1], dim2 = y[2]), dim = c(1,100,100)))
    })
  })
  names(true_covs) <- names(m_true_comp$eigenvals)

  dat_cop <- do.call(rbind, lapply(names(true_covs), function (x) {
    do.call(rbind, lapply(seq_along(auto_cross), function (y) {
      do.call(rbind, lapply(seq_along(eigenvals$mul), function (z) {
        estim_cov <- funData(argvals = list(seq(0,1,length.out = 100),
                                            seq(0,1,length.out = 100)),
                             X = array(covSurv(vals = eigenvals$mul[[z]][[x]],
                                               fcts = eigenfcts$mul[[z]][[x]],
                                               dim1 = auto_cross[[y]][1],
                                               dim2 = auto_cross[[y]][2]),
                                       dim = c(1,100,100)))
        data.frame(it = z,
                   no = factor(paste0("C[",
                                      paste(auto_cross[[y]], collapse = ""),
                                      "]^", x)),
                   y = mrrMSE(fun_true = true_covs[[x]][[y]],
                              fun_estim = estim_cov,
                              flip = FALSE),
                   comp = "Covariance")
      }))
    }))
  }))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Random effects
  load(paste0(folder, "ran_preds", w,".Rdata"))

  dat_ran <- do.call(rbind, lapply(names(ran_preds$mul[[1]]), function (x) {
    do.call(rbind, lapply(seq_along(ran_preds$mul), function (y) {
      randef <- ran_preds$mul[[y]][[x]]
      rantru <- ran_preds$tru[[y]][[x]]
      data.frame(it = y,
                 no = factor(1, labels = x),
                 y = mrrMSE(fun_true = rantru, fun_estim = randef,
                            flip = FALSE),
                 comp = "Fit")
    }))
  }))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Fitted values
  load(paste0(folder, "fitted_cu", w,".Rdata"))

  fit_true <- compute_fitted_sim(fitted_cu = fitted_cu, I = 10, J = 16,
                                 reps = 5)

  dat_fit <- do.call(rbind, lapply(seq_along(fitted_cu$mul), function (x) {
    data.frame(it = x,
               no = factor(1, label = "Y"),
               y = mrrMSE(fun_true = fit_true[[x]],
                          fun_estim = fitted_cu$mul[[x]],
                          flip = FALSE),
               comp = "Fit")
  }))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Covariate effects
  load(paste0(folder, "cov_preds", w,".Rdata"))
  names <- label_cov

  dat_cov <- do.call(rbind, lapply(seq_along(cov_preds$mul), function (x) {
    dat <- data.frame(it = x,
                      no = factor(seq_along(cov_preds$mul[[x]]$fit),
                                  labels = names),
                      y = sapply(seq_along(cov_preds$mul[[x]]$fit),
                                 function (y) {
                        mrrMSE(fun_true = m_true_comp$cov_preds$fit[[y]],
                               fun_estim = cov_preds$mul[[x]]$fit[[y]])
                      }),
                      comp = "Effectfunctions")
  }))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  comb <- c("it", "no", "y", "comp")
  dat <- rbind(dat_fit, dat_ran, dat_err[, comb], dat_cop, dat_val[,comb],
               dat_fun, dat_sco, dat_cov)

  dat

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Evaluate the model components per dimension in the simulation
#------------------------------------------------------------------------------#
sim_eval_dimensions <- function (folder, m_true_comp, label_cov) {

  # Arguments
  # folder      : Folder with saved objects from the simulation
  # m_true_comp : True model components used for the simulation
  # label_cov   : Labels for the covariates

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  load(paste0(folder, "error_var_comp.Rdata"))

  true_sig <- m_true_comp$error_var$modelsig2 /
    m_true_comp$error_var$modelweights

  dat_err_u <- do.call(rbind, lapply(seq_along(error_var$mul),
                                     function (x, true) {
     sigma_hat <- unlist(error_var$uni[[x]])
     data.frame(it = rep(x, times = length(sigma_hat)),
                hat = sigma_hat,
                no = factor(1, labels = c("sigma^2")),
                true = true,
                y = mapply(function (x, y) {
                  rrMSE(theta_true = x, theta_estim = y)
                }, true, sigma_hat),
                comp = factor("sigma^2"),
                dim = factor(c("ACO", "EPG")),
                method = factor("uni"))
  }, true = true_sig))

  dat_err_m <- do.call(rbind, lapply(seq_along(error_var$mul),
                                     function (x, true) {
     sigma_hat <- error_var$mul[[x]]$modelsig2 / error_var$mul[[x]]$modelweights
     data.frame(it = rep(x, times = length(sigma_hat)),
                hat = sigma_hat,
                no = factor(1, labels = c("sigma^2")),
                true = true,
                y = mapply(function (x, y) {
                  rrMSE(theta_true = x, theta_estim = y)
                }, true, sigma_hat),
                comp = factor("sigma^2"),
                dim = factor(c("ACO", "EPG")),
                method = factor("mul"))
   }, true = true_sig))

  dat_err <- rbind(dat_err_m, dat_err_u)
  # dat_err <- dat_err_u
  # dat_err$mul <- dat_err_m$y

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  load(paste0(folder, "ran_preds_comp.Rdata"))

  dat_ran_m <- do.call(rbind, lapply(names(ran_preds$mul[[1]]), function (x) {
    do.call(rbind, lapply(seq_along(ran_preds$mul), function (y) {
      randef <- ran_preds$mul[[y]][[x]]
      rantru <- ran_preds$tru[[y]][[x]]
      data.frame(it = rep(y, times = length(randef)),
                 dim = factor(1:length(randef), labels = c("ACO", "EPG")),
                 y = unlist(urrMSE(fun_true = rantru, fun_estim = randef,
                                   flip = FALSE)),
                 no = factor(x),
                 comp = factor("Fit"),
                 method = factor("mul"))
    }))
  }))

  dat_ran_u <- do.call(rbind, lapply(names(ran_preds$mul[[1]]), function (x) {
    do.call(rbind, lapply(seq_along(ran_preds$mul), function (y) {
      randef <- lapply(ran_preds$uni[[y]], function (z) z[[x]])
      rantru <- ran_preds$tru[[y]][[x]]
      do.call(rbind, lapply(seq_along(c("ACO", "EPG")), function (z) {
        data.frame(it = y,
                   dim = c("ACO", "EPG")[z],
                   no = factor(x),
                   y = mrrMSE(fun_true = rantru[[z]],
                              fun_estim = randef[[z]],
                              flip = FALSE),
                   comp = factor("Fit"),
                   method = factor("uni"))
      }))
    }))
  }))

  dat_ran <- rbind(dat_ran_m, dat_ran_u)
  # dat_ran <- dat_ran_u
  # dat_ran$mul <- dat_ran_m$y

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Fitted values
  load(paste0(folder, "fitted_cu_comp.Rdata"))

  fit_true <- compute_fitted_sim(fitted_cu = fitted_cu, I = 10, J = 16,
                                 reps = 5)

  dat_fit_m <- do.call(rbind, lapply(seq_along(fitted_cu$mul), function (x) {
      data.frame(it = x,
                 dim = c("ACO", "EPG"),
                 no = factor("Y"),
                 y = unlist(urrMSE(fun_true = fit_true[[x]],
                                   fun_estim = fitted_cu$mul[[x]],
                                   flip = FALSE)),
                 comp = "Fit",
                 method = factor("mul"))
  }))

  dat_fit_u <- do.call(rbind, lapply(seq_along(fitted_cu$uni), function (x) {
    do.call(rbind, lapply(seq_along(fitted_cu$uni[[1]]), function (y) {
      data.frame(it = x,
                 dim = c("ACO", "EPG")[y],
                 no = factor("Y"),
                 y = mrrMSE(fun_true = fit_true[[x]][[y]],
                            fun_estim = fitted_cu$uni[[x]][[y]],
                            flip = FALSE),
                 comp = factor("Fit"),
                 method = factor("uni"))
    }))
  }))

  dat_fit <- rbind(dat_fit_m, dat_fit_u)
  # dat_fit <- dat_fit_u
  # dat_fit$mul <- dat_fit_m$y

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Covariance operator
  load(paste0(folder, "eigenvals_comp.Rdata"))
  load(paste0(folder, "eigenfcts_comp.Rdata"))


  auto_cross <- list(c(1,1), c(1,2), c(2,2))
  true_covs <- lapply(names(m_true_comp$eigenvals), function (x) {
    lapply(auto_cross, function (y) {
      funData(argvals = list(seq(0,1,length.out = 100),
                             seq(0,1,length.out = 100)),
              X = array(covSurv(vals = m_true_comp$eigenvals[[x]],
                                fcts = m_true_comp$eigenfcts[[x]],
                                dim1 = y[1], dim2 = y[2]), dim = c(1,100,100)))
    })
  })
  names(true_covs) <- names(m_true_comp$eigenvals)

  dat_cop_m <- do.call(rbind, lapply(names(true_covs), function (x) {
    do.call(rbind, lapply(seq_along(auto_cross), function (y) {
      do.call(rbind, lapply(seq_along(eigenvals$mul), function (z) {
        estim_cov <- funData(argvals = list(seq(0,1,length.out = 100),
                                            seq(0,1,length.out = 100)),
                             X = array(covSurv(vals = eigenvals$mul[[z]][[x]],
                                               fcts = eigenfcts$mul[[z]][[x]],
                                               dim1 = auto_cross[[y]][1],
                                               dim2 = auto_cross[[y]][2]),
                                       dim = c(1,100,100)))
        data.frame(it = z,
                   dim = c("ACO", "cross","EPG")[y],
                   no = factor(paste0("C[",
                                      ifelse(y %% 2 == 0, "cross", "auto"),
                                      #paste(auto_cross[[y]], collapse = ""),
                                      "]^", x)),
                   y = mrrMSE(fun_true = true_covs[[x]][[y]],
                              fun_estim = estim_cov,
                              flip = FALSE),
                   comp = "Covariance",
                   method = "mul")
      }))
    }))
  }))

  dat_cop_u <- do.call(rbind, lapply(names(true_covs), function (x) {
    do.call(rbind, lapply(seq_along(auto_cross), function (y) {
      do.call(rbind, lapply(seq_along(eigenvals$uni), function (z) {

        vals <- lapply(eigenvals$uni[[z]], "[[", x)
        fcts <- lapply(eigenfcts$uni[[z]], "[[", x)

        estim_cov <- funData(argvals = list(seq(0,1,length.out = 100),
                                            seq(0,1,length.out = 100)),
                             X = array(covSurv(vals = vals,
                                               fcts = fcts,
                                               dim1 = auto_cross[[y]][1],
                                               dim2 = auto_cross[[y]][2],
                                               multi = FALSE),
                                       dim = c(1,100,100)))
        data.frame(it = z,
                   dim = c("ACO", "cross", "EPG")[y],
                   no = factor(paste0("C[",
                                      ifelse(y %% 2 == 0, "cross", "auto"),
                                      #paste(auto_cross[[y]], collapse = ""),
                                      "]^", x)),
                   y = mrrMSE(fun_true = true_covs[[x]][[y]],
                              fun_estim = estim_cov,
                              flip = FALSE),
                   comp = "Covariance",
                   method = "uni")
      }))
    }))
  }))

  dat_cop <- rbind(dat_cop_m, dat_cop_u)
  # dat_cop <- dat_cop_u
  # dat_cop$mul <- dat_cop_m$y

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Covariate effects
  load(paste0(folder, "cov_preds_comp.Rdata"))
  names <- label_cov

  dat_cov_m <- do.call(rbind, lapply(seq_along(cov_preds$mul), function (x) {
    do.call(rbind, lapply(seq_along(cov_preds$mul[[x]]$fit), function (y) {
      data.frame(it = x,
                 no = factor(names[y]),
                 dim = c("ACO", "EPG"),
                 y = unlist(urrMSE(fun_true = m_true_comp$cov_preds$fit[[y]],
                                   fun_estim = cov_preds$mul[[x]]$fit[[y]],
                                   flip = FALSE)),
                 comp = factor("Effectfunctions"),
                 method = factor("mul"))
    }))
  }))

  dat_cov_u <- do.call(rbind, lapply(seq_along(cov_preds$mul), function (x) {
    do.call(rbind, lapply(seq_along(cov_preds$uni[[x]]), function (y) {
      do.call(rbind, lapply(seq_along(cov_preds$uni[[x]][[y]]$fit),
                            function (z) {
        aha <- names(cov_preds$uni[[x]][[y]]$fit)
        aha <- sub("famm_cb", "s(t):", aha)
        aha <- gsub("_", ".", aha)
        names(cov_preds$uni[[x]][[y]]$fit) <- sub(".mean|.covariate", "", aha)

        u <- names(cov_preds$uni[[x]][[y]]$fit)[[z]]
        v <- which(names(m_true_comp$cov_preds$fit) == u)

        data.frame(it = x,
                   no = factor(names[v]),
                   dim = c("ACO", "EPG")[y],
                   y = mrrMSE(fun_true = m_true_comp$cov_preds$fit[[u]][[y]],
                              fun_estim = cov_preds$uni[[x]][[y]]$fit[[z]]),
                   comp = factor("Effectfunctions"),
                   method = factor("uni"))
      }))
    }))
  }))

  dat_cov <- rbind(dat_cov_m, dat_cov_u)
  # dat_cov <- dat_cov_u
  # dat_cov$mul <- dat_cov_m$y


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  comb <- c("it", "no", "y", "comp", "dim", "method")
  dat <- rbind(dat_fit, dat_ran, dat_err[, comb], dat_cop, dat_cov)

  dat

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Estimate Covariance Surface
#------------------------------------------------------------------------------#
covSurv <- function (vals, fcts, dim1, dim2, multi = TRUE) {

  # Arguments
  # vals   : Eigenvalues (Vector)
  # fcts   : Eigenfunctions (MultiFunData)
  # dim1   : Choose first dimension (Scalar)
  # dim2   : Choose second dimension (Scalar)
  # multi  : Is the input multivariate (if not then vals and fcts are lists)

  if (multi == TRUE) {

    # Extract the matrix of Eigenfunctions
    mat <- do.call(rbind, lapply(fcts, function (x) t(x@X)))

    # Construct a diagonal matric containing the Eigenvalues
    dia <- diag(vals, nrow = nObs(fcts))

  } else {

    # Create the matrix of Eigenfunctions
    mat <- rbind(cbind(t(fcts[[1]]@X),
                       matrix(0, ncol = nObs(fcts[[2]]), nrow = 100)),
                 cbind(matrix(0, ncol = nObs(fcts[[1]]), nrow = 100),
                       t(fcts[[2]]@X)))

    # Construct a diagonal matric containing the Eigenvalues
    dia <- diag(unlist(vals), nrow = length(unlist(vals)))

  }

  # Compute the Auto- and Cross-covariance
  cov <- mat %*% dia %*% t(mat)

  d1 <- if (dim1 == 1) 1:100 else 101:200
  d2 <- if (dim2 == 1) 1:100 else 101:200

  # Extract the wanted covariance surface
  cov[d1, d2]


}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Plot Covariate Effects Together
#------------------------------------------------------------------------------#
covariate_plot <- function (m_true_comp, effects = c(4, 5, 6), m_fac = 2,
                            limits = limits, dimlabels = c("ACO", "EPG")) {

  # Arguments
  # m_true_comp   : True model components
  # effects       : Number of effect function to be plotted
  # m_fac         : Multiplication factor
  # limits        : Y-limits for the plot

  # Extract covariate information
  covars <- lapply(effects, function (x) {
    funData2DataFrame(fundata = m_true_comp$cov_preds$fit[[x]])
  })

  # Extract information for standard deviation
  covars.se <- lapply(effects, function (x) {
    funData2DataFrame(fundata = m_true_comp$cov_preds$se.fit[[x]])
  })

  # Create the corresponding labels of the plot
  labels <- sapply(effects, function (x) {
    paste0("beta[", x-2, "](t)")
  })

  # Construct the data.frame
  dat <- data.frame(
    t = do.call(c, lapply(covars, "[[", "t")),
    y = do.call(c, lapply(covars, "[[", "y")),
    y_plus = do.call(c, lapply(seq_along(covars), function (x) {
      covars[[x]]$y + m_fac*covars.se[[x]]$y
    })),
    y_minu = do.call(c, lapply(seq_along(covars), function (x) {
      covars[[x]]$y - m_fac*covars.se[[x]]$y
    })),
    effect = factor(rep(effects, times = sapply(covars, nrow)),
                    labels = labels),
    dim = factor(do.call(c, lapply(covars, "[[", "dim")),
                 labels = dimlabels)
  )

  # Plot the data
  ggplot2::ggplot(data = dat, aes(x = t)) +
    geom_line(aes(y = y), size = 0.25) +
    geom_line(aes(y = y_plus), linetype = "dashed", size = 0.25) +
    geom_line(aes(y = y_minu), linetype = "dashed", size = 0.25) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.3) +
    facet_grid(dim ~ effect, labeller = label_parsed) +
    xlab("Normalized Time (t)") +
    ylab(expression("Effect Function (" ~ beta^(d)~"(t) )")) +
    theme_grey(base_size = 8) +
    scale_y_continuous(limits = limits)

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Predict The Mean Function For the FPC Plots
#------------------------------------------------------------------------------#
predict_mean <- function(model, multi = TRUE, dimlabels = c("aco", "epg")) {

  # Arguments
  # model : Multivariate model or list of univariate models
  # multi : Indicating which type of model

  # Handle multivariate models and univariate models differently due to their
  # structure

  if (multi == TRUE) {

    # Use first observation to get the structure of the data.frame
    newdat <- model$model$model[1,]

    # All the covariates have to be set to 0.5
    # BUT only if they are on the same dimension as is computed in the moment
    newdat <- newdat[rep(1, times = 100*length(dimlabels)), ]
    newdat$dim <-factor(rep(dimlabels, each = 100))
    newdat$t <- rep(seq(0, 1, length.out = 100), times = length(dimlabels))
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
      out <- predict.bam(x$fpc_famm_hat_tri_constr$famm_estim,
                         newdata = newdat, type = "terms")
      rowSums(out[, !grepl("^s\\(id\\_", colnames(out))]) +
        x$fpc_famm_hat_tri_constr$intercept
    })
  }

  # Output as multiFunData
  fundata_list <- lapply(seq_along(dimlabels), function (i){
    funData::funData(argvals = seq(0, 1, length.out = 100),
                     X = matrix(out[((i-1)*100+1):(i*100)], ncol = 100,
                                nrow = 1, byrow = TRUE))
  })
  funData::multiFunData(fundata_list)

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Predict The Covariate Effects
#------------------------------------------------------------------------------#
predict_model <- function (model, type = "terms", unconditional = FALSE,
                           grid = seq(0, 1, length.out = 100)) {

  # Arguments:
  # model         : Model with covariates
  # type          : Gam terms to be predicted
  # unconditional : Std conditional on lambda
  # gri           : Grid of evaluation points

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
  out <- predict.bam(model[[grep("^fpc_famm_hat", names(model))]]$famm_estim,
                     newdata = dat, type = type, se.fit = TRUE,
                     unconditional = unconditional)
  out$fit <- cbind(out$fit,
        intercept = rep(model[[grep("^fpc_famm_hat",
                                    names(model))]]$famm_estim$coefficients[1],
                         times = length(grid)))
  out$se.fit <- cbind(out$se.fit,
        intercept = rep(sqrt(model[[grep("^fpc_famm_hat",
                                         names(model))]]$famm_estim$Vp[1,1])))
  out

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Transform Output of Predict Model to Data Frame
#------------------------------------------------------------------------------#
predict2DataFrame <- function (aco_pr, epg_pr, effect,
                               grid = seq(0, 1, length.out = 100)) {

  # Arguments
  # aco_pr  : Output of model aco
  # epg_pr  : Output of model epg
  # effect  : Which effect to extract
  # grid    : Grid of evaluation points

  # Handle intercept differently
  if (effect == 1) {
    y <- c(aco_pr[, effect] + aco_pr[, ncol(aco_pr)],
           epg_pr[, effect] + epg_pr[, ncol(epg_pr)])
  } else {
    y <- c(aco_pr[, effect],
           epg_pr[, effect])
  }

  # Output of data.frame
  data.frame(
    t = rep(grid, times = 2),
    y = y,
    dim = rep(c(1, 2), each = length(grid))
  )

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Plot the Comparison of Univariate and Multivariate Covariate Effects
#------------------------------------------------------------------------------#
cova_comp_plot <- function (m_true_comp, aco_pr, epg_pr,
                            effects = c(4, 5, 6), effects_uni = c(5, 6, 7),
                            m_fac = 2, limits = limits) {

  # Arguments
  # m_true_comp : Model components of the multivariate model
  # aco_pr      : Predictions for aco model
  # epg_pr      : Predictions for epg model
  # effects     : Index number of the effects to be plotted from multivariate m
  # effects_uni : Index number of the effects to be plotted from the univariate
  # m_fac       : Multiplication factor
  # limits      : Y-limits of the Plot


  # List of data.frames of covariate effects for multivariate model
  covars <- lapply(effects, function (x) {
    funData2DataFrame(fundata = m_true_comp$cov_preds$fit[[x]])
  })

  # List of data.frames of covariate effects for univariate models
  covars_u <- lapply(effects_uni, function (x) {
    predict2DataFrame(aco_pr = aco_pr$fit, epg_pr = epg_pr$fit, effect = x)
  })

  # List of covariate credibility intervals for multivariate model
  covars.se <- lapply(effects, function (x) {
    funData2DataFrame(fundata = m_true_comp$cov_preds$se.fit[[x]])
  })

  # List of covariate credibility intervals forunivariate models
  covars_u.se <- lapply(effects_uni, function (x) {
    predict2DataFrame(aco_pr = aco_pr$se.fit, epg_pr = epg_pr$se.fit,
                      effect = x)
  })

  # Create the labels for the plot
  labels <- sapply(effects, function (x) {
    paste0("beta[", x-2, "](t)")
  })

  # Which approach yields "better" results
  difs <- mapply(function (x, y) {
    (x$y - y$y < 0)
  }, covars.se, covars_u.se, SIMPLIFY = FALSE)

  # Construct data.frame
  dat <- data.frame(
    t = c(do.call(c, lapply(covars, "[[", "t")),
          do.call(c, lapply(covars_u, "[[", "t"))),
    y = c(do.call(c, lapply(covars, "[[", "y")),
          do.call(c, lapply(covars_u, "[[", "y"))),
    y_plus = c(do.call(c, lapply(seq_along(covars), function (x) {
      covars[[x]]$y + m_fac*covars.se[[x]]$y
    })), do.call(c, lapply(seq_along(covars_u), function (x) {
      covars_u[[x]]$y + m_fac*covars_u.se[[x]]$y
    }))),
    y_minu = c(do.call(c, lapply(seq_along(covars), function (x) {
      covars[[x]]$y - m_fac*covars.se[[x]]$y
    })), do.call(c, lapply(seq_along(covars_u), function (x) {
      covars_u[[x]]$y - m_fac*covars_u.se[[x]]$y
    }))),
    effect = factor(c(rep(effects, times = sapply(covars, nrow)),
                      rep(effects, times = sapply(covars_u, nrow))),
                    labels = labels),
    dim = factor(c(do.call(c, lapply(covars, "[[", "dim")),
                   do.call(c, lapply(covars_u, "[[", "dim"))),
                 labels = c("ACO", "EPG")),
    method = factor(c(rep(1, times = do.call(sum, lapply(covars, nrow))),
                      rep(rep(c(2, 3), each = 100), times = length(covars_u))),
                    labels = c("multi", "aco", "epg")),
    dif = factor(rep(do.call(c, difs), times = 2),
                 label = c("worse", "better")),
    x_start = rep(rep(c(min(grid), box_grid), times = 2*length(effects)),
                  times = 2),
    x_end = rep(rep(c(box_grid, max(grid)), times = 2*length(effects)),
                times = 2),
    min_val = limits[1],
    max_val = limits[2]
  )

  # Plot the data
  ggplot2::ggplot(data = dat, aes(x = t)) +
    geom_line(aes(y = y, col = method), size = 0.25) +
    geom_line(aes(y = y_plus, col = method), linetype = "dashed", size = 0.25) +
    geom_line(aes(y = y_minu, col = method), linetype = "dashed", size = 0.25) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.3) +
    geom_rect(aes(ymin = min_val, ymax = max_val, xmin = x_start, xmax = x_end,
                  fill = dif), alpha = 0.15) +
    facet_grid(dim ~ effect, labeller = label_parsed) +
    xlab("Normalized Time (t)") +
    ylab(expression("Effect Function (" ~ beta^(d)~"(t) )")) +
    theme_grey(base_size = 8)  +
    theme(legend.position = "none") +
    scale_color_manual(values = c("black", "tomato2", "steelblue3")) +
    scale_fill_manual(values = c("grey80", NA)) +
    scale_y_continuous(expand = c(0, 0))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Convert Univariate Estimates of Eigenfunctions to Data.frame
#------------------------------------------------------------------------------#
fpc2DataFrame <- function(phi_aco, phi_epg, grid = seq(0,1, length.out = 100)) {

  # Arguments
  # phi_aco : Univariate eigenfunction for model aco
  # phi_epg : Univariate eigenfunction for model epg
  # grid    : Grid of evaluation points

  # Construct the data.frame
  dat <- data.frame(
    t = c(rep(grid, times = ncol(phi_aco)), rep(grid, times = ncol(phi_epg))),
    y = c(as.vector(phi_aco), as.vector(phi_epg)),
    dim = c(rep("ACO", length(phi_aco)), rep("EPG", length(phi_epg))),
    obs = c(rep(seq_len(ncol(phi_aco)), each = length(grid)),
            rep(seq_len(ncol(phi_epg)), each = length(grid)))
  )

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Compute the Surface of a Covariance Compoent and Convert it to Data.Frame
#------------------------------------------------------------------------------#
cov2DataFrame <- function(m_aco, m_epg, component) {

  # Arguments
  # m_aco     : Univariate model on dimension aco
  # m_epg     : Univariate model on dimension epg
  # component : Covariance component to extract

  dimnames <- c("ACO", "EPG")

  # Restructure model for easier coding if from univariate source
  model <- list(model_indep = list(m_aco, m_epg))
  names(model$model_indep) <- c("ACO", "EPG")

  # Extract the variance component wanted
  fpc <- lapply(model$model_indep, function (x) {
    y <- list()
    y[[1]] <- x[[grep("^fpc_hat", names(x))]][[paste0("phi_", component,
                                                      "_hat_grid")]]
    y[[2]] <- x[[grep("^fpc_hat", names(x))]][[paste0("nu_", component,
                                                      "_hat")]]
    y
  })

  # Compute the Auto-covariance on each dimension separately
  cov <- lapply(fpc, function (x) {
    x[[1]] %*% diag(x[[2]], nrow = length(x[[2]])) %*% t(x[[1]])
  })

  # Fill in covariance matrix with missing values for the Cross-covariance
  cov <- rbind(cbind(cov[[1]], matrix(NA, ncol = ncol(cov[[1]]),
                                      nrow = nrow(cov[[1]]))),
               cbind(matrix(NA, ncol = ncol(cov[[2]]), nrow = nrow(cov[[2]])),
                     cov[[2]]))

  # Extract the grid on which the fPCs are computed
  grid <- model$model_indep[[1]]$my_grid

  # Name the rows and columns for matrix transformation to data.frame
  rownames(cov) <- colnames(cov) <- rep(grid, times = 2)

  # Create data.frame for plotting
  dat <- data.frame(row_dim = rep(rep(dimnames, each = length(grid)),
                                  times = 2*length(grid)),
                    col_dim = rep(dimnames,
                                  each = 2 * length(grid) * length(grid)),
                    row = as.numeric(rownames(cov)[row(cov)]),
                    col = as.numeric(colnames(cov)[col(cov)]),
                    value = c(cov))

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Correlation Plots for Univariate and Multivariate Scores
#------------------------------------------------------------------------------#
# This code is an adjustment from the package GGally version 1.3.2
# of the function ggcorr
my_corr <- function (data, method = c("pairwise", "pearson"), cor_matrix = NULL,
                     nbreaks = NULL, digits = 2, name = "", low = "#3B9AB2",
                     mid = "#EEEEEE", high = "#F21A00", midpoint = 0,
                     palette = NULL,  geom = "tile", min_size = 2, max_size = 6,
                     label = FALSE, label_alpha = FALSE, label_color = "black",
                     label_round = 1, label_size = 3, limits = c(-1, 1),
                     drop = is.null(limits) || identical(limits, FALSE),
                     layout.exp = 0, legend.position = "right",
                     legend.size = 6, text_size = 2.7,
                     xmin=3.5, xmax=7.5, ymin=0.5, ymax=3.5,
                     xmin_d=5.5, xmax_d=7.5, ymin_d=3.5, ymax_d=5.5, ...) {

  # Arguments
  # See Pacakge Manual
  # Adjusted to plot boxes

  if (is.numeric(limits)) {
    if (length(limits) != 2) {
      stop("'limits' must be of length 2 if numeric")
    }
  }
  if (is.logical(limits)) {
    if (limits) {
      limits <- c(-1, 1)
    }
    else {
      limits <- NULL
    }
  }
  if (length(geom) > 1 || !geom %in% c("blank", "circle", "text",
                                       "tile")) {
    stop("incorrect geom value")
  }
  if (length(method) == 1) {
    method = c(method, "pearson")
  }
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      data = as.data.frame(data)
    }
    x = which(!sapply(data, is.numeric))
    if (length(x) > 0) {
      warning(paste("data in column(s)",
                    paste0(paste0("'", names(data)[x], "'"),
                           collapse = ", "),
                    "are not numeric and were ignored"))
      data = data[, -x]
    }
  }
  if (is.null(cor_matrix)) {
    cor_matrix = cor(data, use = method[1], method = method[2])
  }
  m = cor_matrix
  colnames(m) = rownames(m) = gsub(" ", "_", colnames(m))
  m = data.frame(m * lower.tri(m))
  rownames(m) = names(m)
  m$.ggally_ggcorr_row_names = rownames(m)
  m = reshape::melt(m, id.vars = ".ggally_ggcorr_row_names")
  names(m) = c("x", "y", "coefficient")
  m$coefficient[m$coefficient == 0] = NA
  if (!is.null(nbreaks)) {
    x = seq(-1, 1, length.out = nbreaks + 1)
    if (!nbreaks%%2) {
      x = sort(c(x, 0))
    }
    m$breaks = cut(m$coefficient, breaks = unique(x), include.lowest = TRUE,
                   dig.lab = digits)
  }
  if (is.null(midpoint)) {
    midpoint = median(m$coefficient, na.rm = TRUE)
    message(paste("Color gradient midpoint set at median correlation to",
                  round(midpoint, 2)))
  }
  m$label = round(m$coefficient, label_round)
  p = ggplot2::ggplot(na.omit(m), aes(x, y))
  if (geom == "tile") {
    if (is.null(nbreaks)) {
      p = p + geom_tile(aes(fill = coefficient), color = "grey30")
    }
    else {
      p = p + geom_tile(aes(fill = breaks), color = "white")
    }
    if (is.null(nbreaks) && !is.null(limits)) {
      p = p + scale_fill_gradient2(name, low = low, mid = mid,
                                   high = high, midpoint = midpoint,
                                   limits = limits)
    }
    else if (is.null(nbreaks)) {
      p = p + scale_fill_gradient2(name, low = low, mid = mid,
                                   high = high, midpoint = midpoint)
    }
    else if (is.null(palette)) {
      x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))
      p = p + scale_fill_manual(name, values = x, drop = drop)
    }
    else {
      p = p + scale_fill_brewer(name, palette = palette,
                                drop = drop)
    }
  }
  else if (geom == "circle") {
    p = p + geom_point(aes(size = abs(coefficient) * 1.25),
                       color = "grey50")
    if (is.null(nbreaks)) {
      p = p + geom_point(aes(size = abs(coefficient), color = coefficient))
    }
    else {
      p = p + geom_point(aes(size = abs(coefficient), color = breaks))
    }
    p = p + scale_size_continuous(range = c(min_size, max_size)) +
      guides(size = FALSE)
    r = list(size = (min_size + max_size)/2)
    if (is.null(nbreaks) && !is.null(limits)) {
      p = p + scale_color_gradient2(name, low = low, mid = mid,
                                    high = high, midpoint = midpoint,
                                    limits = limits)
    }
    else if (is.null(nbreaks)) {
      p = p + scale_color_gradient2(name, low = low, mid = mid,
                                    high = high, midpoint = midpoint)
    }
    else if (is.null(palette)) {
      x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))
      p = p + scale_color_manual(name, values = x, drop = drop) +
        guides(color = guide_legend(override.aes = r))
    }
    else {
      p = p + scale_color_brewer(name, palette = palette,
                                 drop = drop) +
        guides(color = guide_legend(override.aes = r))
    }
  }
  else if (geom == "text") {
    if (is.null(nbreaks)) {
      p = p + geom_text(aes(label = label, color = coefficient),
                        size = label_size)
    }
    else {
      p = p + geom_text(aes(label = label, color = breaks),
                        size = label_size)
    }
    if (is.null(nbreaks) && !is.null(limits)) {
      p = p + scale_color_gradient2(name, low = low, mid = mid,
                                    high = high, midpoint = midpoint,
                                    limits = limits)
    }
    else if (is.null(nbreaks)) {
      p = p + scale_color_gradient2(name, low = low, mid = mid,
                                    high = high, midpoint = midpoint)
    }
    else if (is.null(palette)) {
      x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))
      p = p + scale_color_manual(name, values = x, drop = drop)
    }
    else {
      p = p + scale_color_brewer(name, palette = palette,
                                 drop = drop)
    }
  }
  if (label) {
    if (isTRUE(label_alpha)) {
      p = p + geom_text(aes(x, y, label = label, alpha = abs(coefficient)),
                        color = label_color, size = label_size,
                        show.legend = FALSE)
    }
    else if (label_alpha > 0) {
      p = p + geom_text(aes(x, y, label = label), show.legend = FALSE,
                        alpha = label_alpha, color = label_color,
                        size = label_size)
    }
    else {
      p = p + geom_text(aes(x, y, label = label), color = label_color,
                        size = label_size)
    }
  }
  textData <- m[m$x == m$y & is.na(m$coefficient), ]
  xLimits <- levels(textData$y)
  textData$diagLabel <- textData$x
  if (!is.numeric(layout.exp) || layout.exp < 0) {
    stop("incorrect layout.exp value")
  }
  else if (layout.exp > 0) {
    layout.exp <- as.integer(layout.exp)
    textData <- rbind(textData[1:layout.exp, ], textData)
    spacer <- paste(".ggally_ggcorr_spacer_value", 1:layout.exp,
                    sep = "")
    textData$x[1:layout.exp] <- spacer
    textData$diagLabel[1:layout.exp] <- NA
    xLimits <- c(spacer, levels(m$y))
  }
  textData$diagLabel <- names(data)
  p = p + geom_text(data = textData, aes_string(label = "diagLabel"), ...,
                    na.rm = TRUE, parse = TRUE, size = text_size)
  p = p + geom_rect(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                    color="black", fill = "transparent", size = 1)
  p = p + geom_rect(xmin=xmin_d, xmax=xmax_d, ymin=ymin_d, ymax=ymax_d,
                    color="black", fill = "transparent", size = 1,
                    linetype = "dashed")
  p = p + scale_x_discrete(breaks = NULL, limits = xLimits)
  p = p + scale_y_discrete(breaks = NULL, limits = levels(m$y))
  p = p + labs(x = NULL, y = NULL) + coord_equal()
  p = p + theme(panel.background = element_blank(),
                legend.key = element_blank(),
                legend.position = legend.position,
                legend.title = element_text(size = legend.size),
                legend.title.align=0.5,
                legend.text = element_text(size = legend.size), ...)
  p = p + guides(fill = guide_colorbar(barwidth = 6,
                                       title = "Estimated Correlation",
                                       title.position = "top"))

  return(p)
}
#------------------------------------------------------------------------------#
