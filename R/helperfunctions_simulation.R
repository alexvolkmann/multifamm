
# Generate New Data Based on a multiFAMM ----------------------------------

#' Generate Data for Simulation
#'
#' This is an internal function. It creates a new data set based on already
#' given functions.
#'
#' @param I Number of levels for first grouping variable (individuals).
#' @param J Number of levels for second grouping variable (set to NA for
#'   random intercept design).
#' @param nested TRUE: second random effect is nested. FALSE: if it is
#'   a crossed effect. Defaults to FALSE.
#' @param num_dim Number of dimensions.
#' @param lamB Eigenvalues for first grouping variable.
#' @param lamC Eigenvalues for second grouping variable. Can be NULL.
#' @param lamE Eigenvalues for curves-specific deviations.
#' @param normal TRUE: FPC weights are drawn from a normal distribution.
#'   FALSE: they are drawn from a mixture of normals. Defaults to TRUE.
#' @param sigmasq Error variances for each dimension supplied as a list.
#' @param dim_indep TRUE: Use rmvnorm to ensure that error variances are
#'   independent across the dimensions. FALSE: Prolonge Jona's implementation
#'   using . Defaults to TRUE.
#' @param simple_sig TRUE: Random normal numbers weighted with sigma. FALSE:
#'   mvnorm. Defaults to TRUE.
#' @param mu Mean function. Can be NULL if the GAM model is supplied, then
#'   defaults to t + sin(t). If GAM model is supplied, then this argument is
#'   redundant.
#' @param N_B Number of FPCs for first grouping variable. Defaults to 1.
#' @param N_C Number of FPCs for second grouping variable. Can be 0 and defaults
#'   to 1.
#' @param N_E Number of FPCs for curves-specific deviations. Defaults to 1.
#' @param phiB_fun Eigenfunctions of the first grouping variable as a
#'   multiFunData object.
#' @param phiC_fun Eigenfunctions of the second grouping variable as a
#'   multiFunData object. Can be NULL.
#' @param phiE_fun Eigenfunctions of the curve-specific deviations as a
#'   multiFunData object.
#' @param minsize Minimal number of scalar observations per curve. Defaults to
#'   10.
#' @param maxsize Maximal number of observations per curve. Defaults to 20.
#' @param min_grid Minimal value of grid range. Defaults to 0.
#' @param max_grid Maximal value of grid range. Defaults to 1.
#' @param grid_eval Length of the final grid to evaluate the functions. Defaults
#'   to 100.
#' @param min_visit Minimal number of repetitions of each subject-word/session
#'   (first- and second grouping variable) combination (in case of random
#'   intercept design: minimal number of repetitions for first grouping
#'   variable). Defaults to 3.
#' @param max_visit Maximal number of repetitions of each subject-word/session
#'   (first- and second grouping variable) combination (in case of random
#'   intercept design: maximal number of repetitions for first grouping
#'   variable). Defaults to 3.
#' @param use_RI TRUE: data with a random intercept structure are generated.
#'   FALSE: data with crossed/nested random intercepts structure are generated
#'   (Additional layer "words"/"sessions" included). Defaults to FALSE.
#' @param covariate If covariates shall be generated. Defaults to TRUE.
#' @param num_cov Number of covariates if covariate = TRUE. Defaults to 4.
#' @param interaction TRUE if there are interactions between covariates (as in
#'   the phonetic sparseFLMM). Defaults to FALSE.
#' @param which_inter Symmetric matrix specifying the interaction terms (as in
#'   sparseFLMM). Defaults to missing matrix.
#' @param model GAM model from which to extract the covariate functions. Can be
#'   NULL if mu is specified.
#' @param center_scores TRUE: FPC weights are centered. Defaults to FALSE.
#' @param decor_scores TRUE: FPC weights are decorrelated Defaults to FALSE.
gendata <- function(I = 10, J = 10, nested = FALSE, num_dim = 2,
                    lamB, lamC, lamE, normal = TRUE,
                    sigmasq = list(0.05, 4), dim_indep = TRUE,
                    simple_sig = TRUE,
                    mu = function (t) {t + sin(t)},
                    N_B = 1, N_C = 1, N_E = 1,
                    phiB_fun, phiC_fun, phiE_fun,
                    minsize = 10, maxsize = 20,
                    min_grid = 0, max_grid = 1, grid_eval = 100,
                    min_visit = 3, max_visit = 3,
                    use_RI = FALSE, covariate = TRUE, num_cov = 4,
                    interaction = FALSE, which_inter = matrix(NA),
                    model = NULL,
                    center_scores = FALSE, decor_scores = FALSE){

  if(any(!require(data.table) | !require(funData))){
    stop("Packages data.table and funData have to be installed.")
  }

  # Check whether there are enough levels of the words if there are covariates
  # in the crossed design
  # Not really necessary but in this way we make sure that there are no
  # combination of factors that have no observations
  if (covariate) {
    if (!nested) {
      if (!use_RI & J < 2^num_cov) {
        warning("J is fixed to make sure there are all covariate combinations")
        J <- 2^num_cov
      }
    }
  } else {
    num_cov <- 0
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Determine how many curves there are
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Create a vector containing the number of repetitions Hvec
  # Argument prob so that the different sampling numbers are not uniformly
  # distributed
  if (!use_RI) {
    if (min_visit != max_visit) {
      Hvec <- sort(sample(x = (min_visit:max_visit), size = I*J, replace = TRUE,
                          prob = sample(x = 2:5,
                                        size = (max_visit - min_visit + 1),
                                        replace = TRUE)))
    } else {
      Hvec <- rep(min_visit, length = I*J)
    }
  } else {
    if (min_visit != max_visit) {
      Hvec <- sort(sample(x = (min_visit:max_visit), size = I, replace = TRUE,
                          prob = sample(x = 2:5,
                                        size = (max_visit - min_visit + 1),
                                        replace = TRUE)))
    } else {
      Hvec <- rep(min_visit, length = I)
    }
  }

  # Total number of curves
  n <- sum(Hvec)


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Make all combinations of subject/word
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Depending on the Argument use_RI
  if (!use_RI) {
    help <- expand.grid(1:J, 1:I)
    subject <- rep(help$Var2, Hvec)
    word <- rep(help$Var1, Hvec)
    combi_help <- list()
    for (i in 1:(I*J)) {
      combi_help[[i]] <- 1:Hvec[i]
    }
    combi <- unlist(combi_help)

  } else {
    subject <- rep(1:I, times = Hvec)
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Draw number of observed points for each curve
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Different number of observations on each dimension
  if (minsize != maxsize) {
    number <- sample(x = rep(minsize:maxsize), size = num_dim*n, replace = TRUE)
  } else {
    # Number of time points per curve
    number <- rep(minsize, length = num_dim*n)
  }

  # Contains the total amount of observation points of each curve
  number_long <- rep(number, number)

  # Long vector with 1:n each as many times as time points per curve
  n_long <- rep(rep(1:n, num_dim), number)

  # Long vector with subject[1]:subject[n] each as many times as points per
  # curve
  subject_long <- rep(rep(subject, num_dim), number)

  # Long vector with information of dimension
  n_on_dim <- sapply(split(number, rep(1:num_dim, each = n)), sum)
  dim <- factor(rep(paste0("dim", 1:num_dim), n_on_dim))


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Create basis of final data set
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Create curve_info object
  # Also create data set containing true values
  if (!use_RI) {
    word_long <- rep(rep(word, num_dim), number)
    combi_long <- rep(rep(combi, num_dim), number)
    curve_info <- data.table::data.table(dim, n_long, subject_long, word_long,
                                         combi_long, number_long)

    # True values evaluated on regular time points
    curve_true <- data.table::data.table(
      dim = factor(rep(paste0("dim", 1:num_dim), each = n*grid_eval)),
      n_long = rep(rep(1:n, each = grid_eval), times = num_dim),
      subject_long = rep(rep(subject, each = grid_eval), times = num_dim),
      word_long = rep(rep(word, each = grid_eval), times = num_dim),
      t = rep(seq(min_grid, max_grid, length.out = grid_eval),
              times = num_dim*n))

  } else {
    curve_info <- data.table::data.table(dim, n_long, subject_long, number_long)

    # True values evaluated on regular time points
    curve_true <- data.table::data.table(
      dim = factor(rep(paste0("dim", 1:num_dim), each = n*grid_eval)),
      n_long = rep(rep(1:n, each = grid_eval), times = num_dim),
      subject_long = rep(rep(subject, each = grid_eval), times = num_dim),
      t = rep(seq(min_grid, max_grid, length.out = grid_eval),
              times = num_dim*n))

  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Attach covariates and interactions
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Covariates in the crossed design depend upon the words
  if (covariate & !use_RI) {
    if (!nested) {
      covs <- list()
      for (i in 1:num_cov) {
        covs[[paste0("covariate.", i)]] <- c(0, 1)
      }
      covs <- do.call(expand.grid, covs)
      if (J == 2^num_cov) {
        curve_info <- cbind(curve_info,
                            covs[match(curve_info$word_long, rownames(covs)), ])
        curve_true <- cbind(curve_true,
                            covs[match(curve_true$word_long, rownames(covs)), ])
      } else {
        curve_info <- cbind(curve_info,
                            covs[match(curve_info$word_long %% 2^num_cov +1,
                                       rownames(covs)), ])
        curve_true <- cbind(curve_true,
                            covs[match(curve_true$word_long %% 2^num_cov +1,
                                       rownames(covs)), ])
      }
    } else {
      # For the nested scenario, only the Snooker Data scenario is allowed
      if (num_cov != 4) {
        warning(paste0("For the nested scenario, the data have to be similar",
                       " to the Snooker data. num_cov is set to 4."))
        num_cov <- 4
      }

      # Covariate 1: random draw but has to be the same over all observations of
      # the individual.
      cov_1 <- rbinom(n = I, size = 1, prob = 0.5)
      curve_info[, covariate.1 := cov_1[subject_long]]
      curve_true[, covariate.1 := cov_1[subject_long]]

      # Covariate 2: random draw but has to be the same over all observations of
      # the individual.
      cov_2 <- rbinom(n = I, size = 1, prob = 0.5)
      curve_info[, covariate.2 := cov_2[subject_long]]
      curve_true[, covariate.2 := cov_2[subject_long]]

      # Covariate 3: session is the same as word_long.
      curve_info[, covariate.3 := word_long - 1]
      curve_true[, covariate.3 := word_long - 1]

      # Covariate 4: the interaction between group (2) and session (3)
      curve_info[, covariate.4 := covariate.2 * covariate.3]
      curve_true[, covariate.4 := covariate.2 * covariate.3]
    }
  }

  # For the nested scenario, only the Snooker Data scenario is allowed
  if (interaction & nested) {
    warning(paste0("Covariate interactions not implemented for the nested ",
                   "scenario. Interaction set to FALSE. \n"))
    interaction <- FALSE
  }

  # Interaction between the covariates
  # There are only binary covariates
  if (interaction) {
    inters <- list()
    inters_true <- list()
    for (i in 1:num_cov) {
      for (k in 1:num_cov) {
        if (which_inter[i, k] & (i < k)) {
          inters[[paste0("inter.", i, ".", k)]] <-
            curve_info[[paste0("covariate.", i)]] *
            curve_info[[paste0("covariate.", k)]]
          inters_true[[paste0("inter.", i, ".", k)]] <-
            curve_true[[paste0("covariate.", i)]] *
            curve_true[[paste0("covariate.", k)]]
        }
      }
    }
    curve_info <- cbind(curve_info, do.call(cbind,inters))
    curve_true <- cbind(curve_true, do.call(cbind,inters_true))
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Draw observation times for each curve
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Locations are randomly chosen from uniform distribution
  # Locations of measurements for each person as a list, because different
  # numbers
  locs <- list()

  # For each curve i = 1, ..., n on all the dimensions
  for (i in 1:(n*num_dim)) {
    locs[[i]] <- runif(n = number[i], min = min_grid, max = max_grid)
  }

  # Now bring time points in order before functions are constructed
  t <- list()
  o <- lapply(locs, order)
  for (i in 1:(n*num_dim)) {
    t[[i]] <- locs[[i]][o[[i]]]
  }

  # Attach information to the data
  t_vec <- unlist(t)
  curve_info[, t := t_vec]


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Draw scores from normal distributions
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # In normal case: x_i is drawn from N(0, lambda_k)
  if (normal == TRUE) {
    if (N_B > 0) {
      xiB <- mvtnorm::rmvnorm(I, mean = rep(0, N_B), sigma = diag(N_B)*lamB)
    } else {
      xiB <- rep(NA, length = I)
    }

    if (!use_RI) {
      if (N_C > 0) {
        if (!nested) {
          xiC <- mvtnorm::rmvnorm(J, mean = rep(0, N_C), sigma = diag(N_C)*lamC)
        } else {
          xiC <- mvtnorm::rmvnorm(I*J, mean = rep(0, N_C),
                                  sigma = diag(N_C)*lamC)
        }

      } else {
        if(!nested) {
          xiC <- rep(NA, length = J)
        } else {
          xiC <- rep(NA, length = I*J)
        }

      }
    }

    if (N_E > 0) {
      xiE <- mvtnorm::rmvnorm(n, mean = rep(0, N_E), sigma = diag(N_E)*lamE)
    } else {
      xiE <- rep(NA, length = n)
    }
  }

  # Mixture of normals
  if (normal == FALSE) {
    if (nested) {
      stop("Mixture of normals not implemented for nested model.")
    }
    if (N_B > 1) {
      xiB <- matrix(rnorm(I*N_B), I, N_B) %*% diag(sqrt(lamB/2)) +
        matrix(2*rbinom(I*N_B, 1, 0.5) - 1, I, N_B) %*% diag(sqrt(lamB/2))
    } else {
      if (N_B > 0) {
        xiB <- matrix(rnorm(I*N_B), I, N_B) %*% sqrt(lamB/2) +
          matrix(2*rbinom(I*N_B, 1, 0.5) -1, I, N_B) %*% sqrt(lamB/2)
      } else {
        xiB <- rep(NA, length = I)
      }
    }

    if (!use_RI) {
      if (N_C > 1) {
        xiC <- matrix(rnorm(J*N_C), J, N_C) %*% diag(sqrt(lamC/2)) +
          matrix(2*rbinom(J*N_C, 1, 0.5) - 1, J, N_C) %*% diag(sqrt(lamC/2))
      } else {
        if (N_C > 0) {
          xiC <- matrix(rnorm(J*N_C), J, N_C) %*% sqrt(lamC/2) +
            matrix(2*rbinom(J*N_C, 1, 0.5) - 1, J, N_C) %*% sqrt(lamC/2)
        } else {
          xiC <- rep(NA, length = J)
        }
      }
    }

    if (N_E > 1) {
      xiE <- matrix(rnorm(n*N_E), n, N_E) %*% diag(sqrt(lamE/2)) +
        matrix(2*rbinom(n*N_E, 1, 0.5) - 1, n, N_E) %*% diag(sqrt(lamE/2))
    } else {
      if (N_E > 0) {
        xiE <- matrix(rnorm(n*N_E), n, N_E) %*% sqrt(lamE/2) +
          matrix(2*rbinom(n*N_E, 1, 0.5) - 1, n, N_E) %*% sqrt(lamE/2)
      } else {
        xiE <- rep(NA, length = n)
      }
    }
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Center scores
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if (center_scores) {
    if (N_B > 0) {
      for (k in 1:N_B) {
        xiB[, k] <- xiB[, k] - mean(xiB[, k])
      }
    }

    if (!use_RI) {
      if (N_C > 0) {
        for (k in 1:N_C) {
          xiC[, k] <- xiC[, k] - mean(xiC[, k])
        }
      }
    }

    if (N_E > 0) {
      for (k in 1:N_E) {
        xiE[, k] <- xiE[, k] - mean(xiE[, k])
      }
    }
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Decorrelate and scale scores
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Divide by covariance to standardize the score-covariance and then scale with
  # actual covariance
  if (decor_scores) {
    if (N_B > 1)
      xiB <- xiB %*% expm::sqrtm(solve(cov(xiB))) %*% expm::sqrtm(diag(lamB,
                                                                       N_B,
                                                                       N_B))

    if (!use_RI) {
      if (N_C > 1)
        xiC <- xiC %*% expm::sqrtm(solve(cov(xiC))) %*% expm::sqrtm(diag(lamC,
                                                                         N_C,
                                                                         N_C))
    }

    if (N_E > 1)
      xiE <- xiE %*% expm::sqrtm(solve(cov(xiE))) %*% expm::sqrtm(diag(lamE,
                                                                       N_E,
                                                                       N_E))
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Attach scores to data set
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # For subject/first grouping variable
  if (N_B > 0) {
    # Repeat each column according to the data
    xiB_long_help <- list()
    for (k in 1:N_B) {
      xiB_long_help[[k]] <- rep(xiB[, k][rep(subject, num_dim)], number)
    }
    # Retransform the list
    if (N_B > 1) {
      xiB_long <- do.call(cbind, xiB_long_help)
    } else {
      xiB_long <- matrix(xiB_long_help[[1]], ncol = 1)
    }
    # Attach the scores to the data.table
    for (k in 1:N_B) {
      curve_info[, paste0("xiB_long.", k) := xiB_long[, k]]
    }
  }

  # For word/second grouping variable
  if (!use_RI) {
    if (N_C > 0) {
      if (!nested) {
        # Repeat each column according to the data
        xiC_long_help <- list()
        for (k in 1:N_C) {
          xiC_long_help[[k]] <- rep(xiC[, k][rep(word, num_dim)], number)
        }
      } else {
        # In the nested case repeat each score per subject and word
        xiC_long_help <- list()
        for (k in 1:N_C) {
          xiC_long_help[[k]] <- rep(xiC[, k][
            rep(as.numeric(interaction(subject, word, lex.order = TRUE)),
                num_dim)], number)
        }
      }
      # Retransform the list
      if (N_C > 1) {
        xiC_long <- do.call(cbind, xiC_long_help)
      } else {
        xiC_long <- matrix(xiC_long_help[[1]], ncol = 1)
      }
      # Attach the scores to the data.table
      for (k in 1:N_C) {
        curve_info[, paste0("xiC_long.", k) := xiC_long[,k]]
      }
    }
  }

  # For curve-level
  if (N_E > 0) {
    # Repeat each column according to the number of dimensions
    xiE_long_help <- list()
    for (k in 1:N_E) {
      xiE_long_help[[k]] <- rep(xiE[, k][rep(1:n, num_dim)], number)
    }
    # Retransform the list
    xiE_long <- xiE_long_help[[1]]
    if(N_E > 1) {
      xiE_long <- do.call(cbind, xiE_long_help)
    } else {
      xiE_long <- matrix(xiE_long_help[[1]], ncol = 1)
    }
    # Attach the scores to the data.table
    for (k in 1:N_E) {
      curve_info[, paste0("xiE_long.", k) := xiE_long[, k]]
    }
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Compute the Eigenfunctions on each observation time
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # In order to construct the true underlying functions B_i(t), C_j(t) and
  # E_ijh(t) via KL-expansion
  # List of times on each dimension
  t_times <- split(curve_info$t, rep(1:num_dim, table(curve_info$dim)))

  if (N_B > 0) {
    # Separate the dimensions
    phiB_list <- list()
    for (i in 1:num_dim) {
      phiB_list[[i]] <- phiB_fun@.Data[[i]]
    }
    # Use predict function for each funData object separately
    # Order the function values according to the observed values into a matrix
    phiB_list <- mapply(function (x, y) {
      x <- predict_funData(funData = x, new = y)
      t(x@X[, match(y, unlist(argvals(x)))])
    }, phiB_list, t_times, SIMPLIFY = FALSE)
    if (N_B > 1) {
      phiB <- do.call(rbind, phiB_list)
    } else {
      phiB <- matrix(unlist(phiB_list), ncol = 1)
    }

  } else {
    phiB <- rep(NA, length = length(t_vec))
  }

  if (!use_RI) {
    if (N_C > 0) {
      # Separate the dimensions
      phiC_list <- list()
      for (i in 1:num_dim) {
        phiC_list[[i]] <- phiC_fun@.Data[[i]]
      }
      # Use predict function for each funData object separately
      # Order the function values according to the observed values into a matrix
      phiC_list <- mapply(function (x, y) {
        x <- predict_funData(funData = x, new = y)
        t(x@X[, match(y, unlist(argvals(x)))])
      }, phiC_list, t_times, SIMPLIFY = FALSE)
      if (N_C > 1) {
        phiC <- do.call(rbind, phiC_list)
      } else {
        phiC <- matrix(unlist(phiC_list), ncol = 1)
      }
    } else {
      phiC <- rep(NA, length= length(t_vec))
    }
  }

  if (N_E > 0) {
    # Separate the dimensions
    phiE_list <- list()
    for (i in 1:num_dim) {
      phiE_list[[i]] <- phiE_fun@.Data[[i]]
    }
    # Use predict function for each funData object separately
    # Order the function values according to the observed values into a matrix
    phiE_list <- mapply(function (x, y) {
      x <- predict_funData(funData = x, new = y)
      t(x@X[, match(y, unlist(argvals(x)))])
    }, phiE_list, t_times, SIMPLIFY = FALSE)
    if (N_E > 1) {
      phiE <- do.call(rbind, phiE_list)
    } else {
      phiE <- matrix(unlist(phiE_list), ncol = 1)
    }
  } else {
    phiE <- rep(NA, length = length(t_vec))
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Weigh the Eigenfunctions and attach to data
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if (N_B > 0) {
    # matrix with columns the components of each k=1,..,N_B
    B_vec_help <- matrix(NA, ncol = N_B, nrow = length(n_long))
    for (k in 1:N_B) {
      B_vec_help[,k] <- curve_info[[paste("xiB_long.", k, sep = "")]]*phiB[, k]
    }
    B_vec <- rowSums(B_vec_help)
    curve_info[, "B_long" := B_vec]
  } else {
    B_vec <- rep(0, length = length(t_vec))
  }

  if (!use_RI) {
    if (N_C > 0) {
      # matrix with columns the components of each k=1,..,N_B
      C_vec_help <- matrix(NA, ncol = N_C, nrow = length(n_long))
      for (k in 1:N_C) {
        C_vec_help[, k] <- curve_info[[paste("xiC_long.", k, sep = "")]]*
          phiC[, k]
      }
      C_vec <- rowSums(C_vec_help)
      curve_info[, "C_long" := C_vec]
    } else {
      C_vec <- rep(0, length = length(t_vec))
    }
  }

  if (N_E > 0) {
    E_vec_help <- matrix(NA, ncol = N_E, nrow = length(n_long))
    for (k in 1:N_E) {
      E_vec_help[, k] <- curve_info[[paste("xiE_long.", k, sep = "")]]*phiE[, k]
    }
    E_vec <- rowSums(E_vec_help)
    curve_info[, "E_long" := E_vec]
  } else {
    E_vec <- rep(0, length = length(t_vec))
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Draw measurement error
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Weighted standard normal random numbers with simple_sig == TRUE
  # Otherwise use multivariate normal distribution (two implementations)
  if (simple_sig) {
    epsilon_vec <- rnorm(n = nrow(curve_info), mean = 0, sd = 1) * rep(
      sqrt(unlist(sigmasq)), times = table(dim)
    )
  } else {
    # Use lists of length n*num_dim
    epsilon <- list()
    # Jona's original code just prolonged
    if (!dim_indep) {
      for (i in 1:(n*num_dim)) {
        if (number[i] == 1) {
          epsilon[[i]] <- mvtnorm::rmvnorm(1, mean = rep(0, number[i]),
                                           sigma = as.matrix(sigmasq, nrow = 1,
                                                             ncol = 1))
        } else {
          epsilon[[i]] <- mvtnorm::rmvnorm(1, mean = rep(0, number[i]),
                                           sigma = diag(
                                             rep(sigmasq, length = number[i])))
        }
      }
    } else {
      # Own implementation using mvnorm
      numb_dims <- do.call(cbind, split(number, rep(1:num_dim, n)))
      n_dims <- rowSums(numb_dims)
      for (i in 1:n) {

        # Create a list containing the covariance matrices for each subject
        row_list <- mapply(function (x, y) {
          sigma = diag(rep(y, length = x))
        }, numb_dims[i, ], sigmasq, SIMPLIFY = FALSE)

        # If there is only one observation, diag does not work in the same way
        if (1 %in% numb_dims[i, ]) {
          pos <- which(numb_dims[i, ] == 1)
          for (k in pos) {
            row_list[[k]] <- as.matrix(sigmasq[[k]], nrow = 1, ncol = 1)
          }
        }

        # Draw from mvnorm and split the random numbers on
        eps <- mvtnorm::rmvnorm(1, mean = rep(0, n_dims[i]),
                                sigma = as.matrix(do.call(bdiag, row_list)))
        eps <- split(eps, rep(1:num_dim, numb_dims[i, ]))

        # Bring the measurement errors in the right order
        for (k in 0:(num_dim - 1)) {
          epsilon[[i + k*n]] <- eps[[k + 1]]
        }
      }
    }
    epsilon_vec <- unlist(epsilon)
  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Attach the covariate effects
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # If model is NULL then Jona's implementation
  # So far, the same mu function for all dimensions
  if (is.null(model)) {
    mu_list <- list()
    for (i in 1:n) {
      t_use <- subset(curve_info, select = t, n_long == i)
      mu_list[[i]] <- mu(t_use)
    }
    mu_vec <- unlist(mu_list)

  } else {

    # Own implementation that takes the model into account
    # Use model object for prediction
    newdata <- prepare_gam_predict(data = curve_info, num_cov = num_cov,
                                   interaction = interaction,
                                   which_inter = which_inter, model = model)
    newtrue <- prepare_gam_predict(data = curve_true, num_cov = num_cov,
                                   interaction = interaction,
                                   which_inter = which_inter, model = model)

    predictions <- predict(model, newdata = newdata, type = "terms")
    predicttrue <- predict(model, newdata = newtrue, type = "terms")

    # Sum together only the terms without the random effects
    # ignore model terms that contain n_long, subject_long or word_long
    use <- ! grepl(paste(c("subject_long", "word_long", "n_long"),
                         collapse = "|"), colnames(predictions))
    use_true <- ! grepl(paste(c("subject_long", "word_long", "n_long"),
                              collapse = "|"), colnames(predicttrue))
    mu_vec <- rowSums(predictions[, use])
    curve_true[, "mu_true" := rowSums(predicttrue[, use_true])]

  }


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Construct observations y_ij
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if (!use_RI) {
    y_vec <- mu_vec + B_vec + C_vec + E_vec + epsilon_vec
  } else {
    y_vec <- mu_vec + B_vec + E_vec + epsilon_vec
  }

  curve_info[, c("mu_vec", "epsilon_vec", "y_vec") := list(mu_vec, epsilon_vec,
                                                           y_vec)]


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Transform true curves to multiFunData object
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  curve_true <- split(curve_true, curve_true$dim)
  mu <- funData::multiFunData(lapply(curve_true, function (x) {
    funData::funData(argvals = seq(min_grid, max_grid, length.out = grid_eval),
                     X = matrix(x$mu_true, ncol = grid_eval, byrow = TRUE))
  }))


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Compute true functional random effects
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Construct lists for easier handling
  if (!use_RI) {
    scores <- list("E" = xiE, "B" = xiB, "C" = xiC)
    eigfct <- list("E" = phiE_fun, "B" = phiB_fun, "C" = phiC_fun)
  } else {
    scores <- list("E" = xiE, "B" = xiB)
    eigfct <- list("E" = phiE_fun, "B" = phiB_fun)
  }

  # Construct the random effects using the scores and Eigenfunctions
  re <- mapply(function (x, y) {
    lapply(x@.Data, function (z) {
      y %*% z@X
    })

  }, eigfct, scores, SIMPLIFY = FALSE)

  # Convert list to multiFunData object for better comparison
  re <- re[sapply(re, function (x) length(x) > 0)]
  re <- lapply(re, function (x) {
    funData::multiFunData(lapply(x, function (y) {
      funData::funData(argvals = argvals(eigfct[[1]][[1]]),
                       X = y)
    }))
  })


  #++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # Combine to Output
  #++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Initialize output
  res <- list()

  res[["curve_info"]] <- curve_info
  res[["xiB"]] <- xiB
  if (!use_RI) res[["xiC"]] <- xiC
  res[["xiE"]] <- xiE
  res[["mu"]] <- mu
  res[["re"]] <- re

  return(res)

}


# Predict funData Object for Unobserved Values ----------------------------

#' Predict funData Object for Unobserved Values
#'
#' This is an internal function. It takes a funData object and approximates the
#' function values for the supplied new time points using the adapted approxNA()
#' function. The function checks whether the time points in the argument new are
#' already included in the argvals. Predict is probably not the right word but
#' it is used for consistency.
#'
#' @param funData funData object to be evaluated.
#' @param new Vector containing the time points to be approximated.
predict_funData <- function(funData, new) {

  # Are there already observations of this
  new <- new[! new %in% unlist(funData::argvals(funData))]
  new <- unique(new)

  # Reorder the data
  len <- length(new)
  ord <- order(c(unlist(funData::argvals(funData)), new))
  fun_new <- funData::funData(argvals = c(unlist(funData::argvals(funData)),
                                          new)[ord],
                              X = matrix(
                                cbind(funData@X,
                                      matrix(NA, nrow = funData::nObs(funData),
                                             ncol = len))[, ord],
                                nrow = funData::nObs(funData)))

  approxNA(fun_new)

}



# Prepare Prediction of Simulated Data Using MGCV Model Object ------------

#' Prepare Prediction of Simulated Data Using MGCV Model Object
#'
#' This is an internal function. It is used to generate new data sets. It adapts
#' curve_info so that it has the same structure as the original data. Then, it
#' is possible to use the predict function from mgcv.
#'
#' @param data Simulated data set for which the necessary model components are
#'   to be computed.
#' @param num_cov Number of covariates if covariate = TRUE.
#' @param interaction TRUE if there are interactions between covariates (as in
#'   sparseFLMM).
#' @param which_inter Symmetric matrix specifying the interaction terms (as in
#'   sparseFLMM).
#' @param model GAM model from which to extract the covariate functions.
prepare_gam_predict <- function (data, num_cov, interaction, which_inter,
                                 model) {

  # Rename the dimensions to match the model
  levels(data$dim) <- gsub("dim", "", names(model$coefficients)[
    grepl("^dim", names(model$coefficients))])

  # Transform data in order to have factor variables with factors of original
  # model
  facs <- c("n_long", "subject_long", "word_long")
  facs <- facs[facs %in% colnames(model$model)]
  new_values <- data.frame(sapply(model$model[, facs], function(x){
    factor(rep(levels(x)[1], nrow(data)))
  }))
  data.table::setDT(data)[, (facs):= new_values]

  # Interactions of covariates with dimension
  if (num_cov > 0) {
    for (i in 1:num_cov) {

      # Create the necessary interaction variables
      form_temp <- as.formula(paste0("~ 0 + dim:covariate.", i))
      tmp <- model.matrix(form_temp, data = data)
      colnames(tmp) <- sub("\\:covariate", "", colnames(tmp))
      data <- cbind(data, tmp)

    }
  }

  # Include Interaction between covariates
  if (interaction) {
    for (i in 1:num_cov) {
      for (k in 1:num_cov) {
        if (which_inter[i, k] & (i < k)) {

          # Create the necessary interaction variables
          form_temp <- as.formula(paste0("~ 0 + dim:covariate.", i,
                                         ":covariate.", k))
          tmp <- model.matrix(form_temp, data = data)
          colnames(tmp) <- sub("\\:covariate", "\\.inter", colnames(tmp))
          colnames(tmp) <- sub("\\:covariate", "", colnames(tmp))
          data <- cbind(data, tmp)

        }
      }
    }
  }

  # Include weighted principal components for prediction
  wPC <- attr(attr(model$model, which = "terms"), "term.labels")
  wPC <- wPC[grepl("^w[BCE]_", wPC)]
  tmp <- data.frame(matrix(0, nrow = nrow(data), ncol = length(wPC)))
  names(tmp) <- wPC
  data <- cbind(data, tmp)

  data
}



# Create the Coverage Array of Simulated Covariate Effects ----------------

#' Create the Coverage Array of the Simulation
#'
#' This is an internal function. The function takes the index of the covariate
#' to be evaluated and then checks whether the estimated covariate effect of the
#' simulation run covers the true data generating effect function. The output is
#' a logical array where the first dimension gives the dimension of the data,
#' the second dimension gives the time point to be evaluated and the third
#' dimension gives the simulation run.
#'
#' @param sim_curves The large list of simulation results. Use object$mul.
#' @param gen_curves The original data generating curve as part of the output of
#'   multifamm:::extract_components(), so use output$cov_preds.
#' @param effec_index The index position of the effect to be evaluated in the
#'   gen_curves and sim_curves effect lists. If the intercept is to be
#'   evaluated, this can be specified as 1 or 2 (both scalar and functional
#'   intercept are sumed up).
#' @param m_fac Multiplication factor used to create the upper and lower
#'   credibility bounds. Defaults to 1.96 (ca. 95\%).
create_coverage_array <- function (sim_curves, gen_curves, effect_index,
                                   m_fac = 1.96) {

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
    gen <- gen_curves$fit[[1]] + gen_curves$fit[[2]]
  } else {
    sim_up <- lapply(sim_curves, function (it) {
      it$fit[[effect_index]] + m_fac * (it$se.fit[[effect_index]])
    })
    sim_do <- lapply(sim_curves, function (it) {
      it$fit[[effect_index]] - m_fac * (it$se.fit[[effect_index]])
    })
    gen <- gen_curves$fit[[effect_index]]
  }

  # Check if the original data is inside the simulated bounds
  coverage <- mapply(function (sim_up_it, sim_do_it) {
    t(mapply(FUN = function (s_u, s_d, o) {
      o@X > s_d@X & o@X < s_u@X
    }, sim_up_it@.Data, sim_do_it@.Data, gen@.Data))
  }, sim_up, sim_do, SIMPLIFY = FALSE)

  # Order the results in an array
  coverage <- array(unlist(coverage), dim = c(length(gen), 100,
                                              length(sim_curves)))
  coverage

}



# Prepare the Data for the Coverage Plot ----------------------------------

#' Coverage plot helper function
#'
#' This is an internal function. The function takes a list of arrays created by
#' the function create_coverage_array and returns a data.frame ready for
#' plotting.
#'
#' @param cov_list List of arrays containing the evaluations if the estimated
#'   coefficient effect of that simulation run is in the credibility bounds as
#'   given by the function create_coverage_array().
#' @param effect_index Numeric vector of the index positions of the effects to
#'   be plotted.
#' @param dimlabels String vector of labels used in the data set. Defaults to
#'   labels "ACO" and "EPG".
coverage_plot_helper <- function (cov_list, effect_index,
                                  dimlabels = c("ACO", "EPG")) {

  # Get the proportion of covered effects for each index
  prop_list <- list()
  for (ind in seq_along(effect_index)) {
    prop_list[[ind]] <- apply(cov_list[[effect_index[ind]]], MARGIN = c(1, 2),
                              function (it) {sum(it) / length(it)})
  }

  # Rearrange the proportions to a data set
  prop_list <- lapply(prop_list, function (effect) {
    data.frame(t = rep(seq(0, 1, length.out = 100), times = nrow(effect)),
               y = c(t(effect)),
               dim = rep(dimlabels, each = 100))
  })

  # Combine the effects to a data.frame and label them
  dat <- do.call(rbind, prop_list)
  dat$effect <- factor(rep(seq_along(effect_index),
                           each = nrow(prop_list[[1]])),
                       labels = paste0("f[", effect_index - 1, "](t)"))
  dat

}


