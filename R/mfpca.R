

#------------------------------------------------------------------------------#
# Prepare Information Necessary for MFPCA
#------------------------------------------------------------------------------#
prepare_mfpca <- function(model_list, fRI_B){

  # Arguments
  # model_list    : List containing sparseFLMM objects for each dimension
  # fRI_B         : Boolean for including functional random intercept for
  #                   individual (B in Cederbaum)

  # Extract the necessary information from the model list
  model_info <- extract_mfpca_info(model_list = model_list, fRI_B = fRI_B)

  # Create list containing the information on the maximum of FPC
  M_comp <- lapply(model_info, function (x) {
    sum(sapply(x, function(y) y$N))
  })

  # Inflate the fPCs when there are different numbers of fPCs with zero curves
  model_info <- inflate(model_info)

  # Determine which model terms are necessary
  # Look at the scores whether there are missing values
  missing_dim <-lapply(model_info, function(x){
    lapply(lapply(x, "[[", 3), function(y) any(is.na(y)))
  })
  # If the scores are missing for all the dimensions, remove the model term
  keep_terms <- sapply(missing_dim, function(x) !all(unlist(x)))
  model_info <- model_info[keep_terms]
  M_comp <- M_comp[keep_terms]

  # Reconstruct the different functional Random Intercepts
  # Scores %*% Eigcts and saved
  model_info <- lapply(model_info, function(x){
    lapply(x, function(y){
      y$RI <- y[[3]] %*% t(y[[2]])
      y
    })
  })

  # Construct multiFunData objects for each variance component
  multiFun_comp <- lapply(model_info, function(x){
    multiFunData(lapply(x, function(y){
      funData(argvals = list(y$grid), X = y$RI)
    }))
  })

  # Create list containing the information of uniExpansion
  uniExpansions <- lapply(model_info, function(x){
    lapply(x, function(y){
      list(type = "given",
           functions = funData(argvals = list(y$grid),
                               X = t(y[[2]])),
           scores = y[[3]],
           ortho = TRUE)
    })
  })

  out <- mapply(function(x, y, z){
    list(mFData = x,
         uniExpansions = y,
         M_comp = z)
  }, multiFun_comp, uniExpansions, M_comp, SIMPLIFY = FALSE)
  out

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Extract Information Necessary for Preparing the MFPCA
#------------------------------------------------------------------------------#
extract_mfpca_info <- function(model_list, fRI_B){

  # Arguments
  # model_list    : List containing sparseFLMM objects for each dimension
  # fRI_B         : Boolean for including functional random intercept for
  #                   individual (B in Cederbaum)

  # Extract the name necessary to find the results of the decomposition
  method_name <- lapply(model_list,
                        function(x) names(x)[grep("^fpc_hat_", names(x))])

  # Create a list containing the elements of interest per dimension
  out <- mapply(function(x,y){

    # Extract information on grid
    grid <- x$my_grid

    # Extract information on estimated eigenfunctions
    # Correct for unintuitively handling of independent curves of Jona
    E_f <- if(fRI_B) x[[y]]$phi_E_hat_grid else x[[y]]$phi_B_hat_grid
    B_f <- if(fRI_B) x[[y]]$phi_B_hat_grid else x[[y]]$phi_E_hat_grid
    C_f <- x[[y]]$phi_C_hat_grid

    # Extract information on estimated scores
    E_s <- if(fRI_B) x[[y]]$xi_E_hat else x[[y]]$xi_B_hat
    B_s <- if(fRI_B) x[[y]]$xi_B_hat else x[[y]]$xi_E_hat
    C_s <- x[[y]]$xi_C_hat

    # Extract information on number of fPCs
    N_E <- if(fRI_B) x[[y]]$N_E else x[[y]]$N_B
    N_B <- if(fRI_B) x[[y]]$N_B else x[[y]]$N_E
    N_C <- x[[y]]$N_C

    # Also extract number of PCs
    list(E = list(grid = grid, E_f = E_f, E_s = E_s, N = N_E),
         B = list(grid = grid, B_f = B_f, B_s = B_s, N = N_B),
         C = list(grid = grid, C_f = C_f, C_s = C_s, N = N_C))

  }, model_list, method_name, SIMPLIFY = FALSE)

  # Invert the list to combine the results for all the dimensions
  out <- apply(do.call(rbind, out), 2, as.list)
  out

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Inflate the List of fPCs with Zero cCrves
#------------------------------------------------------------------------------#
inflate <- function(model_info){

  # Arguments
  # model_info    : List containing extracted information from function
  #                   extract_mfpca_info

  lapply(model_info, function(x) {

    # For each variance component extract info of number of fpcs
    N <- sapply(x, function(y) {
      y$N
    })

    # If there is the same number on all dimensions then do nothing
    if (length(unique(N)) == 1) {
      return(x)
    } else {

      # If there are different numbers inflate the second and third components
      # Position 2: Eigenfunctions
      # Position 3: Eignescores
      m <- max(N)

      # Unequal and maximum is 0 needs to be handled differently
      # cbind not possible but create new matrix
      if(m == 1) {

        # Number of groups needed for correct form of score matrix
        i <- which.max(N)
        n_groups <- nrow(x[[i]][[3]])

        # Inflation step
        x <- lapply(x, function(y) {
          if(dim(y[[2]])[2] != m) {
            y[[2]] <- matrix(0, ncol = 1, nrow = length(y[[1]]))
            y[[3]] <- matrix(0, ncol = 1, nrow = n_groups)
          }
          y
        })

      } else {

        # Inflation step
        x <- lapply(x, function(y) {

          # Make sure to only change the matrix if necessary
          if(dim(y[[2]])[2] != m) {
            y[[2]] <- cbind(y[[2]], matrix(0, ncol = m - ncol(y[[2]]),
                                           nrow = nrow(y[[2]])))
            y[[3]] <- cbind(y[[3]], matrix(0, ncol = m - ncol(y[[3]]),
                                           nrow = nrow(y[[3]])))
          }
          y
        })
      }
      return(x)
    }
  })


}
#------------------------------------------------------------------------------#

