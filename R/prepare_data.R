

#------------------------------------------------------------------------------#
# Attach Weighted Functional Principal Components to the Data
#------------------------------------------------------------------------------#
attach_wfpc <- function(MFPC, data){

  # Arguments
  # MFPC        : List containing MFPC objects for each coponent of the variance
  # data        : Data where wfPC are attached to


  # Eigenfunctions have to be evaluated on the actual observed observation times
  interpol <- linear_interpol(MFPC, data)

  # Lists of Eigenvalues and argvals per variance component
  values <- lapply(MFPC, function(x) x <- x$values)
  argvals <- lapply(interpol, function(x){
    lapply(x, function(y) y <- unlist(argvals(y)))
  })

  # List of weighted eigenfunctions
  w_interpol <- mapply(function(x, y, z){
    mapply(function(u, v) {
      wE <- t(u@X * sqrt(y))
      rownames(wE) <- v
      wE
    }, x, z, SIMPLIFY = FALSE)
  }, interpol, values, argvals, SIMPLIFY = FALSE)

  # Ensure correct order in data
  # Ordered according to the dimensions
  # Rest of ordering is only for convenience
  data <- data[order(data$dim, data$n_long, data$t), ]

  # Order the weigthed eigenfunctions according to the observations in the data
  # Combine all the dimensions so that the matrices have equal nrow than data
  data_list <- split(data, data$dim)
  w_interpol <- lapply(w_interpol, function(x){
    out <- mapply(function(y, z){
      y[match(z$t, rownames(y)), ]
    }, x, data_list, SIMPLIFY = FALSE)
    if(dim(x[[1]])[2] == 1) out <- lapply(out, matrix, ncol = 1)
    do.call(rbind, out)
  })

  # Rename the elements so they can be distinguished when attached to the data
  w_interpol <- lapply(seq_along(w_interpol), function(i){
    colnames(w_interpol[[i]]) <- paste0("w", names(w_interpol)[[i]], "_",
                                        1:ncol(w_interpol[[i]]))
    w_interpol[[i]]
  })

  # Combine the original data with the weighted Eigenfunctions
  data <- cbind(data, do.call(cbind, w_interpol))
  data

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Evaluate Eigenfunctions on the Observed Data Times
#------------------------------------------------------------------------------#
linear_interpol <- function(MFPC, data){

  # Arguments
  # MFPC        : List containing MFPC objects for each coponent of the variance
  # data        : Data where wfPC are attached to


  # Which observations have to be interpolated
  # Keep track of the order to order the X values appropriately
  obs_t <- lapply(split(data$t, data$dim), unique)
  grid <- unlist(argvals(MFPC[[1]]$functions), recursive = FALSE)
  arg_vals <- mapply(function(x,y){
    unique(c(y, x))
  }, obs_t, grid, SIMPLIFY = FALSE)
  arg_order <- lapply(arg_vals, order)

  # Attach NAs to the observed values to later interpolate them
  # arg_l helps with mapply function as it has the same structure as obs_x
  obs_x <- lapply(MFPC, function(x) X(x$functions))
  arg_l <- lapply(obs_x, function(x) x <- arg_vals)
  obs_x <- mapply(function(x, y){
    mapply(function(u, v){
      u <- cbind(u, matrix(NA, nrow = nrow(u), ncol = length(v)-ncol(u)))
      u
    }, x, y, SIMPLIFY = FALSE)
  }, obs_x, arg_l, SIMPLIFY = FALSE)

  # Reorder the argvals and xvalues
  arg_vals <- mapply(function(x, y){
    x[y]
  }, arg_vals, arg_order, SIMPLIFY = FALSE)
  arg_l <- lapply(obs_x, function(x) x <- arg_order)
  obs_x <- mapply(function(x, y){
    mapply(function(u, v){
      u <- u[, v]
      u
    }, x, y, SIMPLIFY = FALSE)
  }, obs_x, arg_l, SIMPLIFY = FALSE)

  # Convert the data back to funData object to apply the interpolation
  arg_l <- lapply(obs_x, function(x) x <- arg_vals)
  out <- mapply(function(x, y){
    mapply(function(u, v){
      funData(argvals = list(u), X = matrix(v, ncol = length(u)))
    }, x, y, SIMPLIFY = FALSE)
  }, arg_l, obs_x, SIMPLIFY = FALSE)

  # Apply the interpolation
  out <- lapply(out, function(x) lapply(x, approxNA))
  out

}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Adapt funData Function approxNA to Consider x Values
#------------------------------------------------------------------------------#
# Code stems from R package funData Version 1.1
# Only adaptation: argument 'x = unlist(object@argvals)'
approxNA <- function (object) {
  funData::funData(object@argvals, t(zoo::na.approx(t(object@X),
                                                    x = unlist(object@argvals),
                                                    na.rm = FALSE)))
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Prepare Data Set for Modeling with GAM
#------------------------------------------------------------------------------#
prepare_gam <- function(data, bs, bf_covariates, m_mean, covariate,
                        num_covariates, covariate_form, interaction,
                        which_interaction, nested){

  # Arguments
  # data              : Data that contains all the covariates
  # bs                : Spline basis function, only tested for "ps" (as in
  #                       sparseFLMM)
  # bf_covariates     : Basis dimension for all covariates (as in sparseFLMM)
  # m_mean            : Order of penalty for basis function (as in sparseFLMM)
  # covariate         : Covariate effects (as in sparseFLMM)
  # num_covariates    : Number of covariates included in the model (as in
  #                       sparseFLMM)
  # covariate_form    : Vector of strings for type of covariate (as in
  #                       sparseFLMM)
  # interaction       : TRUE if there are interactions between covariates (as in
  #                       sparseFLMM)
  # which_interaction : Symmetric matrix specifying the interaction terms (as in
  #                       sparseFLMM)


  # Transform data in order to have factor variables
  facs <- c("n_long", "subject_long", "word_long")
  facs <- facs[facs %in% colnames(data)]
  setDT(data)[, (facs):= lapply(.SD, as.factor), .SDcols=facs]

  # In the nested model case, the grouping variable word_long has to have I*J
  # different levels instead of only J levels
  if (nested) {
    setnames(data, "word_long", "word_long_orig")
    data.table::set(data, j = "word_long",
                    value = interaction(data$subject_long, data$word_long_orig))
  }

  # Rest of function is only necessary to create the interaction variables for
  # the covariates
  if(!covariate) {
    # No covariates
    return(data)
  }

  # Interactions of covariates with dimension
  terms <- vector()
  for (i in 1:num_covariates) {
    if (covariate_form[i] == "by") {

      # Create the necessary interaction variables
      form_temp <- stats::as.formula(paste0("~ 0 + dim:covariate.", i))
      tmp <- stats::model.matrix(form_temp, data = data)
      colnames(tmp) <- sub("\\:covariate", "", colnames(tmp))
      data <- cbind(data, tmp)

      # Create the necessary gam term
      terms <- c(terms,
                 paste0("s(t, k = ", bf_covariates, ", bs = \"", bs, "\", m = ",
                        paste0("c(", paste(m_mean, collapse = ","), ")"),
                        ", by = ", colnames(tmp), ")"))

    } else {
      cat("Experimental usage of covariable of type smooth\n")
      terms <- c(terms,
                 paste0("ti(t, covariate.", i, ", by = dim, k = ",
                        bf_covariates, ", bs = \"", bs, "\", m = ",
                        paste0("c(", paste(m_mean, collapse = ","), ")"),
                        ", mc = c(0, 1))"))
    }
  }

  # Include Interaction between covariates
  if(interaction) {
    for (i in 1:num_covariates) {
      for (k in 1:num_covariates) {
        if (which_interaction[i, k] & (i < k)) {

          # Create the necessary interaction variables
          form_temp <- stats::as.formula(paste0("~ 0 + dim:covariate.", i,
                                         ":covariate.", k))
          tmp <- stats::model.matrix(form_temp, data = data)
          colnames(tmp) <- sub("\\:covariate", "\\.inter", colnames(tmp))
          colnames(tmp) <- sub("\\:covariate", "", colnames(tmp))
          data <- cbind(data, tmp)

          # Create the necessary gam term
          terms <- c(terms,
                     paste0("s(t, k = ", bf_covariates, ", bs = \"", bs,
                            "\", m = ",
                            paste0("c(", paste(m_mean, collapse = ","), ")"),
                            ", by = ", colnames(tmp), ")"))
        }
      }
    }
  }

  # Output is the transformed data and the covariate terms as attributes
  attr(data, "gam_terms") <- terms
  data


}
#------------------------------------------------------------------------------#

