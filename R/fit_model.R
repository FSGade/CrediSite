#' Fit model from data
#'
#' @param model_type Type of model used to fit. "baseline", "regularised" and "blosum_enhanced" are accepted.
#' @param genotype Genotype dataframe
#' @param phenotype Phenotype vector
#' @param cores STAN: Number of cores
#' @param iter STAN: Number of iterations
#' @param adapt_delta STAN: Adapt delta
#' @param max_treedepth STAN: Max tree depth
#'
#' @return
#' @export
#'
#' @examples fit_model(X, y, model_type="regularised")
fit_model <- function(genotype, phenotype, model_type="regularised", cores=4, iter=2500, adapt_delta=0.95, max_treedepth=15){
  if (is.null(genotype)) {
    stop("Genotype object is NULL.")
  } else if (is.null(phenotype)) {
    stop("Phenotype object is NULL.")
  } else if (!is.data.frame(genotype) || !is.list(genotype)) {
    stop("Genotype object must be a dataframe or list.")
  } else if (!is.vector(phenotype)) {
    stop("Phenotype object must be a vector.")
  }
  if (dim(genotype)[1] != length(phenotype)) {
    stop(paste("Dimension mismatch:", dim(genotype), length(phenotype)))
  }

  aa_list <- strsplit("ARNDCEQGHILKMFPSTWYVX.*", "")[[1]]
  hs_df = 3
  hs_df_global = 1
  hs_df_slab = 4
  hs_scale_global = 1
  hs_scale_slab = 2

  message("Creating model data")
  model_X <- model.matrix(~ . + 0, data=genotype)

  if (model_type=="baseline") {
    model_code = credisite:::baseline_code

    model_data = list(
      N=dim(model_X)[1],
      K=dim(model_X)[2],
      Y=phenotype,
      X=model_X
    )

  } else if (model_type=="regularised") {
    model_code = credisite:::regularised_code

    model_data = list(
      N=dim(model_X)[1],
      K=dim(model_X)[2],
      Y=phenotype,
      X=model_X,
      hs_df=hs_df,
      hs_df_global=hs_df_global,
      hs_df_slab=hs_df_slab,
      hs_scale_global=hs_scale_global,
      hs_scale_slab=hs_scale_slab
    )

  } else if (model_type=="blosum_enhanced") {
    num_pos_aa = rep(0, length(genotype))
    for (i in seq_along(genotype)) {
      num_pos_aa[i] = length(levels(genotype[, i])) - 1
    }

    num_b = sum(num_pos_aa)
    mu_index = rep(0, num_b)
    aa_b = rep("", num_b)

    i = 1
    j = 1
    while (i <= num_b) {
      for (k in 1:num_pos_aa[j]) {
        mu_index[i] = j
        aa_b[i] = levels(genotype[, j])[k+1]
        i = i + 1
      }
      j = j + 1
    }


    COR_MAT = diag(22)
    COR_MAT[1:20, 1:20] = credisite:::BLOSUM62_COR
    colnames(COR_MAT) = c(colnames(credisite:::BLOSUM62_COR), "X", "-")
    rownames(COR_MAT) = c(rownames(credisite:::BLOSUM62_COR), "X", "-")

    cor_mat = diag(num_b)
    cumsum_aa = cumsum(num_pos_aa)
    for (i in seq_along(cumsum_aa)) {
      start_i = cumsum_aa[max(i - 1, 1)]
      if (i != 1) {
        start_i = start_i + 1
      }
      end_i = cumsum_aa[i]
      cor_mat[start_i:end_i,start_i:end_i] =
        COR_MAT[aa_b[start_i:end_i], aa_b[start_i:end_i]]
    }

    model_code = credisite:::blosum_enhanced_code
    model_data = list(
      N=dim(model_X)[1],
      K=dim(model_X)[2],
      L=length(num_pos_aa),
      mu_indices=mu_index,
      cor_mat=cor_mat,
      Y=phenotype,
      X=model_X,
      hs_df=hs_df,
      hs_df_global=hs_df_global,
      hs_df_slab=hs_df_slab,
      hs_scale_global=hs_scale_global,
      hs_scale_slab=hs_scale_slab
    )
  } else {
    stop(paste("Model type, ", model_type, ", was not recognised. Only baseline, regularised and blosum_enhanced are accepted.", sep=""))
  }
  message("Executing RStan model...")
  model = rstan::stan(model_code = model_code,
               data = model_data, cores = cores, iter=iter,
               control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))
  message("Done.")
  #TODO: Rename model variable names
  #TODO: Diagnostics of the resulting model


  return(list(pos_names = colnames(model_X)[-1], model = model))
}
