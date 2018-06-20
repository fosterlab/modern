#' Detect outliers by network reconstruction
#' 
#' \code{modern} is a method to detect outliers in high-dimensional data
#' based on their impact on network reconstruction.
#' The core idea is that the topology of the network reconstructed from a
#' matrix of data should be robust to the inclusion or exclusion of each 
#' individual data point in the matrix.
#' Single points that have a large impact on the global interaction profile
#' of a node (e.g., a gene, protein, or metabolite) compromise the robustness
#' of network inference, and are likely to be outliers. 
#' The degree to which a single point compromises the robustness of the network
#' inference is quantified using autocorrelation. 
#' For each observation of a given node, the correlations between that node and
#' all of its possible neighbors in the network are calculated with and without
#' the inclusion of that observation. 
#' This yields two vectors of correlation coefficients.
#' The correlation between these vectors, or autocorrelation, reflects the 
#' impact of the observation on the global interaction profile of that node, 
#' where a low correlation is indicative of network inference that is strongly
#' dependent on the inclusion or exclusion of that single data point. 
#' This situation is reflective of a likely outlier that compromises the 
#' robustness of network inference. 
#' 
#' @param mat a numeric matrix, with nodes (e.g., analytes such as genes, 
#'   proteins, or metabolites) in columns, and samples in rows
#' @param min_pairs minimum number of paired, non-missing observations 
#'   to calculate a correlation coefficient; correlations between vectors with 
#'   fewer than this number of paired observations will be replaced with
#'   \code{NA}
#' @param method the correlation coefficient to be computed; one of 
#'   \code{"pearson"} (default), \code{"kendall"}, or \code{"spearman"}; 
#'   can be abbreviated  
#' 
#' @return a matrix with identical dimensions to the input matrix, containing
#'   the autocorrelation for each non-missing observation 
#'   
#' @export
detect_outliers = function(mat, min_pairs = 10,
                           method = c("pearson", "kendall", "spearman")) {
  method = match.arg(method)
  
  # first, calculate baseline correlations 
  cor1 = cor(mat, method = method, use = 'pairwise.complete.obs')
  ## censor correlations with too few pairs
  n1 = crossprod(!is.na(mat))
  cor1[n1 < min_pairs] = NA

  # systematically remove each point and calculate autocorrelation
  autocor = mat
  autocor[] = NA
  nodes = colnames(mat)
  pb = progress::progress_bar$new(
    format = "node :what [:bar] :percent eta: :eta",
    clear = F, total = ncol(mat), width = 80)
  for (i in seq_len(ncol(mat))) {
    vec = mat[, i]
    observations = which(!is.na(vec))
    autocors = purrr::map_dbl(observations, function(index) {
      # remove observation 
      vec0 = vec
      vec0[index] = NA
      # recalculate correlations without observations 
      cor2 = cor(mat[, -i], vec0, method = method, 
                 use = 'pairwise.complete.obs')
      # calculate autocorrelation
      autocor = cor(cor2, cor1[-i, i], use = 'pairwise.complete.obs')
    })
    # insert autocorrelations into matrix
    autocor[observations, i] = autocors
    # tick progress bar
    pb$tick(tokens = list(what = sprintf(
      paste0("%-", nchar(ncol(mat)), "s"), i)))
  }
  
  return(autocor)
}
