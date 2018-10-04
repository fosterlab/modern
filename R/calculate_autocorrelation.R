#' Calculate autocorrelation after removal of each observation 
#' 
#' Conventionally, autocorrelation refers to the correlation between
#' a signal and a delayed copy of that signal. 
#' In \code{modern}, autocorrelation is calculated as the correlation 
#' between the global interaction profile of a given node, represented 
#' as a vector of correlation coefficients to all other nodes in the network,
#' before and after the removal of a given data point.
#' The autocorrelation is calculated for each data point in turn, and
#' observations that have a large impact on the global interaction profile are
#' flagged as outliers. 
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
calculate_autocorrelation = function(mat, min_pairs = 20,
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
