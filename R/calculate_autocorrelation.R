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
#' @param n_splits optionally, for very large matrices, split the input into 
#'   subsets of nodes in order to preserve RAM 
#' 
#' @return a matrix with identical dimensions to the input matrix, containing
#'   the autocorrelation for each non-missing observation 
#'   
#' @export
calculate_autocorrelation = function(
  mat,
  min_pairs = 10,
  method = c("pearson", "kendall", "spearman"),
  n_splits = NULL
) {
  method = match.arg(method)
  
  # optionally, split large matrices to fit correlations into memory
  if (is.null(n_splits) || as.integer(n_splits) == 1) {
    n_splits = 1
    matrices = list(mat)
  } else {
    if (is.na(as.integer(n_splits))) {
      stop("n_splits is not an integer: ", n_splits)
    }
    n_splits = as.integer(n_splits)
    nodes = colnames(mat)
    splits = split(nodes, cut(seq_along(nodes), n_splits, labels = FALSE)) 
    matrices = purrr::map(splits, ~ mat[, .])
  }
  
  # create autocorrelation matrix
  autocor = mat
  autocor[] = NA
  
  # process each matrix in sequence
  for (idx in seq_along(matrices)) {
    mat0 = matrices[[idx]]
    if (n_splits > 1) {
      message("working on split ", idx, " of ", n_splits, " ...")
    }
    
    # first, calculate baseline correlations 
    cor1 = cor(mat0, method = method, use = 'pairwise.complete.obs')
    ## censor correlations with too few pairs
    n1 = crossprod(!is.na(mat0))
    cor1[n1 < min_pairs] = NA
    
    # systematically remove each point and calculate autocorrelation
    nodes = colnames(mat0)
    pb = progress::progress_bar$new(
      format = "node :what [:bar] :percent eta: :eta",
      clear = F, total = ncol(mat0), width = 80)
    for (node_idx in seq_len(ncol(mat0))) {
      node = colnames(mat0)[node_idx]
      vec = mat0[, node]
      observations = which(!is.na(vec))
      autocors = purrr::map_dbl(observations, function(index) {
        # remove observation 
        vec0 = vec
        vec0[index] = NA
        # recalculate correlations without observations 
        cor2 = cor(mat0[, -node_idx], vec0, method = method, 
                   use = 'pairwise.complete.obs')
        # calculate autocorrelation
        autocor = cor(cor2, cor1[-node_idx, node_index], 
                      use = 'pairwise.complete.obs')
      })
      # insert autocorrelations into matrix
      autocor[observations, node] = autocors
      # tick progress bar
      pb$tick(tokens = list(what = sprintf(
        paste0("%-", nchar(ncol(mat0)), "s"), node_idx)))
    }
  }
  
  return(autocor)
}
