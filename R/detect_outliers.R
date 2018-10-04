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
#' 
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
#' The matrix of autocorrelations is subsequently converted to a matrix of 
#' Z scores, such that the matrix has a mean of zero and a standard deviation
#' of one. If the matrix contains missing values, this scaling is performed for 
#' each group of columns with equivalent numbers of missing values separately.
#' Optionally, if there are many possible numbers of missing values, the z score
#' can be calculated for approximately equal sized bins of missing value counts
#' using the \code{bins} parameter.
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
#' @param bins optionally, the number of bins into which to group nodes on the
#'   basis of the number of observations
#' @return a matrix with identical dimensions to the input matrix, containing
#'   the autocorrelation Z score assigned to each non-missing observation 
#'   
#' @importFrom purrr %>%
#' @importFrom mltools bin_data
#' @export
detect_outliers = function(mat, min_pairs = 10,
                           method = c("pearson", "kendall", "spearman"),
                           bins = NA) {
  method = match.arg(method)
  
  # calculate autocorrelations
  autocor = calculate_autocorrelation(mat, min_pairs, method)

  # melt to a long data frame
  long = autocor %>% 
    reshape2::melt(varnames = c("sample", "node")) %>%
    dplyr::group_by(node) %>%
    dplyr::mutate(n_obs = sum(!is.na(value))) %>%
    ungroup()
  
  # if there are more bins than missing value counts, bin missing values
  if (!is.na(bins) & dplyr::n_distinct(long$n_obs) > bins) {
    long %<>% mutate(group = mltools::bin_data(
      n_obs, bins = bins, binType = 'quantile'))
  } else {
    long %<>% mutate(group = n_obs)
  }
  
  # convert to Z scores, grouping nodes by # of missing values
  z = long %>% 
    group_by(group) %>%
    dplyr::mutate(z = (value - mean(value, na.rm = T)) / 
                    sd(value, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, node, z) %>%
    dplyr::rename(value = z) %>%
    reshape2::acast(sample ~ node)
  
  return(z)
}
