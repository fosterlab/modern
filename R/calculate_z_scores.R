#' Calculate z scores for outliers
#' 
#' Convert a matrix of autocorrelations to a matrix of 
#' Z scores, such that the matrix has a mean of zero and a standard deviation
#' of one. If the matrix contains missing values, this scaling is performed for 
#' each group of columns with equivalent numbers of missing values separately.
#' Optionally, if there are many possible numbers of missing values, the z score
#' can be calculated for approximately equal sized bins of missing value counts
#' using the \code{bins} parameter.
#' 
#' @param autocor the matrix of autocorrelations
#' @param bins optionally, the number of bins into which to group nodes on the
#'   basis of the number of observations
#' @return a matrix with identical dimensions to the input matrix, containing
#'   the autocorrelation Z score assigned to each non-missing observation 
#'   
#' @importFrom purrr %>%
#' @importFrom dplyr n_distinct
#' @importFrom mltools bin_data
#' @importFrom reshape2 melt acast
#' 
#' @export
calculate_z_scores = function(autocor, bins = 1) {
  # melt to a long data frame
  long = autocor %>% 
    reshape2::melt(varnames = c("sample", "node")) %>%
    dplyr::group_by(node) %>%
    dplyr::mutate(n_obs = sum(!is.na(value))) %>%
    ungroup()
  
  # if there are more bins than missing value counts, bin missing values
  if (!is.na(bins) && !is.null(bins) && n_distinct(long$n_obs) > bins) {
    long %<>% mutate(group = bin_data(
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
