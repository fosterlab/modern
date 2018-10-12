#' Calculate the outlier z score corresponding to a given family-wise error rate
#' 
#' Use Bonferroni correction to calculate the Z score that corresponds to a 
#' family-wide error rate at a given alpha. 
#' 
#' @param mat a numeric matrix, with nodes (e.g., analytes such as genes, 
#'   proteins, or metabolites) in columns, and samples in rows
#' @param alpha type I error associated with the hypothesis test;
#'   defaults to \code{0.05} 
#' 
#' @return the threshold below which z scores can be considered outliers
#'   
#' @export
calculate_fwer = function(mat, alpha = 0.05) {
  qnorm(alpha / sum(!is.na(mat)))
}
