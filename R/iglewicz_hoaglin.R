#' Detect outliers using the modified Z score method of Iglewicz and Hoaglin
#' 
#' Identify outliers within a distribution of numeric values using
#' the modified Z score method. Iglewicz and Hoaglin recommend an absolute
#' Z score threshold of 3.5 to identify potential outliers.
#' 
#' Full details are provided in: 
#' 
#' Boris Iglewicz and David Hoaglin (1993), 
#' "Volume 16: How to Detect and Handle Outliers", 
#' The ASQC Basic References in Quality Control: Statistical Techniques, 
#' Edward F. Mykytka, ed.
#' 
#' @param x distribution to find outliers in
#' @param threshold absolute value of the modified z score threshold above which
#'   to consider a value an outlier; defaults to \code{3.5} on the 
#'   recommendation of Iglewicz and Hoaglin
#' @param return_scores optionally, return the modified z score of each 
#'   observation instead of a masked version of the input vector
#' 
#' @return a vector of the same length as the input, with outliers masked 
#'   (or, if \code{return_scores} is true, the modified z scores of each 
#'   observation)
#' 
#' @export
iglewicz_hoaglin = function(x, threshold = 3.5, return_scores = F) {
  # check input
  if (!is.numeric(x)) 
    stop("could not identify outliers: input is not numeric")  
  # calculate modified Z scores
  med = median(x, na.rm = T)
  MAD = median(abs(x - med), na.rm = T)
  Mi = 0.6745 * (x - med) / MAD
  if (return_scores) {
    return(Mi)
  } else {
    # mask vector
    x[abs(Mi) > threshold] = NA
    return(x)
  }
}
