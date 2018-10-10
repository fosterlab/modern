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
#' 
#' @return a vector of the same length as the input, consisting of modified 
#' Z scores
#' 
#' @export
iglewicz_hoaglin = function(x) {
  # check input
  if (!is.numeric(x)) 
    stop("could not identify outliers: input is not numeric")  
  # calculate modified Z scores
  med = median(x, na.rm = T)
  MAD = median(abs(x - med), na.rm = T)
  Mi = 0.6745 * (x - med) / MAD
}
