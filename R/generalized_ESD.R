#' Detect outliers using the generalized Extreme Studentized Deviate test
#' 
#' Identify outliers within a distribution of numeric values using
#' Rosner's generalized Extreme Studentized Deviate (ESD) test. Unlike the
#' Grubbs or Tietjen-Moore tests, the generalized ESD test does not require
#' the number of outliers to be prespecified. 
#' 
#' 
#' Full details are provided in: 
#' 
#' Rosner, Bernard (May 1983), Percentage Points for a Generalized ESD 
#' Many-Outlier Procedure, \emph{Technometrics}, 25(2), pp. 165-172.
#' 
#' @param x distribution to find outliers in
#' @param max_outliers the prespecified maximum number of outliers; defaults
#'   to \code{5}
#' @param alpha type I error associated with the hypothesis test; 
#'   defaults to \code{0.05}
#' @param return_statistics optionally return the test statistics associated
#'   with the top \code{max_outliers} outliers
#' 
#' @return a vector of the same length as the input, with outliers masked
#'   (or, if \code{return_scores} is true, the test statistics of the top
#'   \code{max_outliers} observations)
#' 
#' @importFrom purrr map_dbl
#' 
#' @export
generalized_ESD = function(x, max_outliers = 5, alpha = 0.05, 
                           return_statistics = T) {
  # check input
  if (!is.numeric(x)) 
    stop("could not identify outliers: input is not numeric")  
  
  # calculate test statistics
  n_outliers = seq_len(max_outliers)
  statistics = rep(NA, max_outliers)
  outliers = rep(NA, max_outliers)
  r = abs(x - mean(x, na.rm = T)) / sd(x, na.rm = T)
  idxs = order(r, decreasing = T)
  x0 = x
  for (i in n_outliers) {
    idx = idxs[1]
    outliers[i] = idx
    statistics[i] = r[idx]
    x0[idx] = NA
    r = abs(x0 - mean(x0, na.rm = T)) / sd(x0, na.rm = T)
    idxs = order(r, decreasing = T)
  }
  
  # calculate critical values for each possible # of outliers
  lambdas = purrr::map_dbl(n_outliers, ~ {
    i = .
    n = sum(!is.na(x))
    df = n - i - 1
    p = 1 - alpha / (2 * (n - i + 1))
    t = suppressWarnings(qt(p, df))
    t * (n - i) / sqrt((n - i - 1 + t**2) * (n - i + 1))
  })
  
  # get largest i with test statistic greater than critical value
  i = suppressWarnings(max(which(statistics > lambdas), na.rm = T))
  
  if (is.finite(i)) {
    if (return_statistics) {
      x0 = rep(NA, length(x))
      x0[outliers[seq_len(i)]] = statistics[seq_len(i)]
      return(x0)
    } else {
      # return masked vector
      x0 = x
      remove = outliers[seq_len(i)]
      x0[remove] = NA
      return(x0)
    }
  } else {
    if (return_statistics) {
      return(rep(NA, length(x)))  
    } else {
      return(x)
    }
  }
}
