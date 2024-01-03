#' @title Obtain posterior estimates 
#'
#' @description This function calculates point and interval estimates of unknown parameters
#' using their posterior samples
#'
#' @param x a vector or matrix of posterior samples. If \code{x} is a matrix, 
#'          each row of the matrix corresponds to posterior samples of an unknown parameter.
#' @param pt a quoted keyword that specifies point estimate.
#'           Supported keywords are \code{"mean"} and \code{"median"}.
#' @param ci a numeric vector of lower and upper probabilities with values in \eqn{[0,1]} 
#'           to generate a credible interval for the unknown parameter(s).
#' @param k an integer indicating the number of decimal places to be used.
#'
#' @return A vector of characters. 
#'         Each element of the vector is in the form: point estimate (lower percentile, upper percentile).
#'         
#' @export
#' 
summarize <- function(x, pt, ci = c(0.025, 0.975), k = 2) {
  UseMethod("summarize", x)
}

#' @title Obtain posterior estimates 
#'
#' @description This function calculates point and interval estimates of unknown parameters
#' using their posterior samples
#'
#' @param x a vector of posterior samples of an unknown parameter.
#' @param pt a quoted keyword that specifies point estimate.
#'           Supported keywords are \code{"mean"} and \code{"median"}.
#' @param ci a numeric vector of lower and upper probabilities with values in \eqn{[0,1]} 
#'           to generate a credible interval for the unknown parameter.
#' @param k an integer indicating the number of decimal places to be used.
#'
#' @return A character in the form: point estimate (lower percentile, upper percentile).
#' 
#' @export
#' 
summarize.numeric <- function(x, 
                              pt = "mean", 
                              ci = c(0.025, 0.975), 
                              k = 2) {
	
  if (pt == "mean"){

    pe <- specRound(mean(x), k)

  } else if (pt == "median"){

    pe <- specRound(median(x), k)

  } else {

    stop("error: point estimate should be mean or median.")

  }

  ci <- specCI(x, ci, k, char = TRUE)
  paste(pe, ci)

}

#' @title Obtain posterior estimates 
#'
#' @description This function calculates point and interval estimates of unknown parameters
#' using their posterior samples
#' 
#' @param x a matrix in which each row corresponds to posterior samples of an unknown parameter.
#' @param pt a quoted keyword that specifies point estimate.
#'           Supported keywords are \code{"mean"} and \code{"median"}.
#' @param ci a numeric vector of lower and upper probabilities with values in \eqn{[0,1]} 
#'           to generate a credible interval for the unknown parameters.
#' @param k an integer indicating the number of decimal places to be used.
#' 
#' @return A vector of characters. 
#'         Each element of the vector is in the form: point estimate (lower percentile, upper percentile). 
#'         
#' @export
#' 
summarize.matrix <- function(x, 
                             pt = "mean", 
                             ci = c(0.025, 0.975), 
                             k = 2) {

  if (pt == "mean"){

    pe <- specRound(rowMeans(x), k)

  } else if (pt == "median"){

    pe <- specRound(apply(x, 1, median), k)

  } else {

    stop("error: point estimate should be mean or median.")

  }

  ci <- apply(x, 1, specCI, ci, k, char = TRUE)
  paste(pe, ci)

}

specRound <- function(x, k = 3) trimws(format(round(x, k), nsmall = k))

specCI <- function(x, ci = c(0.025, 0.975), k = 2, char = TRUE) {

  q <- quantile(x, probs = ci)

  if (char) {

    ci <- paste0('(', specRound(q[1],k), ', ', specRound(q[2],k), ')')

  } else {

    ci <- c(round(q[1], k), round(q[2], k))

  }
  
  ci

}

