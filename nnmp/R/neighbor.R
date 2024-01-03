#' @title Find nearest neighbors 
#' 
#' @description This function finds nearest neighbors for locations in the reference set or in
#' the non-reference set. Neighbors are found in the reference set.
#' 
#' @param neighbor_size neighborhood size (number of nearest-neighbors, say \eqn{L})
#' @param coords an \eqn{n \times 2} matrix of reference set coordinates, e.g., longitude and latitude.
#' @param obs a vector of observations of length \eqn{n} associated with the matrix of reference set coordinates.
#' @param ref if \code{TRUE}, nearest neighbors of reference set are searched. If \code{FALSE}, 
#'            argument \code{nonref} must be specified.
#' @param nonref an \eqn{r \times 2} matrix of non-reference set coordinates, e.g., longitude and latitude.
#'        This argument only needs to be specified when \code{ref = FALSE}.
#' @export
#'
#' @return
#' The function returns a list that comprises:
#' \tabular{ll}{
#' 
#' \code{ne_dist} \tab an \eqn{n \times L} matrix where the \eqn{(i,j)}-th element of the matrix
#' is the Euclidean distance between the \eqn{i}-th reference set location and its \eqn{j}-th nearest neighbor. \cr
#' \tab \cr
#' \code{ne_idx} \tab an \eqn{n \times L} matrix where the \eqn{(i,j)}-th element of the matrix
#' is the row index of the \eqn{i}-th reference set location's \eqn{j}-th nearest neighbor in \code{coords}. \cr
#' \tab \cr
#' \code{ne_obs} \tab an \eqn{n \times L} matrix where the \eqn{(i,j)}-th element of the matrix is the observation
#' associated with the \eqn{i}-th reference set location's \eqn{j}-th nearest neighbor in \code{coords}. \cr
#'}
#'
#' @author Xiaotian Zheng \email{xiaotianzheng42@gmail.com}
#' 
neighbor <- function(neighbor_size, coords, obs, ref = TRUE, nonref = NULL) {

  #----------------------------------------------------
  # nearest neighbor size
  #----------------------------------------------------
  nne <- neighbor_size
  if (nne %% 1 != 0) {
    stop("error: neighbor size must be an integer.")
  }
  
  if (nne < 2) {
    stop("error: neighbor size must be greater than one.")
  }  

  #----------------------------------------------------
  # Checks
  #----------------------------------------------------
  coords <- as.matrix(coords)
  if (nrow(coords) != length(obs)) {
    stop("error: the length of observations is different from 
          the number of rows of reference set coordinates.")
  }
  
  #----------------------------------------------------
  # find nearest neighbors
  #----------------------------------------------------  
  if (ref) {
    
    return(refNe(nne, coords, obs))
    
  } else {

    if (is.null(nonref)) {
      stop("error: when argument ref = FALSE, non-reference set coordinates must be specified.")
    }
    
    return(gridNe(nne, coords, nonref, obs))
    
  }

}

logitGausWeight <- function(nne, ne_dist, zeta, mu, ka, trun) {
  if (trun) {
    return(logitGausWeight2_cpp(nne, ne_dist, zeta, mu, ka))
  } else {
    return(logitGausWeight1_cpp(nne, ne_dist, zeta, mu, ka))
  }
}