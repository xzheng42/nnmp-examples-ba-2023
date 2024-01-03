#' @title Function for fitting a nearest-neighbor mixture process model
#' 
#' @description This function fits a nearest-neighbor mixture process (NNMP) model with continuous 
#' marginal distributions using Markov Chain Monte Carlo simulations.
#' 
#' @param response a vector of complete geostatistical data of length \eqn{n}.
#' @param covars an \eqn{n \times p} matrix of covariates. If an intercept is desired, the first 
#'               column of the matrix should be a vector of ones.
#' @param coords an \eqn{n \times 2} matrix of the observation coordinates, e.g., longitude and latitude.
#' @param neighbor_size neighborhood size (number of nearest neighbors).
#' @param marg_family a quoted keyword that specifies the family of distributions for the marginal.
#'                    Supported keywords are \code{"beta"}, \code{"gamma"}, \code{"gaussian"}, and \code{"sn"}.
#' @param cop_family a quoted keyword that specifies the family of copulas for the bivariate distributions
#'                   that define the first-order conditionals. 
#'                   Supported keywords are \code{"gaussian"}, \code{"gumbel"}, and \code{"clayton"}.
#' @param sp_func a function that introduces spatial dependence for the mixture components. 
#'                The default function is an exponential correlation function now (or its transformation for 
#'                Gumbel and Clayton copulas); users can skip this argument.
#' @param priors a list of priors. Each element of the list corresponds to a combination of an unknown parameter 
#'               and its prior distribution; for example, "\code{phi_invgamma}" and "\code{zeta_invgamma}".
#' @param starting a list of starting values. Each element of the list corresponds to an unknown parameter.
#' @param tuning a list of tuning parameters. Each element of the list corresponds to an unknown parameter.
#'               The value portion of each element defines the variance of the random walk Metropolis sampler 
#'               proposal distribution for the unknown parameter.
#' @param ord an index vector of length \eqn{n} to order the data. The vector defines an ordering based on which
#'            nearest neighbors are searched. If the argument is not specified, the default ordering is random ordering.
#' @param mcmc_settings a list of Markov chain Monte Carlo (MCMC) simulation parameters. 
#'        \code{n_iter}: number of iterations; \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If \code{verbose = TRUE},
#'        \code{n_report} defines the interval to report MCMC progress and Metropolis acceptance rates.
#' @param weights_samples logical; if true, samples of the mixture weights are returned.
#'                        This argument is fixed to be FALSE now; users can skip this argument.
#' @param verbose logical; if true, model specification, MCMC progress, and Metropolis acceptance rates
#'                are printed to the screen.
#' @param neighbor_info logical; if true, a list that comprises each location's neighborhood information is returned.
#' @param model_diag a list of quoted kyewords that specify measurements of model fit to calculate.
#'.            Supported key words are \code{dic} and \code{pplc}. The argument only supports the GNNMP now.
#' @param par_label a vector of numeric values of length \eqn{n}, ranging from 1 to K, 
#'                  where k corresponds to the k-th partition, for k = 1, ..., K.
#'                  The argument is only specified for fitting the extended skew-Gaussian NNMP.
#' 
#' @return
#' An object of class \code{nnmp}. The return object is a list that comprises:
#' 
#' \tabular{ll}{
#' 
#' \code{post_samples} \tab a list of posterior samples. Each element of the list corresponds to an unknown parameter. \cr
#' \tab \cr
#' \code{orig_data} \tab a list that consists of input data 
#'         (\code{response}, \code{covars}, and \code{coords}). \cr
#' \tab \cr
#' \code{ord_data} \tab a list that consists of ordered input data
#'         (\code{ord:} index vector for ordering; 
#'          \code{yy_ord}: ordered responses;
#'          \code{XX_ord}: ordered covariates (if \code{covars = NULL}, \code{XX_ord = NULL});
#'          \code{coords_ord}: ordered coordinates). \cr
#' \tab \cr
#' \code{mod_spec}: \tab a list that consists of model specifications 
#'         (\code{neigbhor_size}: neighborhood size; 
#'          \code{marg_family}: family of marginal distributions;
#'          \code{cop_family}: copula family for the bivariate distributions that define the first-order conditionals;
#'          \code{sp_func}: function that introduces spatial dependence;
#'          \code{priors}: priors for unknown parameters). \cr
#' \tab \cr
#' \code{mcmc_params}: \tab a list that consists of unknown parameters' starting values and Metropolis sampler proposal distribution variances
#'         used in the MCMC simulation (
#'         \code{starting}: starting values;
#'         \code{tuning}: Metropolis sampler proposal distribution variances). \cr
#' \tab \cr
#' \code{mh_data}: \tab a list that consists of each unknown parameter's Metropolis 
#'                      sampler acceptance and rejection indicators.\cr
#' \tab \cr
#' \code{runtime}: \tab running time for MCMC simulation calculated using \code{system.time}.
#'         Running time for searching nearest neighbors is not included. \cr
#' \tab \cr
#' \code{neighbor_info}: \tab returned if \code{neighbor_info = TRUE}. See above the \code{neighbor_info} argument description.\cr
#' \tab \cr
#' \code{par_label}: \tab returned if \code{par_label = TRUE}. See above the \code{par_label} argument description.\cr
#' \tab \cr
#' \code{mod_diag}: \tab returned if \code{mod_diag = TRUE}. See above the \code{mod_diag} argument description.
#' 
#'}
#'
#' 
#' @author Xiaotian Zheng \email{xiaotianzheng42@gmail.com}
#' 
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. 
#' "Nearest-neighbor mixture models for non-Gaussian spatial processes". 
#' Bayesian Anal. 18 (4) 1191 - 1222, December 2023.
#' DOI: \href{https://doi.org/10.1214/23-BA1405}{https://doi.org/10.1214/23-BA1405}.
#' 
#' @export
#'
nnmp <- function(response,
                 covars = NULL,
                 coords,
                 neighbor_size = 10, 
                 marg_family, 
                 cop_family = NULL,
                 sp_func = NULL,
                 priors = NULL,
                 starting = NULL,
                 tuning = NULL,
                 ord = NULL,
                 mcmc_settings,
                 weights_samples = FALSE,
                 verbose = TRUE,
                 neighbor_info = FALSE,
                 model_diag = NULL,
                 par_label = NULL) {
  
  #----------------------------------------------------
  # nearest neighbor size
  #----------------------------------------------------
  nne <- neighbor_size
  if (nne %% 1 != 0) {
    stop("error: neighbor size must be an integer.")
  }
  
  if (nne < 2) {
    stop("error: neighbor size must be greater than one")
  }
  
  
  #----------------------------------------------------
  # family
  #----------------------------------------------------
  if (missing(marg_family)) {
    stop("error: family must be specified")
  }
  
  marg_family <- tolower(marg_family)
  marg_family_names <- c("beta", "gamma", "gaussian", "sn")
  
  if (!marg_family %in% marg_family_names) {
    stop("error: specified family, '", marg_family, "' for the marginal distribution is not a valid option;
         available families are ", paste(marg_family_names, collapse = ", ", sep = "") ,".")
  }  
  
  
  #----------------------------------------------------
  # DAG
  #----------------------------------------------------
  yy <- as.numeric(response)
  coords <- as.matrix(coords)

  orig_data <- list(response = yy, 
                    covars = covars, 
                    coords = coords)  
  
  if (is.null(ord)) {
    ord <- sample(seq_along(response), replace = FALSE)
  }

  if (is.null(covars)) {
    yy_ord <- yy[ord]
    XX_ord <- NULL
    coords_ord <- coords[ord, ]
    ord_data <- list(ord = ord,
                     yy_ord = yy_ord,
                     XX_ord = XX_ord,
                     coords_ord = coords_ord)
  } else {
    yy_ord <- yy[ord]
    XX_ord <- covars[ord, ]
    coords_ord <- coords[ord, ]
    ord_data <- list(ord = ord,
                     yy_ord = yy_ord,
                     XX_ord = XX_ord,
                     coords_ord = coords_ord) 
  }
  
  if (!is.null(par_label)) {
    if (marg_family == "sn") {
      par_label_ord <- par_label[ord]
      ord_data$par_label_ord <- par_label_ord
    } else {
      warning("warning: partition labels are specified, but they are not needed to fit this family of NNMPs.")
    }
  }

  ne_info <- neighbor(nne, coords_ord, yy_ord, ref = TRUE)
  ne_info$ne_size <- nne
  
  
  #----------------------------------------------------
  # function that introduces spatial dependence
  #----------------------------------------------------
  if (is.null(sp_func)) {
    
    if (is.null(cop_family)) {
      
      sp_func <- function(dist, range_param) return(exp(-dist / range_param))
      
    } 
    else if (cop_family == "gaussian") {
      
      sp_func <- function(dist, range_param) return(exp(-dist / range_param))
      
    }
    else if (cop_family == "gumbel") {
      
      sp_func <- function(dist, range_param) 1 / (1 - exp(-dist / range_param))
      
    }
    else if (cop_family == "clayton") {
      
      sp_func <- function(dist, range_param) {
        kk <- exp(-dist / range_param)
        2 * kk / (1 - kk)
      }
      
    }
    
  }
  
  
  #----------------------------------------------------
  # Check some priors
  #----------------------------------------------------
  if (is.null(priors$phi_invgamma)) {
    
    priors$phi_invgamma <- c(3, 1/3)
    
  }
  if (is.null(priors$zeta_invgamma)) {
    
    priors$zeta_invgamma <- c(3, 0.2)
    
  }
  
  if (is.null(priors$ga_gaus)) {
    
    priors$ga_gaus <- list(mean_vec = as.matrix(c(-1.5, 0, 0)),
                           var_mat = 2 * diag(3))
    
  } 
  if (is.null(priors$kasq_invgamma)) {
    
    priors$kasq_invgamma <- c(3, 1)
    
  }
  
  
  #----------------------------------------------------
  # Check some starting values
  #----------------------------------------------------
  if (is.null(starting$ga)) {
    
    starting$ga <- matrix(c(-1.5, 0, 0))
    
  }
  if (is.null(starting$kasq)) {
    
    starting$kasq <- 1
    
  }
  if (is.null(starting$phi)) {
    
    starting$phi <- 1/6
    
  }
  if (is.null(starting$zeta)) {
    
    starting$zeta <- 0.1
    
  }
  
  
  #----------------------------------------------------
  # verbose
  #----------------------------------------------------
  if (verbose) {
    cat("--------------------------------------------------------------\n")
    cat(" << Fitting a nearest-neighbor mxiture process model >>\n")
    cat("--------------------------------------------------------------\n")
    cat("----------------------------------------\n")
    cat("\t  NNMP settings\n");
    cat("----------------------------------------\n")
    
    cat(paste0("Reference set: ", length(response), " observations\n"))
    if (is.null(ord)) {
      cat(paste0("Using ", nne, " nearest neighbors with random ordering\n"))
    } else {
      cat(paste0("Using ", nne, " nearest neighbors with user-specified ordering\n"))
    }
    
    cat(paste0("Marginal distribution family: ", tools::toTitleCase(marg_family), "\n"))
    
    if (!is.null(cop_family)) {
      cat(paste0("Mixture component family: ", tools::toTitleCase(cop_family), " copula\n"))
    } else {
      if (marg_family == "gaussian") {
        cat(paste0("Mixture component family: Bivariate Gaussian\n"))
      } 
      else if (marg_family == "sn") {
        cat(paste0("Mixture component family: Bivariate skew-Gaussian\n"))
      }
    }
    
    if (is.null(sp_func)) {
      cat("Note: model/function for spatial dependence is not specified. Use default functions.\n")
    }
    
    cat("--------------------------------------------\n")
    cat("\t  MCMC settings\n");
    cat("--------------------------------------------\n")
    cat(paste0("Number of iterations: ", mcmc_settings$n_iter, "\n"))
    cat(paste0("Number of burn-in samples: ", mcmc_settings$n_burn, "\n"))
    cat(paste0("Thinning degree: ", mcmc_settings$n_thin, "\n"))
    
  }

  
  #----------------------------------------------------
  # MCMC
  #----------------------------------------------------
  if (marg_family == "gaussian") {

        mcmc_out <- gausNNMP(yy_ord, XX_ord, coords_ord, ne_info, sp_func, priors,
                             starting, tuning, mcmc_settings, verbose)
    
  } 
  else if (marg_family == "sn") {
    
    mcmc_out <- snNNMP(yy_ord, XX_ord, coords_ord, ne_info, sp_func, 
                       priors, starting, tuning, mcmc_settings, verbose, par_label_ord)
    
  }
  else if (marg_family %in% c("gamma", "beta")) {
    
    if (is.null(cop_family)) {
      
      stop(paste0("error: the current package implements ", tools::toTitleCase(marg_family), 
                  " NNMP with Gaussian, Gumbel, or Clayton copulas. Need to specify a copula family"))
      
    } else {
      
      if (cop_family %in% c("gaussian", "gumbel", "clayton")) {
        
        mcmc_out <- copNNMP(yy_ord, XX_ord, coords_ord, ne_info, marg_family, cop_family, 
                            sp_func, priors, starting, tuning, mcmc_settings, verbose)
      
      } else {
        
        stop(paste0("error: the current package implements only ", tools::toTitleCase(marg_family), 
                    " NNMP with Gaussian, Gumbel, or Clayton copulas."))
        
      }
    
    } 
      
  } else {
    
    stop("error: the current package supports families of marginal distributions including
         Gaussian, skew-Gaussian, beta, and gamma.")
  }
  
  post_sams <- mcmc_out$post_sams
  runtime <- mcmc_out$runtime
  mh_data <- mcmc_out$mh_acct_save
  
  #----------------------------------------------------
  # Output
  #----------------------------------------------------  
  if (is.null(cop_family)) {
    cop_family <- NULL
  }
  if (is.null(sp_func)) {
    sp_func <- NULL
  }
  if (is.null(tuning)) {
    tuning <- NULL
  }

  mod_spec <- list(neighbor_size = neighbor_size, 
                   marg_family = marg_family,
                   cop_family = cop_family,
                   sp_func = sp_func,
                   priors = priors)

  mcmc_params <- list(starting = starting, tuning = tuning)

  nnmp_out <- list(
    post_samples = post_sams,
    orig_data = orig_data,
    ord_data = ord_data,
    mod_spec = mod_spec,
    mcmc_params = mcmc_params,
    mcmc_settings = mcmc_settings,
    mh_data = mh_data,
    runtime = runtime
  )
  
  if (neighbor_info) {
    nnmp_out$neighbor_info = ne_info
  }
  
  if (!is.null(par_label)) {
    nnmp_out$par_label <- par_label
  }
  
  
  #----------------------------------------------------
  # Model diagnostic
  #----------------------------------------------------
  if (!is.null(model_diag)) {
    
    # if (verbose) {
    #   cat("--------------------------------------------\n")
    #   cat("\t  Model fit diagnostics\n")
    #   cat("--------------------------------------------\n")
    # }
    
    nnmp_out$mod_diag <- nnmpDiag(nnmp_out, model_diag, yy_ord, XX_ord)
  
  }
  
  class(nnmp_out) <- "nnmp"
  
  nnmp_out
      
  
}
