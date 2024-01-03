#' @title Function for prediction at new locations
#'
#' @description This function produces posterior predictive samples for a set of new locations using an NNMP.
#'
#' @param object an object of class nnmp.
#' @param nonref_covars an \eqn{m \times p} matrix of covariates for \eqn{m} new locations.
#' @param nonref_coords an \eqn{m \times 2} matrix of the corresponding spatial coordinates.
#' @param probs a vector of probabilities with values in \eqn{[0,1]} to generate the corresponding quantiles.
#' @param predict_sam logical; if true, posterior predictive samples for each new location are returned.
#' @param nonref_par_labels a vector of numeric values of length \eqn{m}, indicating which partition
#'                          each new location belongs to.
#' @param compute_control This argument will be used in the future to specific parameters for
#'                        generating samples from a copula (e.g., Gumbel) when numerical method is needed.
#'                        Currently, the desired accuracy (convergence tolerance) is 1e-5, and
#'                        the maximum number of iterations is 1000.
#' @param verbose logical; if true, progress of the prediction is printed to the screen.
#' @param nreport If \code{verbose = TRUE}, \code{n_report} defines the interval to report
#'                the progress of prediction.
#' @param ... additional arguments. No additional arguments are supported currently.
#' 
#' @return
#' The return object is a list that comprises:
#' 
#' @author Xiaotian Zheng \email{xiaotianzheng42@gmail.com}
#' 
#' @references 
#' Zheng, X., Kottas, A., and Sansó, B. 
#' "Nearest-neighbor mixture models for non-Gaussian spatial processes". 
#' Bayesian Anal. 18 (4) 1191 - 1222, December 2023.
#' DOI: \href{https://doi.org/10.1214/23-BA1405}{https://doi.org/10.1214/23-BA1405}.
#' 
#' @exportS3Method 
#'
predict.nnmp <- function(object, 
                         nonref_covars = NULL, 
                         nonref_coords, 
                         probs = c(0.025, 0.5, 0.975), 
                         predict_sam = FALSE,
                         nonref_par_labels = NULL, 
                         compute_control = NULL,
                         verbose = TRUE, 
                         nreport = NULL, 
                         ...) {
  
  if (missing(object)) {
    stop("error: need nnmp object for prediction\n")
    }
  if (!inherits(object, "nnmp")) {
    stop("error: require nnmp object for prediction\n")
  }
  
  
  nne <- object$mod_spec$neighbor_size
  ref_coords <- object$ord_data$coords_ord
  ref_yy <- object$ord_data$yy_ord
  
  
  nonref_ne_info <- neighbor(nne, ref_coords, ref_yy, ref = FALSE, nonref_coords)
  nonref_ne_idx <- nonref_ne_info$ne_index
  nonref_ne_dist <- nonref_ne_info$ne_dist
  
  post_sams <- object$post_samples
  
  DD <- cbind(1, nonref_coords)
  
  
  marg_family <- object$mod_spec$marg_family
  cop_family <- object$mod_spec$cop_family
  if (!is.null(cop_family)) {
    if (cop_family == "gaussian") {
      cop_fam <- 1
    } else if (cop_family == "gumbel") {
      cop_fam <- 2
    } else if (cop_family == "clayton") {
      cop_fam <- 3
    } 
  }
  
  
  if (!is.null(nonref_covars)) {
    XX <- as.matrix(nonref_covars)
  }
  
  
  #----------------------------------------------------
  # verbose
  #----------------------------------------------------
  if (verbose) {
    
    if (is.null(nreport)) {
      nreport <- floor(nrow(nonref_coords) / 10)
    }
    
    nsams <- length(post_sams$zeta)
    
    cat("--------------------------------------------------------------\n")
    cat(" << Nearest-neighbor mixture process model prediction >>\n")
    cat("--------------------------------------------------------------\n")
    cat(paste0("Reference set: ", length(ref_yy), " observations\n"))
    cat(paste0("Prediction over non-reference set: ", nrow(nonref_coords), " locations\n"))
    
    cat("----------------------------------------\n")
    cat("\tNNMP model-fit settings\n");
    cat("----------------------------------------\n")
    cat(paste0("Using ", nne, " nearest neighbors\n")) 
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
    cat(paste0("Number of MCMC samples: ", nsams, "\n"))
    
  }
  
  
  if (marg_family == "gaussian") {
    
    runtime <- system.time(
      pred_out <- predGausNNMP(nne, XX, DD,
                               post_sams$bb, post_sams$zz,
                               post_sams$sigmasq, post_sams$tausq,
                               post_sams$phi, post_sams$zeta,
                               post_sams$ga, post_sams$kasq,
                               nonref_ne_idx, nonref_ne_dist,
                               probs, verbose, nreport, sam = predict_sam)
    )

  }
  else if (marg_family == "sn") {
    
    if (is.null(object$orig_data$covars)) {
      
      runtime <- system.time(
        pred_out <- predSkewGausNNMP_simple(ref_yy, nne, DD, post_sams$la, post_sams$sigmasq, 
                                            post_sams$phi, post_sams$zeta,
                                            post_sams$ga, post_sams$kasq, 
                                            nonref_ne_idx, nonref_ne_dist, probs, 
                                            verbose, nreport, sam = predict_sam)
      )
      
    } else {
      
      if (is.null(nonref_par_labels)) {
        stop("error: partition labels for the non-reference set are not specified.")
      }
      
      ref_XX <- object$ord_data$XX_ord
      ref_par_labels <- object$ord_data$par_label_ord
      
      runtime <- system.time(
        pred_out <- predSkewGausNNMP_covars_lapar(ref_yy, nne, ref_XX, XX, 
                                                  ref_par_labels, nonref_par_labels,
                                                  DD, post_sams$regcoef, post_sams$la,
                                                  post_sams$sigmasq, post_sams$phi,
                                                  post_sams$zeta, post_sams$ga,
                                                  post_sams$kasq, nonref_ne_idx,
                                                  nonref_ne_dist, probs,
                                                  verbose, nreport, sam = predict_sam)
      )
      
    }
    
  }
  else if (marg_family == "beta") {
    
    marg_fam <- 1
    
    if (is.null(object$orig_data$covars)) {
      
      cop_param_list <- sapply(1:length(post_sams$phi), function(x) 
        object$mod_spec$sp_func(nonref_ne_dist, post_sams$phi[x]), simplify = FALSE)
      cop_param <- simplify2array(cop_param_list)
      
      runtime <- system.time(
        pred_out <- predCopNNMP_simple(ref_yy, cop_fam, cop_param, marg_fam,
                                       cbind(post_sams$shape1, post_sams$shape2),
                                       nne, DD, post_sams$zeta, post_sams$ga,
                                       post_sams$kasq, nonref_ne_idx, nonref_ne_dist,
                                       probs, verbose, nreport, sam = predict_sam)
      )
      
    } else {
      
      cop_param_list <- sapply(1:length(post_sams$phi), function(x) 
        object$mod_spec$sp_func(nonref_ne_dist, post_sams$phi[x]), simplify = FALSE)
      cop_param <- simplify2array(cop_param_list)
      
      ref_XX <- object$ord_data$XX_ord
      nonref_XX <- nonref_covars
      
      ref_yy_mu <- apply(post_sams$regcoef, 2, function(x) plogis(ref_XX %*% x))
      ref_yy_shape1 <- t(apply(ref_yy_mu, 1, function(x) x * post_sams$scale))
      ref_yy_shape2 <- t(apply(ref_yy_mu, 1, function(x) (1 - x) * post_sams$scale))
      
      nonref_yy_mu <- apply(post_sams$regcoef, 2, function(x) plogis(nonref_XX %*% x))
      nonref_yy_shape1 <- t(apply(nonref_yy_mu, 1, function(x) x * post_sams$scale))
      nonref_yy_shape2 <- t(apply(nonref_yy_mu, 1, function(x) (1 - x) * post_sams$scale))
      
      runtime <- system.time(
        pred_out <- predCopNNMP_covars(ref_yy, cop_fam, cop_param, marg_fam, 
                                       ref_yy_shape1, ref_yy_shape2, nonref_yy_shape1, nonref_yy_shape2,
                                       nne, DD, post_sams$zeta, post_sams$ga,
                                       post_sams$kasq, nonref_ne_idx, nonref_ne_dist,
                                       probs, verbose, nreport, sam = predict_sam)
      )
      
    }
    
  }
  else if (marg_family == "gamma") {
    
    marg_fam <- 2
    
    if (is.null(object$orig_data$covars)) {
      cop_param_list <- sapply(1:length(post_sams$phi), function(x) 
        object$mod_spec$sp_func(nonref_ne_dist, post_sams$phi[x]), simplify = FALSE)
      cop_param <- simplify2array(cop_param_list)
      
      runtime <- system.time(
        pred_out <- predCopNNMP_simple(ref_yy, cop_fam, cop_param, marg_fam,
                                       cbind(post_sams$shape1, post_sams$shape2),
                                       nne, DD, post_sams$zeta, post_sams$ga,
                                       post_sams$kasq, nonref_ne_idx, nonref_ne_dist,
                                       probs, verbose, nreport, sam = predict_sam)
      )
      
    }
    
  }
  else {
    
    stop("error: this family is not supported currently.")
    
  }
  
  cat("----------------------------------------\n")
  cat(paste0("Prediction running time: ", specRound(runtime[3] / 60, 2), " minutes\n\n"))

  pred_out
  
}


nnmpDiag <- function(object, model_diag, yy_ord, XX_ord) {
  
  post_sams <- object$post_samples
  marg_family <- object$mod_spec$marg_family
  
  model_diag <- tolower(model_diag)
  
  nnmp_diag <- vector("list")
  
  if (marg_family == "gaussian") {
    if ("dic" %in% model_diag) {
      dic <- dicGausNNMP(yy_ord, XX_ord, post_sams$bb, post_sams$zz, post_sams$tausq)
      nnmp_diag$dic <- dic
    }
    
    if ("pplc" %in% model_diag) {
      pplc <- pplcGausNNMP(yy_ord, XX_ord, post_sams$bb, post_sams$zz, post_sams$tausq)
      pplc <- as.numeric(pplc)
      names(pplc) <- c("G", "P", "D")
      nnmp_diag$pplc <- pplc
    }
  }
  
  nnmp_diag
  
}

#' @title Function for the conditional cumulative distribution function at new locations
#' 
#' @description This function produces posterior samples of the conditional cumulative distribution functions of
#' an NNMP for a set of new locations. The function currently supports only stationary beta/gamma NNMPs.
#' 
#' @param object an object of class nnmp.
#' @param nonref_coords an \eqn{m \times 2} matrix of spatial coordinates.
#' @param probs a vector of marginal probabilities corresponding to the NNMP response.
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
pnnmp <- function(object, nonref_coords, probs) {

  nne <- object$mod_spec$neighbor_size
  ref_coords <- object$orig_data$coords
  ref_yy <- object$orig_data$response

  nonref_ne_info <- neighbor(nne, ref_coords, ref_yy, ref = FALSE, nonref_coords)
  nonref_ne_obs <- nonref_ne_info$ne_obs
  nonref_ne_idx <- nonref_ne_info$ne_index
  nonref_ne_dist <- nonref_ne_info$ne_dist
  
  DD <- cbind(1, nonref_coords)
  marg_family <- object$mod_spec$marg_family
  
  cop_family <- object$mod_spec$cop_family
  if (!is.null(cop_family)) {
    if (cop_family == "gaussian") {
      cop_fam <- 1
    } else if (cop_family == "gumbel") {
      cop_fam <- 2
    } else if (cop_family == "clayton") {
      cop_fam <- 3
    }
  }
  
  post_sams <- object$post_samples
  cop_param_list <- sapply(1:length(post_sams$phi), function(x) 
    object$mod_spec$sp_func(nonref_ne_dist, post_sams$phi[x]), simplify = FALSE)
  cop_param <- simplify2array(cop_param_list)
  
  if (!is.null(cop_family)) {
    if (marg_family == "beta") {
      marg_fam <- 1
      cdf_out <- cdfCopNNMP_simple(nonref_ne_obs, nonref_ne_dist,
                                   cop_fam, cop_param, marg_fam,
                                   cbind(post_sams$shape1, post_sams$shape2),
                                   post_sams$zeta, post_sams$ga, post_sams$kasq, 
                                   DD, probs)
    } 
    else if (marg_family == "gamma") {
      marg_fam <- 2
      cdf_out <- cdfCopNNMP_simple(nonref_ne_obs, nonref_ne_dist,
                                   cop_fam, cop_param, marg_fam,
                                   cbind(post_sams$shape1, post_sams$shape2),
                                   post_sams$zeta, post_sams$ga, post_sams$kasq, 
                                   DD, probs)
    } else {
      stop("error: only support beta and gamma marginals now.")
    }
  } else {
    stop("error: only support stationary NNMPs with beta and gamma marginals now.")
  }
  
  cdf_out
  
}

#' summary.nnmp <- function(object, pt = "mean", ci = c(0.025, 0.975), 
#'                          digit = 2, char = TRUE, ...) {
#'   
#'   if (missing(object)) {
#'     stop("error: need nnmp object for summary\n")
#'   }
#'   if (!inherits(object, "nnmp")) {
#'     stop("error: require nnmp object for summary\n")
#'   }
#'   
#'   
#'   post_sams <- object$post_samples
#'   marg_family <- object$mod_spec$marg_family
#'   
#'   
#'   if (marg_family == "gaussian") {
#'     bb_est <- summarize(post_sams$bb, pt = pt, ci = ci, k = digit, char = char)
#'     phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit, char = char)
#'     zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit, char = char)
#'     sigmasq_est <- summarize(post_sams$sigmasq, pt = pt, ci = ci, k = digit, char = char)
#'     tausq_est <- summarize(post_sams$tausq, pt = pt, ci = ci, k = digit, char = char)
#'     ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit, char = char)
#'     kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit, char = char)
#'     tbl_est <- array(NA, dim = c(length(bb_est) + 8, 1))
#'     tbl_est[,1] <- c(bb_est, phi_est, zeta_est, sigmasq_est, tausq_est, ga_est, kasq_est)
#'     rownames(tbl_est) <- c(paste0("beta_", (1:length(bb_est))-1),
#'                            "phi", "zeta", "sigmasq", "tausq",
#'                            "gamma_0", "gamma_1", "gamma_2", "kappasq")
#'     colnames(tbl_est) <- "Gaussian NNMP"
#'   }
#'   
#'   
#'   if (marg_family == "beta") {
#'     if (is.null(object$orig_data$covars)) {
#'       alp_est <- summarize(post_sams$shape1, pt = pt, ci = ci, k = digit, char = char)
#'       be_est <- summarize(post_sams$shape2, pt = pt, ci = ci, k = digit, char = char)
#'       phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit, char = char)
#'       zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit, char = char)
#'       ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit, char = char)
#'       kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit, char = char)
#'       tbl_est <- array(NA, dim = c(8, 1))
#'       tbl_est[,1] <- c(alp_est, be_est, phi_est, zeta_est, ga_est, kasq_est)
#'       rownames(tbl_est) <- c("shape1", "shape2", "phi", "zeta", 
#'                              "gamma_0", "gamma_1", "gamma_2", "kappasq")
#'       colnames(tbl_est) <- "Beta NNMP"
#'     } else {
#'       bb_est <- summarize(post_sams$regcoef, pt = pt, ci = ci, k = digit, char = char)
#'       tau_est <- summarize(post_sams$scale, pt = pt, ci = ci, k = digit, char = char)
#'       phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit, char = char)
#'       zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit, char = char)
#'       ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit, char = char)
#'       kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit, char = char)
#'       tbl_est <- array(NA, dim = c(7 + length(bb_est), 1))
#'       tbl_est[,1] <- c(bb_est, tau_est, phi_est, zeta_est, ga_est, kasq_est)
#'       rownames(tbl_est) <- c(paste0("beta_", 0:(length(bb_est) - 1)),
#'                              "scale", "phi", "zeta", 
#'                              "gamma_0", "gamma_1", "gamma_2", "kappasq")
#'       colnames(tbl_est) <- "Beta NNMP"
#'     }
#'   }
#'   
#'   
#'   if (marg_family == "gamma") {
#'     if (is.null(object$orig_data$covars)) {
#'       alp_est <- summarize(post_sams$shape1, pt = pt, ci = ci, k = digit, char = char)
#'       be_est <- summarize(post_sams$shape2, pt = pt, ci = ci, k = digit, char = char)
#'       phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit, char = char)
#'       zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit, char = char)
#'       ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit, char = char)
#'       kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit, char = char)
#'       tbl_est <- array(NA, dim = c(8, 1))
#'       tbl_est[,1] <- c(alp_est, be_est, phi_est, zeta_est, ga_est, kasq_est)
#'       rownames(tbl_est) <- c("shape1", "shape2", "phi", "zeta", 
#'                              "gamma_0", "gamma_1", "gamma_2", "kappasq")
#'       colnames(tbl_est) <- "Gamma NNMP"
#'     }
#'   }
#'   
#'   
#'   if (marg_family == "sn") {
#'     
#'     if (is.null(object$orig_data$covars)) {
#'       
#'       la_est <- summarize(post_sams$la, pt = pt, ci = ci, k = digit, char = char)
#'       sigmasq_est <- summarize(post_sams$sigmasq, pt = pt, ci = ci, k = digit, char = char)
#'       phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit, char = char)
#'       zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit, char = char)
#'       ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit, char = char)
#'       kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit, char = char)
#'       tbl_est <- array(NA, dim = c(8, 1))
#'       tbl_est[,1] <- c(la_est, sigmasq_est, phi_est, zeta_est, ga_est, kasq_est)
#'       rownames(tbl_est) <- c("lambda", "sigmasq", "phi", "zeta", 
#'                              "gamma_0", "gamma_1", "gamma_2", "kappasq")
#'       colnames(tbl_est) <- "Skew-Gaussian NNMP"
#'       
#'     } else {
#'       
#'       la_est <- summarize(post_sams$la, pt = pt, ci = ci, k = digit, char = char)
#'       bb_est <- summarize(post_sams$regcoef, pt = pt, ci = ci, k = digit, char = char)
#'       sigmasq_est <- summarize(post_sams$sigmasq, pt = pt, ci = ci, k = digit, char = char)
#'       phi_est <- summarize(post_sams$phi, pt = pt, ci = ci, k = digit, char = char)
#'       zeta_est <- summarize(post_sams$zeta, pt = pt, ci = ci, k = digit, char = char)
#'       ga_est <- summarize(post_sams$ga, pt = pt, ci = ci, k = digit, char = char)
#'       kasq_est <- summarize(post_sams$kasq, pt = pt, ci = ci, k = digit, char = char)
#'       
#'       tbl_est <- array(NA, dim = c(length(la_est) + length(bb_est) + 7, 1))
#'       tbl_est[, 1] <- c(la_est, bb_est, sigmasq_est,
#'                         phi_est, zeta_est, ga_est, kasq_est)
#'       rownames(tbl_est) <- c(paste0("lambda_", 1:(length(la_est))), 
#'                              paste0("beta_", 0:(length(bb_est) - 1)), 
#'                              "sigmasq", "phi", "$zeta", 
#'                              "gamma_0", "gamma_1", "gamma_2", "kappasq")
#'       colnames(tbl_est) <- "Extended skew-Gaussian NNMP"
#'       
#'     }
#'   }
#'   
#'   tbl_est
#'   
#' }
