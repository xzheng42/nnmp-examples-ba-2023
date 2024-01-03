copNNMP <- function(yy, XX = NULL, coords, ne_info, marg_family, cop_family, 
                    sp_func, priors, starting, tuning, mcmc_settings, verbose) {
  
  #--------------------------------------------------------------
  # Convert cop_family to integer for faster comparison in loops
  #--------------------------------------------------------------
  if (cop_family == "gaussian") {
    
    cop_fam <- 1
    
  } else if (cop_family == "gumbel") {
    
    cop_fam <- 2
    
  } else if (cop_family == "clayton") {
    
    cop_fam <- 3
    
  } 
  
  #----------------------------------------------------
  # Fit models
  #----------------------------------------------------
  if (marg_family == "beta") 
  {
    marg_fam <- 1
    copBetaNNMP(yy, XX, coords, ne_info, marg_fam, cop_fam, 
                sp_func, priors, starting, tuning, mcmc_settings, verbose)
  }
  else if (marg_family == "gamma") 
  {
    marg_fam <- 2
    copGammaNNMP(yy, XX, coords, ne_info, marg_fam, cop_fam, 
                 sp_func, priors, starting, tuning, mcmc_settings, verbose)
  }
  else {
    stop("error: this family is not avaiable.") 
  }
                        
}


# Update phi of the copula NNMP
updateCopPhi <- function(phi, se_phi, u_phi, v_phi, 
                         dat, dat_ne, dat_rho, data_label, dist_mat, rho_mat,
                         mar_param, mar_family, cop_family, sp_func) {
  

  nn <- length(dat)
  prop_log_phi <- rnorm(1, log(phi), se_phi)
  prop_phi <- exp(prop_log_phi)
  prop_rho_mat <- sp_func(dist_mat, prop_phi)
  prop_dat_rho <- as.numeric(t(sapply(1:nn, function(x) prop_rho_mat[x, data_label[x]])))


  prop_loglik <-
    dgamma(1 / prop_phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / prop_phi) +
    sum(dBiCopMar_cpp2(dat, dat_ne, cop_family, prop_dat_rho, mar_family, mar_param, FALSE, TRUE)) + prop_log_phi
  cur_loglik <-
    dgamma(1 / phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / phi) +
    sum(dBiCopMar_cpp2(dat, dat_ne, cop_family, dat_rho, mar_family, mar_param, FALSE, TRUE)) + log(phi)

  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(phi = prop_phi, rho_mat = prop_rho_mat, dat_rho = prop_dat_rho, accept = 1))
  } else {
    return(list(phi = phi, rho_mat = rho_mat, dat_rho = dat_rho, accept = 0))
  }

}


# Update zeta of the copula NNMP 
updateCopZeta <- function(zeta, se_zeta, nne, dist_mat, mu_t, kasq, data_label, u_zeta, v_zeta){
  
  trun_size <- length(mu_t)
  prop_log_zeta <- rnorm(1, log(zeta), se_zeta)
  prop_zeta <- exp(prop_log_zeta)
  prop_weight_res <- logitGausWeight(nne, dist_mat, prop_zeta, mu_t, rep(sqrt(kasq), trun_size), trun = TRUE)
  prop_weight_mat <- prop_weight_res$weights
  prop_cutoff <- prop_weight_res$cutoff
  weight_res <- logitGausWeight(nne, dist_mat, zeta, mu_t, rep(sqrt(kasq), trun_size), trun = TRUE)
  weight_mat <- weight_res$weights
  cutoff <- weight_res$cutoff
  
  prop_loglik <- 
    dgamma(1 / prop_zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / prop_zeta) + 
    sum(log(sapply(1:trun_size, function(x) prop_weight_mat[x, data_label[x]]))) + prop_log_zeta
  cur_loglik <- 
    dgamma(1 / zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / zeta) + 
    sum(log(sapply(1:trun_size, function(x) weight_mat[x, data_label[x]]))) + log(zeta)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(zeta = prop_zeta, weight_mat = prop_weight_mat, cutoff = prop_cutoff, accept = 1))
  } else {
    return(list(zeta = zeta, weight_mat = weight_mat, cutoff = cutoff, accept = 0))
  }
  
}


# Update gamma of the copula NNMP 
updateCopGa <- function(latent_t, DD, kasq, mu_ga, V_ga) {
  V_ga_inv <- chol2inv(chol(V_ga))
  V_1 <- chol2inv(chol(t(DD) %*% DD / kasq + V_ga_inv))
  mu_ga_1 <- V_1 %*% (t(DD) %*% latent_t / kasq + V_ga_inv %*% mu_ga)
  ga <- t(mvtnorm::rmvnorm(1, mu_ga_1, V_1))
  return(ga)
}


# Update kappa2 of the copula NNMP 
updateCopKasq <- function(latent_t, mu_t, u_kasq, v_kasq) {
  kk <- length(mu_t)
  uu <- u_kasq + kk / 2
  vv <- v_kasq + .5 * sum((latent_t - mu_t)^2)
  kasq <- 1 / rgamma(1, uu, vv)
  return(kasq)
}
