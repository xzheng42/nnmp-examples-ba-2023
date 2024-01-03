copGammaNNMP <- function(yy, XX, coords, ne_info, marg_family, cop_family, 
                         sp_func, priors, starting, tuning, mcmc_settings, verbose) {
  
  
  #-----------------------------------
  # Check priors except for marginals
  #-----------------------------------
  
  # Priors
  u_phi <- priors$phi_invgamma[1]
  v_phi <- priors$phi_invgamma[2]
  u_zeta <- priors$zeta_invgamma[1]
  v_zeta <- priors$zeta_invgamma[2]
  
  mu_ga <- priors$ga_gaus$mean_vec
  V_ga <- priors$ga_gaus$var_mat
  u_kasq <- priors$kasq_invgamma[1]
  v_kasq <- priors$kasq_invgamma[2]
  
  #----------------------------------------------
  # Check tuning parameters except for marginals
  #----------------------------------------------
  se_phi <- tuning$phi
  se_zeta <- tuning$zeta
  
  #----------------------------------------------
  # Check starting values except for marginals
  #----------------------------------------------
  phi <- starting$phi
  zeta <- starting$zeta
  ga <- starting$ga
  kasq <- starting$kasq
  
  #------------------------------------------------------------
  # Check priors, tuining and starting for marginals, and MCMC
  #------------------------------------------------------------
  
  if (is.null(XX)) {
    
    u_alp <- priors$shape1_gamma[1]
    v_alp <- priors$shape1_gamma[2]
    u_be <- priors$shape2_gamma[1]
    v_be <- priors$shape2_gamma[2]
    
    se_alp <- tuning$shape1
    se_be <- tuning$shape2
    
    alp <- starting$shape1
    be <- starting$shape2
    
    runtime <- system.time(
      mcmc_out <- mcmcCopGammaNNMP_simple(yy, coords, ne_info, 
                                          marg_family, cop_family, sp_func, 
                                          mcmc_settings, verbose,
                                           
                                          u_alp, v_alp, u_be, v_be,
                                          u_phi, v_phi, u_zeta, v_zeta,
                                          mu_ga, V_ga, u_kasq,v_kasq,
                                          
                                          se_alp, se_be, se_phi, se_zeta,
                                          
                                          alp, be, phi, zeta, ga, kasq
      )
    )
  }
  
  mcmc_out$runtime <- runtime
  
  mcmc_out
  
  
}




mcmcCopGammaNNMP_simple <- function(yy, coords, ne_info, 
                                    marg_family, cop_family, sp_func, 
                                    mcmc_settings, verbose,
                                   
                                    u_alp, v_alp, u_be, v_be,
                                    u_phi, v_phi, u_zeta, v_zeta,
                                    mu_ga, V_ga, u_kasq,v_kasq,
                                   
                                    se_alp, se_be, se_phi, se_zeta,
                                   
                                    alp, be, phi, zeta, ga, kasq
                                    ) {
  
  
  
  #--------------------------
  # MCMC settings
  #--------------------------
  if (verbose) {
    
    if (is.null(mcmc_settings$n_report)) {
      nreport <- 1000
    } else {
      nreport <- mcmc_settings$n_report
    }
    
    cat("--------------------------------------------\n")
    cat("\t  Running MCMC\n");
    cat("--------------------------------------------\n")
    
  }
  niter <- mcmc_settings$n_iter
  nburn <- mcmc_settings$n_burn
  nthin <- mcmc_settings$n_thin
  
  #--------------------------
  # Initialization
  #--------------------------
  nne <- ne_info$ne_size
  tyy <- yy[-(1:nne)]
  tnn <- length(tyy)
  dist_mat <- ne_info$ne_dist[-(1:nne), ]
  ne_idx_mat <- ne_info$ne_index[-(1:nne), ]
  tyy_ne_mat <- t(sapply(1:tnn, function(x) yy[ne_idx_mat[x,]]))
  DD <- as.matrix(cbind(1, coords))[-(1:nne),]
  mu_t <- DD %*% ga

  weight_res <- logitGausWeight(nne, dist_mat, zeta, mu_t, rep(sqrt(kasq), tnn), TRUE)
  weight_mat <- weight_res$weights
  cutoff <- weight_res$cutoff
  valid_idx <- !is.na(weight_mat)
  if (any(weight_mat[valid_idx]==0)){
    weight_mat[weight_mat==0] <- .Machine$double.eps
  }

  tyy_label <- array(NA, dim = tnn)
  for (i in 1:tnn) {
    tyy_label[i] <- sample(1:nne, size = 1, prob = weight_mat[i,])
  }
  rho_mat <- sp_func(dist_mat, phi)
  tyy_rho <- as.numeric(t(sapply(1:tnn, function(x) rho_mat[x, tyy_label[x]])))
  ne_idx <- as.numeric(t(sapply(1:tnn, function(x) ne_idx_mat[x, tyy_label[x]])))
  tyy_ne <- yy[ne_idx]
  
  # empty stuff to save samples
  # weight_save <- array(NA, dim = c(tnn, nne, niter - nburn))
  alp_save <- array(NA, dim = niter - nburn)
  be_save <- array(NA, dim = niter - nburn)
  phi_save <- array(NA, dim = niter - nburn)
  zeta_save <- array(NA, dim = niter - nburn)
  ga_save <- array(NA, dim = c(ncol(DD), niter - nburn))
  kasq_save <- array(NA, dim = niter - nburn)
  alp_acct_save <- array(0, dim = niter)
  be_acct_save <- array(0, dim = niter)
  phi_acct_save <- array(0, dim = niter)
  zeta_acct_save <- array(0, dim = niter)
  
  
  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0
  
  for (iter in 1:niter) {
    
    start_time <- Sys.time()
    #--------------------------
    # Update alpha
    #--------------------------
    alp_res <- updateCopGammaAlp(alp, se_alp, u_alp, v_alp, 
                                 tyy, tyy_ne, tyy_rho, be, 
                                 marg_family, cop_family)
    alp <- alp_res$alp
    if (alp_res$accept) alp_acct_save[iter] <- 1
    
    #--------------------------
    # Update beta
    #--------------------------
    be_res <- updateCopGammaBe(be, se_be, u_be, v_be, 
                               tyy, tyy_ne, tyy_rho, alp, 
                               marg_family, cop_family)
    be <- be_res$be
    if (be_res$accept) be_acct_save[iter] <- 1
    
    #--------------------------
    # Update phi
    #--------------------------
    mar_param <- cbind(rep(alp, tnn), rep(be, tnn), 
                       rep(alp, tnn), rep(be, tnn))
    phi_res <- updateCopPhi(phi, se_phi, u_phi, v_phi, 
                            tyy, tyy_ne, tyy_rho, tyy_label, dist_mat, rho_mat,
                            mar_param, marg_family, cop_family, sp_func)
    phi <- phi_res$phi
    rho_mat <- phi_res$rho_mat
    tyy_rho <- phi_res$tyy_rho
    if (phi_res$accept) phi_acct_save[iter] <- 1
    
    #--------------------------
    # Update labels
    #--------------------------
    yy_param1 <- rep(alp, tnn)
    yy_param2 <- rep(be, tnn)
    yy_ne_param1 <- matrix(rep(yy_param1, nne), ncol = nne)
    yy_ne_param2 <- matrix(rep(yy_param2, nne), ncol = nne)
    labels <- updateCopLabel(tyy, tyy_ne_mat, nne, weight_mat, rho_mat,
                             cop_family, marg_family, 
                             yy_param1, yy_param2, yy_ne_param1, yy_ne_param2,
                             mu_t, rep(sqrt(kasq), length(mu_t)), cutoff)
    tyy_label <- as.numeric(labels$data_label)
    latent_t <- labels$latent_t
    tyy_rho <- as.numeric(t(sapply(1:tnn, function(x) rho_mat[x, tyy_label[x]])))
    ne_idx <- as.numeric(t(sapply(1:tnn, function(x) ne_idx_mat[x, tyy_label[x]])))
    tyy_ne <- yy[ne_idx]
    
    #--------------------------
    # Update ga
    #--------------------------
    ga <- updateCopGa(latent_t, DD, kasq, mu_ga, V_ga)
    mu_t <- as.numeric(DD %*% ga)
    
    #--------------------------
    # Update kasq
    #--------------------------
    kasq <- updateCopKasq(latent_t, mu_t, u_kasq, v_kasq)
    
    #--------------------------
    # Update zeta and weights
    #--------------------------
    zeta_res <- updateCopZeta(zeta, se_zeta, nne, dist_mat, mu_t, kasq, tyy_label, u_zeta, v_zeta)
    zeta <- zeta_res$zeta
    weight_mat <- zeta_res$weight_mat
    cutoff <- zeta_res$cutoff
    if (zeta_res$accept) zeta_acct_save[iter] <- 1
    
    valid_idx <- !is.na(weight_mat)
    if (any(weight_mat[valid_idx]==0)) weight_mat[weight_mat==0] <- .Machine$double.eps

    #--------------------------
    # Print MCMC progress
    #--------------------------
    end_time <- Sys.time()
    block_runtime <- block_runtime + as.numeric(end_time - start_time)
    
    if (verbose) {
      if (iter %% nreport == 0) {
        cat(paste0("Iterations: ", iter, "/", niter, 
                   "  Percentage: ", specRound(iter / niter * 100, 2), "%\n"))
        ert <- (block_runtime / nreport) * (niter - iter)
        cat(paste0("Estimated remaining time: ", 
                   specRound(ert / 60, 2), " minutes \n"))
        block_runtime <- 0
        cat(paste0("Metropolis Hastings acceptance rates: \n"))
        cat(paste0("  Model parameter           Acceptance rate\n"))
        cat(paste0("  Shape 1                   ", 
                   specRound(sum(alp_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat(paste0("  Shape 2                   ", 
                   specRound(sum(be_acct_save[1:iter]) / iter * 100, 2), "%\n"))        
        cat(paste0("  phi                       ", 
                   specRound(sum(phi_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat(paste0("  zeta                      ", 
                   specRound(sum(zeta_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat("--------------------------------------------\n")
      }
    }  
    
    #--------------------------
    # Save samples
    #--------------------------
    if (iter > nburn){
      # weight_save[, , iter - nburn] <- weight_mat
      alp_save[iter - nburn] <- alp
      be_save[iter - nburn] <- be
      phi_save[iter - nburn] <- phi
      zeta_save[iter - nburn] <- zeta
      ga_save[, iter - nburn] <- ga
      kasq_save[iter - nburn] <- kasq
    }
    
  }
  
  #--------------------------
  # Thinning
  #--------------------------
  selc_index <- seq(1, niter - nburn, by = nthin)
  
  post_sames <- list(#weight = weight_save[, , selc_index],
    shape1 = alp_save[selc_index],
    shape2 = be_save[selc_index],
    phi = phi_save[selc_index],
    zeta = zeta_save[selc_index],
    ga = ga_save[,selc_index],
    kasq = kasq_save[selc_index])
  
  mh_acct_save <- list(shape1_mh = alp_acct_save,
                       shape2_mh = alp_acct_save,
                       phi_mh = phi_acct_save,
                       zeta_mh = zeta_acct_save)
  
  list(post_sams = post_sames, mh_acct_save = mh_acct_save)
  
}


# Update alpha of the copula NNMP with beta marginals
updateCopGammaAlp <- function(alp, se_alp, u_alp, v_alp, 
                              tyy, tyy_ne, tyy_rho, be, 
                              marg_family, cop_family) {
  
  tnn <- length(tyy)
  prop_log_alp <- rnorm(1, log(alp), se_alp)
  prop_alp <- exp(prop_log_alp)
  prop_mar_param <- cbind(rep(prop_alp, tnn), rep(be, tnn),
                          rep(prop_alp, tnn), rep(be, tnn))
  cur_mar_param <- cbind(rep(alp, tnn), rep(be, tnn),
                          rep(alp, tnn), rep(be, tnn))
  
  prop_loglik <- 
    dgamma(prop_alp, u_alp, v_alp, log = TRUE) + 
    sum(dBiCopMar_cpp2(tyy, tyy_ne, cop_family, tyy_rho, marg_family, prop_mar_param, FALSE, TRUE)) + 
    sum(dgamma(tyy, prop_alp, be, log = TRUE)) + prop_log_alp
  cur_loglik <-  
    dgamma(alp, u_alp, v_alp, log = TRUE) + 
    sum(dBiCopMar_cpp2(tyy, tyy_ne, cop_family, tyy_rho, marg_family, cur_mar_param, FALSE, TRUE)) + 
    sum(dgamma(tyy, alp, be, log = TRUE)) + log(alp)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(alp = prop_alp, accept = TRUE))
  } else {
    return(list(alp = alp, accept = FALSE))
  }
  
}

# Update beta of the copula NNMP with beta marginals
updateCopGammaBe <- function(be, se_be, u_be, v_be, 
                             tyy, tyy_ne, tyy_rho, alp, 
                             marg_family, cop_family) {
  
  tnn <- length(tyy)
  prop_log_be <- rnorm(1, log(be), se_be)
  prop_be <- exp(prop_log_be)
  prop_mar_param <- cbind(rep(alp, tnn), rep(prop_be, tnn),
                          rep(alp, tnn), rep(prop_be, tnn))
  cur_mar_param <- cbind(rep(alp, tnn), rep(be, tnn),
                         rep(alp, tnn), rep(be, tnn))
  
  prop_loglik <- 
    dgamma(prop_be, u_be, v_be, log = TRUE) + 
    sum(dBiCopMar_cpp2(tyy, tyy_ne, cop_family, tyy_rho, marg_family, prop_mar_param, FALSE, TRUE)) + 
    sum(dgamma(tyy, alp, prop_be, log = TRUE)) + prop_log_be
  cur_loglik <-  
    dgamma(be, u_be, v_be, log = TRUE) + 
    sum(dBiCopMar_cpp2(tyy, tyy_ne, cop_family, tyy_rho, marg_family, cur_mar_param, FALSE, TRUE)) + 
    sum(dgamma(tyy, alp, be, log = TRUE)) + log(be)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(be = prop_be, accept = TRUE))
  } else {
    return(list(be = be, accept = FALSE))
  }
  
}
