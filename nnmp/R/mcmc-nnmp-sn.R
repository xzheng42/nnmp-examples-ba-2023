snNNMP <- function(yy, XX, coords, ne_info, sp_func, 
                   priors, starting, tuning, mcmc_settings, verbose, par_label) {
  
  #--------------------------
  # Check priors 
  #--------------------------
  u_sigmasq <- priors$sigmasq_invgamma[1]
  v_sigmasq <- priors$sigmasq_invgamma[2]
  
  u_phi <- priors$phi_invgamma[1]
  v_phi <- priors$phi_invgamma[2]
  u_zeta <- priors$zeta_invgamma[1]
  v_zeta <- priors$zeta_invgamma[2]
  
  mu_ga <- priors$ga_gaus$mean_vec
  V_ga <- priors$ga_gaus$var_mat
  u_kasq <- priors$kasq_invgamma[1]
  v_kasq <- priors$kasq_invgamma[2]
  
  #--------------------------
  # Check tuning parameters
  #--------------------------
  se_phi <- tuning$phi
  se_zeta <- tuning$zeta
  se_sigmasq <- tuning$sigmasq
  
  #--------------------------
  # Check starting values
  #--------------------------
  sigmasq <- starting$sigmasq
  ga <- starting$ga
  kasq <- starting$kasq
  phi <- starting$phi
  zeta <- starting$zeta
  
  #------------------------------------------------------------
  # Check priors, tuning and starting for marginals, and MCMC
  #------------------------------------------------------------

  if (is.null(XX)) {

    mu_la <- priors$la_normal[1]
    sigmasq_la <- priors$la_normal[2]

    se_la <- tuning$la

    la <- starting$la
    

    logC1 <- function(yy, yy_ne, tyy_rho, la, sigmasq) {
      
      s2 <- sigmasq * (1 + tyy_rho)    
      tt1 <- la / sqrt(s2) / sqrt(2 * la^2 + s2)
      tt2 <- la / sqrt(sigmasq) / sqrt(la^2 + sigmasq)
      
      numer <- pnorm(tt1 * (yy + yy_ne), log.p = TRUE)
      denom <- pnorm(tt2 * yy_ne, log.p = TRUE)
      
      numer - denom
      
    }
    
    logf1 <- function(yy, yy_ne, tyy_rho, la, sigmasq) {
      
      omega2 <- sigmasq + la^2
      aa <- (tyy_rho * sigmasq + la^2) / omega2
      tomega2 <- omega2 * (1 - aa^2)
      
      log_dens <- dnorm(yy, aa * yy_ne, sqrt(tomega2), log = TRUE)
      
      log_dens
      
    }

    runtime <- system.time(
      mcmc_out <- mcmcSnNNMP_simple(yy, coords, ne_info, sp_func, 
                                    mcmc_settings, verbose,
                                    logC1, logf1,
                                    
                                    mu_la, sigmasq_la, 
                                    u_sigmasq, v_sigmasq, 
                                    u_phi, v_phi, u_zeta, v_zeta,
                                    mu_ga, V_ga, u_kasq,v_kasq,
                                    
                                    se_la, se_sigmasq, se_phi, se_zeta,
                                    
                                    la, sigmasq, phi, zeta, ga, kasq
      ) 
    )
    
  } 
  else 
  {
    
    if (is.null(par_label)) {
      stop("error: partition labels are not specified.")
    }    
    
    mu_la <- priors$la_normal$mu
    sigmasq_la <- priors$la_normal$sigmasq

    se_bb <- tuning$regcoef
    se_la <- tuning$la

    la <- starting$la
    bb <- as.matrix(starting$regcoef)
    
    
    logC2 <- function(yy, yy_ne, dat_rho, xi, xi_ne, la, la_ne, sigmasq) {
      s2 <- (1 + dat_rho) * sigmasq
      denom <- sqrt((1 - dat_rho) * s2) * sqrt((1 - dat_rho) * s2 + la^2 + la_ne^2 - 2 * dat_rho * la * la_ne)
      alp1 <- (la - dat_rho * la_ne) / denom
      alp2 <- (la_ne - dat_rho * la) / denom
      alp_ne <- la_ne / sqrt(sigmasq) 
      numer <- pnorm(alp1 * (yy - xi) + alp2 * (yy_ne - xi_ne), log.p = TRUE)
      denom <- pnorm(alp_ne * (yy_ne - xi_ne) / sqrt(la_ne^2 + sigmasq), log.p = TRUE)
      return(numer - denom)
    }
    
    logf2 <- function(yy, yy_ne, dat_rho, xi, xi_ne, la, la_ne, sigmasq) {
      trho <- dat_rho * sigmasq + la * la_ne
      omega2 <- sigmasq + la^2
      omega2_ne <- sigmasq + la_ne^2
      tga <- trho / omega2_ne
      tomega2 <- omega2 - trho^2 / omega2_ne
      log_dens <- dnorm(yy, xi + tga * (yy_ne - xi_ne), sqrt(tomega2), log = TRUE)
      return(log_dens)
    }
    

    runtime <- system.time(
      mcmc_out <- mcmcSnNNMP_covars_lapar(yy, XX, coords, ne_info, sp_func, 
                                          mcmc_settings, verbose, 
                                          logC2, logf2, par_label,
                                  
                                          mu_la, sigmasq_la, 
                                          u_sigmasq, v_sigmasq, 
                                          u_phi, v_phi, u_zeta, v_zeta,
                                          mu_ga, V_ga, u_kasq,v_kasq,
                                  
                                          se_bb, se_la, se_sigmasq, se_phi, se_zeta,
                                  
                                          bb, la, sigmasq, phi, zeta, ga, kasq) 
    )    
    

  }
  

  
  mcmc_out$runtime <- runtime
  
  mcmc_out
  
  
}

#-------------------------------------------------------------------
# MCMC for stationary skew-Gaussian NNMP
#-------------------------------------------------------------------
mcmcSnNNMP_simple <- function(yy, coords, ne_info, sp_func, 
                              mcmc_settings, verbose, logC, logf,
                              mu_la, sigmasq_la, 
                              u_sigmasq, v_sigmasq, 
                              u_phi, v_phi, u_zeta, v_zeta,
                              mu_ga, V_ga, u_kasq,v_kasq,
                              se_la, se_sigmasq, se_phi, se_zeta,
                              la, sigmasq, phi, zeta, ga, kasq
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
  nn <- length(yy)
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
  
  ## empty stuff to save samples
  # weight_save <- array(NA, dim = c(tnn, nne, niter - nburn))
  la_save <- array(NA, dim = niter - nburn)
  sigmasq_save <- array(NA, dim = niter - nburn)
  phi_save <- array(NA, dim = niter - nburn)
  zeta_save <- array(NA, dim = niter - nburn)
  ga_save <- array(NA, dim = c(ncol(DD), niter - nburn))
  kasq_save <- array(NA, dim = niter - nburn)
  la_acct_save <- array(0, dim = niter)
  sigmasq_acct_save <- array(0, dim = niter)
  phi_acct_save <- array(0, dim = niter)
  zeta_acct_save <- array(0, dim = niter)
  
  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0
  
  for (iter in 1:niter) {
    
    start_time <- Sys.time()
    #--------------------------
    # Update la
    #--------------------------
    la_res <- updateSNnnmpLa(la, se_la, mu_la, sigmasq_la, 
                             tyy, tyy_ne, tyy_rho, sigmasq, logC, logf)
    la <- la_res$la
    la_acct_save[iter] <- la_res$accept
    
    #--------------------------
    # Update sigmasq
    #--------------------------
    sigmasq_res <- updateSNnnmpSigmasq(sigmasq, se_sigmasq, u_sigmasq, v_sigmasq,
                                   tyy, tyy_ne, tyy_rho, la, logC, logf)
    sigmasq <- sigmasq_res$sigmasq
    sigmasq_acct_save[iter] <- sigmasq_res$accept
    
    #--------------------------
    # Update phi
    #--------------------------
    phi_res <- updateSNnnmpPhi(phi, se_phi, u_phi, v_phi,
                               tyy, tyy_ne, tyy_rho, tyy_label, dist_mat, rho_mat,
                               la, sigmasq, sp_func, logC, logf)
    phi <- phi_res$phi
    rho_mat <- phi_res$rho_mat
    tyy_rho <- phi_res$tyy_rho
    phi_acct_save[iter] <- phi_res$accept
    
    #--------------------------
    # Update labels
    #--------------------------
    labels <- updateSNnnmpLabel(tyy, tyy_ne_mat, nne, weight_mat, la, sigmasq, rho_mat,
                                mu_t, rep(sqrt(kasq), length(mu_t)), cutoff)
    tyy_label <- as.numeric(labels$data_label)
    latent_t <- labels$latent_t
    tyy_rho <- as.numeric(t(sapply(1:tnn, function(x) rho_mat[x, tyy_label[x]])))
    ne_idx <- as.numeric(t(sapply(1:tnn, function(x) ne_idx_mat[x, tyy_label[x]])))
    tyy_ne <- yy[ne_idx]
    
    #--------------------------
    # Update ga
    #--------------------------
    ga <- updateSNnnmpGa(latent_t, DD, kasq, mu_ga, V_ga)
    mu_t <- as.numeric(DD %*% ga)
    
    #--------------------------
    # Update kasq
    #--------------------------
    kasq <- updateSNnnmpKasq(latent_t, mu_t, u_kasq, v_kasq)
    
    #--------------------------
    # Update zeta and weights
    #--------------------------
    zeta_res <- updateSNnnmpZeta(zeta, se_zeta, nne, dist_mat, mu_t, kasq, tyy_label, u_zeta, v_zeta)
    zeta <- zeta_res$zeta
    weight_mat <- zeta_res$weight_mat
    cutoff <- zeta_res$cutoff
    zeta_acct_save[iter] <- zeta_res$accept
    
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
        cat(paste0("  la                        ", 
                   specRound(sum(la_acct_save[1:iter]) / iter * 100, 2), "%\n"))
        cat(paste0("  sigmasq                   ", 
                   specRound(sum(sigmasq_acct_save[1:iter]) / iter * 100, 2), "%\n"))        
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
      la_save[iter - nburn] <- la
      sigmasq_save[iter - nburn] <- sigmasq
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
  
  post_sams <- list(#weight = weight_save[, , selc_index],
    la = la_save[selc_index],
    sigmasq = sigmasq_save[selc_index],
    phi = phi_save[selc_index],
    zeta = zeta_save[selc_index],
    ga = ga_save[,selc_index],
    kasq = kasq_save[selc_index])
  
  mh_acct_save <- list(la_mh = la_acct_save,
                       sigmasq_mh = sigmasq_acct_save,
                       phi_mh = phi_acct_save,
                       zeta_mh = zeta_acct_save)
  
  list(post_sams = post_sams, mh_acct_save = mh_acct_save)
  
}

#-------------------------------------------------------------------
# MCMC for skew-Gaussian NNMP with continuous covariates
#-------------------------------------------------------------------
mcmcSnNNMP_covars_lapar <- function(yy, XX, coords, ne_info, sp_func, 
                                    mcmc_settings, verbose, 
                                    logC, logf, par_label,
                                    mu_la, sigmasq_la, 
                                    u_sigmasq, v_sigmasq, 
                                    u_phi, v_phi, u_zeta, v_zeta,
                                    mu_ga, V_ga, u_kasq,v_kasq,
                                    se_bb, se_la, se_sigmasq, se_phi, se_zeta,
                                    bb, la, sigmasq, phi, zeta, ga, kasq) {
  
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
  tyy_par_label <- par_label[-(1:nne)]
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
  yy_nu <- as.vector(XX %*% bb)
  tyy_nu <- yy_nu[-(1:nne)]
  tyy_nu_ne <- yy_nu[ne_idx]
  
  yy_la <- la[par_label]
  tyy_la <- yy_la[-(1:nne)]
  tyy_la_ne <- yy_la[ne_idx]
  tyy_ne_par_label <- par_label[ne_idx]
  
  
  ## empty stuff to save samples
  # weight_save <- array(NA, dim = c(tnn, nne, niter - nburn))
  bb_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  la_save <- array(NA, dim = c(max(tyy_ne_par_label), niter - nburn))
  sigmasq_save <- array(NA, dim = niter - nburn)
  phi_save <- array(NA, dim = niter - nburn)
  zeta_save <- array(NA, dim = niter - nburn)
  ga_save <- array(NA, dim = c(ncol(DD), niter - nburn))
  kasq_save <- array(NA, dim = niter - nburn)
  bb_acct_save <- array(0, dim = c(ncol(XX), niter))
  la_acct_save <- array(0, dim = c(max(tyy_ne_par_label), niter))
  sigmasq_acct_save <- array(0, dim = niter)
  phi_acct_save <- array(0, dim = niter)
  zeta_acct_save <- array(0, dim = niter)
  
  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0
  
  for (iter in 1:niter) {
    
    start_time <- Sys.time()
    #------------------------------------
    # Update egression coefficients beta
    #------------------------------------
    bb_res <- updateSNnnmpCoef_lapar(bb, se_bb, 
                                     tyy, tyy_ne, tyy_rho, XX, 
                                     nne, ne_idx, yy_nu, tyy_nu, tyy_nu_ne, 
                                     tyy_la, tyy_la_ne, sigmasq, logC, logf)
    bb <- bb_res$bb
    bb_acct_save[, iter] <- bb_res$accept
    yy_nu <- bb_res$yy_nu
    tyy_nu <- bb_res$tyy_nu
    tyy_nu_ne <- bb_res$tyy_nu_ne

    #--------------------------
    # Update la
    #--------------------------
    la_res <- updateSNnnmpLa_covars_lapar(la, se_la, mu_la, sigmasq_la, 
                                          tyy, tyy_ne, tyy_rho, tyy_la, tyy_la_ne,
                                          tyy_nu, tyy_nu_ne, sigmasq, logC, logf,
                                          tyy_par_label, tyy_ne_par_label) 
    la <- la_res$la
    la_acct_save[, iter] <- la_res$accept
    yy_la <- la[par_label]
    tyy_la <- la_res$tyy_la
    tyy_la_ne <- la_res$tyy_la_ne
    
    #--------------------------
    # Update sigmasq
    #--------------------------
    sigmasq_res <- updateSNnnmpSigmasq_covars_lapar(sigmasq, se_sigmasq, u_sigmasq, v_sigmasq,
                                                    tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne,
                                                    tyy_la, tyy_la_ne, logC, logf) 
    sigmasq <- sigmasq_res$sigmasq
    sigmasq_acct_save[iter] <- sigmasq_res$accept
    
    #--------------------------
    # Update phi
    #--------------------------
    phi_res <- updateSNnnmpPhi_covars_lapar(phi, se_phi, u_phi, v_phi,
                               tyy, tyy_ne, tyy_rho, tyy_label, dist_mat, rho_mat,
                               tyy_nu, tyy_nu_ne, tyy_la, tyy_la_ne,
                               sigmasq, sp_func, logC, logf)
    
    phi <- phi_res$phi
    rho_mat <- phi_res$rho_mat
    tyy_rho <- phi_res$tyy_rho
    phi_acct_save[iter] <- phi_res$accept
    
    #--------------------------
    # Update labels
    #--------------------------
    tyy_nu_ne_mat <- t(sapply(1:tnn, function(x) yy_nu[ne_idx_mat[x, ]]))
    tyy_la_ne_mat <- t(sapply(1:tnn, function(x) yy_la[ne_idx_mat[x, ]]))
    labels <- updateSNnnmpLabel_covars_lapar(tyy, tyy_ne_mat, nne, weight_mat, 
                                             tyy_nu, tyy_nu_ne_mat, tyy_la, tyy_la_ne_mat,
                                             sigmasq, rho_mat, 
                                             mu_t, rep(sqrt(kasq), length(mu_t)), cutoff)
    tyy_label <- as.numeric(labels$data_label)
    latent_t <- labels$latent_t
    tyy_rho <- as.numeric(t(sapply(1:tnn, function(x) rho_mat[x, tyy_label[x]])))
    ne_idx <- as.numeric(t(sapply(1:tnn, function(x) ne_idx_mat[x, tyy_label[x]])))
    tyy_ne <- yy[ne_idx]
    tyy_nu_ne <- yy_nu[ne_idx]
    tyy_la_ne <- yy_la[ne_idx]
    tyy_ne_par_label <- par_label[ne_idx]
    
    #--------------------------
    # Update ga
    #--------------------------
    ga <- updateSNnnmpGa(latent_t, DD, kasq, mu_ga, V_ga)
    mu_t <- as.numeric(DD %*% ga)
    
    #--------------------------
    # Update kasq
    #--------------------------
    kasq <- updateSNnnmpKasq(latent_t, mu_t, u_kasq, v_kasq)
    
    #--------------------------
    # Update zeta and weights
    #--------------------------
    zeta_res <- updateSNnnmpZeta(zeta, se_zeta, nne, dist_mat, mu_t, kasq, tyy_label, u_zeta, v_zeta)
    zeta <- zeta_res$zeta
    weight_mat <- zeta_res$weight_mat
    cutoff <- zeta_res$cutoff
    zeta_acct_save[iter] <- zeta_res$accept
    
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
        for (j in 1:ncol(XX)) {
          cat(paste0("  beta_", j-1, "                    ", 
                     specRound(sum(bb_acct_save[j,1:iter]) / iter * 100, 2), "%\n")) 
        }
        for (j in 1:max(tyy_ne_par_label)) {
          cat(paste0("  la_", j, "                      ", 
                     specRound(sum(la_acct_save[j,1:iter]) / iter * 100, 2), "%\n")) 
        }         
        cat(paste0("  sigmasq                   ", 
                   specRound(sum(sigmasq_acct_save[1:iter]) / iter * 100, 2), "%\n"))        
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
      bb_save[, iter - nburn] <- bb
      la_save[, iter - nburn] <- la
      sigmasq_save[iter - nburn] <- sigmasq
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
  
  post_sams <- list(#weight = weight_save[, , selc_index],
    regcoef = bb_save[, selc_index],
    la = la_save[, selc_index],
    sigmasq = sigmasq_save[selc_index],
    phi = phi_save[selc_index],
    zeta = zeta_save[selc_index],
    ga = ga_save[,selc_index],
    kasq = kasq_save[selc_index])
  
  mh_acct_save <- list(regcoef_mh = bb_acct_save,
                       la_mh = la_acct_save,
                       sigmasq_mh = sigmasq_acct_save,
                       phi_mh = phi_acct_save,
                       zeta_mh = zeta_acct_save)
  
  list(post_sams = post_sams, mh_acct_save = mh_acct_save)
  
}

# Update lambda of the skew-Gaussian NNMP model
updateSNnnmpLa <- function(la, se_la, mu_la, sigmasq_la, 
                           tyy, tyy_ne, tyy_rho, sigmasq, logC, logf) {

  prop_la <- rnorm(1, la, se_la)

  prop_loglik <- 
    dnorm(prop_la, mu_la, sqrt(sigmasq_la), log = TRUE) + 
    sum(logC(tyy, tyy_ne, tyy_rho, prop_la, sigmasq) + logf(tyy, tyy_ne, tyy_rho, prop_la, sigmasq))

  cur_loglik <- 
    dnorm(la, mu_la, sqrt(sigmasq_la), log = TRUE) + 
    sum(logC(tyy, tyy_ne, tyy_rho, la, sigmasq) + logf(tyy, tyy_ne, tyy_rho, la, sigmasq))

  diff_loglik <- prop_loglik - cur_loglik

  if (diff_loglik > log(runif(1))) {

    list(la = prop_la, accept = 1)

  } else {

    list(la = la, accept = 0)

  }  
}

updateSNnnmpLa_covars_lapar <- function(la, se_la, mu_la, sigmasq_la, 
                                        tyy, tyy_ne, tyy_rho, tyy_la, tyy_la_ne,
                                        tyy_nu, tyy_nu_ne, sigmasq, logC, logf,
                                        tyy_par_label, tyy_ne_par_label) {

  K <- length(la)
  accept <- array(0, dim = K)

  for (k in 1:K) {

      k_idx <- which(tyy_par_label == k)
      ne_k_idx <- which(tyy_ne_par_label == k)
      sub_idx <- union(k_idx, ne_k_idx)

      prop_tyy_la <- tyy_la
      prop_tyy_la_ne <- tyy_la_ne
      prop_la <- rnorm(1, la[k], se_la[k])
      prop_tyy_la[k_idx] <- prop_la
      prop_tyy_la_ne[ne_k_idx] <- prop_la

      sub_tyy <- tyy[sub_idx]
      sub_tyy_ne <- tyy_ne[sub_idx]
      sub_tyy_rho <- tyy_rho[sub_idx]
      sub_tyy_nu <- tyy_nu[sub_idx]
      sub_tyy_nu_ne <- tyy_nu_ne[sub_idx]
      

      prop_loglik <- 
        dnorm(prop_la, mu_la[k], sqrt(sigmasq_la[k]), log = TRUE) + 
        sum(logC(sub_tyy, sub_tyy_ne, sub_tyy_rho, sub_tyy_nu, sub_tyy_nu_ne,
                 prop_tyy_la[sub_idx], prop_tyy_la_ne[sub_idx], sigmasq) + 
            logf(sub_tyy, sub_tyy_ne, sub_tyy_rho, sub_tyy_nu, sub_tyy_nu_ne,
                 prop_tyy_la[sub_idx], prop_tyy_la_ne[sub_idx], sigmasq))
      cur_loglik <- 
        dnorm(la[k], mu_la[k], sqrt(sigmasq_la[k]), log = TRUE) +
        sum(logC(sub_tyy, sub_tyy_ne, sub_tyy_rho, sub_tyy_nu, sub_tyy_nu_ne,
                 tyy_la[sub_idx], tyy_la_ne[sub_idx], sigmasq) +
            logf(sub_tyy, sub_tyy_ne, sub_tyy_rho, sub_tyy_nu, sub_tyy_nu_ne,
                 tyy_la[sub_idx], tyy_la_ne[sub_idx], sigmasq))
      

      diff_loglik <- prop_loglik - cur_loglik
      if (diff_loglik > log(runif(1))) {
        tyy_la <- prop_tyy_la
        tyy_la_ne <- prop_tyy_la_ne
        la[k] <- prop_la
        accept[k] <- 1
      }

  }  

  list(la = la, tyy_la = tyy_la, tyy_la_ne = tyy_la_ne, accept = accept)

}

# Update Sigmasq of the skew-Gaussian NNMP model
updateSNnnmpSigmasq <- function(sigmasq, se_sigmasq, u_sigmasq, v_sigmasq,
                                tyy, tyy_ne, tyy_rho, la, logC, logf) {
  
  prop_log_sigmasq <- rnorm(1, log(sigmasq), se_sigmasq)
  prop_sigmasq <- exp(prop_log_sigmasq)
  
  prop_loglik <- 
    dgamma(1 / prop_sigmasq, u_sigmasq, v_sigmasq, log = TRUE) + 2 * log(1 / prop_sigmasq) + 
    sum(logC(tyy, tyy_ne, tyy_rho, la, prop_sigmasq) + 
        logf(tyy, tyy_ne, tyy_rho, la, prop_sigmasq)) + 
    prop_log_sigmasq

  cur_loglik <- 
    dgamma(1 / sigmasq, u_sigmasq, v_sigmasq, log = TRUE) + 2 * log(1 / sigmasq) + 
    sum(logC(tyy, tyy_ne, tyy_rho, la, sigmasq) + 
        logf(tyy, tyy_ne, tyy_rho, la, sigmasq)) + 
    log(sigmasq)
  
  diff_loglik <- prop_loglik - cur_loglik

  if (diff_loglik > log(runif(1))) {

    list(sigmasq = prop_sigmasq, accept = 1)

  } else {

    list(sigmasq = sigmasq, accept = 0)
  }  
  
}

updateSNnnmpSigmasq_covars_lapar <- function(sigmasq, se_sigmasq, u_sigmasq, v_sigmasq,
                                tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne,
                                tyy_la, tyy_la_ne, logC, logf) {
  
  prop_log_sigmasq <- rnorm(1, log(sigmasq), se_sigmasq)
  prop_sigmasq <- exp(prop_log_sigmasq)
  
  prop_loglik <- 
    dgamma(1 / prop_sigmasq, u_sigmasq, v_sigmasq, log = TRUE) + 2 * log(1 / prop_sigmasq) + 
    sum(logC(tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne, tyy_la, tyy_la_ne, prop_sigmasq) + 
        logf(tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne, tyy_la, tyy_la_ne, prop_sigmasq)) + 
    prop_log_sigmasq

  cur_loglik <- 
    dgamma(1 / sigmasq, u_sigmasq, v_sigmasq, log = TRUE) + 2 * log(1 / sigmasq) + 
    sum(logC(tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne, tyy_la, tyy_la_ne, sigmasq) + 
        logf(tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne, tyy_la, tyy_la_ne, sigmasq)) + 
    log(sigmasq)
  
  diff_loglik <- prop_loglik - cur_loglik

  if (diff_loglik > log(runif(1))) {

    list(sigmasq = prop_sigmasq, accept = 1)

  } else {

    list(sigmasq = sigmasq, accept = 0)

  }  
  
}

# Update phi of the skew-Gaussian NNMP model
updateSNnnmpPhi <- function(phi, se_phi, u_phi, v_phi,
                            dat, dat_ne, dat_rho, data_label, dist_mat, rho_mat,
                            la, sigmasq, sp_func, logC, logf) {
  
  tnn <- length(dat)
  prop_log_phi <- rnorm(1, log(phi), se_phi)
  prop_phi <- exp(prop_log_phi)
  prop_rho_mat <- sp_func(dist_mat, prop_phi)
  prop_dat_rho <- as.numeric(t(sapply(1:tnn, function(x) prop_rho_mat[x, data_label[x]])))
  
  prop_loglik <- 
    dgamma(1 / prop_phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / prop_phi) + 
    sum(logC(dat, dat_ne, prop_dat_rho, la, sigmasq) + logf(dat, dat_ne, prop_dat_rho, la, sigmasq)) + prop_log_phi
  cur_loglik <- 
    dgamma(1 / phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / phi) + 
    sum(logC(dat, dat_ne, dat_rho, la, sigmasq) + logf(dat, dat_ne, dat_rho, la, sigmasq)) + log(phi)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(phi = prop_phi, rho_mat = prop_rho_mat, dat_rho = prop_dat_rho, accept = 1))
  } else {
    return(list(phi = phi, rho_mat = rho_mat, dat_rho = dat_rho, accept = 0))
  }
  
}

updateSNnnmpPhi_covars_lapar <- function(phi, se_phi, u_phi, v_phi,
                            dat, dat_ne, dat_rho, data_label, dist_mat, rho_mat,
                            nu, nu_ne, la, la_ne, sigmasq, sp_func, logC, logf) {
  
  trun_size <- length(dat)
  prop_log_phi <- rnorm(1, log(phi), se_phi)
  prop_phi <- exp(prop_log_phi)
  prop_rho_mat <- sp_func(dist_mat, prop_phi)
  prop_dat_rho <- as.numeric(t(sapply(1:trun_size, function(x) prop_rho_mat[x, data_label[x]])))
  
  prop_loglik <- 
    dgamma(1 / prop_phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / prop_phi) + 
    sum(logC(dat, dat_ne, prop_dat_rho, nu, nu_ne, la, la_ne, sigmasq) + 
          logf(dat, dat_ne, prop_dat_rho, nu, nu_ne, la, la_ne, sigmasq)) + prop_log_phi
  cur_loglik <- 
    dgamma(1 / phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / phi) + 
    sum(logC(dat, dat_ne, dat_rho, nu, nu_ne, la, la_ne, sigmasq) + 
          logf(dat, dat_ne, dat_rho, nu, nu_ne, la, la_ne, sigmasq)) + log(phi)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(phi = prop_phi, rho_mat = prop_rho_mat, dat_rho = prop_dat_rho, accept = 1))
  } else {
    return(list(phi = phi, rho_mat = rho_mat, dat_rho = dat_rho, accept = 0))
  }
  
}

updateSNnnmpCoef_lapar <- function(bb, se_bb, 
                                   tyy, tyy_ne, tyy_rho, XX, 
                                   nne, ne_idx, yy_nu, tyy_nu, tyy_nu_ne, 
                                   tyy_la, tyy_la_ne, sigmasq, logC, logf) {
  
  accept <- array(0, dim = length(bb))

  for(j in seq_along(bb)) {

    prop_bb <- bb
    prop_bb[j] <- rnorm(1, bb[j], se_bb[j])
    prop_yy_nu <- as.vector(XX %*% prop_bb)
    prop_tyy_nu <- prop_yy_nu[-(1:nne)]
    prop_tyy_nu_ne <- prop_yy_nu[ne_idx]

    prop_loglik <-  
      sum(logC(tyy, tyy_ne, tyy_rho, prop_tyy_nu, prop_tyy_nu_ne, tyy_la, tyy_la_ne, sigmasq) + 
          logf(tyy, tyy_ne, tyy_rho, prop_tyy_nu, prop_tyy_nu_ne, tyy_la, tyy_la_ne, sigmasq))
    cur_loglik <- 
      sum(logC(tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne, tyy_la, tyy_la_ne, sigmasq) + 
          logf(tyy, tyy_ne, tyy_rho, tyy_nu, tyy_nu_ne, tyy_la, tyy_la_ne, sigmasq))
    diff_loglik <- prop_loglik - cur_loglik
    
    if (diff_loglik > log(runif(1))) {
      bb <- prop_bb
      yy_nu <- prop_yy_nu
      tyy_nu <- prop_tyy_nu
      tyy_nu_ne <- prop_tyy_nu_ne
      accept[j] <- 1
    }

  }
  
  list(bb = bb, yy_nu = yy_nu, tyy_nu = tyy_nu,  tyy_nu_ne = tyy_nu_ne, accept = accept)
  
}

# Update zeta of the skew-Gaussian NNMP model
updateSNnnmpZeta <- function(zeta, se_zeta, nne, dist_mat, mu_t, kasq, data_label, u_zeta, v_zeta){
  
  tnn <- length(mu_t)
  prop_log_zeta <- rnorm(1, log(zeta), se_zeta)
  prop_zeta <- exp(prop_log_zeta)
  prop_weight_res <- logitGausWeight(nne, dist_mat, prop_zeta, mu_t, rep(sqrt(kasq), tnn), trun = TRUE)
  prop_weight_mat <- prop_weight_res$weights
  prop_cutoff <- prop_weight_res$cutoff
  valid_idx <- !is.na(prop_weight_mat)
  if (any(prop_weight_mat[valid_idx]==0)) prop_weight_mat[prop_weight_mat==0] <- .Machine$double.eps
  weight_res <- logitGausWeight(nne, dist_mat, zeta, mu_t, rep(sqrt(kasq), tnn), trun = TRUE)
  weight_mat <- weight_res$weights
  cutoff <- weight_res$cutoff
  
  prop_loglik <- 
    dgamma(1 / prop_zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / prop_zeta) + 
    sum(log(sapply(1:tnn, function(x) prop_weight_mat[x, data_label[x]]))) + prop_log_zeta

  cur_loglik <- 
    dgamma(1 / zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / zeta) + 
    sum(log(sapply(1:tnn, function(x) weight_mat[x, data_label[x]]))) + log(zeta)
  
  diff_loglik <- prop_loglik - cur_loglik

  if (diff_loglik > log(runif(1))) {

    return(list(zeta = prop_zeta, weight_mat = prop_weight_mat, cutoff = prop_cutoff, accept = 1))

  } else {

    return(list(zeta = zeta, weight_mat = weight_mat, cutoff = cutoff, accept = 0))
    
  }
  
}

# Update gamma of the skew-Gaussian NNMP model
updateSNnnmpGa <- function(latent_t, DD, kasq, mu_ga, V_ga) {
  V_ga_inv <- chol2inv(chol(V_ga))
  V_1 <- chol2inv(chol(t(DD) %*% DD / kasq + V_ga_inv))
  mu_ga_1 <- V_1 %*% (t(DD) %*% latent_t / kasq + V_ga_inv %*% mu_ga)
  ga <- t(mvtnorm::rmvnorm(1, mu_ga_1, V_1))
  return(ga)
}

# Update kappa2 of the skew-Gaussian NNMP model
updateSNnnmpKasq <- function(latent_t, mu_t, u_kasq, v_kasq) {
  kk <- length(mu_t)
  uu <- u_kasq + kk / 2
  vv <- v_kasq + .5 * sum((latent_t - mu_t)^2)
  kasq <- 1 / rgamma(1, uu, vv)
  return(kasq)
}