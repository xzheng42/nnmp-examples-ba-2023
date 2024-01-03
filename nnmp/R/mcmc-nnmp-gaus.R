gausNNMP <- function(yy, XX, coords, ne_info, sp_func, 
                     priors, starting, tuning, mcmc_settings, verbose) {
  
  
  #--------------------------
  # Check priors 
  #--------------------------
  
  # Priors
  if (is.null(priors$sigmasq_invgamma)) {
    u_sigmasq <- 2
    v_sigmasq <- 1
  } else {
    u_sigmasq <- priors$sigmasq_invgamma[1]
    v_sigmasq <- priors$sigmasq_invgamma[2]
  }
  
  if (is.null(priors$tausq_invgamma)) {
    u_tausq <- 2
    v_tausq <- 0.1
  } else {
    u_tausq <- priors$tausq_invgamma[1]
    v_tausq <- priors$tausq_invgamma[2]
  }
  
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
  
  #--------------------------
  # Check starting values
  #--------------------------
  if (is.null(starting$bb)) {
    bb <- c(mean(yy), rep(0, ncol(XX) - 1))
  } else {
    bb <- starting$bb
  }
  ga <- starting$ga
  kasq <- starting$kasq
  phi <- starting$phi
  zeta <- starting$zeta
  
  #--------------------------
  # MCMC
  #--------------------------
  
  runtime <- system.time(
    mcmc_out <- mcmcGausNNMP(yy, XX, coords,
                             ne_info, sp_func, 
                             mcmc_settings, verbose,
                             
                             u_sigmasq,v_sigmasq, u_tausq,v_tausq,
                             u_phi, v_phi, u_zeta, v_zeta,
                             mu_ga, V_ga, u_kasq,v_kasq,
                             
                             se_phi, se_zeta,
                             
                             phi, zeta, bb, ga, kasq
                             )
    )
  
  mcmc_out$runtime <- runtime
  
  mcmc_out
  
  
}


mcmcGausNNMP <- function(yy, XX, coords, 
                         ne_info, sp_func, mcmc_settings, verbose,
                         u_sigmasq,v_sigmasq, u_tausq,v_tausq,
                         u_phi, v_phi, u_zeta, v_zeta,
                         mu_ga, V_ga, u_kasq,v_kasq, 
                         se_phi, se_zeta,
                         phi, zeta, bb, ga, kasq
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
  dist_mat <- ne_info$ne_dist
  ne_idx_mat <- ne_info$ne_index
  DD <- as.matrix(cbind(1, coords))[-(1:2),]
  mu_t <- DD %*% ga
  
  weight_res <- logitGausWeight(nne, dist_mat, zeta, mu_t, rep(sqrt(kasq), nn - 2), FALSE)
  weight_mat <- weight_res$weights
  cutoff <- weight_res$cutoff  
  valid_idx <- !is.na(weight_mat)
  if (any(weight_mat[valid_idx]==0)){
    weight_mat[weight_mat==0] <- .Machine$double.eps
  }
  
  zz_label <- array(NA, dim = nn)
  zz_label[1] <- 0
  zz_label[2] <- 1
  for (i in  3:nne) {
    zz_label[i] <- sample(1:(i-1), size = 1, prob = weight_mat[i, 1:(i-1)])
  }
  for (i in (nne+1):nn) {
    zz_label[i] <- sample(1:nne, size = 1, prob = weight_mat[i,])
  }
  ne_idx <- as.numeric(t(sapply(1:nn, function(x) ne_idx_mat[x, zz_label[x]])))
  zz <- yy - as.numeric(XX %*% bb)
  zz_ne <- c(0, zz[ne_idx[-1]])
  rho_mat <- sp_func(dist_mat, phi)
  zz_rho <- c(0, as.numeric(t(sapply(2:nn, function(x) rho_mat[x, zz_label[x]]))))
  
  ## empty stuff to save samples
  zz_save <- array(NA, dim = c(nn, niter - nburn))
  # weight_save <- array(NA, dim = c(nn, nne, niter - nburn))
  phi_save <- array(NA, dim = niter - nburn)
  sigmasq_save <- array(NA, dim = niter - nburn)
  tausq_save <- array(NA, dim = niter - nburn)
  kasq_save <- array(NA, dim = niter - nburn)
  beta_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  ga_save <- array(NA, dim = c(ncol(DD), niter - nburn))
  zeta_save <- array(NA, dim = niter - nburn)
  phi_acct_save <- array(0, dim = niter)
  zeta_acct_save <- array(0, dim = niter)
  
  
  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0
  
  for (iter in 1:niter) {
    
    start_time <- Sys.time()
    #--------------------------
    # Update tausq
    #--------------------------
    tausq <- updateGnnmpTausq(yy, XX, bb, zz, u_tausq, v_tausq)
    
    #--------------------------
    # Update sigmasq
    #--------------------------
    sigmasq <- updateGnnmpSigmasq(zz, zz_ne, zz_rho, u_sigmasq, v_sigmasq)
    s2 <- sigmasq * (1 - zz_rho^2)
    
    #--------------------------
    # Update phi
    #--------------------------
    phi_res <- updateGnnmpPhi(phi, se_phi, sp_func, dist_mat, rho_mat, zz, zz_ne, 
                              zz_rho, zz_label, sigmasq, s2, u_phi, v_phi)
    phi <- phi_res$phi
    rho_mat <- phi_res$rho_mat
    zz_rho <- phi_res$zz_rho
    if (phi_res$accept) phi_acct_save[iter] <- 1
    
    #--------------------------
    # Update beta
    #--------------------------
    bb <- updateGnnmpBeta(yy, XX, zz, tausq)
    
    #--------------------------
    # Update zz
    #--------------------------
    zz <- updateGnnmpSPE(yy, XX, bb, zz, ne_idx, zz_rho, tausq, sigmasq)
    zz <- as.numeric(zz)
    
    #--------------------------
    # Update labels
    #--------------------------
    zz_ne_mat <- t(sapply(1:nn, function(x) zz[ne_idx_mat[x,]]))
    labels <- updateGnnmpLabel(zz, zz_ne_mat, nne, weight_mat,
                               rep(0, nne), rep(sigmasq, nne), rho_mat, 
                               mu_t, rep(sqrt(kasq), length(mu_t)), cutoff)
    zz_label <- as.numeric(labels$latent_label)
    latent_t <- labels$latent_t[-(1:2)]
    ne_idx <- as.numeric(t(sapply(1:nn, function(x) ne_idx_mat[x, zz_label[x]])))
    zz_ne <- c(0, zz[ne_idx[-1]])
    zz_rho <- c(0, as.numeric(t(sapply(2:nn, function(x) rho_mat[x, zz_label[x]]))))
    
    #--------------------------
    # Update ga
    #--------------------------
    ga <- updateGnnmpGa(latent_t, DD, kasq, mu_ga, V_ga)
    mu_t <- as.numeric(DD %*% ga)
    
    #--------------------------
    # Update kasq
    #--------------------------
    kasq <- updateGnnmpKasq(latent_t, mu_t, u_kasq, v_kasq)
    
    #--------------------------
    # Update zeta and weights
    #--------------------------
    zeta_res <- updateGnnmpZeta(zeta, se_zeta, nne, dist_mat, mu_t, kasq, 
                                zz_label, u_zeta, v_zeta)
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
      zz_save[, iter - nburn] <- zz
      beta_save[, iter - nburn] <- bb
      ga_save[, iter - nburn] <- ga
      sigmasq_save[iter - nburn] <- sigmasq
      tausq_save[iter - nburn] <- tausq
      phi_save[iter - nburn] <- phi
      zeta_save[iter - nburn] <- zeta
      kasq_save[iter - nburn] <- kasq
    }
    
  }
  
  #--------------------------
  # Thinning
  #--------------------------
  selc_index <- seq(1, niter - nburn, by = nthin)
  
  post_sams <- list(#weight = weight_save[, , selc_index],
    zz = zz_save[, selc_index],
    bb = beta_save[,selc_index],
    ga = ga_save[,selc_index],
    sigmasq = sigmasq_save[selc_index],
    tausq = tausq_save[selc_index],
    phi = phi_save[selc_index],
    zeta = zeta_save[selc_index],
    kasq = kasq_save[selc_index])
  
  mh_acct_save <- list(phi_mh = phi_acct_save,
                       zeta_mh = zeta_acct_save)
  
  list(post_sams = post_sams, mh_acct_save = mh_acct_save)
  
}


# Update Tausq of the Gaussian NNMP spatial regression model
updateGnnmpTausq <- function(yy, XX, bb, zz, u_tausq, v_tausq) {
  nn <- length(zz)
  resid <- yy - as.numeric(XX %*% bb) - zz
  uu <- u_tausq + nn / 2
  vv <- v_tausq + .5 * sum(resid^2)
  tausq <- 1 / rgamma(1, uu, vv)
  return(tausq)
}

# Update Sigmasq of the Gaussian NNMP spatial regression model
updateGnnmpSigmasq <- function(zz, zz_ne, zz_rho, u_sigmasq, v_sigmasq) {
  nn <- length(zz)
  uu <- u_sigmasq + nn / 2
  vv <- v_sigmasq + .5 * sum((zz - zz_rho * zz_ne)^2 / (1 - zz_rho^2))
  sigmasq <- 1 / rgamma(1, uu, vv)
  return(sigmasq)
}

# Update phi of the Gaussian NNMP spatial regression model
updateGnnmpPhi <- function(phi, se_phi, sp_func, dist_mat, rho_mat, zz, zz_ne, 
                           zz_rho, zz_label, sigmasq, s2, u_phi, v_phi) {
  
  nn <- length(zz)
  prop_log_phi <- rnorm(1, log(phi), se_phi)
  prop_phi <- exp(prop_log_phi)
  prop_rho_mat <- sp_func(dist_mat, prop_phi)
  prop_zz_rho <- c(0, as.numeric(t(sapply(2:nn, function(x) prop_rho_mat[x, zz_label[x]]))))
  prop_s2 <- sigmasq * (1 - prop_zz_rho^2)
  
  prop_loglik <- 
    dgamma(1/prop_phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / prop_phi) + 
    sum(dnorm(zz, prop_zz_rho * zz_ne, sqrt(prop_s2), log = TRUE)[-1]) + prop_log_phi
  cur_loglik <- 
    dgamma(1/phi, u_phi, v_phi, log = TRUE) + 2 * log(1 / phi) + 
    sum(dnorm(zz, zz_rho * zz_ne, sqrt(s2), log = TRUE)[-1]) + log(phi)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(phi = prop_phi, rho_mat = prop_rho_mat, zz_rho = prop_zz_rho, accept = TRUE))
  } else {
    return(list(phi = phi, rho_mat = rho_mat, zz_rho = zz_rho, accept = FALSE))
  }
  
}

# Update zeta of the Gaussian NNMP spatial regression model
updateGnnmpZeta <- function(zeta, se_zeta, nne, dist_mat, mu_t, kasq, zz_label, u_zeta, v_zeta){
  
  nn <- nrow(dist_mat)
  prop_log_zeta <- rnorm(1, log(zeta), se_zeta)
  prop_zeta <- exp(prop_log_zeta)
  prop_weight_res <- logitGausWeight(nne, dist_mat, prop_zeta, mu_t, rep(sqrt(kasq), nn - 2), FALSE)
  prop_weight_mat <- prop_weight_res$weights
  prop_cutoff <- prop_weight_res$cutoff
  weight_res <- logitGausWeight(nne, dist_mat, zeta, mu_t, rep(sqrt(kasq), nn - 2), FALSE)
  weight_mat <- weight_res$weights
  cutoff <- weight_res$cutoff
  
  prop_loglik <- 
    dgamma(1 / prop_zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / prop_zeta) + 
    sum(log(sapply(3:nn, function(x) prop_weight_mat[x, zz_label[x]]))) + prop_log_zeta
  cur_loglik <- 
    dgamma(1 / zeta, u_zeta, v_zeta, log = TRUE) + 2 * log(1 / zeta) + 
    sum(log(sapply(3:nn, function(x) weight_mat[x, zz_label[x]]))) + log(zeta)
  
  diff_loglik <- prop_loglik - cur_loglik
  if (diff_loglik > log(runif(1))) {
    return(list(zeta = prop_zeta, weight_mat = prop_weight_mat, cutoff = prop_cutoff, accept = TRUE))
  } else {
    return(list(zeta = zeta, weight_mat = weight_mat, cutoff = cutoff, accept = FALSE))
  }
  
}

# Update beta of the Gaussian NNMP spatial regression model
updateGnnmpBeta <- function(yy, XX, zz, tausq) {
  V_1 <- tausq * chol2inv(chol(t(XX) %*% XX))
  mu_beta_1 <- V_1 %*% (t(XX) %*% (yy - zz)) / tausq
  bb <- t(mvtnorm::rmvnorm(1, mu_beta_1, V_1, method = "chol"))
  return(bb)
}

# Update gamma of the Gaussian NNMP spatial regression model
updateGnnmpGa <- function(latent_t, DD, kasq, mu_ga, V_ga) {
  V_ga_inv <- chol2inv(chol(V_ga))
  V_1 <- chol2inv(chol(t(DD) %*% DD / kasq + V_ga_inv))
  mu_ga_1 <- V_1 %*% (t(DD) %*% latent_t / kasq + V_ga_inv %*% mu_ga)
  ga <- t(mvtnorm::rmvnorm(1, mu_ga_1, V_1))
  return(ga)
}

# Update kappa2 of the Gaussian NNMP spatial regression model
updateGnnmpKasq <- function(latent_t, mu_t, u_kasq, v_kasq) {
  kk <- length(mu_t)
  uu <- u_kasq + kk / 2
  vv <- v_kasq + .5 * sum((latent_t - mu_t)^2)
  kasq <- 1 / rgamma(1, uu, vv)
  return(kasq)
}
