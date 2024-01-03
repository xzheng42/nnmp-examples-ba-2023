run_mcmc <- function(all_dat, nne, tuning, mcmc_settings) {
  
  selec_dat <- all_dat$selec_dat
  selec_XX <- all_dat$selec_XX
  selec_par_labels <- all_dat$selec_par_labels
  ord <- all_dat$ord
  
  K <- max(selec_par_labels)
  
  priors <- list("sigmasq_invgamma" = c(3, 1),
                 "la_normal" = list(mu = rep(0, K), sigmasq = rep(5, K)))
  
  starting <- list("la" = rep(0, K), "sigmasq" = 5, 
                   "regcoef" = c(mean(selec_dat[,3]), rep(0, ncol(selec_XX) - 1)))
  
  nnmp_out <- nnmp(response = selec_dat[,3],
                   covars = selec_XX,
                   coords = selec_dat[,1:2],
                   neighbor_size = nne,
                   marg_family = "sn",
                   priors = priors,
                   starting = starting,
                   tuning = tuning,
                   ord = ord,
                   mcmc_settings = mcmc_settings, 
                   verbose = TRUE,
                   par_label = selec_par_labels)
  
  nnmp_out
  
}
