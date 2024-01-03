run_mcmc <- function(all_dat, nne, priors, starting, tuning, mcmc_settings) {
  
  selec_dat <- all_dat$selec_dat
  selec_XX <- all_dat$selec_XX
  ord <- all_dat$ord
  
  nngp_out <- spNNGP(formula = selec_dat[,3] ~ selec_XX[, 2:3],
                     coords = selec_dat[,1:2],
                     starting = starting,
                     tuning = tuning,
                     ord = ord,
                     priors = priors,
                     cov.model = "matern",
                     n.samples = mcmc_settings$niter,
                     n.neighbors = nne,
                     method = "latent",
                     return.neighbor.info = TRUE,
                     n.omp.threads = 1,
                     n.report = 5000)
  
  nngp_out
  
}
