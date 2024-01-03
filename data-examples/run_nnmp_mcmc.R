run_mcmc <- function(all_dat, nne, tuning, mcmc_settings) {
  
  selec_dat <- all_dat$selec_dat
  selec_XX <- all_dat$selec_XX
  ord <- all_dat$ord
  
  nnmp_out <- nnmp(response = selec_dat[,3],
                   covars = selec_XX,
                   coords = selec_dat[,1:2],
                   neighbor_size = nne,
                   marg_family = "gaussian",
                   tuning = tuning,
                   ord = ord,
                   mcmc_settings = mcmc_settings,
                   verbose = TRUE,
                   model_diag = list("dic", "pplc"))
  
  nnmp_out
  
}
