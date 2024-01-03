################################################################################
### SIMULATION STUDY - Model Comparison
### Data generated from a spatial Gaussian copula regression with beta marginals
### Compare NNMP with models: MGP, SPDE, NNGP
################################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(ggplot2)
library(SpatialExtremes)
library(readr)

library(nnmp)

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 10, n_report = 2000)


#-------------------------------------------------------------------------------
# Generate synthetic data
#-------------------------------------------------------------------------------
set.seed(42)

# Build a [0,1] x [0,1] grid at resolution 150 x 150
x_seq <- seq(0, 1, length = 100)
y_seq <- seq(0, 1, length = 100)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
nn <- nrow(grid_loc)

# Generate covariates XX and set true parameters
XX <- as.matrix(rnorm(nrow(grid_loc)))
true_bb <- as.matrix((c(-1.5, -1)))
true_mu <- plogis(cbind(1, XX) %*% true_bb)
true_psi <- 5

# Simulation
range_grid <- seq(0.1, 2, length = 30)
data_list <- vector("list", length = length(range_grid))
names(data_list) <- paste0("data_", seq_along(range_grid))

for (j in seq_along(range_grid)) {
  
  layer_list <- vector("list", length = 9)
  
  # Generate true random field on the grid
  range_val <- range_grid[j]
  gp <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = range_val, smooth = 1))
  yy <- qbeta(pnorm(gp), true_mu * true_psi, (1 - true_mu) * true_psi)
  trf <- cbind(grid_loc, yy)
  colnames(trf) <- c("x", "y", "z")
  layer_list[[1]] <- trf
  
  # Get irregularly spaced data
  obs_num <- 1100
  ref_num <- 800
  test_num <- obs_num - ref_num
  
  layer_list[[2]] <- obs_idx <- sample(1:nn, size = obs_num, replace = FALSE)
  layer_list[[3]] <- ref_idx <- sample(1:obs_num, size = ref_num, replace = FALSE)
  
  layer_list[[4]] <- obs_dat <- trf[obs_idx, ]
  layer_list[[5]] <- obs_XX <- XX[obs_idx, , drop = FALSE]
  
  layer_list[[6]] <- ref_dat <- obs_dat[ref_idx, ]
  layer_list[[7]] <- ref_XX <- obs_XX[ref_idx, , drop = FALSE]
  
  layer_list[[8]] <- test_dat <- obs_dat[-ref_idx, ]
  layer_list[[9]] <- test_XX <- obs_XX[-ref_idx, , drop = FALSE]
  
  names(layer_list) <- c("trf", "obs_idx", "ref_idx", "obs_dat", "obs_XX",
                         "ref_dat", "ref_XX", "test_dat", "test_XX")
  
  data_list[[j]] <- layer_list
  
}


#-------------------------------------------------------------------------------
# Fit NNMP
#-------------------------------------------------------------------------------
priors <- list("scale_gamma" = c(1, 1))

starting <- list("regcoef" = matrix(0, nrow = ncol(cbind(1, XX))), "scale" = 1)

tuning <- list("phi" = 0.2, "zeta" = 0.25,
               "regcoef" = c(0.2, 0.02), "scale" = 0.1)

nne <- 10

set.seed(42)

nnmp_out_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {

  ref_dat <- data_list[[j]]$ref_dat
  ref_XX <- data_list[[j]]$ref_XX

  ord <- sample(1:nrow(ref_XX), replace = FALSE)

  nnmp_out <- nnmp(response = ref_dat[,3],
                   covars = cbind(1, ref_XX),
                   coords = ref_dat[,1:2],
                   neighbor_size = nne,
                   marg_family = "beta",
                   cop_family = "gaussian",
                   priors = priors,
                   starting = starting,
                   tuning = tuning,
                   ord = ord,
                   mcmc_settings = mcmc_settings,
                   verbose = TRUE)

  cat(paste0("Done the ", j, "th run.\n"))

  nnmp_out_list[[j]] <- nnmp_out

  save(nnmp_out, file = paste0(outpath, "nnmp-out/nnmp_out_", j, ".RData"))

}


#-------------------------------------------------------------------------------
# Fit Meshed GP
#-------------------------------------------------------------------------------
library(meshed)

set.seed(42)

mgp_out_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {

  obs_dat <- data_list[[j]]$obs_dat
  ref_idx <- data_list[[j]]$ref_idx
  obs_XX <- data_list[[j]]$obs_XX

  yy <- obs_dat[,3]
  yy[-ref_idx] <- NA
  coords <- obs_dat[, 1:2]
  colnames(coords) <- c("Var1", "Var2")
  mgp_runtime <- system.time(
    mgp_out <- spmeshed(y = yy,
                        x = cbind(1, obs_XX),
                        coords = coords,
                        family = "beta",
                        n_samples = mcmc_settings$n_iter - mcmc_settings$n_burn,
                        n_burn = mcmc_settings$n_burn,
                        n_thin = 1,
                        n_threads = 1,
                        verbose = 5,
                        settings = list(forced_grid = FALSE, low_mem = TRUE),
                        prior = list(phi = c(1, 30))
    ))


  print(mgp_runtime/60)
  cat(paste0("Done the ", j, "th run.\n"))

  mgp_out_list[[j]] <- mgp_out

  save(mgp_out, mgp_runtime, coords, file = paste0(outpath, "mgp-out/mgp_out_", j, ".RData"))

}


#-------------------------------------------------------------------------------
# Fit SPDE-INLA
#-------------------------------------------------------------------------------
library(INLA)

set.seed(42)

spde_out_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {

  ref_dat <- data_list[[j]]$ref_dat
  ref_XX <- data_list[[j]]$ref_XX
  test_dat <- data_list[[j]]$test_dat
  test_XX <- data_list[[j]]$test_XX

  domain <- as.matrix(cbind(c(0,1,1,0,0), c(0, 0, 1, 1,0)))

  coords <- ref_dat[,1:2]
  max_edge <- diff(range(coords[,1])) / (3 * 5) / 2
  mesh <- inla.mesh.2d(coords, max.edge = max_edge, cutoff = 0.01)

  spde <- inla.spde2.matern(mesh)

  A_ref <- inla.spde.make.A(mesh, loc = coords)
  stk_ref <- inla.stack(data = list(y = ref_dat[,3]),
                        A = list(A_ref, 1),
                        tag = "ref",
                        effects = list(list(rf = 1:spde$n.spde),
                                       list(intercept = 1,
                                            x = ref_XX)))

  spde_runtime <- system.time(
    spde_out <- inla(y ~ -1 + intercept + x + f(rf, model = spde),
                     family = "beta",
                     data = inla.stack.data(stk_ref),
                     control.predictor = list(A = inla.stack.A(stk_ref),
                                              compute = TRUE, link = 1),
                     control.compute = list(config =  TRUE))
  )

  cat(paste0("Done the ", j, "th run.\n"))

  spde_out_list[[j]] <- spde_out

  save(spde_out, spde_runtime, mesh, spde, A_ref,
       file = paste0(outpath, "spde-out/spde_out_", j, ".RData"))

}


#-------------------------------------------------------------------------------
# Fit NNGP
#-------------------------------------------------------------------------------
library(spNNGP)

set.seed(42)

starting <- list("phi"=3/0.5, "sigma.sq"=1, "tau.sq"=0.1)
tuning <- list("phi"=0.25)
priors <- list("phi.Unif"=c(1, 30), "sigma.sq.IG"=c(2,1), "tau.sq.IG"=c(2,0.1))

nngp_out_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {
  
  ref_dat <- data_list[[j]]$ref_dat
  ref_XX <- data_list[[j]]$ref_XX
  load(paste0(outpath, "nnmp-out/nnmp_out_", j, ".RData"))
  
  nngp_runtime <- system.time(
    nngp_out <- spNNGP(formula = qlogis(ref_dat[,3]) ~ cbind(1, ref_XX) - 1,
                       coords = ref_dat[,1:2],
                       starting = starting, 
                       tuning = tuning, 
                       ord = nnmp_out$ord_data$ord,
                       priors = priors, 
                       cov.model = "exponential",
                       n.samples = mcmc_settings$n_iter,
                       n.neighbors = nnmp_out$mod_spec$neighbor_size, 
                       method = "latent", 
                       n.omp.threads = 1, 
                       n.report = 2000))
  
  cat(paste0("Done the ", j, "th run.\n"))
  
  nngp_out_list[[j]] <- nngp_out
  
  save(nngp_out, nngp_runtime, file = paste0(outpath, "nngp-out/nngp_out_", j, ".RData"))
  
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Evaluate predictive performance
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# NNMP
#-------------------------------------------------------------------------------
set.seed(42)

nnmp_metric_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {
  
  load(paste0(outpath, "nnmp-out/nnmp_out_", j, ".RData"))
  test_dat <- data_list[[j]]$test_dat
  test_XX <- data_list[[j]]$test_XX
  
  nnmp_test <- predict(nnmp_out, cbind(1, test_XX), test_dat[,1:2], 
                       probs = c(0.025, 0.975), predict_sam = TRUE)
  
  # RMSPE
  nnmp_rmspe <- sqrt(mean((nnmp_test$obs_mu - test_dat[,3])^2))
  
  # MAE
  nnmp_mae <- mean(abs(nnmp_test$obs_mu - test_dat[,3]))
  
  # CRPS
  nnmp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], nnmp_test$sam))
  
  # 95% CI width
  nnmp_qq <- nnmp_test[['obs_qq']]
  nnmp_width <- mean(nnmp_qq[2,] - nnmp_qq[1,])
  
  # 95% CI Coverage rate
  nnmp_cover_boolen <- sapply(1:(test_num), function(x) {
    test_dat[x,3] >= nnmp_qq[1,x] & test_dat[x,3] <= nnmp_qq[2,x]
  })
  nnmp_cover <- sum(nnmp_cover_boolen) / length(nnmp_cover_boolen)
  
  # Put together 
  nnmp_metric <- c(nnmp_rmspe, nnmp_mae, nnmp_crps, nnmp_cover, nnmp_width)
  
  # Save
  nnmp_metric_list[[j]] <- nnmp_metric
  cat(paste0("Done the ", j, "th prediction.\n"))
  
}


#-------------------------------------------------------------------------------
# MGP
#-------------------------------------------------------------------------------
set.seed(42)

mgp_metric_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {
  
  load(file = paste0(outpath, "mgp-out/mgp_out_", j, ".RData"))
  test_dat <- as.data.frame(data_list[[j]]$test_dat)
  colnames(test_dat) <- c("Var1", "Var2", "obs")
  
  ymesh <- data.frame(y_mgp = mgp_out$yhat_mcmc)
  idx <- seq(1, mcmc_settings$n_iter - mcmc_settings$n_burn, by = mcmc_settings$n_thin)
  outdf <- cbind(mgp_out$coordsdata[,1:2], ymesh[, idx])
  test_df <- dplyr::left_join(test_dat, outdf, by = c("Var1", "Var2"))
  
  mgp_test_sam <- as.matrix(test_df[,-(1:3)])
  
  mgp_mean <- rowMeans(mgp_test_sam)
  mgp_q025 <- apply(mgp_test_sam, 1, function(x) quantile(x, probs = 0.025))
  mgp_q975 <- apply(mgp_test_sam, 1, function(x) quantile(x, probs = 0.975))
  
  mgp_rmspe <- sqrt(mean((mgp_mean - test_df$obs)^2))
  
  mgp_mae <- mean(abs(mgp_mean - test_df$obs))
  
  mgp_width <- mean(mgp_q975 - mgp_q025)
  
  mgp_cover_boolen <- sapply(1:test_num, function(x) {
    test_df$obs[x] >= mgp_q025[x] & test_df$obs[x] <= mgp_q975[x]
  })
  mgp_cover <- sum(mgp_cover_boolen) / length(mgp_cover_boolen)
  
  mgp_crps <- mean(scoringRules::crps_sample(test_df$obs, mgp_test_sam))
  
  mgp_metric <- c(mgp_rmspe, mgp_mae, mgp_crps, mgp_cover, mgp_width)
  
  mgp_metric_list[[j]] <- mgp_metric
  cat(paste0("Done the ", j, "th prediction.\n"))
  
}


#-------------------------------------------------------------------------------
# SPDE-INLA
#-------------------------------------------------------------------------------
set.seed(42)

spde_metric_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {
  
  load(file = paste0(outpath, "spde-out/spde_out_", j, ".RData"))
  
  test_dat <- data_list[[j]]$test_dat
  test_XX <- data_list[[j]]$test_XX
  
  nsample <- (mcmc_settings$n_iter - mcmc_settings$n_burn) / mcmc_settings$n_thin
  post_sam <- inla.posterior.sample(nsample, spde_out, add.names = FALSE)
  var_names <- rownames(post_sam[[1]]$latent)
  var_idx <- lapply(c('intercept', 'x', 'rf'), function(nam) 
    which(substr(var_names, 1, nchar(nam)) == nam)) 
  
  mat_sam <- sapply(post_sam, function(iter)
    c(intercept = iter$latent[var_idx[[1]]], 
      x = iter$latent[var_idx[[2]]], 
      rf = iter$latent[var_idx[[3]]]))
  
  hyperparam_sam <- inla.hyperpar.sample(nsample, spde_out)
  phi_sam <- hyperparam_sam[, 1]
  
  A_test <- inla.spde.make.A(mesh, loc = test_dat[,1:2])
  mu_sam <- plogis(as.matrix(cbind(1, test_XX, A_test)) %*% mat_sam)
  spde_test_sam <- t(apply(mu_sam, 1, function(mu) 
    rbeta(test_num, mu * phi_sam, (1 - mu) * phi_sam)))
  
  spde_mean <- rowMeans(spde_test_sam)
  spde_q025 <- apply(spde_test_sam, 1, function(x) quantile(x, probs = 0.025))
  spde_q975 <- apply(spde_test_sam, 1, function(x) quantile(x, probs = 0.975))
  
  spde_rmspe <- sqrt(mean((spde_mean - test_dat[,3])^2))
  
  spde_mae <- mean(abs(spde_mean - test_dat[,3]))
  
  spde_width <- mean(spde_q975 - spde_q025)
  
  spde_cover_boolen <- sapply(1:test_num, function(x) {
    test_dat[x,3] >= spde_q025[x] & test_dat[x,3] <= spde_q975[x]
  })
  spde_cover <- sum(spde_cover_boolen) / length(spde_cover_boolen)
  
  spde_crps <- mean(scoringRules::crps_sample(test_dat[,3], spde_test_sam))
  
  spde_metric <- c(spde_rmspe, spde_mae, spde_crps, spde_cover, spde_width)
  
  spde_metric_list[[j]] <- spde_metric
  cat(paste0("Done the ", j, "th prediction.\n"))
  
}

  
#-------------------------------------------------------------------------------
# NNGP
#-------------------------------------------------------------------------------
set.seed(42)

nngp_metric_list <- vector("list", length = length(range_grid))

for (j in seq_along(range_grid)) {
  
  load(file = paste0(outpath, "nngp-out/nngp_out_", j, ".RData"))
  
  test_dat <- data_list[[j]]$test_dat
  test_XX <- data_list[[j]]$test_XX
  
  nngp_test_out <- predict(nngp_out, X.0 = cbind(1, test_XX), test_dat[,1:2],
                          sub.sample = list(start = mcmc_settings$n_burn, 
                                            end = mcmc_settings$n_iter, 
                                            thin = mcmc_settings$n_thin))
  nngp_test_sam <- plogis(nngp_test_out$p.y.0)
  
  nngp_mean <- rowMeans(nngp_test_sam, 1)
  nngp_q025 <- apply(nngp_test_sam, 1, function(x) quantile(x, probs = 0.025))
  nngp_q975 <- apply(nngp_test_sam, 1, function(x) quantile(x, probs = 0.975))
  
  nngp_rmspe <- sqrt(mean((nngp_mean - test_dat[,3])^2))
  
  nngp_mae <- mean(abs(nngp_mean - test_dat[,3]))
  
  nngp_width <- mean(nngp_q975 - nngp_q025)
  
  nngp_cover_boolen <- sapply(1:test_num, function(x) {
    test_dat[x,3] >= nngp_q025[x] & test_dat[x,3] <= nngp_q975[x]
  })
  nngp_cover <- sum(nngp_cover_boolen) / length(nngp_cover_boolen)
  
  nngp_crps <- mean(scoringRules::crps_sample(test_dat[,3], nngp_test_sam))
  
  nngp_es <- scoringRules::es_sample(test_dat[, 3],  nngp_test_sam)
  
  nngp_vs <- scoringRules::vs_sample(test_dat[, 3], nngp_test_sam, p = 1)
  
  nngp_metric <- c(nngp_rmspe, nngp_mae, nngp_crps, nngp_cover, nngp_width)
  
  nngp_metric_list[[j]] <- nngp_metric
  cat(paste0("Done the ", j, "th prediction.\n"))
  
}


#-------------------------------------------------------------------------------
# Save the prediction results and compare
#-------------------------------------------------------------------------------
save(nnmp_metric_list, mgp_metric_list, spde_metric_list, nngp_metric_list, 
     file = paste0(outpath, "test_metric.RData"))

nnmp_test <- data.frame(nnmp_metric_list)
nnmp_test_avg <- rowMeans(nnmp_test)

mgp_test <- data.frame(mgp_metric_list)
mgp_test_avg <- rowMeans(mgp_test)

spde_test <- data.frame(spde_metric_list)
spde_test_avg <- rowMeans(spde_test)

nngp_test <- data.frame(nngp_metric_list)
nngp_test_avg <- rowMeans(nngp_test)

tbl <- rbind(nnmp_test_avg, mgp_test_avg, spde_test_avg, nngp_test_avg)
colnames(tbl) <- c("RMSPE", "MAE", "CRPS", "95% CI cover", "95% CI width")
rownames(tbl) <- c("NNMP", "MGP", "SPDE", "NNGP")
tbl_latex <- knitr::kable(tbl[,1:5], digits = 3, align = 'c', format = 'latex')
write_file(tbl_latex, paste0(outpath, "beta_s1_pred_comp.txt"))