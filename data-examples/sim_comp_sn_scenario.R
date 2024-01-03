############################################################################
### SIMULATION STUDY: skew-Gaussian NNMP, compared with SG-INLA
### Data generated from a skew-Gaussian random field 
### by Zhang and El-Shaarawi (2010)
############################################################################

rm(list = ls())

outpath <- "/path/to/output/"

# library(RandomFields)
library(SpatialExtremes)
library(ggplot2)
library(readr)

library(nnmp)


#-------------------------------------------------------------------------------
# Generate the true spatial random field
#-------------------------------------------------------------------------------
set.seed(42)

x_seq <- seq(0, 1, length = 100)
y_seq <- seq(0, 1, length = 100)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
zz1 <- abs(t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 0.1, smooth = 1)))
zz2 <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 0.1, smooth = 1))

skewparam <- c(1, 2.5, 5)

dat_save <- vector("list", length = length(skewparam))
for (i in seq_along(skewparam)) {

  obs <- skewparam[i] * zz1  + zz2
  dat <- cbind(grid_loc, obs)
  colnames(dat) <- c("x", "y", "z")
  dat_save[[i]] <- dat
  
}

ref_dat_save <- vector("list", length = length(skewparam))
test_dat_save <- vector("list", length = length(skewparam))

for (i in seq_along(skewparam)) {
  
  dat <- dat_save[[i]]
  
  obs_num <- 1100
  ref_num <- 800
  test_num <- obs_num - ref_num
  
  obs_idx <- sample(1:nrow(grid_loc), size = obs_num, replace = FALSE)
  ref_idx <- sample(1:obs_num, size = ref_num, replace = FALSE)
  
  obs_dat <- dat[obs_idx,]
  ref_dat <- obs_dat[ref_idx,]
  test_dat <- obs_dat[-ref_idx,]
  
  ref_dat_save[[i]] <- ref_dat
  test_dat_save[[i]] <- test_dat
  
}


#-------------------------------------------------------------------------------
# Fit skew-Gaussian NNMP
#-------------------------------------------------------------------------------
set.seed(42)

mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 10, n_report = 2000)             

priors <- list(
  sce1 = list("la_normal" = c(0, 5), "sigmasq_invgamma" = c(2, 1)),
  sce2 = list("la_normal" = c(0, 5), "sigmasq_invgamma" = c(2, 1)),
  sce3 = list("la_normal" = c(0, 5), "sigmasq_invgamma" = c(2, 1))
)                                  

starting <- list(sce1 = list("la" = 0, "sigmasq" = 5), 
                 sce2 = list("la" = 0, "sigmasq" = 5), 
                 sce3 = list("la" = 0, "sigmasq" = 5))  

tuning <- list(
  sce1 = list("phi" = 0.16, "zeta" = 0.3, "la" = 0.7, "sigmasq" = 0.14),
  sce2 = list("phi" = 0.19, "zeta" = 0.35, "la" = 0.7, "sigmasq" = 0.14),
  sce3 = list("phi" = 0.21, "zeta" = 0.5, "la" = 1.6, "sigmasq" = 0.14)
)   

nne <- 10
         
nnmp_out_list <- vector("list", length = length(skewparam))

for (i in seq_along(skewparam)) {
  
  ref_dat <- ref_dat_save[[i]] 
  
  # Ordering
  ord <- sample(1:ref_num, replace = FALSE)
  
  # MCMC
  nnmp_out <- nnmp(response = ref_dat[,3],
                 coords = ref_dat[,1:2],
                 neighbor_size = nne,
                 marg_family = "sn",
                 priors = priors[[i]],
                 starting = starting[[i]],
                 tuning = tuning[[i]],
                 ord = ord,
                 mcmc_settings = mcmc_settings, 
                 verbose = TRUE)

  nnmp_out_list[[i]] <- nnmp_out

}


#-------------------------------------------------------------------------------
# Evaluate skew-Gaussian NNMP predictive performance
#-------------------------------------------------------------------------------
set.seed(42)

coef_skew <- function(yy) {
  yy_mom3 <- mean((yy - mean(yy))^3)
  yy_mom2 <- mean((yy - mean(yy))^2)
  yy_mom3 / yy_mom2^(3/2)
}

nnmp_metric_list <- vector("list", length = length(skewparam))

for (i in seq_along(skewparam)) {
  
  test_dat <- test_dat_save[[i]]
  
  nnmp_out <- nnmp_out_list[[i]]

  nnmp_test <- predict(nnmp_out, nonref_coords = test_dat[,1:2], predict_sam = TRUE)

  # RMSPE
  nnmp_rmspe <- sqrt(mean((nnmp_test$obs_mu - test_dat[,3])^2))
  
  # MAE
  nnmp_mae <- mean(abs(nnmp_test$obs_mu - test_dat[,3]))

  # CRPS
  nnmp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], nnmp_test$sam))
  
  # 95% CI width
  nnmp_qq <- nnmp_test$obs_qq
  nnmp_width <- mean(nnmp_qq[3,] - nnmp_qq[1,])
  
  # 95% CI Coverage rate
  nnmp_cover_boolen <- sapply(1:(test_num), function(x) {
    test_dat[x,3] >= nnmp_qq[1,x] & test_dat[x,3] <= nnmp_qq[3,x]
  })
  nnmp_cover <- sum(nnmp_cover_boolen) / length(nnmp_cover_boolen)
  
  # Skewness
  nnmp_coef_skew_sam <- apply(nnmp_test$sam, 2, coef_skew)
  nnmp_skew_score <- mean(abs(nnmp_coef_skew_sam - coef_skew(test_dat[,3])))
  
  # Put together
  nnmp_metric<- c(nnmp_rmspe, nnmp_mae, nnmp_crps, nnmp_cover, nnmp_width, nnmp_skew_score)
  nnmp_metric_list[[i]] <- nnmp_metric
  
}


#-------------------------------------------------------------------------------
# Fit SG-INLA
#-------------------------------------------------------------------------------
library(INLA)

set.seed(42)

spde_metric_list <- vector("list", length = length(skewparam))
spde_alpha_list <- vector("list", length = length(skewparam))
for (i in seq_along(skewparam)) {
  
  ref_dat <- ref_dat_save[[i]] 
  test_dat <- test_dat_save[[i]] 
  test_num <- nrow(test_dat)
  
  
  coords <- ref_dat[,1:2]
  max_edge <- diff(range(coords[,1])) / (3 * 5) / 3
  mesh <- inla.mesh.2d(coords, max.edge = max_edge, cutoff = 0.01)
  
  spde <- inla.spde2.matern(mesh)
  
  A_ref <- inla.spde.make.A(mesh, loc = coords)
  stk_ref <- inla.stack(data = list(y = ref_dat[,3]), 
                        A = list(A_ref, 1),
                        tag = "ref", 
                        effects = list(list(rf = 1:spde$n.spde),
                                       list(mu = rep(1, nrow(ref_dat)))))
  
  A_test <- inla.spde.make.A(mesh, loc = test_dat[,1:2])
  
  spde_runtime <- system.time(
    spde_out <- inla(y ~ -1 + mu + f(rf, model = spde),
                     family = "sn", 
                     data = inla.stack.data(stk_ref), 
                     control.predictor = list(A = inla.stack.A(stk_ref), compute = TRUE),
                     control.compute = list(config = TRUE))
  )
  
  print(paste0("Done the ", i, "th run."))
  
  nsample <- (mcmc_settings$n_iter - mcmc_settings$n_burn) / mcmc_settings$n_thin
  post_sam <- inla.posterior.sample(nsample, spde_out, add.names = FALSE)
  var_names <- rownames(post_sam[[1]]$latent)

  var_idx <- lapply(c('mu', 'rf'), function(nam)
    which(substr(var_names, 1, nchar(nam)) == nam))
  mat_sam <- sapply(post_sam, function(iter)
    c(mu = iter$latent[var_idx[[1]]],
      rf = iter$latent[var_idx[[2]]]))
  eta_sam <- as.matrix(cbind(1, A_test)) %*% mat_sam
  
  hyperparam_sam <- inla.hyperpar.sample(nsample, spde_out)
  
  prec_sam <- hyperparam_sam[, 1]
  skew_sam <- hyperparam_sam[, 2]
  del_sam <- sqrt(pi/2 *abs(skew_sam)^(2/3) /
                    (abs(skew_sam)^(2/3) + ((4-pi)/2)^(2/3)))
  del_sam = del_sam * sign(skew_sam)
  alp_sam <- del_sam / sqrt(1 - del_sam^2)
  ome_sam <- sqrt(1 / (1 - 2 * del_sam^2 / pi))
  xi_sam <- -ome_sam * del_sam * sqrt(2/pi)
  spde_test_sam <- t(apply(eta_sam, 1, function(x) 
    x + sn::rsn(nsample, xi_sam, ome_sam, alp_sam) / sqrt(prec_sam)))
  
  spde_mean <- rowMeans(spde_test_sam)
  spde_q025 <- apply(spde_test_sam, 1, function(x) quantile(x, probs = 0.025))
  spde_q975 <- apply(spde_test_sam, 1, function(x) quantile(x, probs = 0.975))
  
  spde_rmspe <- sqrt(mean((spde_mean - test_dat[,3])^2))
  
  spde_mae <- mean(abs(spde_mean - test_dat[,3]))
  
  spde_crps <- mean(scoringRules::crps_sample(test_dat[,3], spde_test_sam))
  
  spde_width <- mean(spde_q975 - spde_q025)
  
  spde_cover_boolen <- sapply(1:test_num, function(x) {
    test_dat[x,3] >= spde_q025[x] & test_dat[x,3] <= spde_q975[x]
  })
  spde_cover <- sum(spde_cover_boolen) / length(spde_cover_boolen)
  
  spde_coef_skew_sam <- apply(spde_test_sam, 2, coef_skew)
  spde_skew_score <- mean(abs(spde_coef_skew_sam - coef_skew(test_dat[,3])))
  
  spde_metric <- c(spde_rmspe, spde_mae, spde_crps, spde_cover, spde_width, spde_skew_score)
  
  spde_metric_list[[i]] <- spde_metric
  
  print(paste0("Done the ", i, "th prediction."))
  
}

specRound <- function(x, k = 3) trimws(format(round(x, k), nsmall = k))
nnmp_test <- specRound(as.matrix(data.frame(nnmp_metric_list)))
spde_test <- specRound(as.matrix(data.frame(spde_metric_list)))

tbl <- cbind(nnmp_test[,1], spde_test[,1], rep("", 6),
             nnmp_test[,2], spde_test[,2], rep("", 6),
             nnmp_test[,3], spde_test[,3])

rownames(tbl) <- c("RMSPE", "MAE", "CRPS", "95% CI cover", "95% CI width", "Skewness")
tbl_latex <- knitr::kable(tbl, align = 'c', format = 'latex')
write_file(tbl_latex, file = paste0(outpath, "skew_table_comp.txt"))
