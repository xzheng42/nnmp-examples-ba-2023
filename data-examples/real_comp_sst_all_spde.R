############################################################################
### SST data analysis using skew-Gaussian SPDE-INLA
############################################################################
rm(list = ls())

outpath <- "~/path/to/output/"
datapath <- "/path/to/input/data/"

library(readr)
library(INLA)
library(scoringRules)

#-------------------------------------------------------------------------------
# Import data
#-------------------------------------------------------------------------------
dat <- read_csv(paste0(datapath, "med_sst.csv"))
med_shape <- sf::st_read(paste0(datapath,"med_lme/lme.shp"))
med_shape_sp <- sf::as_Spatial(med_shape)

#-------------------------------------------------------------------------------
# 5-fold cross validation
#-------------------------------------------------------------------------------
set.seed(42)

coef_skew <- function(yy) {
  yy_mom3 <- mean((yy - mean(yy))^3)
  yy_mom2 <- mean((yy - mean(yy))^2)
  yy_mom3 / yy_mom2^(3/2)
}

kk <- 5
folds_idx <- sample(1:kk, nrow(dat), replace = TRUE)

tbl <- array(NA, dim = c(kk, 5))
KK <- max(dat$par_label)
tbl_skew <- array(NA, dim = c(kk, KK))

for (j in 1:kk) {
  
  #-----------------------------------------------------------------------------
  # Fit the model
  #-----------------------------------------------------------------------------
  set.seed(42)
   
  ref_num <- nrow(dat)
  selec_idx <- folds_idx != j
  selec_num <- sum(selec_idx)
  valid_num <- ref_num - selec_num
  
  ref_dat <- as.matrix(dat[, 1:3])
  ref_XX <- cbind(1, ref_dat[, 1:2])
  
  selec_dat <- ref_dat[selec_idx,]
  selec_XX <- ref_XX[selec_idx, ]
  valid_dat <- ref_dat[!selec_idx,]
  valid_XX <- ref_XX[!selec_idx, ]
  
  ref_par_labels <- dat$par_label
  valid_par_labels <- ref_par_labels[!selec_idx]
  
  coords <- selec_dat[, 1:2]
  bound <- inla.sp2segment(med_shape_sp)
  max_edge <- min(c(diff(range(coords[,1])), diff(range(coords[, 2])))) / (3 * 5) / 2
  mesh <- inla.mesh.2d(boundary = bound, loc = coords,max.edge = c(1, 2) * max_edge, cutoff = max_edge / 2)
  spde <- inla.spde2.matern(mesh)
  
  A_ref <- inla.spde.make.A(mesh, loc = coords)
  stk_ref <- inla.stack(data = list(y = selec_dat[,3]), 
                        A = list(A_ref, 1),
                        tag = "ref", 
                        effects = list(list(rf = 1:spde$n.spde),
                                       list(mu = rep(1, selec_num),
                                            lon = coords[, 1],
                                            lat = coords[, 2],
                                            err = 1:nrow(coords))))
  
  spde_runtime <- system.time(
    spde_out <- inla(y ~ -1 + mu + lon + lat + f(rf, model = spde) + f(err, model = "iid"),
                     family = "sn",
                     data = inla.stack.data(stk_ref), 
                     control.predictor = list(A = inla.stack.A(stk_ref), compute = TRUE),
                     control.compute = list(config = TRUE))
  )  
  
  #-----------------------------------------------------------------------------
  # Evaluate predictive performances
  #-----------------------------------------------------------------------------  
  nsample <- 5000
  post_sam <- inla.posterior.sample(nsample, spde_out, add.names = FALSE)
  var_names <- rownames(post_sam[[1]]$latent)
  
  var_idx <- lapply(c('mu', 'lon', 'lat', 'rf'), function(nam)
    which(substr(var_names, 1, nchar(nam)) == nam))
  mat_sam <- sapply(post_sam, function(iter)
    c(mu = iter$latent[var_idx[[1]]],
      lon = iter$latent[var_idx[[2]]],
      lat = iter$latent[var_idx[[3]]],
      rf = iter$latent[var_idx[[4]]]))
  
  A_test <- inla.spde.make.A(mesh, loc = valid_dat[,1:2])
  eta_sam <- as.matrix(cbind(valid_XX, A_test)) %*% mat_sam
  
  hyperparam_sam <- inla.hyperpar.sample(nsample, spde_out)

  prec_sam <- hyperparam_sam[, 1]
  skew_sam <- hyperparam_sam[, 2]
  del_sam <- sqrt(pi/2 *abs(skew_sam)^(2/3) /
                    (abs(skew_sam)^(2/3) + ((4-pi)/2)^(2/3)))
  del_sam = del_sam * sign(skew_sam)
  alp_sam <- del_sam / sqrt(1 - del_sam^2)
  ome_sam <- sqrt(1 / (1 - 2 * del_sam^2 / pi))
  xi_sam <- -ome_sam * del_sam * sqrt(2/pi)
  y_sam <- t(apply(eta_sam, 1, function(x)
    x + sn::rsn(nsample, xi_sam, ome_sam, alp_sam) / sqrt(prec_sam)))
  
  psisq_sam <- 1 / hyperparam_sam[, 5]
  spde_test_sam <- t(apply(y_sam, 1, function(x) rnorm(nsample, x, sqrt(psisq_sam))))
  
  spde_mean <- rowMeans(spde_test_sam)
  spde_q025 <- apply(spde_test_sam, 1, function(x) quantile(x, probs = 0.025))
  spde_q975 <- apply(spde_test_sam, 1, function(x) quantile(x, probs = 0.975))
  
  spde_rmspe <- sqrt(mean((spde_mean - valid_dat[,3])^2))
  
  spde_mae <- mean(abs(spde_mean - valid_dat[,3]))
  
  spde_crps <- mean(scoringRules::crps_sample(valid_dat[,3], spde_test_sam))
  
  spde_width <- mean(spde_q975 - spde_q025)
  
  spde_cover_boolen <- sapply(1:valid_num, function(x) {
    valid_dat[x,3] >= spde_q025[x] & valid_dat[x,3] <= spde_q975[x]
  })
  spde_cover <- sum(spde_cover_boolen) / length(spde_cover_boolen)
  
  spde_metric <- c(spde_rmspe, spde_mae, spde_crps, spde_cover, spde_width)
  
  tbl[j, ] <- spde_metric
  
  cat(paste0("Done the ", j, "-th fold prediction\n\n"))
  
  # Skewness
  dat_coef_skew_par <- array(NA, dim = KK)
  for (l in 1:KK) {
    sst_par_l <- valid_dat[valid_par_labels == l, 3]
    dat_coef_skew_par[l] <- coef_skew(sst_par_l)
  }
  
  spde_skew_score <- array(NA, dim = KK)
  for (l in 1:KK) {
    spde_coef_skew_sam <- apply(spde_test_sam[valid_par_labels==l, ], 2, coef_skew)
    tbl_skew[j, l]  <- mean(abs(spde_coef_skew_sam - dat_coef_skew_par[l]))
  }
  spde_coef_skew_sam <- apply(spde_test_sam, 2, coef_skew)
  
}


tbl_mu <- t(colMeans(tbl))
rownames(tbl_mu) <- "SG-SPDE"
colnames(tbl_mu) <- c("RMSPE", "MAE", "CRPS","95% CI cover", "95% CI width")
tbl_mu_latex <- knitr::kable(tbl_mu, format = "latex", digits = 3, align = "c")
write_file(tbl_mu_latex, paste0(outpath, "sst_spde_cv_pred.txt"))


tbl_skew_mu <- t(colMeans(tbl_skew))
rownames(tbl_skew_mu) <- "SG-SPDE"
colnames(tbl_skew_mu) <- paste0("Sub-basin-", 1:KK)
tbl_skew_mu_latex <- knitr::kable(tbl_skew_mu, format = "latex", digits = 3, align = "c")
write_file(tbl_skew_mu_latex, paste0(outpath, "sst_spde_cv_skew.txt"))

