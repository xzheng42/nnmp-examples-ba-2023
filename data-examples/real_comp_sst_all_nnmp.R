############################################################################
### SST data analysis using extended Skew-GNNMPs
############################################################################
rm(list = ls())

outpath <- "~/path/to/output/"
datapath <- "/path/to/input/data/"

library(readr)
library(parallel)
library(doSNOW)
library(foreach)

library(nnmp)

#-------------------------------------------------------------------------------
# Import data
#-------------------------------------------------------------------------------
dat <- read_csv(paste0(datapath, "med_sst.csv"))

#-------------------------------------------------------------------------------
# Fit extended skew-Gaussian NNMP
#-------------------------------------------------------------------------------
set.seed(42)

kk <- 5
folds_idx <- sample(1:kk, nrow(dat), replace = TRUE)

mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 4, n_report = 2000)
ne_sizes <- seq(5, 50, by = 5)

dat_list <- vector("list", length = kk)
for (j in 1:kk) {
  
  set.seed(42)
  
  # Partition the observations into training and validation sets
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
  selec_par_labels <- ref_par_labels[selec_idx]
  valid_par_labels <- ref_par_labels[!selec_idx]
  
  # Random ordering
  ord <- sample(1:selec_num, replace = FALSE)
  
  dat_list[[j]] <- list(ref_dat = ref_dat, ref_XX = ref_XX, 
                        selec_dat = selec_dat, selec_XX = selec_XX, 
                        valid_dat = valid_dat, valid_XX = valid_XX,
                        ref_par_labels = ref_par_labels,
                        selec_par_labels = selec_par_labels,
                        valid_par_labels = valid_par_labels, ord = ord)
}

# MCMC
tuning <- list("phi" = 0.12, "zeta" = 0.6, "sigmasq" = 0.1,
               "regcoef" = c(0.3, 0.02, 0.008), 
               "la" = c(0.35, 0.6, 1.1, 0.5, 0.7))


closeAllConnections()
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

runtime <- system.time(out_list <- foreach(j = 1:kk, .packages = "nnmp") %:%
                         foreach(k = seq_along(ne_sizes), .packages = "nnmp") %dopar% {
                           set.seed(42)
                           source("path/to/run_ext_sgnnmp_mcmc.R")
                           run_mcmc(dat_list[[j]], ne_sizes[k], tuning, mcmc_settings)
                         }
)

cat("done\n")
cat(paste0("Running time: ", round(runtime[3]/60), " minutes\n"))

stopCluster(cl)

#-------------------------------------------------------------------------------
# Extended skew-Gaussian NNMP prediction 
#-------------------------------------------------------------------------------
coef_skew <- function(yy) {
  yy_mom3 <- mean((yy - mean(yy))^3)
  yy_mom2 <- mean((yy - mean(yy))^2)
  yy_mom3 / yy_mom2^(3/2)
}

tbls <- vector("list", length = kk)
tbl <- array(NA, dim = c(length(ne_sizes), 5))

KK <- max(dat$par_label)
tbls_skew <- array(NA, dim = c(length(ne_sizes), kk, KK))

for (j in 1:kk) {
  
  valid_dat <- dat_list[[j]]$valid_dat
  valid_XX <- dat_list[[j]]$valid_XX
  valid_par_labels <- dat_list[[j]]$valid_par_labels
  
  for (k in seq_along(ne_sizes)) {
    
    set.seed(42)
    
    nnmp_out <- out_list[[j]][[k]]
    nne <- ne_sizes[k]
    
    cat(paste0("Start the ", j, "-th fold prediction with neighborhood size ", k, "\n\n"))
    
    # Prediction
    nnmp_test <- predict(nnmp_out, valid_XX, valid_dat[,1:2], nonref_par_labels = valid_par_labels, predict_sam = TRUE)
    
    # RMSPE
    nnmp_rmspe <- sqrt(mean((nnmp_test$obs_mu - valid_dat[,3])^2))
    
    # MAE
    nnmp_mae <- mean(abs(nnmp_test$obs_mu - valid_dat[,3]))
    
    # CRPS
    nnmp_crps <- mean(scoringRules::crps_sample(valid_dat[,3], nnmp_test$sam))
    
    # 95% CI width
    nnmp_qq <- nnmp_test[['obs_qq']]
    nnmp_width <- mean(nnmp_qq[3,] - nnmp_qq[1,])
    
    # 95% CI Coverage rate
    test_num <- nrow(valid_dat)
    nnmp_cover_boolen <- sapply(1:(test_num), function(x) {
      valid_dat[x,3] >= nnmp_qq[1,x] & valid_dat[x,3] <= nnmp_qq[3,x]
    })
    nnmp_cover <- sum(nnmp_cover_boolen) / length(nnmp_cover_boolen)
    
    # Put together 
    nnmp_metric <- c(nnmp_rmspe, nnmp_mae, nnmp_crps, nnmp_cover, nnmp_width)
    tbl[k, ] <- nnmp_metric
    
    # Skewness
    dat_coef_skew_par <- array(NA, dim = KK)
    for (l in 1:KK) {
      sst_par_l <- valid_dat[valid_par_labels == l, 3]
      dat_coef_skew_par[l] <- coef_skew(sst_par_l)
    }
    
    nnmp_skew_score <- array(NA, dim = KK)
    for (l in 1:KK) {
      nnmp_coef_skew_sam <- apply(nnmp_test$sam[valid_par_labels==l, ], 2, coef_skew)
      tbls_skew[k, j, l]  <- mean(abs(nnmp_coef_skew_sam - dat_coef_skew_par[l]))
    }
    
  }
  
  tbls[[j]] <- tbl
  
}

# Create tables
tbl_rmse <- array(NA, dim = c(length(ne_sizes), kk))
tbl_mae <- array(NA, dim = c(length(ne_sizes), kk))
tbl_crps <- array(NA, dim = c(length(ne_sizes), kk))
tbl_ci <- array(NA, dim = c(length(ne_sizes), kk))
tbl_ciw <- array(NA, dim = c(length(ne_sizes), kk))

for (j in 1:kk) {
  tbl_rmse[,j] <- tbls[[j]][,1]
  tbl_mae[,j] <- tbls[[j]][,2]
  tbl_crps[,j] <- tbls[[j]][,3]
  tbl_ci[,j] <- tbls[[j]][,4]
  tbl_ciw[,j] <- tbls[[j]][,5]
}

tbl_rmse_mu <- rowMeans(tbl_rmse)
tbl_mae_mu <- rowMeans(tbl_mae)
tbl_crps_mu <- rowMeans(tbl_crps)
tbl_ci_mu <- rowMeans(tbl_ci)
tbl_ciw_mu <- rowMeans(tbl_ciw)

tbl_mu <- cbind(tbl_rmse_mu, tbl_mae_mu, tbl_crps_mu, tbl_ci_mu, tbl_ciw_mu)
rownames(tbl_mu) <- paste0("L = ", ne_sizes)
colnames(tbl_mu) <- c("RMSPE", "MAE", "CRPS", "95% CI cover", "95% CI width")
tbl_latex <- knitr::kable(tbl_mu, format = "latex", digits = 3, align = "c")
readr::write_file(tbl_latex, paste0(outpath, "sst_ext_sgnnmp_cv_pred.txt"))

tbl_skew_mu <- array(NA, dim = c(length(ne_sizes), KK))
for (l in 1:KK) {
  tbl_skew_mu[, l] <- rowMeans(tbls_skew[, , l])
}
rownames(tbl_skew_mu) <- paste0("L = ", ne_sizes)
colnames(tbl_skew_mu) <- paste0("Sub-basin-", 1:KK)
tbl_skew_latex <- knitr::kable(tbl_skew_mu, format = "latex", digits = 3, align = "c")
readr::write_file(tbl_skew_latex, paste0(outpath, "sst_ext_sgnnmp_cv_skew.txt"))
