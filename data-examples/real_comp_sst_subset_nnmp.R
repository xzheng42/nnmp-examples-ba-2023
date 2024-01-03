############################################################################
### Regional SST data analysis: Model comparison
############################################################################
rm(list = ls())

outpath <- "/path/to/output/"
datapath <- "/path/to/input/data/"

library(ggplot2)
library(readr)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

library(nnmp)


#-------------------------------------------------------------------------------
# Import and process data
#-------------------------------------------------------------------------------
rawdat <- read_csv(paste0(datapath, "med_sst_raw.csv"))
rawdat <- rawdat %>% 
  arrange(longitude, latitude) %>%
  group_by(longitude, latitude) %>%
  mutate(sst_new = median(sst_celsius))
ind <- duplicated(rawdat[, 1:2])
med_dat <- rawdat[!ind, c(1,2,4)]
colnames(med_dat) <- c("lon", "lat", "sst")

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
sf::sf_use_s2(FALSE)

bbox <- c(minlon = 0, maxlon = 9, minlat = 35.5, maxlat = 44.5)
base <- ggplot(data = world) + geom_sf() + 
  coord_sf(xlim = c(bbox["minlon"], bbox["maxlon"]), 
           ylim = c(bbox["minlat"], bbox["maxlat"]), 
           expand = FALSE)

data_in_box <- med_dat[which(med_dat$lon >= bbox["minlon"] & med_dat$lon <= bbox["maxlon"] & 
                               med_dat$lat >= bbox["minlat"] & med_dat$lat <= bbox["maxlat"]), ]
in_sea_idx <- is.na(maps::map.where("world", data_in_box$lon, data_in_box$lat))
dat <- data_in_box[in_sea_idx, ]

med_shape <- sf::st_read(paste0(datapath,"med_lme/lme.shp"))
med_shape_sp <- sf::as_Spatial(med_shape)
med_map <- maps::map(med_shape_sp, fill = TRUE, plot = FALSE)
in_med_indx <- !is.na(maps::map.where(med_map, dat$lon, dat$lat))
dat <- dat[in_med_indx, ]


#-------------------------------------------------------------------------------
# Fit Gaussian NNMP
#-------------------------------------------------------------------------------
set.seed(42)

kk <- 10
folds_idx <- sample(1:kk, nrow(dat), replace = TRUE)

ne_sizes <- 10:20
dat_list <- vector("list", length = kk)

for (j in 1:kk) {
  
  # Partition the observations into training and validation sets
  ref_num <- nrow(dat)
  selec_idx <- folds_idx != j
  selec_num <- sum(selec_idx)
  valid_num <- ref_num - selec_num
  
  ref_dat <- as.matrix(dat)
  ref_XX <- cbind(1, ref_dat[, 1:2])
  selec_dat <- ref_dat[selec_idx,]
  selec_XX <- ref_XX[selec_idx, ]
  valid_dat <- ref_dat[!selec_idx,]
  valid_XX <- ref_XX[!selec_idx, ]
  
  # Random ordering
  set.seed(42)
  ord <- sample(1:selec_num, replace = FALSE)
  
  dat_list[[j]] <- list(ref_dat = ref_dat, ref_XX = ref_XX, 
                        selec_dat = selec_dat, selec_XX = selec_XX, 
                        valid_dat = valid_dat, valid_XX = valid_XX,
                        ord = ord)
}

# MCMC
mcmc_settings <- list(n_iter = 120000, n_burn = 80000, n_thin = 10, n_report = 5000)
tuning <- list("phi" = 0.22, "zeta" = 0.3)

library(parallel)
library(doSNOW)
library(foreach)

closeAllConnections()
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

cat("Start parallel computing\n\n")
runtime <- system.time(nnmp_out_list <- foreach(j = 1:kk, .packages = "nnmp") %:%
                         foreach(k = seq_along(ne_sizes), .packages = "nnmp") %dopar% {
                           set.seed(42)
                           source("~/path/to/run_nnmp_mcmc.R")
                           run_mcmc(dat_list[[j]], ne_sizes[k], tuning, mcmc_settings)
                         })

cat("done\n")
cat(paste0("Running time: ", round(runtime[3]/60), " minutes\n"))

stopCluster(cl)


#-------------------------------------------------------------------------------
# Evaluate Gaussian NNMP
#-------------------------------------------------------------------------------
tbls <- vector("list", length = kk)
tbl <- array(NA, dim = c(length(ne_sizes), 6))
colnames(tbl) <- c("RMSPE", "CRPS", "95% CI cover", "95% CI width", "PPLC", "DIC")
rownames(tbl) <- paste0("L = ", ne_sizes)

for (j in 1:kk) {
  
  valid_dat <- dat_list[[j]]$valid_dat
  valid_XX <- dat_list[[j]]$valid_XX
  valid_num <- nrow(valid_dat)
  
  for (k in seq_along(ne_sizes)) {
    
    set.seed(42)
    
    nnmp_out <- nnmp_out_list[[j]][[k]]
    nne <- ne_sizes[k]
    
    cat(paste0("Start the ", j, "-th fold prediction with neighborhood size ", k + 9, "\n\n"))
    
    nnmp_test <- predict(nnmp_out, valid_XX, valid_dat[,1:2], predict_sam = TRUE)
    
    # RMSPE
    nnmp_rmspe <- sqrt(mean((as.numeric(nnmp_test$obs_mu) - valid_dat[,3])^2))
    
    # CRPS
    nnmp_crps <- mean(scoringRules::crps_sample(valid_dat[, 3], nnmp_test$obs_sam))
    
    # 95% CI Coverage rate
    nnmp_qq <- nnmp_test[['obs_qq']]
    nnmp_cover_boolen <- sapply(1:(valid_num), function(x) {
      valid_dat[x,3] >= nnmp_qq[1,x] & valid_dat[x,3] <= nnmp_qq[3,x]
    })
    nnmp_cover <- mean(nnmp_cover_boolen)
    
    # 95% CI width
    nnmp_width <- mean(nnmp_qq[3,] - nnmp_qq[1,])
    
    # Put together
    nnmp_metric <- c(nnmp_rmspe, nnmp_crps, nnmp_cover, nnmp_width,
                     nnmp_out$mod_diag$pplc[3], nnmp_out$mod_diag$dic)
    tbl[k, ] <- nnmp_metric
    
    cat(paste0("Done the ", j, "-th fold prediction with neighborhood size ", k + 9, "\n\n"))
    cat("------------------------------------------------------------------------\n\n")
    
  }
  
  tbls[[j]] <- tbl
  
}

# Create tables
tbl_rmse <- array(NA, dim = c(length(ne_sizes), kk))
tbl_crps <- array(NA, dim = c(length(ne_sizes), kk))
tbl_ci <- array(NA, dim = c(length(ne_sizes), kk))
tbl_ciw <- array(NA, dim = c(length(ne_sizes), kk))
tbl_pplc <- array(NA, dim = c(length(ne_sizes), kk))
tbl_dic <- array(NA, dim = c(length(ne_sizes), kk))

for (j in 1:kk) {
  tbl_rmse[,j] <- tbls[[j]][,1]
  tbl_crps[,j] <- tbls[[j]][,2]
  tbl_ci[,j] <- tbls[[j]][,3]
  tbl_ciw[,j] <- tbls[[j]][,4]
  tbl_pplc[,j] <- tbls[[j]][,5]
  tbl_dic[,j] <- tbls[[j]][,6]
}

tbl_rmse_mu <- rowMeans(tbl_rmse)
tbl_ci_mu <- rowMeans(tbl_ci)
tbl_ciw_mu <- rowMeans(tbl_ciw)
tbl_crps_mu <- rowMeans(tbl_crps)
tbl_pplc_mu <- rowMeans(tbl_pplc)
tbl_dic_mu <- rowMeans(tbl_dic)

tbl_mu <- cbind(tbl_rmse_mu, tbl_crps_mu, tbl_ci_mu, tbl_ciw_mu,
                tbl_pplc_mu, tbl_dic_mu)
rownames(tbl_mu) <- paste0("L = ", 10:20)
colnames(tbl_mu) <- c("RMSPE", "CRPS","95% CI cover", "95% CI width", "PPLC", "DIC")
tbl_latex <- knitr::kable(tbl_mu, format = "latex", digits = 3, align = "c")
write_file(tbl_latex, paste0(outpath, "sst_subset_nnmp_cv.txt"))
