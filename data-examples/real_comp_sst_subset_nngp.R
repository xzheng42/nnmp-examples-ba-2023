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
library(spNNGP)


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
# Fit NNGP
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
mcmc_settings <- list(niter = 120000, nburn = 80000, nthin = 10)
starting <- list("phi"= 1, "sigma.sq"= 1, "tau.sq"= 1, "nu" = 0.5)
tuning <- list("phi" = 0.45, "nu" = 0.01)
priors <- list("phi.Unif"= c(3/30, 3/0.1), "nu.Unif" = c(0.1, 2),
               "sigma.sq.IG"= c(2, 1), "tau.sq.IG"= c(2, 1))

library(parallel)
library(doSNOW)
library(foreach)

closeAllConnections()
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

cat("Start parallel computing\n\n")
runtime <- system.time(nngp_out_list <- foreach(j = 1:kk, .packages = "spNNGP") %:%
                         foreach(k = seq_along(ne_sizes), .packages = "spNNGP") %dopar% {
                           set.seed(42)
                           source("~/path/to/run_nngp_mcmc.R")
                           run_mcmc(dat_list[[j]], ne_sizes[k], priors, starting, tuning, mcmc_settings)
                         })

cat("done\n")
cat(paste0("Running time: ", round(runtime[3]/60), " minutes\n"))

stopCluster(cl)


#-------------------------------------------------------------------------------
# Evaluate NNGP
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
    
    nngp_out <- nngp_out_list[[j]][[k]]
    
    cat(paste0("Start the ", j, "-th fold prediction with neighborhood size ", k + 9, "\n\n"))
    
    nngp_test <- predict(nngp_out, valid_XX, valid_dat[,1:2],
                         sub.sample = list(start = mcmc_settings$nburn + 1,
                                           end = mcmc_settings$niter,
                                           thin = mcmc_settings$nthin))
    
    nngp_pred <- rowMeans(nngp_test$p.y.0, 1)
    nngp_q025 <- apply(nngp_test$p.y.0, 1, function(x) quantile(x, probs = 0.025))
    nngp_q975 <- apply(nngp_test$p.y.0, 1, function(x) quantile(x, probs = 0.975))
    
    # RMSPE
    nngp_rmspe <- sqrt(mean((nngp_pred - valid_dat[,3])^2))
    
    # CRPS
    nngp_crps <- mean(scoringRules::crps_sample(valid_dat[,3], nngp_test$p.y.0))
    
    # 95% CI Coverage rate
    nngp_cover_boolen <- sapply(1:valid_num, function(x) {
      valid_dat[x,3] >= nngp_q025[x] & valid_dat[x,3] <= nngp_q975[x]
    })
    nngp_cover <- mean(nngp_cover_boolen)
    
    # 95% CI width
    nngp_width <- mean(nngp_q975 - nngp_q025)
    
    # PPLC and DIC
    nngp_diag <- spDiag(nngp_out, sub.sample = list(start = mcmc_settings$nburn + 1,
                                                    end = mcmc_settings$niter,
                                                    thin = mcmc_settings$nthin))
    
    # Put together
    nngp_metric <- c(nngp_rmspe, nngp_crps, nngp_cover, nngp_width, 
                     unlist(nngp_diag$GP)[3], unlist(nngp_diag$DIC)[1])
    tbl[k, ] <- nngp_metric
    
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
tbl_crps_mu <- rowMeans(tbl_crps)
tbl_ci_mu <- rowMeans(tbl_ci)
tbl_ciw_mu <- rowMeans(tbl_ciw)
tbl_pplc_mu <- rowMeans(tbl_pplc)
tbl_dic_mu <- rowMeans(tbl_dic)

tbl_mu <- cbind(tbl_rmse_mu, tbl_crps_mu, tbl_ci_mu, tbl_ciw_mu, 
                tbl_pplc_mu, tbl_dic_mu)
rownames(tbl_mu) <- paste0("L = ", 10:20)
colnames(tbl_mu) <- c("RMSPE", "CRPS","95% CI cover", "95% CI width", "PPLC", "DIC")
tbl_latex <- knitr::kable(tbl_mu, format = "latex", digits = 3, align = "c")
write_file(tbl_latex, paste0(outpath, "sst_subset_nngp_cv.txt"))

