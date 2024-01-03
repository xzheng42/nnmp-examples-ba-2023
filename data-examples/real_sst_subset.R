############################################################################
### Real data analysis 
### (Section C: regional SST data analysis)
### Produce results with neighborhood size L = 13
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

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

#-------------------------------------------------------------------------------
# Import, process, and plot data
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

bbox_med <- c(minlon = -6.5, maxlon = 37, minlat = 29.5, maxlat = 46.5)
data_in_box <- med_dat[which(med_dat$lon >= bbox_med["minlon"] & med_dat$lon <= bbox_med["maxlon"] & 
                               med_dat$lat >= bbox_med["minlat"] & med_dat$lat <= bbox_med["maxlat"]), ]
in_sea_idx <- is.na(maps::map.where("world", data_in_box$lon, data_in_box$lat))
med_dat <- data_in_box[in_sea_idx, ]

med_shape <- sf::st_read(paste0(datapath,"med_lme/lme.shp"))
med_shape_sp <- sf::as_Spatial(med_shape)
med_map <- maps::map(med_shape_sp, fill = TRUE, plot = FALSE)
in_med_indx <- !is.na(maps::map.where(med_map, med_dat$lon, med_dat$lat))
med_dat <- med_dat[in_med_indx, ]

bbox <- c(minlon = 0, maxlon = 9, minlat = 35.5, maxlat = 44.5)
data_in_box <- med_dat[which(med_dat$lon >= bbox["minlon"] & med_dat$lon <= bbox["maxlon"] & 
                               med_dat$lat >= bbox["minlat"] & med_dat$lat <= bbox["maxlat"]), ]
in_sea_idx <- is.na(maps::map.where("world", data_in_box$lon, data_in_box$lat))
dat <- data_in_box[in_sea_idx, ]

base <- ggplot() + geom_sf(data = med_shape) + 
  coord_sf(xlim = c(bbox["minlon"], bbox["maxlon"]), 
           ylim = c(bbox["minlat"], bbox["maxlat"]), 
           expand = FALSE)


plot_dat_save <- base + theme_bw() +
  geom_point(aes(lon, lat, color = sst), data = dat, size = 3) + 
  labs(color = "", x = "Longitude", y = "Latitude") + 
  scale_colour_gradientn(colours = myPalette(100), limits = c(min(med_dat[,3]), max(med_dat[,3]))) + 
  theme(
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
    axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
    legend.key.height= unit(2, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"), 
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent")) 
png(paste0(outpath, "selec_sst", ".png"), width=600, height=600, pointsize=20, bg = "transparent")  
print(plot_dat_save)
dev.off()

#-------------------------------------------------------------------------------
# Fit Gaussian NNMP
#-------------------------------------------------------------------------------
set.seed(42)

ref_num <- nrow(dat)
ref_idx <- sample(1:ref_num, ref_num, FALSE)
ref_dat <- as.matrix(dat[ref_idx, ])
ref_XX <- cbind(1, ref_dat[, 1:2])

# Random ordering
ord <- sample(1:ref_num, replace = FALSE)

# MCMC
tuning <- list("phi" = 0.22, "zeta" = 0.4)

mcmc_settings <- list(n_iter = 120000, n_burn = 80000, n_thin = 10, n_report = 2000)

nne <- 13

nnmp_out <- nnmp(response = ref_dat[,3],
                 covars = ref_XX,
                 coords = ref_dat[,1:2],
                 neighbor_size = nne,
                 marg_family = "gaussian",
                 tuning = tuning,
                 ord = ord,
                 mcmc_settings = mcmc_settings, 
                 verbose = TRUE)

#-------------------------------------------------------------------------------
# Gaussian NNMP prediction on a grid
#-------------------------------------------------------------------------------
lon_grid <- seq(bbox["minlon"], bbox["maxlon"], length = 150)
lat_grid <- seq(bbox["minlat"], bbox["maxlat"], length = 150)
whole_grid <- expand.grid(lon_grid, lat_grid)
colnames(whole_grid) <- c("lon", "lat")
grid_in_box <- whole_grid[which(whole_grid$lon >= bbox["minlon"] & whole_grid$lon <= bbox["maxlon"] & 
                          whole_grid$lat >= bbox["minlat"] & whole_grid$lat <= bbox["maxlat"]), ]
grid_in_sea_idx <- is.na(maps::map.where("world", grid_in_box$lon, grid_in_box$lat))                          
grid_loc <- grid_in_box[grid_in_sea_idx, ]
in_med_indx <- !is.na(maps::map.where(med_map, grid_loc$lon, grid_loc$lat))
grid_loc <- as.matrix(grid_loc[in_med_indx, ])
grid_XX <- cbind(1, grid_loc[, 1:2])

nnmp_pred <- predict(nnmp_out, grid_XX, grid_loc, predict_sam = FALSE)
nnmp_surf <- as.data.frame(cbind(grid_loc, nnmp_pred$obs_mu, t(nnmp_pred$obs_qq)))
colnames(nnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")
nnmp_lat_surf <- as.data.frame(cbind(grid_loc, nnmp_pred$zz_mu, t(nnmp_pred$zz_qq)))
colnames(nnmp_lat_surf) <- c("x", "y", "mu", "q025", "q50", "q975")

plot_save_gnnmp <- base + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = q50), nnmp_surf, interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = c(min(med_dat[,3]), max(med_dat[,3]))) + 
  labs(fill = "", x = "Longitude", y = "Latitude") + 
  theme(
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "sst_subset_nnmp_obs_median_", nne, ".png"), width=600, height=600, pointsize=20, bg = "transparent")  
print(plot_save_gnnmp)
dev.off()

plot_save_gnnmp <- base + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = q50), nnmp_lat_surf, interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = c(-4, 4)) + 
  labs(fill = "", x = "Longitude", y = "Latitude") + 
  theme(
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "sst_subset_nnmp_spe_median_", nne, ".png"), width=600, height=600, pointsize=20, bg = "transparent")  
print(plot_save_gnnmp)
dev.off()  

#-------------------------------------------------------------------------------
# Fit NNGP with matern covariance function
#-------------------------------------------------------------------------------
library(spNNGP)

set.seed(42)

starting <- list("phi" = 1, "sigma.sq" = 1, "tau.sq" = 1, "nu" = 0.5)
tuning <- list("phi" = 0.5, "nu" = 0.01)
priors <- list("phi.Unif"= c(3/30, 3/0.1), "nu.Unif" = c(0.1, 2),
               "sigma.sq.IG"= c(2, 1), "tau.sq.IG"= c(2, 1))

nne <- 13

nngp_out <- spNNGP(formula = ref_dat[,3] ~ ref_XX[, 2:3], 
                   coords = ref_dat[,1:2],
                   starting = starting, 
                   tuning = tuning, 
                   ord = ord,
                   priors = priors, 
                   cov.model = "matern",
                   n.samples = mcmc_settings$n_iter,
                   n.neighbors = nne,
                   method = "latent", 
                   n.omp.threads = 1, 
                   n.report = 5000)

save(nngp_out, file = paste0(outpath, "nngp_matern_out_nne_13.RData"))

#-------------------------------------------------------------------------------
# NNGP prediction on a grid
#-------------------------------------------------------------------------------
nngp_pred <- predict(nngp_out, grid_XX, grid_loc, 
                     sub.sample = list(start = mcmc_settings$n_burn + 1, 
                                       end = mcmc_settings$n_iter, 
                                       thin = mcmc_settings$n_thin),
                     n.report = 1000)
nngp_surf <- as.data.frame(cbind(grid_loc, apply(nngp_pred$p.y.0, 1, median)))
nngp_lat_surf <- as.data.frame(cbind(grid_loc, apply(nngp_pred$p.w.0, 1, median)))
colnames(nngp_surf) <- c("x", "y", "z")
colnames(nngp_lat_surf) <- c("x", "y", "z")

plot_save_ngpp <- base + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = z), nngp_surf, interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = c(min(med_dat[,3]), max(med_dat[,3]))) + 
  labs(fill = "", x = "Longitude", y = "Latitude") + 
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "sst_subset_nngp_matern_obs_median_", nne, ".png"), width=600, height=600, pointsize=20, bg = "transparent")  
print(plot_save_ngpp)
dev.off()

plot_save_ngpp <- base + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = z), nngp_lat_surf, interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = c(-4, 4)) + 
  labs(fill = "", x = "Longitude", y = "Latitude") + 
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "sst_subset_nngp_matern_spe_median_", nne, ".png"), width=600, height=600, pointsize=20, bg = "transparent")  
print(plot_save_ngpp)
dev.off()


