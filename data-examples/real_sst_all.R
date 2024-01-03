############################################################################
### SST data analysis using the extended Skew-GNNMP
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
# Import, process and plot data
#-------------------------------------------------------------------------------
# Import and process data
rawdat <- read_csv(paste0(datapath, "med_sst_raw.csv"))
rawdat <- rawdat %>% 
    arrange(longitude, latitude) %>%
    group_by(longitude, latitude) %>%
    mutate(sst_new = median(sst_celsius)) 
ind <- duplicated(rawdat[, 1:2])
med_dat <- rawdat[!ind, c(1,2,4)]
colnames(med_dat) <- c("lon", "lat", "sst")
med_dat <- med_dat %>% arrange(lon, lat)

# Mediterranean sea
world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
sf::sf_use_s2(FALSE)

bbox <- c(minlon = -6.5, maxlon = 37, minlat = 29.5, maxlat = 46.5)
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

# Partition the area according to basins
ICES_shape <- sf::st_read(paste0(datapath, "ICES_ecoregions/ICES_ecoregions_20171207_erase_ESRI.shp"))
med_area_shape <- ICES_shape %>% filter(OBJECTID %in% c(2, 4, 5, 6, 7, 8, 17))
med_trhh_shape <- sf::st_read(paste0(datapath,"med_trhh/iho.shp"))

med_west_shape <- ICES_shape %>% filter(Ecoregion == "Western Mediterranean Sea" )
med_adr_shape <- ICES_shape %>% filter(Ecoregion == "Adriatic Sea" )
med_ionian_shape <- ICES_shape %>% filter(Ecoregion == "Ionian Sea and the Central Mediterranean Sea" )

med_trhh_shape_sp <- sf::as_Spatial(med_trhh_shape)
med_west_shape_sp <- sf::as_Spatial(med_west_shape)
med_adr_shape_sp <- sf::as_Spatial(med_adr_shape)
med_ionian_shape_sp <- sf::as_Spatial(med_ionian_shape)

med_trhh_map <- maps::map(med_trhh_shape_sp, fill = TRUE, plot = FALSE)
med_west_map <- maps::map(med_west_shape_sp, fill = TRUE, plot = FALSE)
med_adr_map <- maps::map(med_adr_shape_sp, fill = TRUE, plot = FALSE)
med_ionian_map <- maps::map(med_ionian_shape_sp, fill = TRUE, plot = FALSE)

dat$in_trhh_indx <- !is.na(maps::map.where(med_trhh_map, dat$lon, dat$lat))
dat$in_west_indx <- !is.na(maps::map.where(med_west_map, dat$lon, dat$lat))
dat$in_west_most_indx <- dat$in_west_indx & !dat$in_trhh_indx
dat$in_adr_index <- !is.na(maps::map.where(med_adr_map, dat$lon, dat$lat))
dat$in_ionian_index <- !is.na(maps::map.where(med_ionian_map, dat$lon, dat$lat)) | dat$lon < 21

dat <- dat %>% mutate(area = ifelse(in_trhh_indx, "Tyrrhenian Sea", 
                             ifelse(in_west_most_indx, "Westernmost",
                                    ifelse(in_adr_index, "Adriatic Sea", 
                                           ifelse(in_ionian_index, "Ionian Sea", "Aegean-Levantine Sea")))))
dat <- dat %>% mutate(par_label = ifelse(in_trhh_indx, 2,
                                    ifelse(in_west_most_indx, 1,
                                           ifelse(in_adr_index, 3,
                                                  ifelse(in_ionian_index, 4, 5)))))

# Plot histograms of residuals
ss <- lm(sst ~ lon + lat, data = dat)
dat$resid <- ss$residuals

png(paste0(outpath, "med_partition_resid_hist", ".png"), width=1000, height=600, pointsize=20, bg = "transparent")
par(mfrow = c(2, 3))
hist(dat$resid[dat$area == "Westernmost"], xlab = "Residuals", main = "Westernmost Mediterranean Sea", probability = TRUE)
hist(dat$resid[dat$area == "Tyrrhenian Sea"], xlab = "Residuals", main = "Tyrrhenian Sea", probability = TRUE)
hist(dat$resid[dat$area == "Adriatic Sea"], xlab = "Residuals", main = "Adriatic Sea", probability = TRUE)
hist(dat$resid[dat$area == "Ionian Sea"], xlab = "Residuals", main = "Ionian Sea", probability = TRUE)
hist(dat$resid[dat$area == "Aegean-Levantine Sea"], xlab = "Residuals", main = "Aegean-Levantine Sea", probability = TRUE)
dev.off()

# Plot observations
plot_dat_save <- 
  ggplot() + geom_sf(data = med_shape) + theme_bw() +
  coord_sf(xlim = c(bbox["minlon"], bbox["maxlon"]), ylim = c(bbox["minlat"], bbox["maxlat"]), expand = FALSE) + 
  geom_point(aes(lon, lat, color = sst), data = dat, size = 2) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(color = "", x = "Longitude", y = "Latitude") + 
  scale_colour_gradientn(colours = myPalette(100)) +
  theme(axis.text = element_text(size = 20),
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
png(paste0(outpath, "med_sst", ".png"), width=1000, height=550, pointsize=20, bg = "transparent")
print(plot_dat_save)
dev.off()


#-------------------------------------------------------------------------------
# Fit extended skew-Gaussian NNMP
#-------------------------------------------------------------------------------
set.seed(42)

# Use all observations as referenece set
ref_dat <- dat %>% select(lon, lat, sst) %>% as.matrix
ref_XX <- cbind(1, ref_dat[, 1:2])

# Random ordering
ord <- sample(1:nrow(ref_dat), replace = FALSE)

# MCMC
mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 4, n_report = 2000)

par_labels <- as.numeric(dat$par_label)

K <- max(par_labels)

priors <- list("sigmasq_invgamma" = c(3, 1),
               "la_normal" = list(mu = rep(0, K), sigmasq = rep(5, K)))

starting <- list("la" = rep(0, K), "sigmasq" = 5, 
                 "regcoef" = c(mean(ref_dat[,3]), rep(0, ncol(ref_XX) - 1)))

tuning <- list(
  sce1 = list("phi" = 0.12, "zeta" = 0.6, "sigmasq" = 0.1,
              "regcoef" = c(0.3, 0.02, 0.008),
              "la" = c(0.35, 0.7, 1.2, 0.5, 0.7)),
  sce2 = list("phi" = 0.12, "zeta" = 0.6, "sigmasq" = 0.1,
              "regcoef" = c(0.3, 0.02, 0.008), 
              "la" = c(0.35, 0.6, 1.1, 0.5, 0.7)),
  sce3 = list("phi" = 0.12, "zeta" = 0.7, "sigmasq" = 0.1,
              "regcoef" = c(0.3, 0.02, 0.008),
              "la" = c(0.35, 0.55, 0.9, 0.5, 0.6)))

ne_sizes <- c(10, 15, 20)
nnmp_out_list <- vector("list", length = length(ne_sizes))


for (i in seq_along(ne_sizes)) {
  
  nne <- ne_sizes[i]
  nnmp_out <- nnmp(response = ref_dat[,3],
                   covars = ref_XX,
                   coords = ref_dat[,1:2],
                   neighbor_size = nne,
                   marg_family = "sn",
                   priors = priors,
                   starting = starting,
                   tuning = tuning[[i]],
                   ord = ord,
                   mcmc_settings = mcmc_settings, 
                   verbose = TRUE,
                   par_label = par_labels)

  nnmp_out_list[[i]] <- nnmp_out

}

save(nnmp_out_list, file = paste0(outpath, "nnmp_out_la_10_15_20", ".RData"))

#-------------------------------------------------------------------------------
# Extended skew-Gaussian NNMP prediction 
#-------------------------------------------------------------------------------
set.seed(42)

# Create a grid
lon_grid <- seq(bbox["minlon"], bbox["maxlon"], length = 150)
lat_grid <- seq(bbox["minlat"], bbox["maxlat"], length = 150)
grid_loc <- expand.grid(lon_grid, lat_grid)
colnames(grid_loc) <- c("lon", "lat")

# Filter points and add partition labels
med_area_shape_sp <- sf::as_Spatial(med_area_shape)
med_area_map <- maps::map(med_area_shape_sp, fill = TRUE, plot = FALSE)
in_area_idx <- !is.na(maps::map.where(med_area_map, grid_loc$lon, grid_loc$lat))
grid_loc <- grid_loc[in_area_idx, ]
in_med_indx <- !is.na(maps::map.where(med_map, grid_loc$lon, grid_loc$lat))
grid_loc <- grid_loc[in_med_indx, ]

grid_loc$in_trhh_indx <- !is.na(maps::map.where(med_trhh_map, grid_loc$lon, grid_loc$lat)) 
grid_loc$in_west_indx <- !is.na(maps::map.where(med_west_map, grid_loc$lon, grid_loc$lat))
grid_loc$in_west_most_indx <- grid_loc$in_west_indx & !grid_loc$in_trhh_indx
grid_loc$in_adr_index <- !is.na(maps::map.where(med_adr_map, grid_loc$lon, grid_loc$lat))
grid_loc$in_ionian_index <- !is.na(maps::map.where(med_ionian_map, grid_loc$lon, grid_loc$lat)) 

grid_loc <- grid_loc %>% mutate(par_label = ifelse(in_trhh_indx, 2,
                                                   ifelse(in_west_most_indx, 1,
                                                          ifelse(in_adr_index, 3,
                                                                 ifelse(in_ionian_index, 4, 5)))))
grid_loc$par_label <- as.factor(grid_loc$par_label)

# Filter points in Med sea 
in_sea_idx <- is.na(maps::map.where("world", grid_loc$lon, grid_loc$lat))
grid_loc <- grid_loc[in_sea_idx, ]

# Covariate
grid_par_labels <- as.numeric(grid_loc$par_label)
grid_loc <- as.matrix(grid_loc[,1:2])
grid_XX <- cbind(1, grid_loc[, 1:2])

# Prediction
nnmp_pred_list <- vector("list", length = length(ne_sizes))

for (i in seq_along(ne_sizes)) {
  
  nnmp_pred <- predict(nnmp_out, grid_XX, grid_loc, nonref_par_labels = grid_par_labels)
  nnmp_pred_list[[i]] <- nnmp_pred
  
}
  
save(nnmp_pred_list, file = paste0(outpath, "nnmp_pred.RData"))

for (i in seq_along(ne_sizes)) {
  
  nne <- ne_sizes[i]
  nnmp_pred <- nnmp_pred_list[[i]]
  
  nnmp_surf <- as.data.frame(cbind(grid_loc, nnmp_pred$obs_mu, t(nnmp_pred$obs_qq)))
  colnames(nnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")
  
  plot_save_gnnmp <- ggplot() + geom_sf(data = med_shape) + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = q50), nnmp_surf, interpolate = T) +
    scale_fill_gradientn(colours = myPalette(100), limits = range(dat$sst)) + 
    labs(fill = "", x = "Longitude", y = "Latitude") + 
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.key.height= unit(2, 'cm'),
          axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"))
  png(paste0(outpath, "sst_sgnnmp_median_", nne, ".png"), width=1000, height=600, pointsize=20, bg = "transparent")
  print(plot_save_gnnmp)
  dev.off()
  
  plot_save_gnnmp <- ggplot() + geom_sf(data = med_shape) + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = q975 - q025), nnmp_surf, interpolate = T) +
    scale_fill_gradientn(colours = myPalette(100), limits = c(0, 12)) + 
    labs(fill = "", x = "Longitude", y = "Latitude") + 
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.key.height= unit(2, 'cm'),
          axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"))
  png(paste0(outpath, "sst_sgnnmp_95ci_width_", nne, ".png"), width=1000, height=600, pointsize=20, bg = "transparent")
  print(plot_save_gnnmp)
  dev.off()
  
}
 
  
