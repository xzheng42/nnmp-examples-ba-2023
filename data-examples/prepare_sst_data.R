################################################################################
### Mediterranean Sea partitions
################################################################################
rm(list = ls())

outpath <- "/path/to/output/"
datapath <- "/path/to/input/data/"

library(ggplot2)
library(readr)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

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

# Write data
write_csv(dat[, c("lon", "lat", "sst", "area", "par_label")], file = paste0(datapath, "med_sst.csv"))

#-------------------------------------------------------------------------------
# Figure 5(a) of the main paper
#-------------------------------------------------------------------------------
# Create a grid
lon_grid <- seq(bbox["minlon"], bbox["maxlon"], length = 200)
lat_grid <- seq(bbox["minlat"], bbox["maxlat"], length = 200)
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

# Plot partitions
plot_dat_save <- 
  ggplot() + geom_sf(data = med_shape) + 
  geom_point(aes(lon, lat, color = par_label), data = grid_loc, size = 2) + theme_bw() +
  annotate(geom = "text", x = 3, y = 38, label = "Westernmost Mediterranean Sea", fontface = "italic", color = "grey22", size = 7) + 
  annotate(geom = "text", x = 12.5, y = 40, label = "Tyrrhenian Sea", fontface = "italic", color = "grey22", size = 7) + 
  annotate(geom = "text", x = 15.5, y = 43, label = "Adriatic Sea", fontface = "italic", color = "grey22", size = 7) + 
  annotate(geom = "text", x = 17, y = 35, label = "Ionian Sea", fontface = "italic", color = "grey22", size = 7) + 
  annotate(geom = "text", x = 29.5, y = 33.5, label = "Aegean-Levantine Sea", fontface = "italic", color = "grey22", size = 7) + 
  coord_sf(xlim = c(bbox["minlon"], bbox["maxlon"]), ylim = c(bbox["minlat"], bbox["maxlat"]), expand = FALSE) + 
  labs(color = "", x = "Longitude", y = "Latitude") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        legend.key.height= unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "med_grid_parition", ".png"), width=1000, height=600, pointsize=20, bg = "transparent")
print(plot_dat_save)
dev.off()
