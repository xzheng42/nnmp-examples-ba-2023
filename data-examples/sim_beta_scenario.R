############################################################################
### SIMULATION STUDY: copula NNMP with beta marignals
### Data generated from a spatial Gaussian copula with beta marginals
############################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(SpatialExtremes)
library(ggplot2)

library(nnmp)

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))


#-------------------------------------------------------------------------------
# Generate the true spatial random field
#-------------------------------------------------------------------------------
set.seed(32)

# Simulate from a Gaussian process on a grid
x_seq <- seq(0, 1, length = 100)
y_seq <- seq(0, 1, length = 100)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
zz <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 0.1, smooth = 1))

# Simulate proportions 
a_0 <- 3
b_0 <- 6
yy <- qbeta(pnorm(zz), a_0, b_0)

# Plot the true field
dat <- cbind(grid_loc, yy)
colnames(dat) <- c("lon", "lat", "yy")

plot_save <- ggplot() + theme_bw() + 
  geom_raster(aes(x = lon, y = lat, fill = yy), data = as.data.frame(dat), interpolate = T) + 
  scale_fill_gradientn(colours = myPalette(100), limits = c(0, 1)) + 
  labs(fill = "", x = "Easting", y = "Northing") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "beta_true_field", ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()


#-------------------------------------------------------------------------------
# Fit beta NNMP
#-------------------------------------------------------------------------------
set.seed(42)

# Randomly select part of the gridded data as observations
obs_num <- 2500
obs_idx <- sample(1:nrow(grid_loc), obs_num, FALSE)
obs_dat <- dat[obs_idx,]

# Partition the observations into training and testing sets
# Use the whole training set as reference set 
ref_num <- 2000
test_num <- obs_num - ref_num
ref_idx <- sample(1:obs_num, ref_num, FALSE)
ref_dat <- obs_dat[ref_idx,]
test_dat <- obs_dat[-ref_idx,]

# Random ordering
ord <- sample(1:ref_num, replace = FALSE)

# MCMC
priors <- list("shape1_gamma" = c(1, 1), "shape2_gamma" = c(1, 1))

starting <- list("shape1" = 1, "shape2" = 1)

tuning <- list("phi" = 0.11, "zeta" = 0.25, "shape1" = 0.1, "shape2" = 0.1) 

mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 10, n_report = 2000)

nne <- 10

nnmp_out <- nnmp(response = ref_dat[,3],
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


#-------------------------------------------------------------------------------
# Plot stationary marginal distribution
#-------------------------------------------------------------------------------
post_sams <- nnmp_out$post_samples
data_grid <- seq(0, 1, length = 2000) 
true_yy <- dbeta(data_grid, a_0, b_0)
dens_gaus_sam <- sapply(data_grid, function(x)  dbeta(x, post_sams$shape1, post_sams$shape2))
dens_gaus_mean <- colMeans(dens_gaus_sam)
dens_gaus_qq <- apply(dens_gaus_sam, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
dens <- cbind(data_grid, ref_dat[,3], true_yy, dens_gaus_mean, t(dens_gaus_qq))
colnames(dens) <- c("grid", "dat", "true", "dens_gaus_mu", "dens_gaus_qq_025", "dens_gaus_qq_975")
plot_hist_save <- ggplot(as.data.frame(dens), aes(x = dat)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.1, fill="white", color="darkgrey", size = 1) + 
  geom_line(aes(x = grid, y = true, color = "True"), linetype = "dashed", size = 1.5) +
  geom_ribbon(aes(x = grid, ymin = dens_gaus_qq_025, ymax = dens_gaus_qq_975, fill = "Gaus 95% CI"), alpha = 0.5) + 
  geom_line(aes(x = grid, y = dens_gaus_mean, color = "Gaus mean"), linetype = "dashed", size = 1.5) + 
  scale_fill_manual("", values = c("brown1")) +
  scale_colour_manual("", values = c("brown1", "green")) +
  labs(fill = "", x = "Data", y = "Density") +  
  theme_bw() + ylim(0, 3.5) +
  theme(plot.margin = unit(c(0.3, 1.5, 0.5, 1.5), "cm")) +
  theme(legend.position="none",
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "beta_marginals.png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_hist_save)
dev.off()


#-------------------------------------------------------------------------------
# Beta NNMP prediction on a grid
#-------------------------------------------------------------------------------
nonref_grid <- as.matrix(expand.grid(seq(0.005, 1, by = 0.01), seq(0.005, 1, by = 0.01)))
nnmp_pred <- predict(nnmp_out, nonref_coords = nonref_grid, predict_sam = FALSE)

nnmp_surf <- as.data.frame(cbind(nonref_grid, nnmp_pred[[1]], t(nnmp_pred[[2]])))
colnames(nnmp_surf) <- c("x", "y", "mu", "q005", "q50", "q95")
plot_save_nnmp <- ggplot(as.data.frame(nnmp_surf), aes(x = x, y = y, fill = q50)) +
  theme_bw() + geom_raster(interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = range(0,1)) + 
  labs(fill = "", x = "Easting", y = "Northing") + 
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "beta_pred_", "q50_", nne, ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save_nnmp)
dev.off()