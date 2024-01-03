############################################################################
### SIMULATION STUDY: NNMP with skew-Gaussian marginals
### Data generated from a skew-Gaussian random field
### by Zhang and El-Shaarawi (2010)
############################################################################

rm(list = ls())

outpath <- "/path/to/output/"

library(SpatialExtremes)
library(ggplot2)

library(nnmp)


#-------------------------------------------------------------------------------
# Generate the true spatial random field
#-------------------------------------------------------------------------------
set.seed(42)

# Simulate from Gaussian processes on a grid
x_seq <- seq(0, 1, length = 100)
y_seq <- seq(0, 1, length = 100)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
zz1 <- abs(t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 0.1, smooth = 1)))
zz2 <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 0.1, smooth = 1))

skewparam <- c(1, 2.5, 5)
dat_save <- vector("list", length = length(skewparam))
for (i in seq_along(skewparam)) {
  
  # Simulate from a skew-Gaussian process with the i-th skewness parameter
  yy <- skewparam[i] * zz1  + zz2
  
  # Plot the true field
  dat <- cbind(grid_loc, yy)
  colnames(dat) <- c("lon", "lat", "yy")
  dat_save[[i]] <- dat

  plot_save <- ggplot(as.data.frame(dat), aes(x = lon, y = lat, fill = yy)) +
  geom_raster(interpolate = T) + theme_bw() +
  scale_fill_gradient2(high = "blue", low = "red", mid = "white", limits = range(yy), midpoint = 0,
                       space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
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
  png(paste0(outpath, "sn_field_la_", skewparam[i], ".png"),
             width=600, height=500, pointsize=20, bg = "transparent")
  print(plot_save)
  dev.off()

}


#-------------------------------------------------------------------------------
# Fit skew-Gaussian NNMP
#-------------------------------------------------------------------------------
mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 10, n_report = 2000)             

priors <- list(
  sce1 = list("la_normal" = c(0, 5), "sigmasq_invgamma" = c(3, 1)),
  sce2 = list("la_normal" = c(1, 5), "sigmasq_invgamma" = c(3, 1)),
  sce3 = list("la_normal" = c(1, 5), "sigmasq_invgamma" = c(3, 1))
)                                  

starting <- list(sce1 = list("la" = 0, "sigmasq" = 5), 
                 sce2 = list("la" = 1, "sigmasq" = 5), 
                 sce3 = list("la" = 1, "sigmasq" = 5))  

tuning <- list(
  sce1 = list("phi" = 0.12, "zeta" = 0.3, "la" = 0.35, "sigmasq" = 0.1),
  sce1 = list("phi" = 0.15, "zeta" = 0.5, "la" = 0.55, "sigmasq" = 0.1),
  sce3 = list("phi" = 0.16, "zeta" = 0.7, "la" = 0.95, "sigmasq" = 0.1)
)   

nne <- 10

nnmp_out_list <- vector("list", length = length(skewparam))
test_dat_save <- vector("list", length = length(skewparam))

for (i in seq_along(skewparam)) {

  dat <- dat_save[[i]]

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
  test_dat_save[[i]] <- obs_dat[-ref_idx,]

  # Random ordering
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
  
save(nnmp_out_list, file = paste0(outpath, "sim_sn.RData"))


#-------------------------------------------------------------------------------
# Plot stationary marginal distribution
#-------------------------------------------------------------------------------
bandwidth <- c(1, 1.5, 2.5)
grids <- rbind(seq(-5, 10, length = 2000), 
               seq(-5, 10, length = 2000), 
               seq(-5, 20, length = 2000))

for (i in seq_along(skewparam)) {
  
  post_sams <- nnmp_out_list[[i]]$post_samples
  
  
  data_grid <- grids[i, ]
  true_dens <- sn::dsn(data_grid, 0, sqrt(1 + skewparam[[i]]^2), skewparam[[i]])
  dens_sam <- sapply(data_grid, function(x) 
       sn::dsn(x, 0, sqrt(post_sams$la^2 + post_sams$sigmasq), post_sams$la/sqrt(post_sams$sigmasq)))
  dens_mean <- colMeans(dens_sam)
  dens_qq <- apply(dens_sam, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  dens <- cbind(data_grid, true_dens, nnmp_out_list[[i]]$orig_data$response, dens_mean, t(dens_qq))
  colnames(dens) <- c("grid", "true_dens", "dat", "dens_mu", "dens_qq_025", "dens_qq_975")
  plot_hist_save <- ggplot(as.data.frame(dens), aes(x = dat)) + theme_bw() + 
      ylim(c(0, 0.4)) + 
      geom_histogram(aes(y=..density..), binwidth = bandwidth[i], fill="white", color="darkgrey", size = 1) +
      geom_line(aes(x = grid, y = true_dens, color = "TRUE"), linetype = "dashed", size = 1.5) +  
      geom_ribbon(aes(x = grid, ymin = dens_qq_025, ymax = dens_qq_975, fill = "Pointwise 95% CI"), alpha = 0.5) + 
      geom_line(aes(x = grid, y = dens_mean, color = "Poseterior mean"), linetype = "dashed", size = 1.5) + 
      scale_colour_manual("", values = c( "brown1", "green")) +
      scale_fill_manual("", values = "brown1") +
      labs(fill = "", x = "Data", y = "Density") + 
      theme(legend.position="none",
            axis.text = element_text(size = 25),
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
  png(paste0(outpath, "sn_hist_la_", skewparam[i], ".png"), width=600, height=500, pointsize=20, bg = "transparent")
  print(plot_hist_save)
  dev.off()
}


#-------------------------------------------------------------------------------
# Skew-Gaussian NNMP prediction on a grid
#-------------------------------------------------------------------------------
for (i in seq_along(skewparam)) {
  
  nnmp_out <- nnmp_out_list[[i]]
  nonref_grid <- as.matrix(expand.grid(seq(0.005, 1, by = 0.01), seq(0.005, 1, by = 0.01)))
  nnmp_pred <- predict(nnmp_out, nonref_coords = nonref_grid, predict_sam = FALSE)

  nnmp_surf <- as.data.frame(cbind(nonref_grid, nnmp_pred[[1]], t(nnmp_pred[[2]])))
  colnames(nnmp_surf) <- c("x", "y", "mu", "q025", "q50", "q975")
  
  plot_save_nnmp <- ggplot(as.data.frame(nnmp_surf), aes(x = x, y = y, fill = q50)) +
    geom_raster(interpolate = T) + theme_bw() +
    scale_fill_gradient2(high = "blue", low = "red", mid = "white", 
                         limits = range(nnmp_out$orig_data$response), space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
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
  png(paste0(outpath, "sn_pred_la_", skewparam[i], ".png"), 
             width=600, height=500, pointsize=20, bg = "transparent")
  print(plot_save_nnmp)
  dev.off()

}
