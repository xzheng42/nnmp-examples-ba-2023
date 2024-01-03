############################################################################
### SIMULATION STUDY: copula NNMP with gamma marignals
### Data generated from a spatial Student-t copula with gamma marginals
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
set.seed(42)

# Simulate from a Student-t process on a grid
x_seq <- seq(0, 1, length = 100)
y_seq <- seq(0, 1, length = 100)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
gp <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget =  0, sill = 1, range = 1/12, smooth = 1))
dof <- 10
zz <- gp * sqrt(dof / rchisq(1, dof)) 

# Simulate proportions 
a_0 <- 2
b_0 <- 2
yy <- qgamma(pt(zz, dof), a_0, b_0)

# Plot the true field
dat <- cbind(grid_loc, yy)
colnames(dat) <- c("lon", "lat", "yy")

plot_save <- ggplot() + theme_bw() + 
  geom_raster(aes(x = lon, y = lat, fill = yy), data = as.data.frame(dat), interpolate = T) + 
  scale_fill_gradientn(colours = myPalette(100)) +
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
png(paste0(outpath, "true_st_gamma", ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()


#-------------------------------------------------------------------------------
# Fit gamma NNMP
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

tuning <- list(gaus = list("phi" = 0.1, "zeta" = 0.18, "shape1" = 0.13, "shape2" = 0.11), 
               gum = list("phi" = 0.08, "zeta" = 0.16, "shape1" = 0.1, "shape2" = 0.07))

mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 10, n_report = 2000)

nne <- 10

nnmp_gaus_out <- nnmp(response = ref_dat[,3],
                      coords = ref_dat[,1:2],
                      neighbor_size = nne,
                      marg_family = "gamma",
                      cop_family = "gaussian",
                      priors = priors,
                      starting = starting,
                      tuning = tuning[["gaus"]],
                      ord = ord,
                      mcmc_settings = mcmc_settings, 
                      verbose = TRUE)

nnmp_gum_out <- nnmp(response = ref_dat[,3],
                     coords = ref_dat[,1:2],
                     neighbor_size = nne,
                     marg_family = "gamma",
                     cop_family = "gumbel",
                     priors = priors,
                     starting = starting,
                     tuning = tuning[["gum"]],
                     ord = ord,
                     mcmc_settings = mcmc_settings, 
                     verbose = TRUE)


#-------------------------------------------------------------------------------
# Plot stationary marginal distribution
#-------------------------------------------------------------------------------
gaus_post_sams <- nnmp_gaus_out$post_samples
gum_post_sams <- nnmp_gum_out$post_samples
data_grid <- seq(0, range(ref_dat[,3])[2]+1, length = 2000) 
true_yy <- dgamma(data_grid, 2, 2)
dens_gaus_sam <- sapply(data_grid, function(x)  dgamma(x, gaus_post_sams$shape1, gaus_post_sams$shape2))
dens_gum_sam <- sapply(data_grid, function(x)  dgamma(x, gum_post_sams$shape1, gum_post_sams$shape2))
dens_gaus_mean <- colMeans(dens_gaus_sam)
dens_gum_mean <- colMeans(dens_gum_sam)
dens_gaus_qq <- apply(dens_gaus_sam, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
dens_gum_qq <- apply(dens_gum_sam, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
dens <- cbind(data_grid, ref_dat[,3], true_yy, dens_gaus_mean, t(dens_gaus_qq), dens_gum_mean, t(dens_gum_qq))

colnames(dens) <- c("grid", "dat", "true", "dens_gaus_mu", "dens_gaus_qq_025", "dens_gaus_qq_975", 
                    "dens_gum_mu", "dens_gum_qq_025", "dens_gum_qq_975")

plot_hist_save <- ggplot(as.data.frame(dens), aes(x = dat)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.25, fill="white", color="darkgrey", size = 1) + 
  geom_ribbon(aes(x = grid, ymin = dens_gaus_qq_025, ymax = dens_gaus_qq_975, fill = "Gaus 95% CI"), alpha = 0.5) + 
  geom_ribbon(aes(x = grid, ymin = dens_gum_qq_025, ymax = dens_gum_qq_975, fill = "Gum 95% CI"), alpha = 0.5) +
  geom_line(aes(x = grid, y = dens_gaus_mean, color = "Gaus mean"), linetype = "dashed", size = 1.5) + 
  geom_line(aes(x = grid, y = dens_gum_mean, color = "Gum mean"), linetype = "dashed", size = 1.5) +
  geom_line(aes(x = grid, y = true, color = "True"), linetype = "dashed", size = 1.5) + 
  scale_fill_manual("", values = c("brown1", "blue3")) +
  scale_colour_manual("", values = c("brown1", "blue3", "green")) +
  labs(fill = "", x = "Data", y = "Density") +  
  theme_bw() + ylim(0, 1.1) + 
  theme(plot.margin = unit(c(0, 1, 0.5, 1.5), "cm")) + 
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
png(paste0(outpath, "gamma_marginals.png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_hist_save)
dev.off()


#-------------------------------------------------------------------------------
# Gamma NNMP prediction on a grid
#-------------------------------------------------------------------------------
nonref_grid <- as.matrix(expand.grid(seq(0.005, 1, by = 1/100), seq(0.005, 1, by = 1/100)))

nnmp_gaus_pred <- predict(nnmp_gaus_out, nonref_coords = nonref_grid, predict_sam = FALSE)
nnmp_gum_pred <- predict(nnmp_gum_out, nonref_coords = nonref_grid, predict_sam = FALSE)

nnmp_gaus_surf <- as.data.frame(cbind(nonref_grid, nnmp_gaus_pred[[1]], t(nnmp_gaus_pred[[2]])))
nnmp_gum_surf <- as.data.frame(cbind(nonref_grid, nnmp_gum_pred[[1]], t(nnmp_gum_pred[[2]])))

colnames(nnmp_gaus_surf) <- c("x", "y", "mu", "q005", "q50", "q95")
colnames(nnmp_gum_surf) <- c("x", "y", "mu", "q005", "q50", "q95")

plot_save_gnnmp <- ggplot(as.data.frame(nnmp_gaus_surf), aes(x = x, y = y, fill = q50)) +
  theme_bw() + geom_raster(interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = range(yy)) + 
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
png(paste0(outpath,"gamma_gaus_q50_", nne, ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save_gnnmp)
dev.off()


plot_save_gnnmp <- ggplot(as.data.frame(nnmp_gum_surf), aes(x = x, y = y, fill = q50)) +
  theme_bw() + geom_raster(interpolate = T) +
  scale_fill_gradientn(colours = myPalette(100), limits = range(yy)) + 
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
png(paste0(outpath,"gamma_gum_q50_", nne, ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save_gnnmp)
dev.off()


#-------------------------------------------------------------------------------
# Gamma NNMP conditional survival probabilities
#-------------------------------------------------------------------------------
# Plot two selected locations
selec_sites <- rbind(c(0.5448276, 0.6034483), c(0.7, 0.55))
selec_sites_df <- cbind(1:nrow(selec_sites), selec_sites)

colnames(selec_sites_df) <- c("idx", "x", "y")
plot_save <- ggplot() + theme_bw() + 
  geom_raster(aes(x = lon, y = lat, fill = yy), data = as.data.frame(dat), interpolate = T) + 
  geom_point(aes(x = x, y = y), data = as.data.frame(selec_sites_df), size = 8, shape = 4, stroke = 2) + 
  geom_text(aes(label = idx, x = x, y = y), data = as.data.frame(selec_sites_df), vjust = 2, hjust = 0.8, size = 10) + 
  scale_fill_gradientn(colours = myPalette(100)) + 
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
png(paste0(outpath, "true_st_gamma_cross", ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()

# Neighbor information of the selected sites
selec_sites_ne_info <- neighbor(nne, ref_dat[,1:2], ref_dat[,3], ref = FALSE, selec_sites)
ne_index <- selec_sites_ne_info$ne_index
ne_dist <- selec_sites_ne_info$ne_dist
ne_obs <- selec_sites_ne_info$ne_obs

# Compute true conditional survival probabilities 
probs <- seq(0.5, 1, length = 100)

true_sur_prob <- array(NA, c(length(probs), nrow(selec_sites)))
expCovMat <- function(coords, phi){
  dd <- as.matrix(dist(coords))
  covmat <- exp(-dd / 0.1)
  covmat
}

response <- qt(probs, dof)
for (j in 1:nrow(selec_sites)){
  coord_j <- rbind(selec_sites[j,], ref_dat[ne_index[j,1:nne], 1:2])
  covmat_j <- expCovMat(coord_j, phi_t)
  lat_st_j <- as.matrix(qt(pgamma(ne_obs[j,1:nne], a_0, b_0), dof))
  for (i in seq_along(probs)) {
    cov_vec <- as.matrix(covmat_j[1,-1])
    con_mean <- t(cov_vec) %*% solve(covmat_j[-1,-1], lat_st_j)
    dd1 <- as.numeric(t(lat_st_j) %*% solve(covmat_j[-1,-1], lat_st_j))
    con_var <- (dof + dd1) / (dof + nne) * (1 - t(cov_vec) %*% solve(covmat_j[-1,-1], cov_vec))
    true_sur_prob[i,j] <- pt((response[i] - con_mean) / sqrt(con_var), dof + nne, lower.tail = FALSE)
  }
}

# Compute conditional survival probabilities 
nnmp_gaus_cdf <- pnnmp(nnmp_gaus_out, selec_sites, probs)
nnmp_gum_cdf <- pnnmp(nnmp_gum_out, selec_sites, probs)

# Plot conditional survival probabilities 
for (j in 1:nrow(selec_sites)) {
  gaus_sur_prob_mean <- rowMeans(1 - nnmp_gaus_cdf[j,,])
  gum_sur_prob_mean <- rowMeans(1 - nnmp_gum_cdf[j,,])
  gaus_sur_prob_qq <- apply(1 - nnmp_gaus_cdf[j,,], 1, function(x) quantile(x, probs = c(0.025, 0.975))) 
  gum_sur_prob_qq <- apply(1 - nnmp_gum_cdf[j,,], 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  sur_probs <- cbind(probs, true_sur_prob[,j], gaus_sur_prob_mean, t(gaus_sur_prob_qq), gum_sur_prob_mean, t(gum_sur_prob_qq))
  colnames(sur_probs) <- c("xx", "true", "gaus_mu", "gaus_q025", "gaus_q975", "gum_mu", "gum_q025", "gum_975")
  plot_hist_save <- ggplot(as.data.frame(sur_probs), aes(x = xx)) + 
    geom_ribbon(aes(x = xx, ymin = gaus_q025, ymax = gaus_q975, fill = "Gaus 95% CI"), alpha = 0.5) + 
    geom_ribbon(aes(x = xx, ymin = gum_q025, ymax = gum_975, fill = "Gum 95% CI"), alpha = 0.5) + 
    geom_line(aes(x = xx, y = gaus_mu, color = "Gaus mean"), linetype = "dashed", size = 1.5) + 
    geom_line(aes(x = xx, y = gum_mu, color = "Gum mean"), linetype = "dashed", size = 1.5) + 
    geom_line(aes(x = xx, y = true, color = "true"), size = 1.5, linetype = "dashed") + 
    scale_fill_manual("", values = c("brown1", "blue3")) + 
    scale_colour_manual("", values = c("brown1", "blue3", "green")) + 
    labs(fill = "", x = "Marginal probability", y = "Conditional survival probability") + #, 
    theme_bw() + ylim(0, 1) +  xlim(0.5, 1) + 
    theme(plot.margin = unit(c(0, 1, 0.5, 1.5), "cm")) + 
    theme(legend.position="none",
          axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent"))
  png(paste0(outpath, "gamma_sur_prob_model_", j, ".png"), width=600, height=500, pointsize=20, bg = "transparent")
  print(plot_hist_save)
  dev.off()
}