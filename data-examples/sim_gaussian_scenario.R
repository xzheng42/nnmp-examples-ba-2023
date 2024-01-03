############################################################################
### SIMULATION STUDY: Comparison with the NNGP 
### Data generated from a full GP
### Simulation setting follows Datta et al. (2016)
############################################################################
rm(list = ls())

outpath <- "/path/to/output/"

library(SpatialExtremes)
library(ggplot2)
library(readr)

library(nnmp)


#-------------------------------------------------------------------------------
# Generate the true spatial random field
#-------------------------------------------------------------------------------
set.seed(42)

# Simulate from a Gaussian process on a grid
x_seq <- seq(0, 1, length = 100)
y_seq <- seq(0, 1, length = 100)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
zz <- t(rgp(1, grid_loc, cov.mod = "powexp", nugget = 0, sill = 1, range = 1/12, smooth = 1))

# Simulate data using a spatially varying regression
XX <- cbind(1, rnorm(length(zz)))
bb <- as.matrix((c(1, 5)))
tau2 <- 0.1
yy <- XX %*% bb + zz + rnorm(length(zz), 0, sqrt(tau2))

# Plot the true field
dat <- cbind(grid_loc, yy)
colnames(dat) <- c("lon", "lat", "yy")

plot_save <- ggplot() + theme_bw() + 
  geom_raster(aes(x = lon, y = lat, fill = zz), data = as.data.frame(dat), interpolate = T) + 
  scale_fill_gradient2(low = "red", mid = "white",  high = "blue", 
                       midpoint = 0, space = "Lab", limits = range(zz),
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
png(paste0(outpath, "true_gp", ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()


#-------------------------------------------------------------------------------
# Fit Gaussian NNMP
#-------------------------------------------------------------------------------
set.seed(42)

# Randomly select part of the gridded data as observations 
obs_num <- 2500
obs_idx <- sample(1:nrow(grid_loc), obs_num, FALSE)
obs_dat <- dat[obs_idx,]
obs_XX <- XX[obs_idx,]

# Partition the observations into training and testing sets
# Use the whole training set as reference set 
ref_num <- 2000
test_num <- obs_num - ref_num
ref_idx <- sample(1:obs_num, ref_num, FALSE)
ref_dat <- obs_dat[ref_idx,]
ref_XX <- obs_XX[ref_idx, ]
test_dat <- obs_dat[-ref_idx,]
test_XX <- obs_XX[-ref_idx, ]

# Random ordering
ord <- sample(1:ref_num, replace = FALSE)

# MCMC
tuning <- list("phi" = 0.1, "zeta" = 0.15)

mcmc_settings <- list(n_iter = 30000, n_burn = 10000, n_thin = 10, n_report = 2000)

nne <- 10

nnmp_out <- nnmp(response = ref_dat[,3],
                 covars = ref_XX,
                 coords = ref_dat[,1:2],
                 neighbor_size = nne,
                 marg_family = "gaussian",
                 tuning = tuning,
                 ord = ord,
                 mcmc_settings = mcmc_settings, 
                 verbose = TRUE, 
                 model_diag = list("dic", "pplc"))


#-------------------------------------------------------------------------------
# Evaluate Gaussian NNMP
#-------------------------------------------------------------------------------
# Prediction
nnmp_test <- predict(nnmp_out, test_XX, test_dat[,1:2], predict_sam = TRUE)

# RMSPE
nnmp_rmspe <- sqrt(mean((nnmp_test$obs_mu - test_dat[,3])^2))

# 95% CI width
nnmp_qq <- nnmp_test[['obs_qq']]
nnmp_width <- mean(nnmp_qq[3,] - nnmp_qq[1,])

# 95% CI Coverage rate
nnmp_cover_boolen <- sapply(1:(test_num), function(x) {
  test_dat[x,3] >= nnmp_qq[1,x] & test_dat[x,3] <= nnmp_qq[3,x]
})
nnmp_cover <- sum(nnmp_cover_boolen) / length(nnmp_cover_boolen)

# CRPS
nnmp_crps <- mean(scoringRules::crps_sample(test_dat[, 3], nnmp_test$obs_sam))

# Put together 
nnmp_metric <- c(nnmp_rmspe, nnmp_crps, nnmp_cover, nnmp_width, 
                 nnmp_out$mod_diag$pplc[3], nnmp_out$mod_diag$dic)


#-------------------------------------------------------------------------------
# Gaussian NNMP prediction on a grid
#-------------------------------------------------------------------------------
x_seq <- seq(0.005, 1, by = 0.01)
y_seq <- seq(0.005, 1, by = 0.01)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
XX <- cbind(1, rnorm(nrow(grid_loc)))

nnmp_pred <- predict(nnmp_out, XX, grid_loc, predict_sam = FALSE)
nnmp_surf <- cbind(grid_loc, nnmp_pred$zz_mu, t(nnmp_pred$zz_qq))
colnames(nnmp_surf) <- c("lon", "lat", "mu", "q025", "q50", "q975")

plot_save_nnmp <- ggplot(as.data.frame(nnmp_surf), aes(x = lon, y = lat, fill = q50)) +
  theme_bw() + geom_raster(interpolate = T) +
  scale_fill_gradient2(low = "red", mid = "white", limits = range(zz),
                     high = "blue", space = "Lab",
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
 
png(paste0(outpath, "sim_gnnmp_median_", nne, ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save_nnmp)
dev.off()


#-------------------------------------------------------------------------------
# Fit NNGP
#-------------------------------------------------------------------------------
library(spNNGP)

set.seed(42)

nsamples <- mcmc_settings$n_iter
starting <- list("phi"=3/0.5, "sigma.sq"=1, "tau.sq"=1)
tuning <- list("phi"=0.2,  "sigma.sq"=0.5, "tau.sq"=0.5)
priors <- list("phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2,1), "tau.sq.IG"=c(2,0.1))

nngp_out <- spNNGP(formula = ref_dat[,3] ~ ref_XX[, 2], 
                 coords = ref_dat[,1:2],
                 starting = starting, 
                 tuning = tuning, 
                 ord = ord,
                 priors = priors, 
                 cov.model = "exponential",
                 n.samples = nsamples,
                 n.neighbors = nne,
                 method = "latent", 
                 n.omp.threads = 1, 
                 n.report = 5000)


#-------------------------------------------------------------------------------
# Evaluate NNGP
#-------------------------------------------------------------------------------
# Prediction
nngp_test <- predict(nngp_out, test_XX, test_dat[,1:2],
                     sub.sample = list(start = mcmc_settings$n_burn+1, 
                                       end = nsamples, 
                                       thin = mcmc_settings$n_thin))

# Obtain posterior predictive samples and compute summaries
nngp_test_sam <- nngp_test$p.y.0
nngp_mean <- rowMeans(nngp_test_sam, 1)
nnmp_q025 <- apply(nngp_test_sam, 1, function(x) quantile(x, probs = 0.025))
nnmp_q975 <- apply(nngp_test_sam, 1, function(x) quantile(x, probs = 0.975))

# RMSPE
nngp_rmspe <- sqrt(mean((nngp_mean - test_dat[,3])^2))

# 95% CI width
nngp_width <- mean(nnmp_q975 - nnmp_q025)

# 95$ CI coverage rate
nngp_cover_boolen <- sapply(1:test_num, function(x) {
  test_dat[x,3] >= nnmp_q025[x] & test_dat[x,3] <= nnmp_q975[x]
})
nngp_cover <- sum(nngp_cover_boolen) / length(nngp_cover_boolen)

# CRPS
nngp_crps <- mean(scoringRules::crps_sample(test_dat[,3], nngp_test_sam))

# PPLC and DIC
nngp_diag <- spDiag(nngp_out, sub.sample = list(start = mcmc_settings$n_burn+1, 
                                                end = nsamples, 
                                                thin = mcmc_settings$n_thin))

# Put together
nngp_metric <- c(nngp_rmspe, nngp_crps, nngp_cover, nngp_width, 
                 unlist(nngp_diag$GP)[3], unlist(nngp_diag$DIC)[1])


#-------------------------------------------------------------------------------
# NNGP prediction on a grid
#-------------------------------------------------------------------------------
x_seq <- seq(0.005, 1, by = 0.01)
y_seq <- seq(0.005, 1, by = 0.01)
grid_loc <- as.matrix(expand.grid(x_seq, y_seq))
XX <- cbind(1, rnorm(nrow(grid_loc)))
nngp_pred <- predict(nngp_out, XX, grid_loc, 
                     sub.sample = list(start = mcmc_settings$n_burn+1, 
                                       end = nsamples, 
                                       thin = mcmc_settings$n_thin),
                     n.report = 1000)
nngp_pred_sam <- nngp_pred$p.w.0
nngp_pred_mean <- rowMeans(nngp_pred_sam, 1)
nngp_pred_qq <- apply(nngp_pred_sam, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
nngp_surf <- cbind(grid_loc, nngp_pred_mean, t(nngp_pred_qq))
colnames(nngp_surf) <- c("lon", "lat", "mu", "q025", "q50", "q975")

plot_save_ngpp <- ggplot(as.data.frame(nngp_surf), aes(x = lon, y = lat, fill = q50)) +
  theme_bw() +  geom_raster(interpolate = T) +
  scale_fill_gradient2(low = "red", mid = "white",limits = range(zz),
                       high = "blue", space = "Lab",
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
png(paste0(outpath, "sim_nngp_median_", 10, ".png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save_ngpp)
dev.off()


#-------------------------------------------------------------------------------
# Model comparison (Table 1 of the supplementary material)
#-------------------------------------------------------------------------------
tbl_pred <- cbind(nngp_metric, nnmp_metric)
colnames(tbl_pred) <- c("NNGP", "NNMP")
rownames(tbl_pred) <- c("RMSPE", "CRPS","95% CI cover", "95% CI width", "PPLC", "DIC")

tbl_pred_latex <- knitr::kable(t(tbl_pred), digits = 3, align = 'c', format = 'latex')
write_file(tbl_pred_latex, paste0(outpath, "sim_gnnmp_nngp_pred_comp.txt"))
