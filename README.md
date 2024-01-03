## NNMP: Nearest-neighbor mixture process

This repository contains the R package [**nnmp**](https://github.com/xzheng42/nnmp-examples-ba-2023/tree/main/nnmp) (currently developer's version) 
and R scripts to reproduce the numerical results in

Xiaotian Zheng. Athanasios Kottas. Bruno Sansó. 
"Nearest-Neighbor Mixture Models for Non-Gaussian Spatial Processes." Bayesian Anal. 18 (4) 1191 - 1222, December 2023. https://doi.org/10.1214/23-BA1405

### Installing and using the **nnmp** package

You can install the package with **devtools**
```
devtools::install_github("xzheng42/nnmp-examples-ba-2023", subdir = "nnmp")
library(nnmp)
```

Main functions of the package are `nnmp` and `predict.nnmp`:

- `nnmp` fits an NNMP model via Markov chain Monte Carlo (MCMC).
- `predict.nnmp` (or simply `predict`) generates posterior predictive samples for a set of new locations given an nnmp object.

Notes: The current version was tested on macOS 10.15.7 under R version 4.2.2 and on Fedora Linux 38 under R versions 4.1.3 and 4.3.2.

### Workflow to reproduce numerical results

R scripts to reproduce results in the paper are available in 
[*data-examples/*](https://github.com/xzheng42/nnmp-examples-ba-2023/tree/main/data-examples),
and [*data/*](https://github.com/xzheng42/nnmp-examples-ba-2023/tree/main/data) contains sea surface temperature (SST) data of 
the Mediterranean Sea and relevant shape files.

- Run all simulation experiments: `run_all_sim_rscripts.R`.
- Run all SST data examples: `run_all_real_rscripts.R`.
- Run simulation experiments in Section 5.1: 
  `sim_gamma_scenario.R` (gamma NNMP); 
  `sim_beta_scenario.R` (beta NNMP).
- Run simulation experiments in Section B.1 of the supplementary material: 
  `sim_gaussian_scenario.R` (Gaussian NNMP vs. NNGP);
  `sim_comp_beta_first_scenario.R`and `sim_comp_beta_second_scenario.R` (beta NNMP vs. MGP, SPDE-INLA, NNGP);
  `sim_sn_scenario.R` (skew-Gaussian NNMP);
  `sim_comp_sn_scenario.R` (skew-Gaussian NNMP vs. SPDE-INLA).

- Run SST data examples in Section 5.1: `prepare_sst_data.R` and `real_sst_all.R`.
- Run SST regional data analyses in Section B.1: 
  `real_comp_sst_subset_nnmp.R` and `real_comp_sst_subset_nngp.R` (Gaussian NNMP vs. NNGP); 
  `real_sst_subset.R`.
- Run SST global data analyses in Section B.1: 
  `real_comp_sst_all_nnmp.R` and `real_comp_sst_all_spde.R` (extended skew-Gaussian NNMP vs. SPDE-INLA).
