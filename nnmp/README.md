## NNMP: Nearest-neighbor mixture process

This is the R package **nnmp** (currently developer's version) for the paper:

Xiaotian Zheng. Athanasios Kottas. Bruno Sansó. 
"Nearest-Neighbor Mixture Models for Non-Gaussian Spatial Processes." Bayesian Anal. 18 (4) 1191 - 1222, December 2023. https://doi.org/10.1214/23-BA1405

You can install the package with **devtools**
```
devtools::install_github("xzheng42/nnmp-examples-ba-2023", subdir = "nnmp")
```

Main functions of the package are `nnmp` and `predict.nnmp`:

- `nnmp` fits an NNMP model via Markov chain Monte Carlo (MCMC).
- `predict.nnmp` (or simply `predict`) generates posterior predictive samples for a set of new locations given an nnmp object.

Detailed guidelines for using the functions are referred to their help pages in R. 
R scripts to reproduce results in the paper are available in [*data-examples/*](https://github.com/xzheng42/nnmp-examples-ba-2023/tree/main/data-examples)
with instructions available in [*nnmp-examples-ba-2023*](https://github.com/xzheng42/nnmp-examples-ba-2023/).

Notes: The current version was tested on macOS 10.15.7 under R version 4.2.2 and on Fedora Linux 38 under R versions 4.1.3 and 4.3.2.
