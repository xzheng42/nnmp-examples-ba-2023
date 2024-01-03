################################################################################
### Master function: run all R scripts for simulation experiments
################################################################################
setwd("set/directory/to/data-examples/")


#-------------------------------------------------------------------------------
# First simulation experiment of the main paper (Section 5.1)
# Produce Figure 2 of the main paper
#-------------------------------------------------------------------------------
source("sim_gamma_scenario.R")


#-------------------------------------------------------------------------------
# Second simulation experiment of the main paper (Section 5.1)
# Produce Figure 3 of the main paper
#-------------------------------------------------------------------------------
source("sim_beta_scenario.R")


#-------------------------------------------------------------------------------
# First simulation experiment of the supplementary material (Section B.1)
# Produce Figure 1 and Table 1 of the supplementary material
#-------------------------------------------------------------------------------
source("sim_gaussian_scenario.R")


#-------------------------------------------------------------------------------
# Second simulation experiment of the supplementary material (Section B.2)
# Produce Tables 2 and 3 of the supplementary material
#-------------------------------------------------------------------------------
# First scenario
source("sim_comp_beta_first_scenario.R")
# Second scenario
source("sim_comp_beta_second_scenario.R")


#-------------------------------------------------------------------------------
# Third simulation experiment of the supplementary material (Section B.3)
#-------------------------------------------------------------------------------
# Produce Figures 2 and 3 of the supplementary material
source("sim_sn_scenario.R")
# Produce Table 4 of the supplementary material
source("sim_comp_sn_scenario.R")