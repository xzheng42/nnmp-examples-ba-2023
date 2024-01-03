################################################################################
### Master function: run all R scripts for SST data examples
################################################################################
setwd("set/directory/to/data-examples/")


#-------------------------------------------------------------------------------
# SST data example of the main paper
#-------------------------------------------------------------------------------
# Produce Figure 5(a) of the main paper
source("prepare_sst_data.R")

# Produce Figures 4 and 5(b)-(d) of the main paper
source("real_sst_all.R")


#-------------------------------------------------------------------------------
# SST regional data analysis of the supplementary material (Section C.1)
# Produce Tables 5 and 6 of the supplementary material
#-------------------------------------------------------------------------------
source("real_comp_sst_subset_nnmp.R")
source("real_comp_sst_subset_nngp.R")


#-------------------------------------------------------------------------------
# SST regional data analysis of the supplementary material (Section C.1)
# Produce Figure 4 of the supplementary material
#-------------------------------------------------------------------------------
source("real_sst_subset.R")


#-------------------------------------------------------------------------------
# SST global data analysis of the supplementary material (Section C.2)
# Produce Tables 7 and 8 of the supplementary material
#-------------------------------------------------------------------------------
source("real_comp_sst_all_nnmp.R")
source("real_comp_sst_all_spde.R")


