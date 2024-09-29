############### take a look at spatial domain subcustering of habenula ###########

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(rearrr)
  library(reticulate)
  library(anndata)
  library(plyr)
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

############ extract metadata from anndata file for 687997 ################

sis <- read_h5ad("/data/merscope_638850_687997_mouseadult_each_dataset_staligner_harmony_integration_fullcell_spatial_subdomain/mouse_687997_fullcell_spatial_subdommain.h5ad")
metadata_sis <- sis$obs
