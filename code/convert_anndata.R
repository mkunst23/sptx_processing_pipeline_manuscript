############ extract data tables from anndata files ##################

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

sis <- read_h5ad("/data/merscope_638850_mouseadult_processed/whole_dataset/mouse_687997_filtered.h5ad")
metadata_sis <- sis$obs
save(metadata_sis, file = "/scratch/metadata_sis.rda")

# extract coordinates and and add sections
coordinates_sis <- as.matrix(sis$obsm$spatial_rotate)
rownames(coordinates_sis) <- rownames(sis$X)
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_coordinate","y_coordinate")
save(coordinates_sis, file = "/scratch/coordinates_sis.rda")


# extract count matrix
data_sis <- sis$X
data_sis <- as.data.frame(data_sis)
# Convert row names into a column
data_sis <- rownames_to_column(data_sis, var = "sample_name")
save(data_sis, file = "/scratch/data_sis.rda")
