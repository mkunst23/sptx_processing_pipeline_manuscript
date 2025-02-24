
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

######################## extract anndata for 638850 #######################

sis <- read_h5ad("/data/merscope_638850_mouseadult_registered/whole_dataset/mouse_638850_registered.h5ad")
metadata_sis <- sis$obs
save(metadata_sis, file = "/scratch/638850_metadata_sis.rda")

# extract coordinates and and add sections
coordinates_sis <- as.matrix(sis$obsm$reconstructed_coords)
rownames(coordinates_sis) <- rownames(sis$X)
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_reconstructed","y_reconstructed","z_reconstructed")
save(coordinates_sis, file = "/scratch/638850_reconstructed_coordinates_sis.rda")

# extract reconstructed coordinates
coordinates_sis <- as.matrix(sis$obsm$spatial_rotate)
rownames(coordinates_sis) <- rownames(sis$X)
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_coordinate","y_coordinate")
save(coordinates_sis, file = "/scratch/638850_coordinates_sis.rda")

# extract count matrix
data_sis <- sis$X
data_sis <- as.data.frame(data_sis)
# Convert row names into a column
data_sis <- rownames_to_column(data_sis, var = "sample_name")
save(data_sis, file = "/scratch/638850_data_sis.rda")

# extract spatial domains
sd_domains <- read_h5ad("/data/merscope_638850_mouseadult_processed_clust_QC_filt/whole_dataset/mouse_638850_processed.h5ad")
metadata_sd_domains <- sd_domains$obs
metadata_sd_domains <- metadata_sd_domains %>% 
  select(production_cell_id,
         leiden_res_0.6_knn_8,
         leiden_res_0.8_knn_8,
         leiden_res_1.0_knn_8,
         leiden_res_1.2_knn_8,
         leiden_res_1.4_knn_8,
         leiden_res_1.6_knn_8)

fwrite(metadata_sd_domains,
       '/scratch/mouse_638850_sd.csv',
       row.names = F)

################# extract vpt anndata for 638850 ################
vpt <- read_h5ad("/data/merscope_638850_mouseadult_processed_VPT/whole_dataset/mouse_638850_filtered.h5ad")
metadata_vpt <- vpt$obs
save(metadata_vpt, file = "/scratch/metadata_vpt.rda")



################### extract anndata for 687997 #####################
sis <- read_h5ad("/data/merscope_687997_mouseadult_registered/whole_dataset/mouse_687997_registered.h5ad")
metadata_sis <- sis$obs
save(metadata_sis, file = "/scratch/687997_metadata_sis.rda")

# extract coordinates and and add sections
coordinates_sis <- as.matrix(sis$obsm$reconstructed_coords)
rownames(coordinates_sis) <- rownames(sis$X)
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_reconstructed","y_reconstructed","z_reconstructed")
save(coordinates_sis, file = "/scratch/687997_reconstructed_coordinates_sis.rda")

# extract reconstructed coordinates
coordinates_sis <- as.matrix(sis$obsm$spatial_rotate)
rownames(coordinates_sis) <- rownames(sis$X)
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_coordinate","y_coordinate")
save(coordinates_sis, file = "/scratch/687997_coordinates_sis.rda")

# extract count matrix
data_sis <- sis$X
data_sis <- as.data.frame(data_sis)
# Convert row names into a column
data_sis <- rownames_to_column(data_sis, var = "sample_name")
save(data_sis, file = "/scratch/687997_data_sis.rda")
