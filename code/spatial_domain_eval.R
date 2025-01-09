#################### spatial_domain_evaluation ######################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)   
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(anndata)
  library(googlesheets4)
  library(data.table)
  library(reshape)
  library(forcats)
  library(entropy)
  library(aplot)
  library(vegan)
  library(ggridges)
  library(scales)
  library(DescTools)
  library(RColorBrewer)
  library(philentropy)
  library(igraph)
  library(paletteer)
  library(randomcoloR)
  library(stringr)
})

############### setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

example_sections <- c("1199651069",
                      "1199651057",
                      "1199651036",
                      "1199651012",
                      "1199650953",
                      "1199650941")

neuronal_classes <- c("01 IT-ET Glut",
                      "02 NP-CT-L6b Glut",
                      "03 OB-CR Glut",
                      "06 CTX-CGE GABA",
                      "07 CTX-MGE GABA",
                      "08 CNU-MGE GABA",
                      "09 CNU-LGE GABA",
                      "10 LSX GABA",
                      "11 CNU-HYa GABA",
                      "12 HY GABA",        
                      "13 CNU-HYa Glut",
                      "14 HY Glut",
                      "15 HY Gnrh1 Glut",
                      "16 HY MM Glut",
                      "17 MH-LH Glut",
                      "18 TH Glut",
                      "19 MB Glut",
                      "20 MB GABA",
                      "21 MB Dopa",
                      "22 MB-HB Sero",
                      "23 P Glut",
                      "24 MY Glut",
                      "25 Pineal Glut",
                      "26 P GABA",
                      "27 MY GABA",
                      "28 CB GABA",
                      "29 CB Glut")

white_matter_classes <- "31 OPC-Oligo"

ventricles_classes <- c("325 CHOR NN",
                        "323 Ependymal NN")

wm_1_0 <- c("0","5","11","25")
ventricles_1_0 <- c("3","8") 

wm_1_2 <- c("5","7","11","18")
ventricles_1_2 <- c("1","2")   

wm_1_4 <- c("1","7","9","30")
ventricles_1_4 <- c("3","27")   
    
# extract clustering results from anndata file
ad <- read_h5ad("/data/merscope_638850_mouseadult_processed_filt_z_1.0_1.2_1.4/mouse_638850_filt_10_12_14.h5ad")
# extract metadata
metadata <- ad$obs
# save as rda file
metadata_spatial_clustering_eval <- metadata %>% 
  select(leiden_res_1.0_knn_8_filt,
         leiden_res_1.2_knn_8_filt,
         leiden_res_1.4_knn_8_filt,
         qc_passed,
         hrc_mmc_subclass_name,
         hrc_mmc_class_name,
         production_cell_id)

save(metadata_spatial_clustering_eval, file = "/scratch/638850_spatial_clustering_eval.rda")
# extract relevant columns

# load metadata as rda file
load("/scratch/638850_spatial_clustering_eval.rda")
# read in coordinates
load("/scratch/638850_coordinates_sis.rda")

metadata_spatial_clustering_eval <- merge(metadata_spatial_clustering_eval,
                                          coordinates_sis,
                                          by.x = "production_cell_id",
                                          by.y = 0)

metadata_spatial_clustering_eval$leiden_res_1.0_knn_8_filt <- as.character(metadata_spatial_clustering_eval$leiden_res_1.0_knn_8_filt)
metadata_spatial_clustering_eval$leiden_res_1.2_knn_8_filt <- as.character(metadata_spatial_clustering_eval$leiden_res_1.2_knn_8_filt)
metadata_spatial_clustering_eval$leiden_res_1.4_knn_8_filt <- as.character(metadata_spatial_clustering_eval$leiden_res_1.4_knn_8_filt)


metadata_spatial_clustering_eval <- metadata_spatial_clustering_eval %>% 
  filter(qc_passed == TRUE)

# classify sections
metadata_spatial_clustering_eval <- metadata_spatial_clustering_eval %>% 
  mutate(N_G_V_ratio = case_when(hrc_mmc_class_name %in% neuronal_classes ~ "neuronal",
                                 hrc_mmc_class_name %in% white_matter_classes ~ "glial",
                                 hrc_mmc_subclass_name %in% ventricles_classes ~ "ventricles",
                                 TRUE ~ NA_character_))


# subset only spatial domains and sections I need for each resolution
metadata_section_plot <- metadata_spatial_clustering_eval %>% 
  filter(section %in% example_sections) %>% 
  filter


# plot example sections


