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
  library(randomcoloR)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

############ extract metadata from anndata file for 687997 ################

sis <- read_h5ad("/data//merscope_687997_mouseadult_staligner_harmony_integration_fullcell_spatial_subdomain/mouse_687997_fullcell_spatial_subdomain.h5ad")
metadata_sis <- sis$obs

coordinates_sis <- as.matrix(sis$obsm$spatial_rotate)
rownames(coordinates_sis) <- metadata_sis$production_cell_id
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_coordinate","y_coordinate")

metadata_sis <- merge(metadata_sis,
                      coordinates_sis,
                      by.x = "production_cell_id",
                      by.y = 0)

metadata_subset <- metadata_sis %>% 
  filter(final_filter == F) %>% 
  select(production_cell_id,
         section,
         x_coordinate,
         y_coordinate,
         flat_CDM_class_name,
         flat_CDM_subclass_name,
         flat_CDM_supertype_name,
         flat_CDM_cluster_name,
         leiden_res_1.2_knn_8,
         leiden_STAligner_Harmony_1.2_sub0.2_knn_8,   
         leiden_STAligner_Harmony_1.2_sub0.4_knn_8,   
         leiden_STAligner_Harmony_1.2_sub0.5_knn_8,   
         leiden_STAligner_Harmony_1.2_sub0.6_knn_8,   
         leiden_STAligner_Harmony_1.2_sub0.8_knn_8,   
         leiden_STAligner_Harmony_1.2_sub1.0_knn_8,   
         leiden_STAligner_Harmony_1.2_sub1.2_knn_8,   
         leiden_STAligner_Harmony_1.2_sub1.4_knn_8,   
         leiden_STAligner_Harmony_1.2_sub1.5_knn_8,   
         leiden_STAligner_Harmony_1.2_sub1.6_knn_8,   
         leiden_STAligner_Harmony_1.2_sub1.8_knn_8,   
         leiden_STAligner_Harmony_1.2_sub2.0_knn_8)


# set of sections
section_set <- c("1288760867","1288760871","1288760876","1288760880","1288760883","1288760887")

example_sis <- metadata_subset %>% 
  filter(section %in% section_set) %>% 
  filter(leiden_res_1.2_knn_8 == "47")

# plot sptial domains
leiden_cluster_colors <- randomColor(length(unique(metadata_sis$leiden_res_1.2_knn_8)), luminosity="light")
leiden_cluster_palette <- setNames(leiden_cluster_colors, unique(metadata_sis$leiden_res_1.2_knn_8))

for (i in section_set) {
  
  plot_data <- metadata_subset %>% 
    filter(section == i) %>% 
    filter(leiden_res_1.2_knn_8 == "47") %>% 
    filter(x_coordinate > -750) %>% 
    filter(x_coordinate < 850) %>% 
    filter(y_coordinate > -1200) %>% 
    filter(y_coordinate < 0)
  
  plot <- ggplot(plot_data,
         aes(x=x_coordinate,
             y=y_coordinate,
             color=leiden_res_1.2_knn_8
         )) +
    geom_point(size=2,
               stroke=0,
               shape=19) +
    coord_fixed() +
    scale_color_manual(values=leiden_cluster_palette) +
    scale_y_reverse() +
    theme(legend.position = "none") +
    guides(color = "none") +
    theme_void()
  
  ggsave(filename = paste0("/results/habenula_",i,".png"), 
         plot = plot, 
         width = 6,
         height = 6, 
         dpi = 300)
}


# plot sptial domains
leiden_cluster_colors <- randomColor(length(unique(metadata_sis$leiden_STAligner_Harmony_1.2_sub1.8_knn_8)), luminosity="light")
leiden_cluster_palette <- setNames(leiden_cluster_colors, unique(metadata_sis$leiden_STAligner_Harmony_1.2_sub1.8_knn_8))

for (i in section_set) {
  
  plot_data <- metadata_subset %>% 
    filter(section == i) %>% 
    filter(leiden_res_1.2_knn_8 == "47") %>% 
    filter(x_coordinate > -750) %>% 
    filter(x_coordinate < 850) %>% 
    filter(y_coordinate > -1200) %>% 
    filter(y_coordinate < 0)
  
  plot <- ggplot(plot_data,
         aes(x=x_coordinate,
             y=y_coordinate,
             color=leiden_STAligner_Harmony_1.2_sub1.8_knn_8
         )) +
    geom_point(size=1.5,
               stroke=0,
               shape=19) +
    coord_fixed() +
    scale_color_manual(values=leiden_cluster_palette) +
    scale_y_reverse() +
    theme(legend.position = "none") +
    guides(color = "none") +

    theme_void()
  
  ggsave(filename = paste0("/results/habenula_subdomain",i,".png"), 
         plot = plot, 
         width = 6,
         height = 6, 
         dpi = 300)
}

ggplot(example_sis,
       aes(x=x_coordinate,
           y=y_coordinate,
           color=leiden_STAligner_Harmony_1.2_sub1.8_knn_8
       )) +
  geom_point(size=.5,
             stroke=0,
             shape=19) +
  coord_fixed() +
  #ggtitle("SIS Segmentation") +
  scale_color_manual(values=leiden_cluster_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~section) +
  theme_void()

