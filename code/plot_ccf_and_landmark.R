suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(cowplot)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(data.table)
  library(rearrr)
  library(googlesheets4)
  library(randomcoloR)
  library(forcats)
})

# setup environment
options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()

# cluster colors
ss <- "https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?usp=sharing"
cluster_colors <- read_sheet(ss, sheet="clusters")
cluster_color_palette <- setNames(cluster_colors$cluster_color, cluster_colors$cluster_label)

# ccf colors
ccf <- "https://docs.google.com/spreadsheets/d/1QOhsYhlsk2KE2pZuSnEwUhEIG4mh782PcAgHrtJyRYI/edit?gid=0#gid=0"
ccf_anno <- read_sheet(ccf, sheet = "CCFv3")
ccf_anno$color_hex_triplet <- gsub("^", "#", ccf_anno$color_hex_triplet)
ccf_color_palette <- setNames(ccf_anno$color_hex_triplet,ccf_anno$acronym)

sd <- "https://docs.google.com/spreadsheets/d/1UpfDYDbq3s0Jeh8cLg9Y0yhPCJxH_H8WG4Y0y0_mNHg/edit?gid=70495470#gid=70495470"
spatial_domains <- read_sheet(sd, sheet = "knn_8_res_1-2_30_687997_1")
spatial_domains <- spatial_domains[,1:4]

ccf_broad <- read_sheet(sd, sheet = "CCF_broad_region_color")
ccf_broad_palette <- setNames(ccf_broad$color_hex_tripet,ccf_broad$ccf_broad)

anat <- "https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=2094634713#gid=2094634713"
ccf_landmarks <- read_sheet(anat, sheet = "CCF_landmark_color")
ccf_landmark_palette <- setNames(ccf_landmarks$Color,ccf_landmarks$CCF_acronym)

anat_df <- read_sheet(anat, sheet = "cluster_region_annotation")
anat_df <- anat_df %>% 
  select(cluster_id_label,
         registration_landmark)

# pick sections
example_sections <- c("1288760916",
                      "1288760891",
                      "1288760867",
                      "1288760822")

# load metadata
load("/scratch/687997_metadata_sis.rda")

# select relevant columns
metadata_subset <- metadata_sis %>% 
  filter(final_filter == F) %>% 
  select(production_cell_id,
         section,
         leiden_res_1.2_knn_8,
         flat_CDM_class_name,
         flat_CDM_subclass_name,
         flat_CDM_supertype_name,
         flat_CDM_cluster_name,
         volume_x,
         volume_y,
         volume_z) %>% 
  filter(section %in% example_sections)

metadata_subset <- merge(metadata_subset,
                         spatial_domains,
                         by.x = "leiden_res_1.2_knn_8",
                         by.y = "Cluster_id",
                         all.x = T,
                         all.y = F)

metadata_subset <- merge(metadata_subset,
                         anat_df,
                         by.x = "flat_CDM_cluster_name",
                         by.y = "cluster_id_label",
                         all.x = T,
                         all.y = F)

# create colormap for leiden clustering
leiden_cluster_colors <- randomColor(length(unique(metadata_sis$leiden_res_1.2_knn_8)), luminosity="dark")
leiden_cluster_palette <- setNames(leiden_cluster_colors, unique(metadata_sis$leiden_res_1.2_knn_8))

# plot sections
plot <- ggplot(subset(metadata_subset, spatial_domain_level_1 != "LQ"),
               aes(x=volume_x,
                   y=volume_y,
                   color = leiden_res_1.2_knn_8
               )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Spatial Domains") +
  scale_color_manual(values=leiden_cluster_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/spatial_domains.png", 
       plot = plot, 
       width = 5,
       height = 12, 
       dpi = 300)

plot <- ggplot(subset(metadata_subset, spatial_domain_level_1 != "LQ"),
       aes(x=volume_x,
           y=volume_y,
           color = ccf_broad
       )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Broad Regions") +
  scale_color_manual(values=ccf_broad_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/ccf_broad.png", 
       plot = plot, 
       width = 5,
       height = 12, 
       dpi = 300)

plot <- ggplot(subset(metadata_subset, registration_landmark!= "NA"),
       aes(x=volume_x,
           y=volume_y,
           color = flat_CDM_cluster_name
       )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Cell Types") +
  scale_color_manual(values=cluster_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/clusters.png", 
       plot = plot, 
       width = 5,
       height = 12, 
       dpi = 300)

plot <- ggplot(subset(metadata_subset, registration_landmark!= "NA"),
       aes(x=volume_x,
           y=volume_y,
           color = registration_landmark
       )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("CCF landmarks") +
  scale_color_manual(values=ccf_landmark_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/registration_landmarks.png", 
       plot = plot, 
       width = 5,
       height = 12, 
       dpi = 300)