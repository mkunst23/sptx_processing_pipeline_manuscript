################ plot clusters and landmarks associated with ##################

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
  library(cowplot)
})


############### setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()


cl.df <- read_sheet("https://docs.google.com/spreadsheets/d/1T86EppILF-Q97OjJyd4R2FjUmPUrdYBUEbbiXqgWEoE/edit#gid=1639101745",sheet = "cl.df.v9_230722")
# select only necessary columns 
cl.df <- cl.df %>% 
  dplyr::select(cluster_id_label,
                nt_type_label,
                nt_type_combo_label)


# load landmark annotation per cluster
cl.anat.df <- read_sheet("https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=612105623#gid=612105623", sheet = "cluster_region_annotation")
cl.anat.df <- cl.anat.df %>% 
  select(cluster_id_label,
         broad_region,
         registration_landmark)

broad_ccf_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=1244943357#gid=1244943357", sheet = "broad_landmarks_anno")
broad_ccf_color_palette <- setNames(broad_ccf_color$broad_region_color, broad_ccf_color$broad_region)
# order by graph order
broad_ccf_color <- broad_ccf_color %>% 
  arrange(desc(graph_order)) 


# generate color palette for ccf landmark color
landmark_color <- read_sheet("https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=2094634713#gid=2094634713", sheet = "CCF_landmark_color")
landmark_color_palette <- setNames(landmark_color$registration_landmark_color, landmark_color$registration_landmark)

# order by graph order
landmark_color <- landmark_color %>% 
  arrange(graph_order) 

# generate color palette for cluster labels 
# generate color palette for ccf landmark color
cluster_color <- read_sheet("https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?gid=1882477354#gid=1882477354", sheet = "clusters")
cluster_color_palette <- setNames(cluster_color$cluster_pal, cluster_color$cluster_label)

# order by graph order
cluster_color <- cluster_color %>% 
  arrange(cluster_id) 

# read in metadata file
load("/scratch/638850_metadata_sis.rda")
# load reconstructed coordinates
load("/scratch/638850_reconstructed_coordinates_sis.rda")

metadata_sis <- merge(metadata_sis,
                      coordinates_sis,
                      by = 0)


metadata_sis <- merge(metadata_sis,
                      cl.anat.df,
                      by.x = "flat_CDM_cluster_name",
                      by.y = "cluster_id_label",
                      all.x = T,
                      all.y = F)

##################### plot example sections #####################

metadata_subset <- metadata_sis %>% 
  filter(final_filter == F) %>%
  filter(!is.na(broad_region)) %>% 
  select(cell_id,
         section,
         flat_CDM_cluster_name,
         broad_region,
         registration_landmark,
         x_reconstructed,
         y_reconstructed,
         z_reconstructed)

# only plot hemisections
example_sections <- c("1199650941",
                      "1199650950",
                      "1199650965",
                      "1199650975",
                      "1199650999",
                      "1199651012",
                      "1199651024",
                      "1199651039",
                      "1199651057",
                      "1199651072",
                      "1199651084",
                      "1199651103")

plot_data <- metadata_subset %>% 
  filter(section %in% example_sections) %>% 
  filter(x_reconstructed < 5.6)


plot <- ggplot(plot_data,
               aes(x=x_reconstructed,
                   y=y_reconstructed,
                   color = flat_CDM_cluster_name
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=cluster_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_example_sections_clusters_broad.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)

plot_data <- metadata_subset %>% 
  filter(section %in% example_sections) %>% 
  filter(x_reconstructed > 5.6)

plot <- ggplot(plot_data,
               aes(x=x_reconstructed,
                   y=y_reconstructed,
                   color = broad_region
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=broad_ccf_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_example_sections_ccf_broad.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)

################## plot sections for landmarks ###################

metadata_subset <- metadata_sis %>% 
  filter(final_filter == F) %>%
  filter(!is.na(registration_landmark)) %>% 
  select(cell_id,
         section,
         flat_CDM_cluster_name,
         broad_region,
         registration_landmark,
         x_reconstructed,
         y_reconstructed,
         z_reconstructed)

# only plot hemisections
example_sections <- c("1199650941",
                      "1199650950",
                      "1199650965",
                      "1199650975",
                      "1199650999",
                      "1199651012",
                      "1199651024",
                      "1199651039",
                      "1199651057",
                      "1199651072",
                      "1199651084",
                      "1199651103")

plot_data <- metadata_subset %>% 
  filter(section %in% example_sections) %>% 
  filter(x_reconstructed < 5.6)


plot <- ggplot(plot_data,
               aes(x=x_reconstructed,
                   y=y_reconstructed,
                   color = flat_CDM_cluster_name
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=cluster_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_example_sections_clusters_landmarks.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)

plot_data <- metadata_subset %>% 
  filter(section %in% example_sections) %>% 
  filter(x_reconstructed > 5.6)

plot <- ggplot(plot_data,
               aes(x=x_reconstructed,
                   y=y_reconstructed,
                   color = registration_landmark
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=landmark_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 4, 
                                                  shape = 15))) +
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

# Extract the legend
legend <- get_legend(plot)

# Plot only the legend
legend_plot <- plot_grid(legend)

# Save the legend to a PDF
ggsave("/results/legend_ccf_landmarks.pdf", 
       legend_plot, 
       width = 6, 
       height = 5)

plot <- ggplot(plot_data,
               aes(x=x_reconstructed,
                   y=y_reconstructed,
                   color = registration_landmark
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=landmark_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_example_sections_ccf_landmark.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)
