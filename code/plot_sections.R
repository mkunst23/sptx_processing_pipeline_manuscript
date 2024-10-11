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
  library(reticulate)
  library(anndata)
  library(googlesheets4)
  library(randomcoloR)
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()

sd <- "https://docs.google.com/spreadsheets/d/1UpfDYDbq3s0Jeh8cLg9Y0yhPCJxH_H8WG4Y0y0_mNHg/edit?gid=70495470#gid=70495470"
spatial_domains <- read_sheet(sd, sheet = "knn_8_res_1-2_30_687997_1")
spatial_domains <- spatial_domains[,1:4]

ccf_broad <- read_sheet(sd, sheet = "CCF_broad_region_color")
ccf_broad_palette <- setNames(ccf_broad$color_hex_tripet,ccf_broad$ccf_broad)

anat <- "https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=2094634713#gid=2094634713"
ccf_landmarks <- read_sheet(anat, sheet = "CCF_landmark_color")
ccf_landmark_palette <- setNames(ccf_landmarks$Color,ccf_landmarks$CCF_acronym)

# read in old segmentation
load("/scratch/687997_metadata_sis.rda")

# read in coordinates
load("/scratch/687997_coordinates_sis.rda")

metadata_sis <- merge(metadata_sis,
                      coordinates_sis,
                      by.x = "production_cell_id",
                      by.y = 0)


# create colormap for leiden clustering
leiden_cluster_colors <- randomColor(length(unique(metadata_sis$leiden_res_1.2_knn_8)), luminosity="dark")
leiden_cluster_palette <- setNames(leiden_cluster_colors, unique(metadata_sis$leiden_res_1.2_knn_8))


metadata_subset <- metadata_sis %>% 
  filter(final_filter == F) %>%
  select(production_cell_id,
         section,
         x_coordinate,
         y_coordinate,
         leiden_res_1.2_knn_8,
         CCF_level1,
         CCF_level2)

  
metadata_subset <- merge(metadata_subset,
                           spatial_domains,
                           by.x = "leiden_res_1.2_knn_8",
                           by.y = "Cluster_id",
                           all.x = T,
                           all.y = F)
    
metadata_subset$section <- as.character(metadata_subset$section)

sections <- unique(metadata_subset$section)


for (i in sections) {
  # subset dataset to sections
  plot_data <- metadata_subset %>% 
    filter(section == i)
  
  plot <- ggplot(plot_data,
                 aes(x=x_coordinate,
                     y=y_coordinate,
                     color = leiden_res_1.2_knn_8
                 )) +
    geom_point(size=.5,
               stroke=0,
               shape=19,) +
    coord_fixed() +
    scale_color_manual(values=leiden_cluster_palette) +
    scale_y_reverse() +
    theme(legend.position = "none") +
    guides(color = "none") +
    theme_void() +
    theme(strip.text = element_blank())
  
  ggsave(filename = paste0("/results/spatial_domain_",i,".png"), 
         plot = plot, 
         width = 6,
         height = 6, 
         dpi = 300)
}




# plot sections
plot1 <- ggplot(example_vpt,
       aes(x=corrected_x,
           y=corrected_y,
           )) +
  geom_point(size=.1,stroke=0,shape=19, color = "grey") +
  coord_fixed() +
  ggtitle("VPT Segmentation") +
  #scale_color_manual(values=CCF_regions_color) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~section, nrow = 1) +
  theme_void()

example_sis <-coordinates_sis %>% 
  filter(section %in% section_set)

# plot sections
plot2 <- ggplot(example_sis,
       aes(x=x_coordinate,
           y=y_coordinate,
       )) +
  geom_point(size=.1,stroke=0,shape=19, color = "grey") +
  coord_fixed() +
  ggtitle("SIS Segmentation") +
  #scale_color_manual(values=CCF_regions_color) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~section, nrow = 1) +
  theme_void()

combined_plot <- plot_grid(plot1, plot2, nrow = 2)
combined_plot

ggsave(filename = "/results/example_sections.pdf", 
       plot = combined_plot, 
       width = 10,
       height = 4, 
       dpi = 160)



# Define the coordinates for the rectangle
xmin <- 500
xmax <- 2000
ymin <- 1500
ymax <- 3000

plot <- ggplot(subset(example_vpt, section == "1199650944"),
       aes(x=corrected_x,
           y=corrected_y,
       )) +
  geom_point(size=.5,stroke=0,shape=19, color = "grey") +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = NA, color = "red", linewidth = 1) +  # Add open rectangle
  coord_fixed() +
  #scale_color_manual(values=CCF_regions_color) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void()

ggsave(filename = "/results/vpt_example_section_1199650944.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)

high_zoom <- example_vpt %>% 
  filter(section == "1199650944") %>% 
  filter(corrected_x > xmin) %>% 
  filter(corrected_x < xmax)%>% 
  filter(corrected_y > ymin) %>% 
  filter(corrected_y < ymax) 
  

plot <- ggplot(high_zoom,
       aes(x=corrected_x,
           y=corrected_y,
       )) +
  geom_point(size=1,stroke=0,shape=19, color = "grey") +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = NA, color = "red", linewidth = 1) +  # Add open rectangle
  coord_fixed() +
  xlim(xmin,xmax) +
  #scale_color_manual(values=CCF_regions_color) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void()

ggsave(filename = "/results/vpt_example_section_1199650944_zoom.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)



# Define the coordinates for the rectangle
xmin <- -3000
xmax <- -1500
ymin <- -1200
ymax <- 300

plot <- ggplot(subset(coordinates_sis, section == "1199650944"),
               aes(x=x_coordinate,
                   y=y_coordinate,
               )) +
  geom_point(size=.5,stroke=0,shape=19, color = "grey") +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = NA, color = "red", linewidth = 1) +  # Add open rectangle
  coord_fixed() +
  #scale_color_manual(values=CCF_regions_color) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void()

ggsave(filename = "/results/sis_example_section_1199650944.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)

high_zoom <- coordinates_sis %>% 
  filter(section == "1199650944") %>% 
  filter(x_coordinate > xmin) %>% 
  filter(x_coordinate < xmax)%>% 
  filter(y_coordinate > ymin) %>% 
  filter(y_coordinate < ymax) 


plot <- ggplot(high_zoom,
               aes(x=x_coordinate,
                   y=y_coordinate,
               )) +
  geom_point(size=1,stroke=0,shape=19, color = "grey") +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = NA, color = "red", linewidth = 1) +  # Add open rectangle
  coord_fixed() +
  xlim(xmin,xmax) +
  #scale_color_manual(values=CCF_regions_color) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void()

ggsave(filename = "/results/sis_example_section_1199650944_zoom.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)
