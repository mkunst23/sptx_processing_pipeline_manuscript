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
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=20)

# read in old segmentation
#vpt <- read_h5ad("/data/merscope_638850_mouseadult_processed_VPT/whole_dataset/mouse_638850_filtered.h5ad")
#metadata_vpt <- vpt$obs
#save(metadata_vpt, file = "/scratch/metadata_vpt.rda")
load("/scratch/metadata_vpt.rda")

# read in new segmentation
sis <- read_h5ad("/data/merscope_638850_mouseadult_processed_SIS/whole_dataset/mouse_638850_filtered.h5ad")

# extract coordinates and and add sections
coordinates_sis <- as.matrix(sis$obsm$spatial_rotate)
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_coordinate","y_coordinate")
coordinates_sis$section <- sis$obs$section


# set of sections
section_set <- c("1199650950","1199650944","1199650938","1199651036","1199651048")

example_vpt <- metadata_vpt %>% 
  filter(section %in% section_set)

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
