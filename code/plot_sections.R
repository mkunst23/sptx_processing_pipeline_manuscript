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
vpt <- read_h5ad("/data/VPT_processed_data/mouse_638850_processed_VPT.h5ad")

metadata_vpt <- vpt$obs

# read in new segmentation
sis <- read_h5ad("/data/SIS_processed_data/mouse_638850_processed_SIS.h5ad")

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
