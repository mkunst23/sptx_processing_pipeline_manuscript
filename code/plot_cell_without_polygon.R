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

sis <- read_h5ad("/data/SIS_processed_data/mouse_638850_processed_SIS.h5ad")

# extract coordinates and and add sections
coordinates_sis <- as.matrix(sis$obsm$spatial_rotate)
coordinates_sis <- as.data.frame(coordinates_sis)
colnames(coordinates_sis) <- c("x_coordinate","y_coordinate")
coordinates_sis$section <- sis$obs$section
coordinates_sis$volume <- sis$obs$volume
coordinates_sis$final_filter <- sis$obs$final_filter

final_coordinates <- coordinates_sis %>%
  filter(final_filter == FALSE) %>% 
  filter(is.na(volume))

ggplot(final_coordinates,
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
  facet_wrap(~section, nrow = 6) +
  theme_void()

