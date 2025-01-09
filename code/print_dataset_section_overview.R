################### plot section overview ###################

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

# read in section info
section_link <- "https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=1482302853#gid=1482302853"
section_info <- read_sheet(section_link, sheet = "section_metadata")

# extract section list
section_list <- section_info$barcode
section_list <- na.omit(section_list)

# 

# for loop for each barcode
i = 5

ggplot(plot_data,
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
