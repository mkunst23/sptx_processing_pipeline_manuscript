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
  library(forcats)
})


options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()

############# plot section overview #############

# read in section info
section_link <- "https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=1482302853#gid=1482302853"
section_info <- read_sheet(section_link, sheet = "section_metadata")

# extract section list
section_list <- section_info$barcode
section_list <- na.omit(section_list)

# load metadata
load("/scratch/638850_metadata_sis.rda")

# select relevant columns
metadata_plot <- metadata_sis %>% 
  filter(final_filter == FALSE) %>% 
  select(production_cell_id,
         section,
         volume_x,
         volume_y)

load("/scratch/638850_coordinates_sis.rda")

metadata_plot <- merge(metadata_plot,
                       coordinates_sis,
                       by.x = "production_cell_id",
                       by.y = 0,
                       all.x = T,
                       all.y = F)

for (i in section_list) {
  # subset dataset to sections
  plot_data <- metadata_plot %>% 
    filter(section == i)
  
  plot <- ggplot(plot_data,
                 aes(x=x_coordinate,
                     y=y_coordinate
                 )) +
    geom_point(size=.5,
               stroke=0,
               shape=19,) +
    coord_fixed() +
    scale_y_reverse() +
    theme(legend.position = "none") +
    guides(color = "none") +
    theme_void() +
    theme(strip.text = element_blank())
  
  ggsave(filename = paste0("/results/section_",i,".png"), 
         plot = plot, 
         width = 10,
         height = 10, 
         dpi = 300)
}
