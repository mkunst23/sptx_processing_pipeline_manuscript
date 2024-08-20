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
final_vpt <- metadata_vpt %>%
  filter(final_filter == FALSE) %>% 
  select(volume,
         n_genes_by_counts,
         total_counts,
         total_counts_per_cell_volume)

final_vpt$tool <- "VPT"

ggplot(final_vpt, aes(x = factor(1), y = volume)) + 
  geom_violin(fill="lightblue") +
  coord_trans(y = "log10") +
#ylim(0,1000) +
labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

# read in new segmentation
sis <- read_h5ad("/data/SIS_processed_data/mouse_638850_processed_SIS.h5ad")
metadata_sis <- sis$obs
final_sis <- metadata_sis %>% 
  filter(final_filter == FALSE) %>% 
  select(volume,
         n_genes_by_counts,
         total_counts,
         total_counts_per_cell_volume)

final_sis$tool <- "SIS"

ggplot(final_sis, aes(x = factor(1), y = volume)) + 
  geom_violin(fill="lightblue") +
  coord_trans(y = "log10") +
  #ylim(0,1000) +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() # Remove x-axis elements
