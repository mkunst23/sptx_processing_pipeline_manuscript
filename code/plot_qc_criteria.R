############# script to plot basic qc measures ################

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
  library(ttr)
  library(devtools)
  library(RColorBrewer)
  library(googlesheets4)
  library(plyr)
  library(forcats)
})

# load data segmented with in-house cellpose model
load("/scratch/metadata_sis.rda")

# filter out relvevant columns
filtered_metadata_sis <- metadata_sis %>% 
  select()


