######### script to describe the result of doublet detection ###############

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

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

# load corrected cell coordinates
load("/scratch/coordinates_sis.rda")
# load data segmented with in-house cellpose model
load("/scratch/metadata_sis.rda")

# add rotated coordinates to metadata
metadata_sis <- merge(metadata_sis,
                      coordinates_sis,
                      by = 0)

