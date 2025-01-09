
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
  library(car)
})


options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

load("/scratch/638850_metadata_sis.rda")
load("/scratch/metadata_vpt.rda")