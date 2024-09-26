########## region analysis plotting #############

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)   
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(anndata)
  library(googlesheets4)
  library(data.table)
  library(reshape)
  library(forcats)
  library(entropy)
  library(aplot)
  library(vegan)
  library(ggridges)
  library(scales)
  library(DescTools)
  library(RColorBrewer)
  library(philentropy)
  library(igraph)
  library(paletteer)
})

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

# function to calculate gini coefficient for each row of a data frame
calcGini <- function(row) {
  Gini(row, unbiased = FALSE)
}

# read in metadata file
load("/scratch/638850_metadata_sis.rda")

# general data structure maintenance
metadata_sis$subclass_id_label <- as.character(meta$subclass_id_label)
metadata_sis$supertype_id_label <- as.character(meta$supertype_id_label)
metadata_sis$cluster_id_label <- as.character(meta$cluster_id_label)
