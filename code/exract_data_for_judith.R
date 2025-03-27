##################### extract data for judith ##########################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(anndata)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(googlesheets4)
  library(data.table)
  library(reshape)
  library(forcats)
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()


example_sections <- c("1199650981",
                      "1199650978",
                      "1199650975",
                      "1199650972",
                      "1199650968",
                      "1199650965",
                      "1199650962",
                      "1199650959",
                      "1199650956",
                      "1199650953",
                      "1199650950",
                      "1199650944",
                      "1199650941",
                      "1199650938",
                      "1199650935",
                      "1199650932",
                      "1199650929")


# read in metadata file
metadata <- fread("/data/merscope_638850_mouseadult_registered_v2/whole_dataset/mouse_638850_registered.csv")

metadata$section <- as.character(metadata$section)

metadata <- metadata %>% 
  filter(final_qc_passed == T) %>%
  filter(section %in% example_sections) %>% 
  select(production_cell_id,
         section,
         hrc_mmc_class_name,
         hrc_mmc_subclass_name,
         hrc_mmc_supertype_name,
         hrc_mmc_cluster_name,
         center_x,
         center_y)

# save new metadata file
fwrite(metadata,
       "/scratch/medulla_section_for_registration.csv",
       row.names = F)

