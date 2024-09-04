# script to plot expression of genes in cells

suppressPackageStartupMessages({
  library(anndata)
  library(scrattch.vis)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ttr)
})

# select marker genes:
# exict:Slc17a7,Slc17a6,Nrn1,
# inhib:Slc32a1,Gad2,
# olig:Sox10,Mog
# astro:Aqp4,Gfap,Glis3,
# immune:Arhgap15,Ctss
# vasc: Igf2,Tbx3
marker_genes <- c("Slc17a7","Slc17a6","Nrn1","Slc32a1","Gad2","Sox10","Mog","Aqp4","Gfap","Glis3","Arhgap15","Ctss","Igf2","Tbx3")

vpt <- read_h5ad("/data/merscope_638850_mouseadult_processed_VPT/whole_dataset/mouse_638850_filtered.h5ad")
metadata_vpt <- vpt$obs
filtered_metadata_vpt <- metadata_vpt %>%
  filter(final_filter == FALSE) %>%
  select(CDM_cluster_name,
         CDM_supertype_name,
         CDM_subclass_name,
         CDM_class_name)
rm(vpt)
rm(metadata_vpt)
data <- vpt$X
data <- as.data.frame(data)

# Extract row names from df2
row_names_to_keep <- rownames(filtered_metadata_vpt)

# Subset df1 using the row names from df2
subset_data <- data[row_names_to_keep, ]

marker_genes_data <- data %>% 
  select(all_of(marker_genes))

vpt_data_combined <- merge(filtered_metadata_vpt,
                           marker_genes_data,
                           by = 0)

vpt_data_combined$Row.names <- NULL

# Convert from wide to long format using column indices
plot_data_vpt <- vpt_data_combined %>%
  pivot_longer(
    cols = 5:18,           # Specify the columns to pivot by their indices
    names_to = "gene",    # Name of the new key column
    values_to = "expression"   # Name of the new value column
  )

plot_data_vpt$tool <- "SIS"
