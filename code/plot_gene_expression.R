# script to plot expression of genes in cells

suppressPackageStartupMessages({
  library(anndata)
  library(scrattch.vis)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ttr)
  library(data.table)
  library(scrattch.vis)
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()

# select marker genes:
# exict:Slc17a7,Slc17a6,Nrn1,
# inhib:Slc32a1,Gad2,
# olig:Sox10,Mog
# astro:Aqp4,Gfap,Glis3,
# immune:Arhgap15,Ctss
# vasc: Igf2,Tbx3
marker_genes <- c("Slc17a7","Slc17a6","Nrn1","Slc32a1","Gad2","Sox10","Mog","Aqp4","Gfap","Glis3","Arhgap15","Ctss","Igf2","Tbx3")

anno <- read_sheet("https://docs.google.com/spreadsheets/d/1T86EppILF-Q97OjJyd4R2FjUmPUrdYBUEbbiXqgWEoE/edit?gid=674473196#gid=674473196",sheet = "cl.df.v9_230722")

anno <- anno %>% 
  select(cl,
         supertype_label,
         subclass_label,
         class_label)

ColorPalCluster <- read_sheet("https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?gid=1882477354#gid=1882477354", sheet = "clusters")
ColorPalCluster <- ColorPalCluster %>% 
  select(cl,
         cluster_label,
         cluster_id,
         cluster_color)

ColorPalSupert <- read_sheet("https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?gid=1882477354#gid=1882477354", sheet = "supertypes")
ColorPalSupert <- ColorPalSupert %>% 
  select(supertype_label,
         supertype_id,
         supertype_color)

ColorPalSubclass <- read_sheet("https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?gid=1882477354#gid=1882477354", sheet = "subclasses")
ColorPalSubclass <- ColorPalSubclass %>% 
  select(subclass_label,
         subclass_id,
         subclass_color)

ColorPalClass <- read_sheet("https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?gid=1882477354#gid=1882477354", sheet = "classes")
ColorPalClass <- ColorPalClass %>% 
  select(class_label,
         class_id,
         class_color)

####################### process SIS data ################

# load anndata file for vpt segmentaiton 
vpt <- read_h5ad("/data/merscope_638850_mouseadult_processed_VPT/whole_dataset/mouse_638850_filtered.h5ad")

# extract metadata and keep only cell type assignment
metadata_vpt <- vpt$obs

filtered_metadata_vpt <- metadata_vpt %>%
  filter(final_filter == FALSE) %>%
  select(flat_CDM_cluster_alias) 

# Convert row names into a column
filtered_metadata_vpt <- rownames_to_column(filtered_metadata_vpt, var = "sample_name")

filtered_metadata_vpt <- merge(filtered_metadata_vpt,
                               anno,
                               by.x = "flat_CDM_cluster_alias",
                               by.y = "cl",
                               all.x = T,
                               all.y = F)

filtered_metadata_vpt <- merge(filtered_metadata_vpt,
                               ColorPalCluster,
                               by.x = "flat_CDM_cluster_alias",
                               by.y = "cl",
                               all.x = T,
                               all.y = F)

filtered_metadata_vpt <- merge(filtered_metadata_vpt,
                               ColorPalSupert,
                               by = "supertype_label",
                               all.x = T,
                               all.y = F)

filtered_metadata_vpt <- merge(filtered_metadata_vpt,
                               ColorPalSubclass,
                               by = "subclass_label",
                               all.x = T,
                               all.y = F)

filtered_metadata_vpt <- merge(filtered_metadata_vpt,
                               ColorPalClass,
                               by = "class_label",
                               all.x = T,
                               all.y = F)



# extract count matrix and convert to data frame
data <- vpt$X
data <- as.data.frame(data)
# Convert row names into a column
data <- rownames_to_column(data, var = "sample_name")

# remove unnecessay datasets to save space
rm(vpt)
#rm(metadata_vpt)

# Extract row names from metadata
cells_to_keep <- filtered_metadata_vpt$sample_name

# Subset count matrix to high quality cells
subset_data <- data %>% 
  filter(sample_name %in% cells_to_keep)

# Subset count martrix to marker genes
marker_genes_data <- subset_data %>% 
  select(all_of(marker_genes))

group_violin_plot(t(marker_genes_data), 
                  filtered_metadata_vpt, 
                  genes = c("Slc17a7","Slc32a1","Mog"), 
                  grouping = "class", 
                  log_scale = FALSE,
                  font_size = 5,
                  rotate_counts = TRUE)


# Combine datasets
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

plot_data_vpt$tool <- "VPT"
rm(vpt_data_combined)


################## process SIS data ################

# load anndata file for vpt segmentaiton 
sis <- read_h5ad("/data/merscope_638850_mouseadult_processed_SIS/whole_dataset/mouse_638850_filtered.h5ad")

# extract metadata and keep only cell type assignment
metadata_sis <- sis$obs
filtered_metadata_sis <- metadata_sis %>%
  filter(final_filter == FALSE) %>%
  select(CDM_cluster_name,
         CDM_supertype_name,
         CDM_subclass_name,
         CDM_class_name)

# extract count matrix and convert to data frame
data <- sis$X
data <- as.data.frame(data)

# remove unnecessay datasets to save space
rm(sis)
rm(metadata_sis)

# Extract row names from metadata
row_names_to_keep <- rownames(filtered_metadata_sis)

# Subset count matrix to high quality cells
subset_data <- data[row_names_to_keep, ]

# Subset count martrix to marker genes
marker_genes_data <- subset_data %>% 
  select(all_of(marker_genes))

# Combine datasets
sis_data_combined <- merge(filtered_metadata_sis,
                           marker_genes_data,
                           by = 0)

sis_data_combined$Row.names <- NULL

# Convert from wide to long format using column indices
plot_data_sis <- sis_data_combined %>%
  pivot_longer(
    cols = 5:18,           # Specify the columns to pivot by their indices
    names_to = "gene",    # Name of the new key column
    values_to = "expression"   # Name of the new value column
  )

plot_data_sis$tool <- "SIS"

################# plot combined dataset ###################

# merge datasets
list_of_dfs <- list(plot_data_vpt, plot_data_sis)
plot_data <- rbindlist(list_of_dfs)

ggplot(plot_data, 
       aes(x = CDM_class_name, 
           y = expression, 
           fill = tool)) + 
  geom_split_violin(alpha = .4) +
  # geom_boxplot(width = .2, 
  #              alpha = .6, 
  #              fatten = NULL, 
  #              show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  facet_grid(gene ~ .) +
  labs(x = "", y = "# of total mRNA molecules per cell") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
