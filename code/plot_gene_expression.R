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
  library(googlesheets4)
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

# list of classes 
class_list <- c("IT-ET Glut","CTX-CGE GABA","CTX-MGE GABA")

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

####################### process VPT data ################

# load anndata file for vpt segmentaiton 
#vpt <- read_h5ad("/data/merscope_638850_mouseadult_processed_VPT/whole_dataset/mouse_638850_filtered.h5ad")

# extract metadata and keep only cell type assignment
#metadata_vpt <- vpt$obs
load("/scratch/metadata_vpt.rda")

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
# data_vpt <- vpt$X
# data_vpt <- as.data.frame(data_vpt)
# Convert row names into a column
# data_vpt <- rownames_to_column(data_vpt, var = "sample_name")
# save(data_vpt, file = "/scratch/data_vpt.rda")
load("/scratch/data_vpt.rda")

# remove unnecessay datasets to save space
#rm(vpt)

#rm(metadata_vpt)

# Extract row names from metadata
cells_to_keep <- filtered_metadata_vpt$sample_name

# Subset count matrix to high quality cells
subset_data_vpt <- data_vpt %>% 
  filter(sample_name %in% cells_to_keep)

# Subset count martrix to marker genes
marker_genes_data <- subset_data_vpt %>% 
  select(sample_name,
         all_of(marker_genes))

# list of subclasses

plot_anno_vpt <-  filtered_metadata_vpt %>% 
  filter(subclass_label == "Sncg Gaba")

plot_data_vpt <- subset_data_vpt %>% 
  filter(sample_name %in% plot_anno_vpt$sample_name)


group_violin_plot(plot_data_vpt, 
                  plot_anno_vpt, 
                  genes = c("Slc17a7","Slc17a6","Nrn1",
                            "Slc32a1","Gad2",
                            "Sox10","Mog",
                            "Aqp4","Gfap","Glis3",
                            "Igf2","Tbx3",
                            "Arhgap15","Ctss"), 
                  grouping = "class", 
                  log_scale = FALSE,
                  font_size = 8,
                  rotate_counts = TRUE)


group_heatmap_plot(subset_data_vpt, 
                   filtered_metadata_vpt, 
                   genes = c("Slc17a7","Slc17a6","Nrn1",
                             "Slc32a1","Gad2",
                             "Sox10","Mog",
                             "Aqp4","Gfap","Glis3",
                             "Igf2","Tbx3",
                             "Arhgap15","Ctss"), 
                   grouping = "class", 
                   stat = "tmean",
                   log_scale = TRUE,
                   font_size = 8,
                   rotate_counts = TRUE)

################### prepare datasets for combining ######################

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


####################### process SIS data ##############################

# load anndata file for vpt segmentaiton 
#sis <- read_h5ad("/data/merscope_638850_mouseadult_processed_SIS/whole_dataset/mouse_638850_filtered.h5ad")
# extract metadata and keep only cell type assignment
#metadata_sis <- sis$obs
load("/scratch/metadata_sis.rda")

filtered_metadata_sis <- metadata_sis %>%
  filter(final_filter == FALSE) %>%
  select(flat_CDM_cluster_alias) 

# Convert row names into a column
filtered_metadata_sis <- rownames_to_column(filtered_metadata_sis, var = "sample_name")

filtered_metadata_sis <- merge(filtered_metadata_sis,
                               anno,
                               by.x = "flat_CDM_cluster_alias",
                               by.y = "cl",
                               all.x = T,
                               all.y = F)

filtered_metadata_sis <- merge(filtered_metadata_sis,
                               ColorPalCluster,
                               by.x = "flat_CDM_cluster_alias",
                               by.y = "cl",
                               all.x = T,
                               all.y = F)

filtered_metadata_sis <- merge(filtered_metadata_sis,
                               ColorPalSupert,
                               by = "supertype_label",
                               all.x = T,
                               all.y = F)

filtered_metadata_sis <- merge(filtered_metadata_sis,
                               ColorPalSubclass,
                               by = "subclass_label",
                               all.x = T,
                               all.y = F)

filtered_metadata_sis <- merge(filtered_metadata_sis,
                               ColorPalClass,
                               by = "class_label",
                               all.x = T,
                               all.y = F)



# extract count matrix and convert to data frame
# data_sis <- sis$X
# data_sis <- as.data.frame(data_sis)
# Convert row names into a column
# data_sis <- rownames_to_column(data_sis, var = "sample_name")
# save(data_sis, file = "/scratch/data_sis.rda")
load("/scratch/data_sis.rda")

# remove unnecessay datasets to save space
#rm(sis)

# Extract row names from metadata
cells_to_keep <- filtered_metadata_sis$sample_name

# Subset count matrix to high quality cells
subset_data_sis <- data_sis %>% 
  filter(sample_name %in% cells_to_keep)

# Subset count martrix to marker genes
marker_genes_data <- subset_data_sis %>% 
  select(sample_name,
         all_of(marker_genes))

plot_anno_sis <-  filtered_metadata_sis %>% 
  filter(class_label %in% class_list)

plot_data_sis <- subset_data_sis %>% 
  filter(sample_name %in% plot_anno_sis$sample_name)


group_violin_plot(plot_data_sis, 
                  plot_anno_sis, 
                  genes = c("Slc17a7","Slc17a6","Nrn1",
                            "Slc32a1","Gad2",
                            "Sox10","Mog",
                            "Aqp4","Gfap","Glis3",
                            "Igf2","Tbx3",
                            "Arhgap15","Ctss"), 
                  grouping = "subclass", 
                  log_scale = FALSE,
                  font_size = 8,
                  rotate_counts = TRUE)


################### Combine datasets ####################
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
