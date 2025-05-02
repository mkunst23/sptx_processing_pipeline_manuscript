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

# rename cell label column to sample name
metadata_vpt <- metadata_vpt %>% 
  rename(sample_name = cell_label)

filtered_metadata_vpt <- metadata_vpt %>%
  filter(final_filter == FALSE) %>%
  select(CDM_cluster_alias, sample_name) 


filtered_metadata_vpt <- merge(filtered_metadata_vpt,
                               anno,
                               by.x = "CDM_cluster_alias",
                               by.y = "cl",
                               all.x = T,
                               all.y = F)

filtered_metadata_vpt <- merge(filtered_metadata_vpt,
                               ColorPalCluster,
                               by.x = "CDM_cluster_alias",
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

# add vpt to end of each class name
filtered_metadata_vpt <- filtered_metadata_vpt %>%
  mutate(class_label = paste0(class_label, " vpt"))

group_dot_plot(subset_data_vpt, 
               filtered_metadata_vpt, 
               genes = c("Slc17a7","Slc17a6","Nrn1",
                         "Slc32a1","Gad2",
                         "Sox10","Mog",
                         "Aqp4","Gfap","Glis3",
                         "Igf2","Tbx3",
                         "Arhgap15","Ctss"), 
               grouping = "class", 
               log_scale = TRUE,
               font_size = 12,
               max_size = 10,
               rotate_counts = T)


rm(data_vpt, metadata_vpt)
gc()

####################### process SIS data ##############################

# load anndata file for vpt segmentaiton 
#sis <- read_h5ad("/data/merscope_638850_mouseadult_processed_SIS/whole_dataset/mouse_638850_filtered.h5ad")
# extract metadata and keep only cell type assignment
#metadata_sis <- sis$obs
load("/scratch/metadata_sis.rda")

metadata_sis <- metadata_sis %>% 
  rename(sample_name = production_cell_id)

filtered_metadata_sis <- metadata_sis %>%
  filter(final_filter == FALSE) %>%
  select(sample_name,CDM_cluster_alias) 

filtered_metadata_sis <- merge(filtered_metadata_sis,
                               anno,
                               by.x = "CDM_cluster_alias",
                               by.y = "cl",
                               all.x = T,
                               all.y = F)

filtered_metadata_sis <- merge(filtered_metadata_sis,
                               ColorPalCluster,
                               by.x = "CDM_cluster_alias",
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

# add vpt to end of each class name
filtered_metadata_sis <- filtered_metadata_sis %>%
  mutate(class_label = paste0(class_label, " sis"))

filtered_metadata_sis$class_id <- filtered_metadata_sis$class_id + 45

group_dot_plot(subset_data_sis, 
               filtered_metadata_sis, 
               genes = c("Slc17a7","Slc17a6","Nrn1",
                         "Slc32a1","Gad2",
                         "Sox10","Mog",
                         "Aqp4","Gfap","Glis3",
                         "Igf2","Tbx3",
                         "Arhgap15","Ctss"), 
               grouping = "class", 
               log_scale = TRUE,
               font_size = 18,
               max_size = 15,
               rotate_counts = T)

rm(data_sis, metadata_sis)
gc()

################### Combine datasets ####################

data_combined <- rbind(subset_data_vpt,
                       subset_data_sis)

metadata_combined <- rbind(filtered_metadata_vpt,
                           filtered_metadata_sis)

# extract class label as data frame
class_labels_order <- unique(as.data.frame(metadata_combined$class_label))
colnames(class_labels_order) <- "class_label"
class_labels_order$class_id <- 1:68
class_labels_order$class_label <- factor(class_labels_order$class_label,
                             levels = c("IT-ET Glut sis",
                                        "IT-ET Glut vpt",
                                        "NP-CT-L6b Glut sis",
                                        "NP-CT-L6b Glut vpt",
                                        "OB-CR Glut sis",
                                        "OB-CR Glut vpt",
                                        "DG-IMN Glut sis",
                                        "DG-IMN Glut vpt",
                                        "CNU-HYa Glut sis",
                                        "CNU-HYa Glut vpt",
                                        "HY Glut sis",
                                        "HY Glut vpt",
                                        "HY MM Glut sis",
                                        "HY MM Glut vpt",
                                        "HY Gnrh1 Glut sis",
                                        "HY Gnrh1 Glut vpt",
                                        "TH Glut sis",
                                        "TH Glut vpt",
                                        "MH-LH Glut sis",
                                        "MH-LH Glut vpt",
                                        "Pineal Glut sis",
                                        "Pineal Glut vpt",
                                        "MB Glut sis",
                                        "MB Glut vpt",
                                        "MB Dopa sis",
                                        "MB Dopa vpt",
                                        "MB-HB Sero sis",
                                        "MB-HB Sero vpt",
                                        "P Glut sis",
                                        "P Glut vpt",
                                        "MY Glut sis",
                                        "MY Glut vpt",
                                        "CB Glut sis",
                                        "CB Glut vpt",
                                        "OB-IMN GABA sis",
                                        "OB-IMN GABA vpt",
                                        "CTX-CGE GABA sis",
                                        "CTX-CGE GABA vpt",
                                        "CTX-MGE GABA sis",
                                        "CTX-MGE GABA vpt",
                                        "LSX GABA sis",
                                        "LSX GABA vpt",
                                        "CNU-LGE GABA sis",
                                        "CNU-LGE GABA vpt",
                                        "CNU-MGE GABA sis",
                                        "CNU-MGE GABA vpt",
                                        "CNU-HYa GABA sis",
                                        "CNU-HYa GABA vpt",
                                        "HY GABA sis",
                                        "HY GABA vpt",
                                        "MB GABA sis",
                                        "MB GABA vpt",
                                        "P GABA sis",
                                        "P GABA vpt",
                                        "MY GABA sis",
                                        "MY GABA vpt",
                                        "CB GABA sis",
                                        "CB GABA vpt",
                                        "Astro-Epen sis",
                                        "Astro-Epen vpt",
                                        "OPC-Oligo sis",
                                        "OPC-Oligo vpt",
                                        "OEC sis",
                                        "OEC vpt",
                                        "Immune sis",
                                        "Immune vpt",
                                        "Vascular sis",
                                        "Vascular vpt"))

class_labels_order <- class_labels_order %>%
  arrange(class_label)

class_labels_order$class_id <- 1:68

metadata_combined$class_id <- NULL

metadata_combined <- merge(metadata_combined,
                           class_labels_order,
                           by = "class_label",
                           all.x = T,
                           all.y = F)

plot <- group_dot_plot(data_combined, 
               metadata_combined, 
               genes = c("Slc17a7","Slc17a6","Nrn1","Adcyap1","Nhlh2",
                         "Slc32a1","Gad2",
                         "Sox10","Mog",
                         "Aqp4","Gfap","Glis3",
                         "Ighm","Lsp1","Arhgap15","Ctss",
                         "Igf2","Tbx3"), 
               grouping = "class", 
               log_scale = TRUE,
               font_size = 14,
               max_size = 15,
               rotate_counts = T)

ggsave(filename = "/results/gene_expression_comparison.pdf", 
       plot = plot, 
       width = 20,
       height = 10, 
       dpi = 300)


########## subset to only 1 gene and two classes
  
# list of classes 
class_list <- c("IT-ET Glut sis",
                "IT-ET Glut vpt",
                "CTX-CGE GABA sis",
                "CTX-CGE GABA vpt",
                "CTX-MGE GABA sis",
                "CTX-MGE GABA vpt")

metadata_combined_subset <- metadata_combined %>% 
  filter(class_label %in% class_list)

cells_to_keep2 <- metadata_combined_subset %>% 
  pull(sample_name)

data_combined_subset <- data_combined %>% 
  filter(sample_name %in% cells_to_keep2)

plot <- group_dot_plot(data_combined_subset, 
                       metadata_combined_subset, 
                       genes = c("Slc17a7","Nrn1"), 
                       grouping = "class", 
                       log_scale = TRUE,
                       font_size = 14,
                       max_size = 30,
                       rotate_counts = T)

ggsave(filename = "/results/gene_expression_comparison_subset.pdf", 
       plot = plot, 
       width = 20,
       height = 10, 
       dpi = 300)

########## legacy code for split violin plot ##############  
  
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
