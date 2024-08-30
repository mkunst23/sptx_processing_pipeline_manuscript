#TODO wait for incongruent genes
#TODO add different sizes cell types

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
  library(ggrastr)
  library(forcats)
  library(scrattch.vis)
})


options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()

ss <- "https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?usp=sharing"
class_colors <- read_sheet(ss,sheet="classes")
class_color_palette <- setNames(class_colors$class_color, class_colors$class_id_label)
subclass_colors <- read_sheet(ss, sheet="subclasses")
subclass_color_palette <- setNames(subclass_colors$subclass_color, subclass_colors$subclass_id_label)
supertype_colors <- read_sheet(ss, sheet="supertypes")
supertype_color_palette <- setNames(supertype_colors$supertype_color_new, supertype_colors$supertype_id_label)

#TODO add in doublet detection

# read in old segmentation
#vpt <- read_h5ad("/data/merscope_638850_mouseadult_processed_VPT/whole_dataset/mouse_638850_filtered.h5ad")
#metadata_vpt <- vpt$obs
#save(metadata_vpt, file = "/scratch/metadata_vpt.rda")
load("/scratch/metadata_vpt.rda")
final_vpt <- metadata_vpt %>%
  filter(final_filter == FALSE) %>%
  select(volume,
         n_genes_by_counts,
         total_counts,
         total_counts_per_cell_volume,
         flat_CDM_class_name,
         flat_CDM_subclass_name,
         flat_CDM_supertype_name,
         flat_CDM_cluster_name,
         flat_CDM_cluster_avg_correlation,
         section,
         incongruous_genes_pct)

final_vpt$tool <- "VPT"

# read in new segmentation
#sis <- read_h5ad("/data/SIS_processed_data/mouse_638850_processed_SIS.h5ad")
#metadata_sis <- sis$obs
#save(metadata_sis, file = "/scratch/metadata_sis.rda")
load("/scratch/metadata_sis.rda")
final_sis <- metadata_sis %>% 
  filter(final_filter == FALSE) %>%
  filter(!is.na(volume)) %>% 
  select(volume,
         n_genes_by_counts,
         total_counts,
         total_counts_per_cell_volume,
         flat_CDM_class_name,
         flat_CDM_subclass_name,
         flat_CDM_supertype_name,
         flat_CDM_cluster_name,
         flat_CDM_cluster_avg_correlation,
         section,
         incongruous_genes_pct)

final_sis$tool <- "SIS"

# merge datasets
list_of_dfs <- list(final_vpt, final_sis)
combined_meta <- rbindlist(list_of_dfs)

volume_meds <- ddply(combined_meta, .(tool), summarise, med = median(volume))

# plot comparison between volumes
plot <- ggplot(combined_meta, aes(x = factor(1), 
                          y = volume, 
                          fill = tool)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  #ylim(0,1000) +
  labs(x = "", y = "cell volume [µm^3]") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/volume_qc.pdf", 
       plot = plot, 
       width = 5,
       height = 10, 
       dpi = 300)

genes_meds <- ddply(combined_meta, 
                    .(tool), summarise, 
                    med = median(n_genes_by_counts))

# plot comparison between genes detected
plot <- ggplot(combined_meta, aes(x = factor(1), y = n_genes_by_counts, fill = tool)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  #ylim(0,1000) +
  labs(x = "", y = "#genes per cell") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/number_of_genes_qc.pdf", 
       plot = plot, 
       width = 5,
       height = 10, 
       dpi = 300)

counts_meds <- ddply(combined_meta, .(tool), summarise, med = median(total_counts))

# plot comparison between total counts detected
plot <- ggplot(combined_meta, aes(x = factor(1), y = total_counts, fill = tool)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  #ylim(0,1000) +
  labs(x = "", y = "# of total mRNA molecules per cell") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/number_of_spots_qc.png", 
       plot = plot,
       width = 5,
       height = 10, 
       dpi = 72)

density_meds <- ddply(combined_meta, .(tool), summarise, med = median(total_counts_per_cell_volume))

# plot comparison between total counts detected
plot <- ggplot(combined_meta, aes(x = factor(1), y = total_counts_per_cell_volume, fill = tool)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  #ylim(0,1000) +
  labs(x = "", y = "# of total mRNA molecules per cell") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/density_qc.pdf", 
       plot = plot, 
       width = 5,
       height = 10, 
       dpi = 300)

incongruent_meds <- ddply(combined_meta, .(tool), summarise, med = median(incongruous_genes_pct))

# plot comparison between total counts detected
plot <- ggplot(combined_meta, aes(x = factor(1), y = log2(incongruous_genes_pct+1), fill = tool)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  #ylim(0,1000) +
  labs(x = "", y = "% of incongruent genes per cell [log2]") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/incongruent_qc.pdf", 
       plot = plot, 
       width = 5,
       height = 10, 
       dpi = 300)


#TODO plot avg.coorelation split up by cluster
# plot comparison between total counts detected
plot <- ggplot(combined_meta, 
               aes(x = flat_CDM_class_name, 
                   y = flat_CDM_cluster_avg_correlation, 
                   fill = tool)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  labs(x = "", y = "# of total mRNA molecules per cell") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "/results/correlation_by_class_qc.pdf", 
       plot = plot, 
       width = 20,
       height = 8, 
       dpi = 300)

# example small cell 046 Vip Gaba
#                    037 DG Glut
#                    314 CB Granule Glut
# example big cell 022 L5 ET CTX Glut 
#                  058 PAL-STR Gaba-Chol
#                  215 SNc-VTA-RAmb Foxa1 Dopa

# subset dataset to two subclasses
vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name %in% c("046 Vip Gaba","022 L5 ET CTX Glut"))

plot1 <- ggplot(vol_comp_1, 
                aes(x = tool, 
                    y = volume, 
                    fill = flat_CDM_subclass_name)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Subclass") +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name %in% c("046 Vip Gaba","058 PAL-STR Gaba-Chol"))

plot2 <- ggplot(vol_comp_1, 
                aes(x = tool, 
                    y = volume, 
                    fill = flat_CDM_subclass_name)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Subclass") +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name %in% c("046 Vip Gaba","215 SNc-VTA-RAmb Foxa1 Dopa"))

plot3 <- ggplot(vol_comp_1, 
       aes(x = tool, 
           y = volume, 
           fill = flat_CDM_subclass_name)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Subclass") +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plots <- plot_grid(plot1, plot2, plot3, nrow = 1)
ggsave(filename = "/results/small_vs_large_1_qc.pdf", 
       plot = combined_plots, 
       width = 12,
       height = 7, 
       dpi = 300)


# look at random cortical cell types

vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name %in% c("046 Vip Gaba","047 Sncg Gaba"))

plot <- ggplot(vol_comp_1, 
                aes(x = tool, 
                    y = volume, 
                    fill = flat_CDM_subclass_name)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Subclass") +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# look at volume distribution
ggsave(filename = "/results/random_cortical_cell_types_qc.pdf", 
       plot = plot, 
       width = 12,
       height = 7, 
       dpi = 300)

cluster_names <- sort(unique(metadata_sis$flat_CDM_subclass_name))
cluster_names <- as.character(cluster_names)
 
# look at distribution within subclass 047
vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name == cluster_names[8])

plot <- ggplot(vol_comp_1, 
               aes(x = flat_CDM_supertype_name, 
                   y = volume,
                   fill = flat_CDM_supertype_name)) + 
  geom_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F) +
  #coord_trans(y = "log10") +
  ggtitle(cluster_names[8]) +
  scale_fill_manual(values=supertype_color_palette) +
  facet_wrap(~tool) +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fixed_name <-  gsub("/", "-", cluster_names[8])

# look at volume distribution
ggsave(filename = paste0("/results/volume_",fixed_name,".pdf"), 
       plot = plot, 
       width = 12,
       height = 7, 
       dpi = 300)


subclass_name <- cluster_names[25]

vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name == subclass_name)

plot <- ggplot(vol_comp_1, 
               aes(x = flat_CDM_supertype_name, 
                   y = volume,
                   fill = flat_CDM_supertype_name)) + 
  geom_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F) +
  #coord_trans(y = "log10") +
  ggtitle(subclass_name) +
  scale_fill_manual(values=supertype_color_palette) +
  facet_wrap(~tool) +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fixed_name <-  gsub("/", "-", subclass_name)

# look at volume distribution
ggsave(filename = paste0("/results/volume_",fixed_name,".pdf"), 
       plot = plot, 
       width = 12,
       height = 7, 
       dpi = 300)


subclass_name <- cluster_names[69]

vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name == subclass_name)

plot <- ggplot(vol_comp_1, 
               aes(x = flat_CDM_supertype_name, 
                   y = volume,
                   fill = flat_CDM_supertype_name)) + 
  geom_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F) +
  #coord_trans(y = "log10") +
  ggtitle(subclass_name) +
  scale_fill_manual(values=supertype_color_palette) +
  facet_wrap(~tool) +
  labs(x = "", y = "cell volume") + # Customize axis labels
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fixed_name <-  gsub("/", "-", subclass_name)

# look at volume distribution
ggsave(filename = paste0("/results/volume_",fixed_name,".pdf"), 
       plot = plot, 
       width = 12,
       height = 7, 
       dpi = 300)

# plot number of cells per section for sis or vpt
summary_data <- combined_meta %>%
  group_by(section, tool) %>%
  dplyr::summarise(count = n()) %>%
  ungroup()

# Create the bar graph
plot <- ggplot(summary_data, 
       aes(x = fct_rev(section), 
           y = count, 
           fill = tool)) +
  geom_bar(stat = "identity", position = position_dodge2()) +
  labs(title = "Number of cells per section",
       x = "Section",
       y = "# of high quality cells") +
  scale_fill_brewer(palette = "Dark2", name = "Segmentation Tool") +
  theme_minimal() + # Remove x-axis elements
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "/results/cell_per_section.pdf", 
       plot = plot, 
       width = 12,
       height = 5, 
       dpi = 300)


# count subclasses and find the one with the biggest difference

summary_data_sis <- combined_meta %>%
  filter(tool == "SIS") %>% 
  group_by(flat_CDM_subclass_name) %>%
  dplyr::summarise(count = n()) %>%
  ungroup()

# make pie chart
plot <- ggplot(summary_data_sis, aes(x = "", y = count, fill = flat_CDM_subclass_name)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = subclass_color_palette) +
  labs(title = "SIS Segmentation",
       x = NULL,
       y = NULL) +
  theme_void() +
  theme(legend.position = "none")

ggsave(filename = "/results/subclass_distribution_sis.pdf", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 300)

# plot distribution for VPT
summary_data_vpt <- combined_meta %>%
  filter(tool == "VPT") %>% 
  group_by(flat_CDM_subclass_name) %>%
  dplyr::summarise(count = n()) %>%
  ungroup()

plot <- ggplot(summary_data_vpt, aes(x = "", y = count, fill = flat_CDM_subclass_name)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = subclass_color_palette) +
  labs(title = "VPT Segmentation",
       x = NULL,
       y = NULL) +
  theme_void() +
  theme(legend.position = "none")

ggsave(filename = "/results/subclass_distribution_vpt.pdf", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 300)

# merge and calculate ratio
subclass_merge <- merge(summary_data_sis,
                        summary_data_vpt,
                        by = "flat_CDM_subclass_name")

subclass_merge$subclass_merge_ratio <- subclass_merge$count.x/subclass_merge$count.y
subclass_merge$absolute_change <- subclass_merge$count.x - subclass_merge$count.y
