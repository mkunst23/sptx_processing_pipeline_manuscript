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
  library(plyr)
})


options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

#TODO add in class label and cluster avg correlation
#TODO add in doublet detection

# read in old segmentation
vpt <- read_h5ad("/data/VPT_processed_data/mouse_638850_processed_VPT.h5ad")
metadata_vpt <- vpt$obs
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
         flat_CDM_cluster_avg_correlation)

final_vpt$tool <- "VPT"

# read in new segmentation
sis <- read_h5ad("/data/SIS_processed_data/mouse_638850_processed_SIS.h5ad")
metadata_sis <- sis$obs
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
         flat_CDM_cluster_avg_correlation)

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

ggsave(filename = "/results/number_of_spots_qc.pdf", 
       plot = plot, 
       width = 5,
       height = 10, 
       dpi = 300)

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
  filter(flat_CDM_subclass_name %in% c("046 Vip Gaba","047 Sncg Gaba"))

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




vol_comp_1 <- combined_meta %>% 
  filter(flat_CDM_subclass_name %in% c("037 DG Glut","022 L5 ET CTX Glut"))

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
  filter(flat_CDM_subclass_name %in% c("037 DG Glut","058 PAL-STR Gaba-Chol"))

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
  filter(flat_CDM_subclass_name %in% c("037 DG Glut","215 SNc-VTA-RAmb Foxa1 Dopa"))

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

ggsave(filename = "/results/small_vs_large_2_qc.pdf", 
       plot = combined_plots, 
       width = 12,
       height = 7, 
       dpi = 300)

 