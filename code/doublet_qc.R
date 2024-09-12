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

example_sections = c("1199650972",
                     "1199651042",
                     "1199651078",
                     "1199651097")

filter_palette <- c("grey","red")

# load corrected cell coordinates
load("/scratch/coordinates_sis.rda")
# load data segmented with in-house cellpose model
load("/scratch/metadata_sis.rda")

# add rotated coordinates to metadata
metadata_sis <- merge(metadata_sis,
                      coordinates_sis,
                      by = 0)

metadata_sis <- metadata_sis %>% 
  rownames_to_column( var = "cell_label")

filtered_metadata_sis <- metadata_sis %>%
  filter(genes_filter == F) %>% 
  filter(mrna_filter == F) %>% 
  filter(blanks_filter == F) %>% 
  select(cell_label,
         x_coordinate,
         y_coordinate,
         volume,
         section,
         n_genes_by_counts,
         total_counts,
         total_counts_Blank,
         pct_counts_Blank,
         genes_filter,
         mrna_filter,
         blanks_filter,
         basic_qc_filter,
         dif,
         doublets_thr,
         doublets_filter,
         incongruous_genes_pct,
         total_counts_per_cell_volume,
         flat_CDM_class_name,
         flat_CDM_class_avg_correlation,
         CDM_class_color,
         flat_CDM_subclass_name,
         flat_CDM_subclass_avg_correlation,
         CDM_subclass_color,
         flat_CDM_supertype_name,
         flat_CDM_supertype_avg_correlation,
         CDM_supertype_color,
         flat_CDM_cluster_name,
         flat_CDM_cluster_avg_correlation,
         CDM_cluster_color)

example_sections_plot <- filtered_metadata_sis %>% 
  filter(section %in% example_sections)

# print diff distribution and set threshold

doublet_threshold <- unique(filtered_metadata_sis$doublets_thr)[1]

plot <- ggplot(filtered_metadata_sis, 
       aes(x = dif)) + 
  geom_histogram(binwidth = .1, 
                 fill = "lightblue") +
  labs(title = "Difference histogram", 
       x = "Dif", 
       y = "Frequency") +
  geom_vline(xintercept = doublet_threshold, 
             linetype = "dashed", 
             color = "black") +
  theme_minimal()

# save plot
ggsave(filename = "/results/doublets_dif_distribution.png", 
       plot = plot, 
       width = 10,
       height = 5, 
       dpi = 300)

# supplemental figure for location of doublets

# plot sections
plot <- ggplot(example_sections_plot,
               aes(x=x_coordinate,
                   y=y_coordinate,
                   color = doublets_filter
               )) +
  geom_point(size=.1,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Doublet Filter") +
  scale_color_manual(values=filter_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sections_doublets_filter.png", 
       plot = plot, 
       width = 5,
       height = 15, 
       dpi = 300)


# have split violin plots for gene, volumes, 

# plot comparison between genes detected
plot <- ggplot(filtered_metadata_sis, aes(x = factor(1), 
                                  y = n_genes_by_counts, 
                                  fill = doublets_filter)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", name = "Cell classification") +
  #ylim(0,1000) +
  labs(x = "", y = "#genes per cell") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/doublet_number_of_genes_qc.pdf", 
       plot = plot, 
       width = 5,
       height = 10, 
       dpi = 300)

# plot comparison between total counts detected
plot <- ggplot(filtered_metadata_sis, 
               aes(x = factor(1), 
                   y = total_counts, 
                   fill = doublets_filter)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", 
                    name = "Cell classification") +
  #ylim(0,1000) +
  labs(x = "", y = "# of total mRNA molecules per cell") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/number_of_spots_qc.png", 
       plot = plot,
       width = 5,
       height = 10, 
       dpi = 300)

# plot comparison between volumes
plot <- ggplot(filtered_metadata_sis, 
       aes(x = factor(1), 
           y = volume, 
           fill = doublets_filter)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", 
                    name = "Cell classification") +
  #ylim(0,1000) +
  labs(x = "", y = "# of total mRNA molecules per cell") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/doublet_volumes_qc.png", 
       plot = plot,
       width = 5,
       height = 10, 
       dpi = 300)

plot <- ggplot(filtered_metadata_sis, 
       aes(x = factor(1), 
           y = total_counts_per_cell_volume, 
           fill = doublets_filter)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .2, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", 
                    name = "Cell classification") +
  #ylim(0,1000) +
  labs(x = "", y = "counts per µm3") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/doublet_cell_densities.png", 
       plot = plot,
       width = 5,
       height = 10, 
       dpi = 300)

# plot incongruent genes

plot <- ggplot(filtered_metadata_sis, 
       aes(x = factor(1), 
           y = incongruous_genes_pct, 
           fill = doublets_filter)) + 
  geom_split_violin(alpha = .4) +
  geom_boxplot(width = .1, 
               alpha = .6, 
               fatten = NULL, 
               show.legend = F) +
  stat_summary(fun = "median", 
               show.legend = F,
               position = position_dodge(.2)) +
  #coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Dark2", 
                    name = "Cell classification") +
  #ylim(0,1000) +
  labs(x = "", y = "counts per µm3") + # Customize axis labels
  theme_minimal() # Remove x-axis elements

ggsave(filename = "/results/doublet_incongruent_genes.png", 
       plot = plot,
       width = 5,
       height = 10, 
       dpi = 300)
