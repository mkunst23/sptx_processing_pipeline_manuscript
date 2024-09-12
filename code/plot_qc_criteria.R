############# script to plot basic qc measures ################

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

filter_palette <- c("grey","red")

example_sections = c("1199651009",
                     "1199650956",
                     "1199651045",
                     "1199651112")

# load corrected cell coordinates
load("/scratch/coordinates_sis.rda")
# load data segmented with in-house cellpose model
load("/scratch/metadata_sis.rda")

# add rotated coordinates to metadata
metadata_sis <- merge(metadata_sis,
                      coordinates_sis,
                      by = 0)

# filter out relvevant columns
filtered_metadata_sis <- metadata_sis %>% 
  select(cell_id,
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

################# plot genes detected #################

genes_filter_threshold <- 6

plot <- ggplot(filtered_metadata_sis,
       aes(x = factor(1),
           y = n_genes_by_counts)) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  geom_hline(yintercept = genes_filter_threshold, 
             linetype = "dashed", 
             color = "black") +
  labs(x = "",
       y = "# of detected genes per cell") + 
  theme_minimal()

ggsave(filename = "/results/distribution_genes_filter.pdf", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)


# plot sections
plot <- ggplot(example_sections_plot,
                aes(x=x_coordinate,
                    y=y_coordinate,
                    color = genes_filter
                )) +
  geom_point(size=.1,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Gene Filter") +
  scale_color_manual(values=filter_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sections_genes_filter.pdf", 
       plot = plot, 
       width = 5,
       height = 15, 
       dpi = 300)

################ plot spots detected ###################

spots_filter_threshold <- 30

plot <- ggplot(filtered_metadata_sis,
               aes(x = factor(1),
                   y = total_counts)) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  geom_hline(yintercept = spots_filter_threshold, 
             linetype = "dashed", 
             color = "black") +
  coord_trans(y = "log10") +
  labs(x = "",
       y = "# of detected spots per cell") + 
  theme_minimal()


ggsave(filename = "/results/distribution_spots_filter.pdf", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)



# plot sections
plot <- ggplot(example_sections_plot,
               aes(x=x_coordinate,
                   y=y_coordinate,
                   color = mrna_filter
               )) +
  geom_point(size=.1,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Gene Filter") +
  scale_color_manual(values=filter_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())


ggsave(filename = "/results/sections_spots_filter.pdf", 
       plot = plot, 
       width = 5,
       height = 12, 
       dpi = 300)

################# plot blank percentage ###############

filtered_metadata_sis_2 <- filtered_metadata_sis %>% 
  filter(genes_filter == F) %>% 
  filter(mrna_filter == F)

blank_filter_threshold <- 2

plot <- ggplot(filtered_metadata_sis_2,
               aes(x = factor(1),
                   y = pct_counts_Blank)) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  geom_hline(yintercept = blank_filter_threshold, 
             linetype = "dashed", 
             color = "black") +
  #coord_trans(y = "log2") +
  labs(x = "",
       y = "# of detected spots per cell") + 
  theme_minimal()


ggsave(filename = "/results/distribution_blank_filter.pdf", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)


example_sections_plot_2 <- filtered_metadata_sis_2 %>% 
  filter(section %in% example_sections)

example_sections_plot_2 <- example_sections_plot_2 %>%
  mutate(blank_fixed_filter = pct_counts_Blank > 2)

# plot sections
plot <- ggplot(example_sections_plot_2,
               aes(x=x_coordinate,
                   y=y_coordinate,
                   color = blank_fixed_filter
               )) +
  geom_point(size=.1,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Blank Filter") +
  scale_color_manual(values=filter_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())


ggsave(filename = "/results/sections_blank_filter.pdf", 
       plot = plot, 
       width = 5,
       height = 12, 
       dpi = 300)

############### plot volume distribution ####################

plot <- ggplot(filtered_metadata_sis,
               aes(x = factor(1),
                   y = volume)) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  # geom_hline(yintercept = spots_filter_threshold, 
  #            linetype = "dashed", 
  #            color = "black") +
  #coord_trans(y = "log10") +
  labs(x = "",
       y = "# of detected spots per cell") + 
  theme_minimal()


ggsave(filename = "/results/distribution_volume.pdf", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)

############### plot distributions of qc criteria in vpt segmented cells ##########

# load data segmented with vpt cellpose model
load("/scratch/metadata_vpt.rda")

filtered_metadata_vpt <- metadata_vpt %>% 
  select(cell_label,
         volume,
         section,
         n_genes_by_counts,
         total_counts,
         total_counts_Blank,
         pct_counts_Blank,
         genes_filter,
         mrna_lower_filter,
         mrna_upper_filter,
         blanks_filter,
         basic_qc_filter)

################### plot gene filter #####################

genes_filter_threshold <- 15

plot <- ggplot(filtered_metadata_vpt,
               aes(x = factor(1),
                   y = n_genes_by_counts)) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .05,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  geom_hline(yintercept = genes_filter_threshold, 
             linetype = "dashed", 
             color = "black") +
  labs(x = "",
       y = "# of detected genes per cell") + 
  theme_minimal()

ggsave(filename = "/results/distribution_genes_filter_vpt.pdf", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)

################### plot spot filter ##################

spots_filter_threshold <- 40

plot <- ggplot(filtered_metadata_vpt,
               aes(x = factor(1),
                   y = (total_counts + 1))) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  geom_hline(yintercept = spots_filter_threshold, 
             linetype = "dashed", 
             color = "black") +
  coord_trans(y = "log2") +
  labs(x = "",
       y = "# of detected spots per cell") + 
  theme_minimal()


ggsave(filename = "/results/distribution_spots_filter_vpt.pdf", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)

###################### plot volume filter #####################

gene_volume_threshold <- 100

plot <- ggplot(filtered_metadata_vpt,
               aes(x = factor(1),
                   y = volume)) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  geom_hline(yintercept = gene_volume_threshold,
             linetype = "dashed",
             color = "black") +
  #coord_trans(y = "log10") +
  labs(x = "",
       y = "cell volumes") + 
  theme_minimal()


ggsave(filename = "/results/distribution_volume_vpt.pdf", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)
