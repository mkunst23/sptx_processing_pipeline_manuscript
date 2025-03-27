############## script to plot spatial domain sections #####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(cowplot)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(data.table)
  library(rearrr)
  library(reticulate)
  library(anndata)
  library(googlesheets4)
})

############### setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

# load spatial domains
sd.df <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domains")
# select only necessary columns 
sd.df <- sd.df %>% 
  dplyr::select(Cluster_id,
                spatial_domain_level_1,
                spatial_domain_level_2,
                graph_order,
                ccf_broad)

spatial_domain_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domains")
spatial_domain_palette <- setNames(spatial_domain_color$spatial_domain_level_2_color, spatial_domain_color$spatial_domain_level_2)

# order by graph order
spatial_domain_color <- spatial_domain_color %>% 
  arrange(desc(graph_order))

spatial_domain_1_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domain_level_1_color")
spatial_domain_1_palette <- setNames(spatial_domain_1_color$spatial_domain_level_1_color, spatial_domain_1_color$spatial_domain_level_1)

# order by graph order
spatial_domain_1_color <- spatial_domain_1_color %>% 
  arrange(desc(graph_order))


ss <- "https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?usp=sharing"
class_colors <- read_sheet(ss,sheet="classes")
class_color_palette <- setNames(class_colors$class_color, class_colors$class_id_label)

# section to plot
# LQ_1: 1199651045
# Borders: 1199651018 1199651029 1199651002

# load in metadata file
metadata <- fread("/data/merscope_638850_mouseadult_registered_v2/whole_dataset/mouse_638850_registered.csv")
metadata$section <- as.character(metadata$section)

# filter out only relevant columns
metadata <- metadata %>% 
  select(production_cell_id,
         n_genes_by_counts,
         total_counts,
         section,
         final_qc_passed,
         hrc_mmc_subclass_name,
         hrc_mmc_class_name,
         volume_x,
         volume_y,
         volume_z)

sd_domains <- fread('/scratch/mouse_638850_sd.csv')

metadata <- merge(metadata,
                  sd_domains,
                  by = 'production_cell_id',
                  all.x = T,
                  all.y = F)

metadata <- merge(metadata,
                  sd.df,
                  by.x = "leiden_res_1.4_knn_8",
                  by.y = "Cluster_id",
                  all.x = T,
                  all.y = F)

metadata_subset <- metadata %>% 
  filter(final_qc_passed == T) %>%
  filter(!is.na(spatial_domain_level_1))

################## plot low quality spatial domains ########################

# select sections and domains

example_sections <- c("1199651045")
example_domains <- c("LQ_1")

metadata_plot <- metadata_subset %>% 
  filter(section %in% example_sections)

# Add a new column for color and transparency
metadata_plot <- metadata_plot %>%
  mutate(
    plot_color = ifelse(spatial_domain_level_2 %in% example_domains, 
                        spatial_domain_palette[spatial_domain_level_2], 
                        "grey"),
    plot_alpha = ifelse(spatial_domain_level_2 %in% example_domains, 1, 0.3)
  )

plot <- ggplot(metadata_plot,
               aes(x = volume_x,
                   y = volume_y,
                   color = plot_color,
                   alpha = plot_alpha)) +
  geom_point(size = 1,
             stroke = 0,
             shape = 19) +
  coord_fixed() +
  scale_color_identity() +  # Use the exact colors specified in the data
  scale_alpha_identity() +  # Use the exact alpha values specified in the data
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none", alpha = "none") +
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_LQ_1_sections.png", 
       plot = plot, 
       width = 6,
       height = 4, 
       dpi = 160)

qc_plot <- metadata_subset %>% 
  filter(spatial_domain_level_2 %in% example_domains)

plot <- ggplot(qc_plot,
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
  labs(x = "",
       y = "# of detected genes per cell") + 
  theme_minimal()

ggsave(filename = "/results/lq_genes_per_cell.png", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)

median_gene_lq <- median(qc_plot$n_genes_by_counts)

plot <- ggplot(qc_plot,
               aes(x = factor(1),
                   y = log(total_counts))) +
  geom_violin(fill = "red") +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  labs(x = "",
       y = "# of detected genes per cell") + 
  theme_minimal()

ggsave(filename = "/results/lq_counts_per_cell.png", 
       plot = plot, 
       width = 5,
       height = 8, 
       dpi = 300)

median_spot_lq <- median(qc_plot$total_counts)


################# plot border spatial domains #################

example_sections <- c("1199651018","1199651039","1199651002")
example_domains <- c("Borders_1","Borders_2")

metadata_plot <- metadata_subset %>% 
  filter(section %in% example_sections)

# Add a new column for color and transparency
metadata_plot <- metadata_plot %>%
  mutate(
    plot_color = ifelse(spatial_domain_level_2 %in% example_domains, 
                        spatial_domain_palette[spatial_domain_level_2], 
                        "grey"),
    plot_alpha = ifelse(spatial_domain_level_2 %in% example_domains, 1, 0.3)
  )

plot <- ggplot(metadata_plot,
               aes(x = volume_x,
                   y = volume_y,
                   color = plot_color,
                   alpha = plot_alpha)) +
  geom_point(size = 0.5,
             stroke = 0,
             shape = 19) +
  coord_fixed() +
  scale_color_identity() +  # Use the exact colors specified in the data
  scale_alpha_identity() +  # Use the exact alpha values specified in the data
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none", alpha = "none") +
  facet_wrap(~fct_rev(section), nrow = 3) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_borders_sections.png", 
       plot = plot, 
       width = 6,
       height = 4, 
       dpi = 160)

################# plot stacked bar graph for classes ##################

class_distribution <- metadata_subset %>% 
  filter(spatial_domain_level_2 %in% example_domains)

# Create a summary of counts for each class_label
class_counts <- class_distribution %>%
  group_by(hrc_mmc_class_name) %>%
  summarise(count = n()) %>%
  ungroup()


# Calculate positions for labels
class_counts <- class_counts %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    label_position = (ymax + ymin) / 2,  # Position for the label
    label = paste0(hrc_mmc_class_name, " (", round(fraction * 100, 1), "%)")  # Label text
  )


# Plot the stacked bar graph
plot <- ggplot(class_counts, aes(x = "", 
                         y = count, 
                         fill = hrc_mmc_class_name)) +
  geom_bar(stat = "identity", 
           position = "stack") +
  scale_fill_manual(values = class_color_palette) +
  coord_polar(theta = "y") +  # Optional: Use this for a pie chart-like appearance
  labs(x = NULL, y = "Count", 
       fill = "Class Label") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        panel.grid = element_blank())   # Remove grid lines

ggsave(filename = "/results/borders_class_distribution.png", 
       plot = plot, 
       width = 8,
       height = 6, 
       dpi = 160)
