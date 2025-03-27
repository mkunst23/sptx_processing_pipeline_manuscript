suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)   
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(anndata)
  library(googlesheets4)
  library(data.table)
  library(reshape)
  library(forcats)
  library(entropy)
  library(aplot)
  library(vegan)
  library(ggridges)
  library(scales)
  library(DescTools)
  library(RColorBrewer)
  library(philentropy)
  library(igraph)
  library(paletteer)
  library(anndata)
  library(cowplot)
})

###################### setup environment ########################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

######################### read in anno files ####################

sd.df <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domains")
# select only necessary columns 
sd.df <- sd.df %>% 
  dplyr::select(Cluster_id,
                spatial_domain_level_1,
                spatial_domain_level_2,
                graph_order,
                ccf_broad)

# load spatial domain color palette
spatial_domain_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domains")
spatial_domain_palette <- setNames(spatial_domain_color$spatial_domain_level_2_color, spatial_domain_color$spatial_domain_level_2)

# order by graph order
spatial_domain_color <- spatial_domain_color %>% 
  arrange(desc(graph_order))

spatial_domain_1_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domain_level_1_color")
spatial_domain_1_palette <- setNames(spatial_domain_1_color$spatial_domain_level_1_color, spatial_domain_1_color$spatial_domain_level_1)

spatial_domain_1_color <- spatial_domain_1_color %>% 
  arrange(graph_order)

broad_ccf_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "CCF_broad_region_color")
broad_ccf_color_palette <- setNames(broad_ccf_color$color_hex_tripet, broad_ccf_color$ccf_broad)

broad_ccf_color <- broad_ccf_color %>% 
  arrange(desc(graph_order))

# load color palette for section
section_cluster_colors <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=1482302853#gid=1482302853", sheet = "section_metadata")
section_cluster_palette <- setNames(section_cluster_colors$color,section_cluster_colors$barcode)


################## read in anndata file and filter/enhance metadata ##########

# load anndata file of downsampled data
gridded_ad <- read_h5ad("/data/merscope_638850_mouseadult_clustered_sections_QC_filt/whole_dataset/mouse_638850_res_1.4_clustered.h5ad")

# extract umap metadata
umap <- as.data.frame(gridded_ad$obsm["X_umap"])
rownames(umap) <- rownames(gridded_ad$X)

# extract metadata
gridded_metadata <- gridded_ad$obs

# merge data
gridded_metadata <- merge(gridded_metadata,
                          umap,
                          by = 0)

gridded_metadata <- merge(gridded_metadata,
                      sd.df,
                      by.x = "leiden_res_1.4_knn_8",
                      by.y = "Cluster_id",
                      all.x = T,
                      all.y = F)

gridded_metadata <- gridded_metadata %>% 
  filter(spatial_domain_level_1 != "LQ") %>% 
  filter(spatial_domain_level_1 != "Borders")


################# plot umap with no color code ####################

#plot umap for spatial domains with color legend
plot <- ggplot(gridded_metadata,
               aes(x=X_umap.1,
                   y=X_umap.2
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_y_reverse() +
  theme(legend.position = "none") +
  theme_void() +
  theme(strip.text = element_blank()) 

ggsave(filename = "/results/umap.png", 
       plot = plot, 
       width = 10,
       height = 10, 
       dpi = 160)


################## plot umap for spatial domain sd2 ######################

# reorder spatial domains
gridded_metadata$spatial_domain_level_2 <- factor(gridded_metadata$spatial_domain_level_2, 
                                                  levels = spatial_domain_color$spatial_domain_level_2)

#plot umap for spatial domains with color legend
plot <- ggplot(gridded_metadata,
       aes(x=X_umap.1,
           y=X_umap.2,
           color = spatial_domain_level_2
       )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=spatial_domain_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  theme_void() +
  theme(strip.text = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4, 
                                                  shape = 15)))
# Extract the legend
legend <- get_legend(plot)

# Plot only the legend
legend_plot <- plot_grid(legend)

# Save the legend to a PDF
ggsave("/results/umap_legend_sd2.pdf", 
       legend_plot, 
       width = 6, 
       height = 5)

plot <- ggplot(gridded_metadata,
               aes(x=X_umap.1,
                   y=X_umap.2,
                   color = spatial_domain_level_2
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=spatial_domain_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void() +
  theme(strip.text = element_blank()) 


ggsave(filename = "/results/umap_sd_domains.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)


################## plot umap for spatial domain sd1 ######################

gridded_metadata$spatial_domain_level_1 <- factor(gridded_metadata$spatial_domain_level_1, 
                                                  levels = spatial_domain_1_color$spatial_domain_level_1)


#plot umap for spatial domains with color legend
plot <- ggplot(gridded_metadata,
               aes(x=X_umap.1,
                   y=X_umap.2,
                   color = spatial_domain_level_1
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=spatial_domain_1_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  theme_void() +
  theme(strip.text = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4, 
                                                  shape = 15)))
# Extract the legend
legend <- get_legend(plot)

# Plot only the legend
legend_plot <- plot_grid(legend)

# Save the legend to a PDF
ggsave("/results/legend_sd1.pdf", 
       legend_plot, 
       width = 6, 
       height = 5)

plot <- ggplot(gridded_metadata,
               aes(x=X_umap.1,
                   y=X_umap.2,
                   color = spatial_domain_level_1
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=spatial_domain_1_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void() +
  theme(strip.text = element_blank()) 


ggsave(filename = "/results/umap_sd1_domains.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)


###################### plot umaps for ccf_broad ###################

# reorder spatial domains
gridded_metadata$ccf_broad <- factor(gridded_metadata$ccf_broad, 
                                                  levels = broad_ccf_color$ccf_broad)


# plot legend
plot <- ggplot(gridded_metadata,
               aes(x=X_umap.1,
                   y=X_umap.2,
                   color = fct_rev(ccf_broad)
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=broad_ccf_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void() +
  theme(strip.text = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4, 
                                                  shape = 15)))
# Extract the legend
legend <- get_legend(plot)

# Plot only the legend
legend_plot <- plot_grid(legend)

# Save the legend to a PDF
ggsave("/results/legend_ccf_broad.pdf", 
       legend_plot, 
       width = 6, 
       height = 5)

# plot umap without legend
plot <- ggplot(gridded_metadata,
               aes(x=X_umap.1,
                   y=X_umap.2,
                   color = ccf_broad
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=broad_ccf_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/umap_ccf_broad.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)


###################### plot umaps for sections ###################


# plot legend
plot <- ggplot(gridded_metadata,
               aes(x=X_umap.1,
                   y=X_umap.2,
                   color = fct_rev(section)
               )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=section_cluster_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void() +
  theme(strip.text = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4, 
                                                  shape = 15)))
# Extract the legend
legend <- get_legend(plot)

# Plot only the legend
legend_plot <- plot_grid(legend)

# Save the legend to a PDF
ggsave("/results/legend_sections.pdf", 
       legend_plot, 
       width = 6, 
       height = 5)

# plot umap without legend
plot <- ggplot(gridded_metadata,
       aes(x=X_umap.1,
           y=X_umap.2,
           color = section
       )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=section_cluster_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/umap_sections.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)
