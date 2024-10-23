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
  library(randomcoloR)
  library(cowplot)
})

###################### setup environment ########################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

######################### read in anno files ####################

sd.df <- read_sheet("https://docs.google.com/spreadsheets/d/1UpfDYDbq3s0Jeh8cLg9Y0yhPCJxH_H8WG4Y0y0_mNHg/edit?gid=70495470#gid=70495470", sheet = "638850_v3")
# select only necessary columns 
sd.df <- sd.df %>% 
  dplyr::select(Cluster_id,
                spatial_domain_level_1,
                spatial_domain_level_2,
                graph_order,
                ccf_broad)

# load spatial domain color palette
spatial_domain_color <- read_sheet("https://docs.google.com/spreadsheets/d/1UpfDYDbq3s0Jeh8cLg9Y0yhPCJxH_H8WG4Y0y0_mNHg/edit?gid=70495470#gid=70495470", sheet = "638850_v3")
# reorder by graph order
spatial_domain_color <- spatial_domain_color %>%
  arrange(graph_order) 
spatial_domain_palette <- setNames(spatial_domain_color$spatial_domain_level_2_color, spatial_domain_color$spatial_domain_level_2)

################## read in anndata file and filter/enhance metadata ##########

# load anndata file of downsampled data
gridded_ad <- read_h5ad("/data/merfish_638850_mouseadult_spatial_domain_30/mouse_638850_spacel.h5ad")

# ectract umap metadata
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
                      by.x = "leiden_res_1.2_knn_8",
                      by.y = "Cluster_id",
                      all.x = T,
                      all.y = F)

gridded_metadata <- gridded_metadata %>% 
  filter(spatial_domain_level_1 != "LQ")

# reorder spatial domains
gridded_metadata$spatial_domain_level_2 <- factor(gridded_metadata$spatial_domain_level_2, 
                                                  levels = spatial_domain_color$spatial_domain_level_2)

################## plot umap for spatial domain ######################

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
ggsave("/results/legend_sd2.pdf", 
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


###################### plot umaps for sections ###################

# create colormap for sections

section_cluster_colors <- randomColor(length(unique(gridded_metadata$section)), luminosity="dark")
section_cluster_palette <- setNames(section_cluster_colors, unique(gridded_metadata$section))

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

ggsave(filename = "/results/umap.pdf", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)
