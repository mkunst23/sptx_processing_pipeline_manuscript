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
})

############### setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

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

# create colormap for leiden clustering
leiden_cluster_colors <- randomColor(length(unique(gridded_metadata$leiden_res_1.2_knn_8)), luminosity="dark")
leiden_cluster_palette <- setNames(leiden_cluster_colors, unique(gridded_metadata$leiden_res_1.2_knn_8))


# plot umap
ggplot(gridded_metadata,
       aes(x=X_umap.1,
           y=X_umap.2,
           color = leiden_res_1.2_knn_8
       )) +
  geom_point(size=.1,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=leiden_cluster_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  theme_void() +
  theme(strip.text = element_blank())

# create colormap for leiden clustering
section_cluster_colors <- randomColor(length(unique(gridded_metadata$section)), luminosity="dark")
section_cluster_palette <- setNames(section_cluster_colors, unique(gridded_metadata$section))

ggplot(gridded_metadata,
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
