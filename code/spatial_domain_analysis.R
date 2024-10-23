#################### spatial_domain_analysis ######################

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
  library(randomcoloR)
})

############### setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

proportion_colors <- c('#005DFFFF',
                       "#0080FFFF",
                       "#00A2FFFF", 
                       '#00C3FFFF', 
                       "#00E5FFFF",
                       "#00FF4DFF",
                       "#00FF2BFF",
                       "#00FF08FF",
                       '#1AFF00FF', 
                       "#3CFF00FF", 
                       "#5DFF00FF", 
                       "#80FF00FF", 
                       "#A2FF00FF",
                       "#C3FF00FF",
                       "#E6FF00FF", 
                       "#FFFF00FF",
                       "#FFF514FF", 
                       "#FFEC28FF", 
                       "#FFE53CFF", 
                       "#FFE04FFF", 
                       "#FFDC63FF",
                       "#FFD100FF",
                       "#FFC500FF",
                       "#FFB900FF",
                       "#FFAE00FF",
                       "#FFA200FF",
                       "#FF9700FF",
                       "#FF8B00FF",
                       "#FF8000FF",
                       "#FF7400FF",
                       "#FF6800FF",
                       "#FF5D00FF",
                       "#FF5100FF",
                       "#FF4600FF",
                       "#FF3A00FF",
                       "#FF2E00FF",
                       "#FF2300FF",
                       "#FF1700FF")

cl.df <- read_sheet("https://docs.google.com/spreadsheets/d/1T86EppILF-Q97OjJyd4R2FjUmPUrdYBUEbbiXqgWEoE/edit#gid=1639101745",sheet = "cl.df.v9_230722")
# select only necessary columns 
cl.df <- cl.df %>% 
  dplyr::select(cluster_id_label,
                nt_type_label,
                nt_type_combo_label)


sd.df <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domains")
# select only necessary columns 
sd.df <- sd.df %>% 
  dplyr::select(Cluster_id,
                spatial_domain_level_1,
                spatial_domain_level_2,
                graph_order,
                ccf_broad)

broad_ccf_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "CCF_broad_region_color")
broad_ccf_color_palette <- setNames(broad_ccf_color$color_hex_tripet, broad_ccf_color$ccf_broad)

spatial_domain_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "spatial_domains")
spatial_domain_palette <- setNames(spatial_domain_color$spatial_domain_level_2_color, spatial_domain_color$spatial_domain_level_2)

# read in metadata file
load("/scratch/638850_metadata_sis.rda")

metadata_sis <- metadata_sis %>%
  mutate(temp = flat_CDM_class_name) %>%
  separate(temp, 
           into = c("flat_CDM_class_id", "flat_CDM_class_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata_sis <- metadata_sis %>%
  mutate(temp = flat_CDM_subclass_name) %>%
  separate(temp, 
           into = c("flat_CDM_subclass_id", "flat_CDM_subclass_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata_sis <- metadata_sis %>%
  mutate(temp = flat_CDM_supertype_name) %>%
  separate(temp, 
           into = c("flat_CDM_supertype_id", "flat_CDM_supertype_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata_sis <- metadata_sis %>%
  mutate(temp = flat_CDM_cluster_name) %>%
  separate(temp, 
           into = c("flat_CDM_cluster_id", "flat_CDM_cluster_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata_sis <- merge(metadata_sis,
                      cl.df,
                      by.x = "flat_CDM_cluster_name",
                      by.y = "cluster_id_label",
                      all.x = T,
                      all.y = F)


metadata_sis <- merge(metadata_sis,
                      sd.df,
                      by.x = "leiden_res_1.2_knn_8",
                      by.y = "Cluster_id",
                      all.x = T,
                      all.y = F)


metadata_subset <- metadata_sis %>% 
  filter(final_filter == F) %>%
  filter(!is.na(spatial_domain_level_1)) %>% 
  filter(spatial_domain_level_1 != "LQ") %>% 
  select(cell_id,
         section,
         leiden_res_1.2_knn_8,
         flat_CDM_class_name,
         flat_CDM_class_id,
         flat_CDM_subclass_name,
         flat_CDM_subclass_id,
         flat_CDM_supertype_name,
         flat_CDM_supertype_id,
         flat_CDM_cluster_label,
         flat_CDM_cluster_id,
         nt_type_label,
         nt_type_combo_label,
         structure_acronym,
         structure_id,
         spatial_domain_level_1,
         spatial_domain_level_2,
         graph_order,
         ccf_broad,
         CCF_level1,
         CCF_level2,
         volume_x,
         volume_y,
         volume_z,
         mapped)

# grab colors for the different levels of cell types. Thsoe will be used later for plotting
ss <- "https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?usp=sharing"
class_colors <- read_sheet(ss,sheet="classes")
subclass_colors <- read_sheet(ss, sheet="subclasses")
supertype_colors <- read_sheet(ss, sheet="supertypes")
rcol <- read_sheet(ss, sheet = "Region")
nt_type_col <- read_sheet(ss, sheet = "nt_type")


############ create helper data frames for plotting ####################
region.mapping <- metadata_subset %>% 
  select(ccf_broad,spatial_domain_level_2,graph_order) %>% 
  group_by(spatial_domain_level_2) %>% 
  slice(1) %>% 
  arrange(desc(graph_order)) %>% 
  drop_na()

add.meta <- metadata_subset %>% 
  select(flat_CDM_subclass_name, flat_CDM_subclass_id, flat_CDM_class_name) %>% 
  group_by_all() %>% 
  unique %>% 
  arrange(flat_CDM_class_name, flat_CDM_subclass_name)

example_sections <- c("1199650953","1199651021","1199651042","1199651084")

plot_data <- metadata_subset %>% 
  filter(section %in% example_sections)

plot <- ggplot(plot_data,
               aes(x=volume_x,
                   y=volume_y,
                   color = spatial_domain_level_2
               )) +
  geom_point(size=.5,
             stroke=0,
             shape=19,) +
  coord_fixed() +
  scale_color_manual(values=spatial_domain_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_example_sections.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)

####################### create base heatmap #####################

dominance_scores <- metadata_subset %>%
  count(spatial_domain_level_2, flat_CDM_subclass_name) %>%
  group_by(spatial_domain_level_2) %>%
  mutate(Dominance = n/max(n))

dominance_scores$spatial_domain_level_2 <- factor(dominance_scores$spatial_domain_level_2, 
                                      levels = region.mapping$spatial_domain_level_2)

subclass_dominance <- ggplot(dominance_scores, aes(x = flat_CDM_subclass_name, 
                                                   y = spatial_domain_level_2, 
                                                   fill = Dominance)) +
  geom_tile() +
  scale_fill_gradientn(colors = proportion_colors) +
  theme(panel.background = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        panel.border = element_rect(colour = "black", fill = NA)
  )

# plot subclass tiles
cols <- setNames(subclass_colors$subclass_color, 
                 subclass_colors$subclass_id_label)

l2plot <- ggplot(data=add.meta) + 
  geom_tile(aes(x=flat_CDM_subclass_name, 
                y=1, 
                fill=flat_CDM_subclass_name)) +
  scale_fill_manual(values=cols, 
                    name = "Subclass", 
                    limits = force) +
  theme_void() +
  ylab( "Subclass") +
  theme(axis.title.y = element_text(size=8, 
                                    hjust=1, 
                                    vjust=0.5),
        legend.position = "none")

# plot class tiles
cols <- setNames(class_colors$class_color, 
                 class_colors$class_id_label)

clplot <- ggplot(data=add.meta) + 
  geom_tile(aes(x=flat_CDM_subclass_name, 
                y=1, 
                fill=flat_CDM_class_name)) +
  scale_fill_manual(values=cols, 
                    name = "Class", 
                    limits = force) +
  theme_void() +
  ylab( "Class") +
  theme(axis.title.y = element_text(size=8, 
                                    hjust=1, 
                                    vjust=0.5),
        #axis.text.x = element_text(angle = 45, 
        #                           hjust = 1),
        legend.position = "right")


sd2_plot <- ggplot(data=region.mapping) + 
  geom_tile(aes(y=spatial_domain_level_2, 
                x=1, 
                fill=spatial_domain_level_2)) +
  scale_fill_manual(values=spatial_domain_palette, 
                    name = "spatial domain", 
                    limits = force) +
  #theme_void() +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

#plot CCF color tiles
ccfplot <- ggplot(data=region.mapping) + 
  geom_tile(aes(y=spatial_domain_level_2, 
                x=1, 
                fill=ccf_broad)) +
  scale_fill_manual(values=broad_ccf_color_palette, 
                    name = "CCF", 
                    limits = force) +
  #theme_void() +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")


final_subclass_dominance <- subclass_dominance %>%
  insert_bottom(l2plot, height = 0.03) %>%
  insert_bottom(clplot, height = 0.03)  %>% 
  insert_left(sd2_plot, width=0.01) %>% 
  insert_left(ccfplot, width=0.01)

ggsave(filename = "/results/sd_heatmap.pdf", 
       plot = final_subclass_dominance, 
       width = 14,
       height = 6, 
       dpi = 160)

