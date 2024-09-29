########## region analysis plotting #############

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
})

############### setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

############# setup helper functions ##################

# function to calculate gini coefficient for each row of a data frame
calcGini <- function(row) {
  Gini(row, unbiased = FALSE)
}

################ setup annotation files ###############

# define grey and white matter
grey_matter <- c("CB",
                 "CTXsp",
                 "HIP",
                 "HY",
                 "Isocortex",
                 "LSX",
                 "MB",
                 "MY",
                 "OLF",
                 "P",
                 "PAL",
                 "RHP",
                 "sAMY",
                 "STR",
                 "TH")
white_matter <- c("fiber tracts",
                  "VS")

# define glial subclasses
glia_sc <- c("316 Bergmann NN",
             "317 Astro-CB NN",
             "318 Astro-NT NN",
             "319 Astro-TE NN",
             "320 Astro-OLF NN",
             "321 Astroependymal NN",
             "326 OPC NN",
             "327 Oligo NN",
             "328 OEC NN")

# define neuronal classes
neuronal_classes <- c("01 IT-ET Glut",
                      "02 NP-CT-L6b Glut",
                      "03 OB-CR Glut",
                      "04 DG-IMN Glut",
                      "05 OB-IMN GABA",
                      "06 CTX-CGE GABA",
                      "07 CTX-MGE GABA",
                      "08 CNU-MGE GABA",
                      "09 CNU-LGE GABA",
                      "10 LSX GABA",
                      "11 CNU-HYa GABA",
                      "12 HY GABA",
                      "13 CNU-HYa Glut",
                      "14 HY Glut",
                      "15 HY Gnrh1 Glut",
                      "16 HY MM Glut",
                      "17 MH-LH Glut",
                      "18 TH Glut",
                      "19 MB Glut",
                      "20 MB GABA",
                      "21 MB Dopa",
                      "22 MB-HB Sero",
                      "23 P Glut",
                      "24 MY Glut",
                      "26 P GABA",
                      "27 MY GABA",
                      "28 CB GABA",
                      "29 CB Glut")

cl.df <- read_sheet("https://docs.google.com/spreadsheets/d/1T86EppILF-Q97OjJyd4R2FjUmPUrdYBUEbbiXqgWEoE/edit#gid=1639101745",sheet = "cl.df.v9_230722")
# select only necessary columns (labels plus max region and ratio)
cl.df <- cl.df %>% 
  dplyr::select(cluster_id_label,
                nt_type_label,
                nt_type_combo_label)

# create CCF region colors
ccf_anno <- read_sheet("https://docs.google.com/spreadsheets/d/1QOhsYhlsk2KE2pZuSnEwUhEIG4mh782PcAgHrtJyRYI/edit#gid=0",
                       sheet="CCFv3")
# create color scheme
CCF_regions_color <- ccf_anno %>% 
  select(acronym,color_hex_triplet) 
CCF_regions_color$color_hex_triplet <- gsub("^", "#", CCF_regions_color$color_hex_triplet)
CCF_regions_color <- setNames(CCF_regions_color$color_hex_triplet, CCF_regions_color$acronym)

ccf_anno <- ccf_anno %>% 
  select(acronym,
         graph_order)

# grab colors for the different levels of cell types. Thsoe will be used later for plotting
ss <- "https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?usp=sharing"
class_colors <- read_sheet(ss,sheet="classes")
subclass_colors <- read_sheet(ss, sheet="subclasses")
supertype_colors <- read_sheet(ss, sheet="supertypes")
rcol <- read_sheet(ss, sheet = "Region")
nt_type_col <- read_sheet(ss, sheet = "nt_type")

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


############## process metadata file ####################

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

metadata_subset <- metadata_sis %>% 
  filter(final_filter == F) %>% 
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
         CCF_level1,
         CCF_level2,
         volume_x,
         volume_y,
         volume_z,
         mapped)

metadata_subset <- merge(metadata_subset,
                         ccf_anno,
                         by.x = "structure_acronym",
                         by.y = "acronym",
                         all.x = T,
                         all.y = F)

# assign neuron/glia metadata_subset to cells
metadata_subset <- metadata_subset %>% 
  mutate(N_G_ratio = case_when(flat_CDM_class_name %in% neuronal_classes ~ "neuronal",
                               flat_CDM_subclass_name %in% glia_sc ~ "glial",
                               TRUE ~ NA_character_))

############################### filter metadata ################################
metadata_subset <- metadata_subset %>% 
  filter(CCF_level1 %in% grey_matter)

# remove cells not assigned to one of the major regions
metadata_subset <- metadata_subset %>% 
  filter(CCF_level2 != "NA")

############ create helper data frames for plotting ####################
region.mapping <- metadata_subset %>% 
  select(CCF_level1,CCF_level2,graph_order) %>% 
  group_by(CCF_level2) %>% 
  slice(1) %>% 
  arrange(desc(graph_order)) %>% 
  drop_na()

region.mapping_broad <- metadata_subset %>% 
  select(CCF_level1,graph_order) %>% 
  group_by(CCF_level1) %>% 
  slice(1) %>% 
  arrange(desc(graph_order)) %>% 
  drop_na()

add.meta <- metadata_subset %>% 
  select(flat_CDM_subclass_name, flat_CDM_subclass_id, flat_CDM_class_name) %>% 
  group_by_all() %>% 
  unique %>% 
  arrange(flat_CDM_class_name, flat_CDM_subclass_name)

add.meta_st <- metadata_subset %>% 
  select(flat_CDM_supertype_name, flat_CDM_supertype_id, flat_CDM_class_name) %>% 
  group_by_all() %>% 
  unique %>% 
  arrange(flat_CDM_class_name, flat_CDM_supertype_name)

############ calculate cell types per region at different levels #############
# calculates subclass per region
subclass_per_region <- dplyr::count(metadata_subset, flat_CDM_subclass_name, CCF_level2)
subclass_per_region <- spread(subclass_per_region, key = CCF_level2, value = n)
subclass_per_region <- subclass_per_region %>% 
  replace(is.na(.), 0)

subclass_per_region_perc <- subclass_per_region %>%
  mutate(across(-flat_CDM_subclass_name, ~ ./rowSums(subclass_per_region[, -1])))

subclass_per_region_norm <- subclass_per_region %>% 
  tibble::column_to_rownames(var = "flat_CDM_subclass_name")

subclass_per_region_norm <- t(t(subclass_per_region_norm) / colSums(subclass_per_region_norm))

subclass_per_region_norm_max <- subclass_per_region %>% 
  tibble::column_to_rownames(var = "flat_CDM_subclass_name")

subclass_per_region_norm_max <- t(t(subclass_per_region_norm_max) / apply(subclass_per_region_norm_max, 1 ,max))

subclass_per_region_norm_max <- as.data.frame(subclass_per_region_norm_max)

subclass_per_region_norm_max <- subclass_per_region_norm_max %>% 
  tibble::rownames_to_column(var = "flat_CDM_subclass_name")

subclass_per_region_norm_max_long <- melt(subclass_per_region_norm_max, 
                                          id.vars = "flat_CDM_subclass_name")

# calculates supertype per region
supertype_per_region <- dplyr::count(metadata_subset, flat_CDM_supertype_name, CCF_level2)
supertype_per_region <- spread(supertype_per_region, key = CCF_level2, value = n)
supertype_per_region <- supertype_per_region %>% 
  replace(is.na(.), 0)

supertype_per_region_perc <- supertype_per_region %>%
  mutate(across(-flat_CDM_supertype_name, ~ ./rowSums(supertype_per_region[, -1])))


supertype_per_region_norm <- supertype_per_region %>% 
  tibble::column_to_rownames(var = "flat_CDM_supertype_name")

supertype_per_region_norm <- t(t(supertype_per_region_norm) / colSums(supertype_per_region_norm))

####################### plot nt type distribution ##############################
metadata_nt <- metadata_subset %>%
  filter(nt_type_label != "NA") %>% 
  filter(flat_CDM_class_name %in% neuronal_classes) %>%
  select(CCF_level1,
         CCF_level2,
         nt_type_label)

metadata_nt$nt_type_label <- factor(metadata_nt$nt_type_label, 
                                    levels = c("Glut",
                                               "GABA",
                                               "Glut-GABA",
                                               "GABA-Glyc",
                                               "Dopa",
                                               "Sero",
                                               "Chol",
                                               "Nora",
                                               "Hist"))

metadata_nt$CCF_level1 <- factor(metadata_nt$CCF_level1, 
                                 levels = region.mapping_broad$CCF_level1)

cols <- setNames(nt_type_col$nt_type_color, 
                 nt_type_col$nt_type_label)

E_I_level_1 <- ggplot(data = metadata_nt, aes(x = fct_rev(CCF_level1), fill = nt_type_label)) +
  geom_bar(position = "fill") +
  labs(x = "Broad CCF region", y = "%", title = "Transmitter ratio per region") +
  scale_fill_manual(values = cols) +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

cols <- setNames(nt_type_col$nt_type_color, 
                 nt_type_col$nt_type_label)

metadata_nt$CCF_level2 <- factor(metadata_nt$CCF_level2, 
                                 levels = region.mapping$CCF_level2)

E_I_level_2 <- ggplot(data = na.omit(metadata_nt), 
                      aes(x = CCF_level2, 
                          fill = nt_type_label)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols) +
  ylab( "NT type distribution") +
  theme_void() +
  theme(legend.position = "right") +
  coord_flip()

###################### plot neuron to glia ratio ##########################
metadata_ng <- metadata_subset %>% 
  drop_na(N_G_ratio) %>%
  drop_na(CCF_level2) %>% 
  select(CCF_level1,
         CCF_level2,
         N_G_ratio)

metadata_ng$N_G_ratio <- factor(metadata_ng$N_G_ratio, 
                                levels = c("neuronal",
                                           "glial"))

metadata_ng$CCF_level1 <- factor(metadata_ng$CCF_level1, 
                                 levels = region.mapping_broad$CCF_level1)

color_ng <- c("#00203FFF",
              "#ADEFD1FF")

names(color_ng) <- c("neuronal",
                     "glial")

N_G <- ggplot(data = metadata_ng, 
              aes(x = fct_rev(CCF_level1), 
                  fill = N_G_ratio)) +
  geom_bar(position = "fill") +
  labs(x = "Broad CCF region",
       y = "%", 
       title = "Neuron to Glia ratio") +
  scale_fill_manual(values = color_ng) +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,
                               hjust = 1)
  )

metadata_ng$CCF_level2 <- factor(metadata_ng$CCF_level2, 
                                 levels = region.mapping$CCF_level2)

N_G_level_2 <- ggplot(data = metadata_ng, 
                      aes(x = CCF_level2, 
                          fill = N_G_ratio)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = color_ng) +
  ylab( "Neuron to Glia ratio") +
  theme_void() +
  theme(legend.position = "right") +
  coord_flip()

# calculate number of cells per region
cell_count_per_region <- colSums(subclass_per_region[,-1])
cell_count_per_region <- as.data.frame(cell_count_per_region)
cell_count_per_region <- cell_count_per_region %>% 
  tibble::rownames_to_column("CCF_level2")

cell_count_per_region$CCF_level2 <- factor(cell_count_per_region$CCF_level2, 
                                           levels = region.mapping$CCF_level2)
region_cell_count <- ggplot(cell_count_per_region, 
                            aes(x = CCF_level2, 
                                y = log10(cell_count_per_region))) +
  geom_bar(stat = "identity", 
           fill = "Grey") +
  ylab( "cells per region") +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.ticks.y=element_blank(),
    axis.text.y = element_blank(),
    axis.title.y=element_blank(),
    legend.position = "none"
  )

####################### create base heatmap #####################

dominance_scores <- metadata_subset %>%
  count(CCF_level2, flat_CDM_subclass_name) %>%
  group_by(CCF_level2) %>%
  mutate(Dominance = n/max(n))

dominance_scores$CCF_level2 <- factor(dominance_scores$CCF_level2, 
                                      levels = region.mapping$CCF_level2)

dominance_scores$subclass_id_label <- factor(dominance_scores$flat_CDM_subclass_name, 
                                             levels = add.meta$flat_CDM_subclass_name)

subclass_dominance <- ggplot(dominance_scores, aes(x = flat_CDM_subclass_name, 
                                                   y = CCF_level2, 
                                                   fill = Dominance)) +
  geom_tile() +
  scale_fill_gradientn(colors = proportion_colors) +
  theme(panel.background = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        panel.border = element_rect(colour = "black", fill = NA)
  )

# change data fram from wide to long for plotting of heatmap
for_heatmap <- melt(subclass_per_region_perc)
# merge heatmap information (on more refined regions)
for_heatmap <- for_heatmap %>% 
  left_join(region.mapping, 
            by=c("variable"="CCF_level2"))

for_heatmap <- for_heatmap %>% 
  drop_na()

colnames(for_heatmap) <- c("flat_CDM_subclass_name",
                           "CCF_level2",
                           "perc_subclass", 
                           "CCF_level1",
                           "graph_order")

# for broad labels
for_heatmap$CCF_broad <- factor(for_heatmap$CCF_level1, 
                                levels = region.mapping_broad$CCF_level1)
# for minor regions
for_heatmap$CCF_level2 <- factor(for_heatmap$CCF_level2, 
                                 levels = region.mapping$CCF_level2)

# retain order of supertype_id_label 
hm_sc <- data.frame(flat_CDM_subclass_name = unique(for_heatmap$flat_CDM_subclass_name))
hm_sc <- hm_sc %>% 
  left_join(add.meta) %>% 
  drop_na()

hm_sc$flat_CDM_subclass_name <- factor(hm_sc$flat_CDM_subclass_name, 
                                  levels = add.meta$flat_CDM_subclass_name)

# subclass tiles
cols <- setNames(subclass_colors$subclass_color, 
                 subclass_colors$subclass_id_label)

l2plot <- ggplot(data=hm_sc) + 
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

cols <- setNames(class_colors$class_color, 
                 class_colors$class_id_label)

clplot <- ggplot(data=hm_sc) + 
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
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = "right")

