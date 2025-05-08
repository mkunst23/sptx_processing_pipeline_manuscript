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

folder_path <- "/results/638850"

# Check if the folder exists
if (!dir.exists(folder_path)) {
  # Create the folder if it does not exist
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created:", folder_path, "\n")
} else {
  cat("Folder already exists:", folder_path, "\n")
}

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
                 "STRv",
                 "STRd",
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
# select only necessary columns 
rnaseq_anno <- cl.df %>% 
  select(subclass_id_label,
         class_id_label)

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
metadata <- fread("/data/mouse_638850_registration//whole_dataset/mouse_638850_registered.csv")

metadata <- metadata %>%
  mutate_at(vars(CCF_level1, CCF_level2), ~ na_if(., ""))

metadata <- metadata %>%
  mutate(temp = hrc_mmc_class_name) %>%
  separate(temp, 
           into = c("hrc_mmc_class_id", "hrc_mmc_class_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata <- metadata %>%
  mutate(temp = hrc_mmc_subclass_name) %>%
  separate(temp, 
           into = c("hrc_mmc_subclass_id", "hrc_mmc_subclass_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata <- metadata %>%
  mutate(temp = hrc_mmc_supertype_name) %>%
  separate(temp, 
           into = c("hrc_mmc_supertype_id", "hrc_mmc_supertype_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata <- metadata %>%
  mutate(temp = hrc_mmc_cluster_name) %>%
  separate(temp, 
           into = c("hrc_mmc_cluster_id", "hrc_mmc_cluster_string"), 
           sep = " ", 
           extra = "merge", 
           fill = "right")

metadata <- merge(metadata,
                      cl.df,
                      by.x = "hrc_mmc_cluster_name",
                      by.y = "cluster_id_label",
                      all.x = T,
                      all.y = F)

metadata_subset <- metadata %>% 
  filter(final_qc_passed == T) %>% 
  select(production_cell_id,
         section,
         hrc_mmc_class_name,
         hrc_mmc_class_id,
         hrc_mmc_subclass_name,
         hrc_mmc_subclass_id,
         hrc_mmc_supertype_name,
         hrc_mmc_supertype_id,
         hrc_mmc_cluster_label,
         hrc_mmc_cluster_id,
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
  mutate(N_G_ratio = case_when(hrc_mmc_class_name %in% neuronal_classes ~ "neuronal",
                               hrc_mmc_subclass_name %in% glia_sc ~ "glial",
                               TRUE ~ NA_character_))

############################### filter metadata ################################
metadata_subset <- metadata_subset %>% 
  filter(CCF_level1 %in% grey_matter)

# remove cells not assigned to one of the major regions
metadata_subset <- metadata_subset %>% 
  filter(CCF_level1 != "NA") %>% 
  filter(CCF_level1 != "" )

################### code to clean out bad alignment cells #################
# calculate number of cells per broad region
subclass_per_region_level1 <- dplyr::count(metadata_subset, hrc_mmc_subclass_name, CCF_level1)
subclass_per_region_level1 <- spread(subclass_per_region_level1, key = CCF_level1, value = n)
subclass_per_region_level1 <- subclass_per_region_level1 %>% 
  replace(is.na(.), 0)

# calculate percentage of cells per subclass per broad region
subclass_per_region_level1_perc <- subclass_per_region_level1 %>%
  mutate_at(vars(-hrc_mmc_subclass_name), list(~ ./rowSums(subclass_per_region_level1[, -1])))

# convert subclass_id_label to row names
subclass_per_region_level1_perc <- subclass_per_region_level1_perc %>% 
  tibble::column_to_rownames("hrc_mmc_subclass_name")

# filter out cells if their percentage in a region is less than 5%
indices <- which(subclass_per_region_level1_perc > 0 & subclass_per_region_level1_perc < 0.05, arr.ind = TRUE)
row_names <- rownames(subclass_per_region_level1_perc)[indices[,1]]
col_names <- colnames(subclass_per_region_level1_perc)[indices[,2]]

# make data frame with subclass name and region that have less than 5%
bad_alignment <- cbind(row_names,
                       col_names)
bad_alignment <- as.data.frame(bad_alignment)
colnames(bad_alignment) <- c("hrc_mmc_subclass_name","CCF_level1")

# extend to cluster_id_label
# merge data frames
bad_alignment <- merge(bad_alignment,
                       rnaseq_anno,
                       by.x = "hrc_mmc_subclass_name",
                       by.y = "subclass_id_label",
                       all = T)

# only keep neuronal classes
bad_alignment <- bad_alignment %>% 
  filter(class_id_label %in% neuronal_classes) %>% 
  select(hrc_mmc_subclass_name,
         CCF_level1) %>% 
  unique()

# extract cells with bad alignment
bad_cells <- metadata_subset %>% 
  inner_join(bad_alignment, by = c("hrc_mmc_subclass_name", "CCF_level1")) %>% 
  pull(production_cell_id)

metadata_subset <- metadata_subset %>% 
  filter(!production_cell_id %in% bad_cells)

# save cleaned data
fwrite(metadata_subset,
       "/scratch/638850_metadata_cleaned.csv",
       row.names = F)

############ create helper data frames for plotting ####################
region.mapping <- metadata_subset %>% 
  select(CCF_level1,CCF_level2,graph_order) %>% 
  group_by(CCF_level2) %>% 
  slice(1) %>% 
  arrange(desc(graph_order)) %>% 
  drop_na() %>% 
  filter(CCF_level2 != "")

region.mapping_broad <- metadata_subset %>% 
  select(CCF_level1,graph_order) %>% 
  group_by(CCF_level1) %>% 
  slice(1) %>% 
  arrange(desc(graph_order)) %>% 
  drop_na()

add.meta <- metadata_subset %>% 
  select(hrc_mmc_subclass_name, hrc_mmc_subclass_id, hrc_mmc_class_name) %>% 
  group_by_all() %>% 
  unique %>% 
  arrange(hrc_mmc_class_name, hrc_mmc_subclass_name)

add.meta_st <- metadata_subset %>% 
  select(hrc_mmc_supertype_name, hrc_mmc_supertype_id, hrc_mmc_class_name) %>% 
  group_by_all() %>% 
  unique %>% 
  arrange(hrc_mmc_class_name, hrc_mmc_supertype_name)

############ calculate cell types per region at different levels #############
# calculates subclass per region
subclass_per_region <- dplyr::count(metadata_subset, hrc_mmc_subclass_name, CCF_level2)
subclass_per_region <- spread(subclass_per_region, key = CCF_level2, value = n)
subclass_per_region <- subclass_per_region %>% 
  replace(is.na(.), 0)

subclass_per_region$`<NA>` <- NULL

subclass_per_region_perc <- subclass_per_region %>%
  mutate(across(-hrc_mmc_subclass_name, ~ ./rowSums(subclass_per_region[, -1])))

subclass_per_region_norm <- subclass_per_region %>% 
  tibble::column_to_rownames(var = "hrc_mmc_subclass_name")

subclass_per_region_norm <- t(t(subclass_per_region_norm) / colSums(subclass_per_region_norm))

subclass_per_region_norm_max <- subclass_per_region %>% 
  tibble::column_to_rownames(var = "hrc_mmc_subclass_name")

subclass_per_region_norm_max <- t(t(subclass_per_region_norm_max) / apply(subclass_per_region_norm_max, 1 ,max))

subclass_per_region_norm_max <- as.data.frame(subclass_per_region_norm_max)

subclass_per_region_norm_max <- subclass_per_region_norm_max %>% 
  tibble::rownames_to_column(var = "hrc_mmc_subclass_name")

subclass_per_region_norm_max_long <- melt(subclass_per_region_norm_max, 
                                          id.vars = "hrc_mmc_subclass_name")

# calculates supertype per region
supertype_per_region <- dplyr::count(metadata_subset, hrc_mmc_supertype_name, CCF_level2)
supertype_per_region <- spread(supertype_per_region, key = CCF_level2, value = n)
supertype_per_region <- supertype_per_region %>% 
  replace(is.na(.), 0)

supertype_per_region$`<NA>` <- NULL

supertype_per_region_perc <- supertype_per_region %>%
  mutate(across(-hrc_mmc_supertype_name, ~ ./rowSums(supertype_per_region[, -1])))


supertype_per_region_norm <- supertype_per_region %>% 
  tibble::column_to_rownames(var = "hrc_mmc_supertype_name")

supertype_per_region_norm <- t(t(supertype_per_region_norm) / colSums(supertype_per_region_norm))

####################### plot nt type distribution ##############################
metadata_nt <- metadata_subset %>%
  filter(nt_type_label != "NA") %>% 
  filter(hrc_mmc_class_name %in% neuronal_classes) %>%
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

ggsave(filename = "/results/EI_ratio.pdf", 
       plot = E_I_level_1, 
       width = 10,
       height = 5, 
       dpi = 300)

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

ggsave(filename = "/results/NG_ratio.pdf", 
       plot = N_G, 
       width = 10,
       height = 5, 
       dpi = 300)

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
  dplyr::count(CCF_level2, hrc_mmc_subclass_name) %>%
  group_by(CCF_level2) %>%
  mutate(Dominance = n/max(n)) %>% 
  drop_na()

dominance_scores$CCF_level2 <- factor(dominance_scores$CCF_level2, 
                                      levels = region.mapping$CCF_level2)

dominance_scores$subclass_id_label <- factor(dominance_scores$hrc_mmc_subclass_name, 
                                             levels = add.meta$hrc_mmc_subclass_name)

subclass_dominance <- ggplot(dominance_scores, aes(x = hrc_mmc_subclass_name, 
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

colnames(for_heatmap) <- c("hrc_mmc_subclass_name",
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
hm_sc <- data.frame(hrc_mmc_subclass_name = unique(for_heatmap$hrc_mmc_subclass_name))
hm_sc <- hm_sc %>% 
  left_join(add.meta) %>% 
  drop_na()

hm_sc$hrc_mmc_subclass_name <- factor(hm_sc$hrc_mmc_subclass_name, 
                                  levels = add.meta$hrc_mmc_subclass_name)

# subclass tiles
cols <- setNames(subclass_colors$subclass_color, 
                 subclass_colors$subclass_id_label)

l2plot <- ggplot(data=hm_sc) + 
  geom_tile(aes(x=hrc_mmc_subclass_name, 
                y=1, 
                fill=hrc_mmc_subclass_name)) +
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
  geom_tile(aes(x=hrc_mmc_subclass_name, 
                y=1, 
                fill=hrc_mmc_class_name)) +
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

cols <- setNames(class_colors$class_color, 
                 class_colors$class_id_label)

clplot <- ggplot(data=hm_sc) + 
  geom_tile(aes(x=hrc_mmc_subclass_name, y=1, fill=hrc_mmc_class_name)) +
  scale_fill_manual(values=cols, name = "Class", limits = force) +
  theme_void() +
  ylab( "Class") +
  theme(axis.title.y = element_text(size=8, 
                                    hjust=1, 
                                    vjust=0.5),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = "right")

# reorder region mapping
region.mapping$CCF_level2 <- factor(region.mapping$CCF_level2, 
                                    levels = fct_rev(region.mapping$CCF_level2))

# code to plot bar graphs for brain regions
ccfplot <- ggplot(data=region.mapping) + 
  geom_tile(aes(y=fct_rev(CCF_level2), 
                x=1, fill=CCF_level1)) +
  scale_fill_manual(values=CCF_regions_color, 
                    name = "CCF", 
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

ylabs <- unique(for_heatmap[,c("CCF_level1",
                               "CCF_level2")])
ylabs$y <- 1:nrow(ylabs)

ylabs.2 = ylabs %>%
  group_by(grp = rleid(CCF_level1), 
           CCF_level1) %>%    
  summarise(y = round((min(y)+ max(y))/2 )) %>% 
  left_join(select(ylabs, y, CCF_level2))

# Broad regions
lab2 <- ggplot(data=ylabs.2, aes(y=CCF_level2, 
                                 x=-50, 
                                 label = CCF_level1)) +   
  geom_text(hjust=1.0,
            vjust=0.5,
            size=3) +
  theme_void()+
  coord_cartesian( clip = "off") 

###################### calculate gini coefficient for subclass ########################
gini_coefficient <- apply(subclass_per_region[,-1], 1, calcGini)
names(gini_coefficient) <- subclass_per_region$hrc_mmc_subclass_name
gini_coefficient <- as.data.frame(gini_coefficient)

cell_count <- rowSums(subclass_per_region[,-1])
names(cell_count) <- subclass_per_region$hrc_mmc_subclass_name
cell_count <- as.data.frame(cell_count)

subclass_meta <- merge(gini_coefficient,
                       cell_count,
                       by=0)

colnames(subclass_meta)[1] <- "hrc_mmc_subclass_name"

subclass_meta <- merge(subclass_meta,
                       add.meta,
                       by = "hrc_mmc_subclass_name",
                       all.x = T,
                       all.y = F)

subclass_meta$name <- "Name"

subclass_meta$hrc_mmc_subclass_name <- factor(subclass_meta$hrc_mmc_subclass_name, 
                                          levels = add.meta$hrc_mmc_subclass_name)

subclass_cell_count <- ggplot(subclass_meta, 
                              aes(x = hrc_mmc_subclass_name, 
                                  y = log10(cell_count))) +
  geom_bar(stat = "identity", fill = "Grey") +
  ylab( "cells per subclass") +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    #axis.ticks.y=element_blank(),
    #axis.text.y = element_blank(),
    axis.title.y=element_blank(),
    legend.position = "none"
  )

# calculate normalized gini coefficient
gini_coefficient_norm <- apply(subclass_per_region_norm, 1, calcGini)
gini_coefficient_norm <- as.data.frame(gini_coefficient_norm)
colnames(gini_coefficient_norm)[1] <- "gini_coefficient_norm"
gini_coefficient_norm <- gini_coefficient_norm %>% 
  tibble::rownames_to_column(var = "hrc_mmc_subclass_name")


subclass_meta <- merge(subclass_meta,
                       gini_coefficient_norm,
                       by = "hrc_mmc_subclass_name",
                       all.x = T,
                       all.y = T)

subclass_meta$hrc_mmc_subclass_name <- factor(subclass_meta$hrc_mmc_subclass_name, 
                                          levels = add.meta$hrc_mmc_subclass_name)

subclass_meta_long <- tidyr::gather(subclass_meta, 
                                    key = "hrc_mmc_subclass_name", 
                                    value = "value", 
                                    -gini_coefficient_norm, 
                                    -cell_count)

###################### calculate gini coefficient for supertypes ########################
gini_coefficient_st <- apply(supertype_per_region[,-1], 1, calcGini)
names(gini_coefficient_st) <- supertype_per_region$hrc_mmc_supertype_name
gini_coefficient_st <- as.data.frame(gini_coefficient_st)

cell_count_st <- rowSums(supertype_per_region[,-1])
names(cell_count_st) <- supertype_per_region$hrc_mmc_supertype_name
cell_count_st <- as.data.frame(cell_count_st)

supertype_meta <- merge(gini_coefficient_st,
                        cell_count_st,
                        by=0)

colnames(supertype_meta)[1] <- "hrc_mmc_supertype_name"

supertype_meta <- merge(supertype_meta,
                        add.meta_st,
                        by = "hrc_mmc_supertype_name",
                        all.x = T,
                        all.y = F)


supertype_meta$name <- "Name"

# calculate normalized gini coefficient
gini_coefficient_norm_st <- apply(supertype_per_region_norm, 1, calcGini)
gini_coefficient_norm_st <- as.data.frame(gini_coefficient_norm_st)
colnames(gini_coefficient_norm_st)[1] <- "gini_coefficient_norm"
gini_coefficient_norm_st <- gini_coefficient_norm_st %>% 
  tibble::rownames_to_column(var = "hrc_mmc_supertype_name")


supertype_meta <- merge(supertype_meta,
                        gini_coefficient_norm_st,
                        by = "hrc_mmc_supertype_name",
                        all.x = T,
                        all.y = T)

supertype_meta$hrc_mmc_supertype_name <- factor(supertype_meta$hrc_mmc_supertype_name, 
                                                 levels = add.meta_st$hrc_mmc_supertype_name)

################# plot distribution for a gini coefficient ###################
gini_colors <- paletteer_c("grDevices::terrain.colors", 30)

gini_distribution <- ggplot(subclass_meta, 
                            aes(x = gini_coefficient_norm, 
                                y = name, 
                                fill = after_stat(x))) +
  geom_density_ridges_gradient() +
  theme_minimal() +
  scale_fill_gradientn(colors = gini_colors)

ggsave(filename = "/results/gini_sc.pdf", 
       plot = gini_distribution, 
       width = 10,
       height = 5, 
       dpi = 300)

gini_distribution_st <- ggplot(supertype_meta, 
                               aes(x = gini_coefficient_norm, 
                                   y = name, 
                                   fill = after_stat(x))) +
  geom_density_ridges_gradient() +
  theme_minimal() +
  scale_fill_gradientn(colors = gini_colors)


ggsave(filename = "/results/gini_st.pdf", 
       plot = gini_distribution, 
       width = 10,
       height = 5, 
       dpi = 300)

gini_plot <- ggplot(data=subclass_meta_long) + 
  geom_tile(aes(x=value, 
                y=hrc_mmc_subclass_name, 
                fill=gini_coefficient_norm)) +
  scale_fill_gradientn(colors = gini_colors) +
  theme_void() +
  ylab( "Gini coefficient") +
  theme(axis.title.y = element_text(size=8, 
                                    hjust=1, 
                                    vjust=0.5),
        legend.position = "right")

cols <- setNames(class_colors$class_color, 
                 class_colors$class_id_label)

gini_distribution_class <- ggplot(subclass_meta, 
                                  aes(x = gini_coefficient_norm, 
                                      y = hrc_mmc_class_name, 
                                      fill = hrc_mmc_class_name)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_manual(values=cols) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = "/results/gini_sc_by_class.pdf", 
       plot = gini_distribution_class, 
       width = 10,
       height = 5, 
       dpi = 300)

cols <- setNames(class_colors$class_color, 
                 class_colors$class_id_label)

st_gini_distribution_class <- ggplot(supertype_meta, 
                                     aes(x = gini_coefficient_norm, 
                                         y = hrc_mmc_class_name, 
                                         fill = hrc_mmc_class_name)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_manual(values=cols) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = "/results/gini_st_by_class.pdf", 
       plot = gini_distribution_class, 
       width = 10,
       height = 5, 
       dpi = 300)

############# calculate shannon diversity for subclasses ###################
subclass_abundance <- subclass_per_region[, 2:ncol(subclass_per_region)]

shann_d_values <- subclass_per_region %>%
  column_to_rownames(var = "hrc_mmc_subclass_name") %>%
  summarise(across(everything(),
                   vegan::diversity,
                   index = "shannon")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "CCF_level2")

colnames(shann_d_values)[2] <- "sd_index"

shann_d_values$normalized_shannon <- shann_d_values$sd_index / log(ncol(subclass_abundance))

shann_d_values <- merge(shann_d_values,
                        region.mapping,
                        by = "CCF_level2",
                        all.x = T,
                        all.y = F)

shann_d_values$CCF_level2 <- factor(shann_d_values$CCF_level2, 
                                    levels = fct_rev(region.mapping$CCF_level2))

shann_d_values$name <- "Name"

############# calculate shannon diversity for supertypes ###################
supertype_abundance <- supertype_per_region[, 2:ncol(subclass_per_region)]

shann_d_values_st <- supertype_per_region %>%
  column_to_rownames(var = "hrc_mmc_supertype_name") %>%
  summarise(across(everything(),
                   vegan::diversity,
                   index = "shannon")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "CCF_level2")

colnames(shann_d_values_st)[2] <- "sd_index"

shann_d_values_st$normalized_shannon <- shann_d_values_st$sd_index / log(ncol(supertype_abundance))

shann_d_values_st <- merge(shann_d_values_st,
                           region.mapping,
                           by = "CCF_level2",
                           all.x = T,
                           all.y = F)

shann_d_values_st$CCF_level2 <- factor(shann_d_values_st$CCF_level2, 
                                       levels = fct_rev(region.mapping$CCF_level2))

shann_d_values_st$name <- "Name"

#################### plot shannon diverity ##########################

shannon_colors <- paletteer_c("ggthemes::Orange-Blue Diverging", 30, direction = -1)

shann_d_distribution <- ggplot(shann_d_values, 
                               aes(x = sd_index, 
                                   y = name, 
                                   fill = after_stat(x))) +
  geom_density_ridges_gradient() +
  theme_minimal() +
  scale_fill_gradientn(colors = shannon_colors)

ggsave(filename = "/results/shannon_dist.pdf", 
       plot = shann_d_distribution, 
       width = 10,
       height = 5, 
       dpi = 300)

shann_d_values$CCF_level1 <- factor(shann_d_values$CCF_level1, 
                                    levels = region.mapping_broad$CCF_level1)

# split up by CCF_level1 (change order)
shann_d_distribution_ccf <- ggplot(shann_d_values, 
                                   aes(x = sd_index, 
                                       y = CCF_level1, 
                                       fill = CCF_level1)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_manual(values = CCF_regions_color) +
  xlim(1, 5) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = "/results/shannon_dist_by_region.pdf", 
       plot = shann_d_distribution_ccf, 
       width = 10,
       height = 5, 
       dpi = 300)

shann_d_values_st$CCF_level1 <- factor(shann_d_values_st$CCF_level1, 
                                       levels = region.mapping_broad$CCF_level1)

# split up by CCF_level1 (change order)
shann_d_distribution_ccf_st <- ggplot(shann_d_values_st, 
                                      aes(x = sd_index, 
                                          y = CCF_level1, 
                                          fill = CCF_level1)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_manual(values = CCF_regions_color) +
  xlim(1, 5) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = "/results/shannon_dist_by_region_st.pdf", 
       plot = shann_d_distribution_ccf_st, 
       width = 10,
       height = 5, 
       dpi = 300)


shannon_plot <- ggplot(data=shann_d_values) +
  geom_tile(aes(y=fct_rev(CCF_level2), 
                x=1, 
                fill=sd_index)) +
  scale_fill_gradientn(colors = shannon_colors) +
  theme_void() +
  theme(legend.position = "right")

shannon_plot_st <- ggplot(data=shann_d_values_st) +
  geom_tile(aes(y=fct_rev(CCF_level2), 
                x=1, 
                fill=sd_index)) +
  scale_fill_gradientn(colors = shannon_colors) +
  theme_void() +
  theme(legend.position = "right")

################## make final plot ###################

final_subclass_dominance <- subclass_dominance %>%
  insert_bottom(l2plot, height = 0.03) %>%
  insert_bottom(clplot, height = 0.03) %>%
  insert_left(E_I_level_2, width=0.06) %>% 
  insert_left(N_G_level_2, width=0.04) %>% 
  insert_left(ccfplot, width=0.01) %>%
  #insert_left(lab2, width = 0.2) %>% 
  insert_top(gini_plot,height = 0.03) %>% 
  insert_top(subclass_cell_count,height = 0.03) %>% 
  insert_right(shannon_plot, width=0.02) %>%
  insert_right(shannon_plot_st, width=0.02) %>%
  insert_right(region_cell_count, width=0.06)

ggsave(filename = "/results/638850_heatmap.pdf", 
       plot = final_subclass_dominance, 
       width = 10,
       height = 7, 
       dpi = 300)


