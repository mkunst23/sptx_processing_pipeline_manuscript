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


################# setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()


######### load anno files ###########

# load cluster taxonomy
#cl.df <- read_sheet("https://docs.google.com/spreadsheets/d/1T86EppILF-Q97OjJyd4R2FjUmPUrdYBUEbbiXqgWEoE/edit#gid=1639101745",sheet = "cl.df.v9_230722")
# select only necessary columns 
#cl.df <- cl.df %>% 
#  dplyr::select(cluster_id_label,
#                nt_type_label,
#                nt_type_combo_label)

# load landmark annotation per cluster
cl.anat.df <- read_sheet("https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=612105623#gid=612105623", sheet = "cluster_region_annotation")
cl.anat.df <- cl.anat.df %>% 
  select(cluster_id_label,
         broad_region,
         registration_landmark)

###### load color palettes

broad_ccf_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "broad_landmarks_anno")
broad_ccf_color_palette <- setNames(broad_ccf_color$broad_region_color, broad_ccf_color$broad_region)
# order by graph order
broad_ccf_color <- broad_ccf_color %>% 
  arrange(desc(graph_order)) 

# load ccf color
ccf_color <- read_sheet("https://docs.google.com/spreadsheets/d/1QOhsYhlsk2KE2pZuSnEwUhEIG4mh782PcAgHrtJyRYI/edit?gid=0#gid=0", sheet = "CCFv3")
ccf_color <- ccf_color %>% 
  select(acronym,
         color_hex_triplet,
         graph_order) %>% 
  mutate(color_hex_triplet = paste0("#", color_hex_triplet))

ccf_color_palette <- setNames(ccf_color$color_hex_triplet, ccf_color$acronym)

ccf_color <- ccf_color %>% 
  arrange(desc(graph_order))

# generate color palette for ccf landmark color
landmark_color <- read_sheet("https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=2094634713#gid=2094634713", sheet = "CCF_landmark_color")
landmark_color_palette <- setNames(landmark_color$registration_landmark_color, landmark_color$registration_landmark)

# order by graph order
landmark_color <- landmark_color %>% 
  arrange(graph_order) 


# read in metadata file
metadata <- fread("/data/merscope_638850_mouseadult_registered_v2/whole_dataset/mouse_638850_registered.csv")

metadata <- metadata %>%
  mutate_at(vars(CCF_level1, CCF_level2), ~ na_if(., ""))


# add region information for cell types
metadata <- merge(metadata,
                  cl.anat.df,
                  by.x = "hrc_mmc_cluster_name",
                  by.y = "cluster_id_label",
                  all.x = T,
                  all.y = F)


metadata_subset <- metadata %>% 
  filter(final_qc_passed == T) %>%
  filter(!is.na(structure_acronym)) %>%
  select(production_cell_id,
         section,
         hrc_mmc_class_name,
         hrc_mmc_subclass_name,
         hrc_mmc_supertype_name,
         hrc_mmc_cluster_name,
         structure_acronym,
         structure_id,
         broad_region,
         registration_landmark,
         CCF_level1,
         CCF_level2,
         volume_x,
         volume_y,
         volume_z)

############## plot Jaccard overlay with broad landmark clusters ############

# subset data to relevant features
jaccard_df <- metadata_subset %>%
  filter(!is.na(CCF_level1)) %>% 
  filter(!is.na(broad_region)) %>%
  select(production_cell_id,
         CCF_level1,
         broad_region)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$broad_region
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$CCF_level1
names(ref.cl) <- rownames(jaccard_df)
ref.cl <- as.factor(ref.cl)

# create a table of the number of cells in each cluster
tb <- table(cl, ref.cl)

# set minimum threshold for number of cells in a cluster
min.th = 1

# remove rows with less than min.th cells
tb.df <- as.data.frame(tb)
tb.df <- tb.df[tb.df$Freq >= min.th,]

# create a table of the number of cells in each cluster
cl.size <- table(cl)
ref.cl.size <- table(ref.cl)

# Compute Jaccard statistics for each pair of clusters
tb.df$jaccard <- as.vector(tb.df$Freq / (cl.size[as.character(tb.df[,1])] + ref.cl.size[as.character(tb.df[,2])] - tb.df$Freq))


# reorder ref.cl and cl
tb.df$cl <- factor(tb.df$cl,
                   levels = ccf_color$acronym) 

tb.df$ref.cl <- factor(tb.df$ref.cl,
                       levels = broad_ccf_color$broad_region) 

tb.df <- tb.df %>% 
  drop_na()

# make a dot plot where the size of the dot is proportional to tb.df$Freq and the color is proportional to tb.df$jaccard
plot <- ggplot(tb.df, 
               aes(x = fct_rev(cl), 
                   y = ref.cl)) + 
  geom_point(aes(size = sqrt(Freq),
                 color = jaccard)) + 
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size = 8),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust = 1, 
                                 size = 8),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        panel.grid.major = element_line(color = "grey80"), # Major grid lines
        panel.grid.minor = element_line(color = "grey90"),  # Minor grid lines
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Spatial domaind level 1", y = "Region specific clusters") +
  scale_color_gradient(low = "yellow", 
                       high = "darkblue") + 
  scale_size(range = c(0, 5)) 

ggsave(filename = "/results/jaccard_ccf_level1.pdf", 
       plot = plot, 
       width = 7,
       height = 5, 
       dpi = 160)



############## plot Jaccard overlay with fine landmark clusters ############

# subset data to relevant features
jaccard_df <- metadata_subset %>%
  filter(!is.na(CCF_level2)) %>% 
  filter(!is.na(registration_landmark)) %>%
  select(production_cell_id,
         CCF_level2,
         registration_landmark)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$registration_landmark
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$CCF_level2
names(ref.cl) <- rownames(jaccard_df)
ref.cl <- as.factor(ref.cl)

# create a table of the number of cells in each cluster
tb <- table(cl, ref.cl)

# set minimum threshold for number of cells in a cluster
min.th = 1

# remove rows with less than min.th cells
tb.df <- as.data.frame(tb)
tb.df <- tb.df[tb.df$Freq >= min.th,]

# create a table of the number of cells in each cluster
cl.size <- table(cl)
ref.cl.size <- table(ref.cl)

# Compute Jaccard statistics for each pair of clusters
tb.df$jaccard <- as.vector(tb.df$Freq / (cl.size[as.character(tb.df[,1])] + ref.cl.size[as.character(tb.df[,2])] - tb.df$Freq))


# reorder ref.cl and cl
tb.df$cl <- factor(tb.df$cl,
                   levels = ccf_color$acronym) 

tb.df$ref.cl <- factor(tb.df$ref.cl,
                       levels = landmark_color$registration_landmark) 

tb.df <- tb.df %>% 
  drop_na()


# make a dot plot where the size of the dot is proportional to tb.df$Freq and the color is proportional to tb.df$jaccard
plot <- ggplot(tb.df, 
               aes(x = fct_rev(cl), 
                   y = ref.cl)) + 
  geom_point(aes(size = sqrt(Freq),
                 color = jaccard)) + 
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size = 8),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust = 1, 
                                 size = 8),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        panel.grid.major = element_line(color = "grey80"), # Major grid lines
        panel.grid.minor = element_line(color = "grey90"),  # Minor grid lines
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Spatial domaind level 1", y = "Region specific clusters") +
  scale_color_gradient(low = "yellow", 
                       high = "darkblue") + 
  scale_size(range = c(0, 5)) 

ggsave(filename = "/results/jaccard_ccf_level2.pdf", 
       plot = plot, 
       width = 7,
       height = 5, 
       dpi = 160)





