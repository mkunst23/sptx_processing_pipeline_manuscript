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
  library(mclust)
})


################# setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()


################################# load anno files ############################

# load landmark annotation per cluster
cl.anat.df <- read_sheet("https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=612105623#gid=612105623", sheet = "cluster_region_annotation")
cl.anat.df <- cl.anat.df %>% 
  select(cluster_id_label,
         broad_region,
         registration_landmark)

# load ccf annotation
ccf_anno <- read_sheet("https://docs.google.com/spreadsheets/d/1QOhsYhlsk2KE2pZuSnEwUhEIG4mh782PcAgHrtJyRYI/edit?gid=0#gid=0",sheet = "CCFv3")
ccf_anno <- ccf_anno %>% 
  select(acronym,
         CCF_registration)

sections_to_keep <- c("1199650929", 
                      "1199650932", 
                      "1199650935", 
                      "1199650938", 
                      "1199650941", 
                      "1199650944", 
                      "1199650950", 
                      "1199650953", 
                      "1199650956", 
                      "1199650959",
                      "1199650962", 
                      "1199650965", 
                      "1199650968", 
                      "1199650972", 
                      "1199650975", 
                      "1199650978",
                      "1199650981", 
                      "1199650984",
                      "1199650999", 
                      "1199651002",
                      "1199651005", 
                      "1199651009", 
                      "1199651012", 
                      "1199651015", 
                      "1199651018",
                      "1199651021",
                      "1199651024", 
                      "1199651027", 
                      "1199651033", 
                      "1199651036",
                      "1199651039", 
                      "1199651042", 
                      "1199651045", 
                      "1199651048", 
                      "1199651054", 
                      "1199651057", 
                      "1199651060", 
                      "1199651063", 
                      "1199651066", 
                      "1199651069",
                      "1199651072", 
                      "1199651075", 
                      "1199651078", 
                      "1199651081", 
                      "1199651084", 
                      "1199651091", 
                      "1199651094",
                      "1199651097", 
                      "1199651100", 
                      "1199651103",
                      "1199651106", 
                      "1199651109", 
                      "1199651112", 
                      "1199651115", 
                      "1199651121",
                      "1199651127", 
                      "1199651130",
                      "1199651133", 
                      "1199651136")

########################### load color palettes ###########################

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
  arrange(desc(graph_order)) 

###################################### analyze new registration ###################################

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

metadata <- merge(metadata,
                  ccf_anno,
                  by.x = "structure_acronym",
                  by.y = "acronym",
                  all.x = T,
                  all.y = F)

metadata_subset <- metadata %>% 
  filter(final_qc_passed == T) %>%
  filter(!is.na(structure_acronym)) %>%
  filter(section %in% sections_to_keep) %>% 
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
         CCF_registration,
         CCF_level1,
         volume_x,
         volume_y,
         volume_z)

############## plot Jaccard overlay with broad landmark clusters ############

# subset data to relevant features
jaccard_df <- metadata_subset %>%
  filter(!is.na(CCF_level1)) %>% 
  filter(!is.na(broad_region)) %>%
  filter(CCF_level1 != "fiber tracts") %>% 
  filter(broad_region != "fiber tracts") %>% 
  select(production_cell_id,
         CCF_level1,
         broad_region)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# calculate the ari
ari <- adjustedRandIndex(jaccard_df$broad_region, jaccard_df$CCF_level1)

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
  scale_size(range = c(0, 5)) + 
  annotate("text", 
           x = Inf, 
           y = Inf, 
           label = paste("ARI:", round(ari, 2)), 
           hjust = 1.1, 
           vjust = 1.1, 
           size = 5, 
           color = "black")
  

ggsave(filename = "/results/jaccard_CCFbroad_CCFlevel1.pdf", 
       plot = plot, 
       width = 7,
       height = 5, 
       dpi = 160)



############## plot Jaccard overlay with fine landmark clusters ############

# subset data to relevant features
jaccard_df <- metadata_subset %>%
  filter(!is.na(CCF_registration)) %>% 
  filter(!is.na(registration_landmark)) %>% 
  select(production_cell_id,
         CCF_registration,
         registration_landmark)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# calculate the ari
ari <- adjustedRandIndex(jaccard_df$registration_landmark, jaccard_df$CCF_registration)

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$registration_landmark
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$CCF_registration
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
                   levels = landmark_color$registration_landmark) 

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
  scale_size(range = c(0, 5))  + 
  annotate("text", 
           x = Inf, 
           y = Inf, 
           label = paste("ARI:", round(ari, 2)), 
           hjust = 1.1, 
           vjust = 1.1, 
           size = 5, 
           color = "black")

ggsave(filename = "/results/jaccard_ccflandmarks.pdf", 
       plot = plot, 
       width = 9,
       height = 7, 
       dpi = 160)


###################################### analyze old registration ###################################


registration_old <- fread('/data/merscope_6538850_old_reg_SIS_seg/ccf_cells_with_structures.csv')

# remove the first two columns from the data frame
registration_old <- registration_old %>% 
  select(-1,-2)

registration_old <- registration_old %>%
  mutate(cell_id = substr(cell_id, 3, nchar(cell_id) - 1))

# add cluster information (temporary)
cluster_info <- metadata %>% 
  select(production_cell_id,
         final_qc_passed,
         hrc_mmc_class_name,
         hrc_mmc_subclass_name,
         hrc_mmc_supertype_name,
         hrc_mmc_cluster_name)

registration_old <- merge(registration_old,
                          cluster_info,
                          by.x = "cell_id",
                          by.y = "production_cell_id",
                          all.x = T,
                          all.y = F)

# convert empty cell to NA
registration_old <- registration_old %>%
  mutate_at(vars(CCF_level1, CCF_level2), ~ na_if(., ""))


# add region information for cell types
registration_old <- merge(registration_old,
                          cl.anat.df,
                          by.x = "hrc_mmc_cluster_name",
                          by.y = "cluster_id_label",
                          all.x = T,
                          all.y = F)


registration_old <- merge(registration_old,
                          ccf_anno,
                          by.x = "structure_acronym",
                          by.y = "acronym",
                          all.x = T,
                          all.y = F)

registration_old_subset <- registration_old %>% 
  filter(final_qc_passed == T) %>%
  filter(!is.na(structure_acronym)) %>%
  select(cell_id,
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
         CCF_registration,
         ccf_x,
         ccf_y,
         ccf_z)

############## plot Jaccard overlay with broad landmark clusters ############

# subset data to relevant features
jaccard_df <- registration_old_subset %>%
  filter(!is.na(CCF_level1)) %>% 
  filter(!is.na(broad_region)) %>%
  filter(CCF_level1 != "fiber tracts") %>% 
  filter(broad_region != "fiber tracts") %>% 
  select(cell_id,
         CCF_level1,
         broad_region)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("cell_id")

# calculate the ari
ari <- adjustedRandIndex(jaccard_df$broad_region, jaccard_df$CCF_level1)

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
  scale_size(range = c(0, 5))  + 
  annotate("text", 
           x = Inf, 
           y = Inf, 
           label = paste("ARI:", round(ari, 2)), 
           hjust = 1.1, 
           vjust = 1.1, 
           size = 5, 
           color = "black")

ggsave(filename = "/results/jaccard_ccf_level1_old.pdf", 
       plot = plot, 
       width = 7,
       height = 5, 
       dpi = 160)

############## plot Jaccard overlay with fine landmark clusters ############

# subset data to relevant features
jaccard_df <- registration_old_subset %>%
  filter(!is.na(CCF_registration)) %>% 
  filter(!is.na(registration_landmark)) %>%
  select(cell_id,
         CCF_registration,
         registration_landmark)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("cell_id")

# calculate the ari
ari <- adjustedRandIndex(jaccard_df$registration_landmark, jaccard_df$CCF_registration)

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$registration_landmark
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$CCF_registration
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
                   levels = landmark_color$registration_landmark) 

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
  scale_size(range = c(0, 5))  + 
  annotate("text", 
           x = Inf, 
           y = Inf, 
           label = paste("ARI:", round(ari, 2)), 
           hjust = 1.1, 
           vjust = 1.1, 
           size = 5, 
           color = "black")

ggsave(filename = "/results/jaccard_ccf_level2_old.pdf", 
       plot = plot, 
       width = 9,
       height = 7, 
       dpi = 160)

