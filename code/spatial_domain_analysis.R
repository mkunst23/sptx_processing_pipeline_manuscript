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
  library(stringr)
})

############### setup environment ######################

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores=30)

gs4_deauth()

# setup proportion colors for dominance scores
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

# read in cluster information
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

# load landmark annotation per cluster
cl.anat.df <- read_sheet("https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=612105623#gid=612105623", sheet = "cluster_region_annotation")
cl.anat.df <- cl.anat.df %>% 
  select(cluster_id_label,
         broad_region,
         registration_landmark)


# ccf color palette
ccf_color <- read_sheet("https://docs.google.com/spreadsheets/d/1QOhsYhlsk2KE2pZuSnEwUhEIG4mh782PcAgHrtJyRYI/edit?gid=0#gid=0", sheet = "CCFv3")
ccf_color <- ccf_color %>% 
  select(acronym,
         color_hex_triplet,
         graph_order) %>% 
  mutate(color_hex_triplet = paste0("#",color_hex_triplet))

ccf_color_palette <- setNames(ccf_color$color_hex_triplet, ccf_color$acronym)
# order by graph order
ccf_color <- ccf_color %>% 
  arrange(desc(graph_order)) 


broad_ccf_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "CCF_broad_region_color")
broad_ccf_color_palette <- setNames(broad_ccf_color$color_hex_tripet, broad_ccf_color$ccf_broad)
# order by graph order
broad_ccf_color <- broad_ccf_color %>% 
  arrange(desc(graph_order)) 

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

# generate color palette for ccf landmark color
landmark_color <- read_sheet("https://docs.google.com/spreadsheets/d/1v7PDfLc_9vOcuz_WKk5ZiSUKjbu3xCux9VanVOA7QOE/edit?gid=2094634713#gid=2094634713", sheet = "CCF_landmark_color")
landmark_color_palette <- setNames(landmark_color$registration_landmark_color, landmark_color$registration_landmark)

# order by graph order
landmark_color <- landmark_color %>% 
  arrange(desc(graph_order)) 

broad_landmarks_color <- read_sheet("https://docs.google.com/spreadsheets/d/1SdmQooCJtqq__n0D7INn12yTmIHwsJ6KeO5mGhR2OBU/edit?gid=0#gid=0", sheet = "broad_landmarks_anno")
broad_landmarks_color_palette <- setNames(broad_landmarks_color$broad_region_color, broad_landmarks_color$broad_region)

# order by graph order
broad_landmarks_color <- broad_landmarks_color %>% 
  arrange(desc(graph_order))


# load in metadata file
metadata <- fread("/data/merscope_638850_mouseadult_registered_v2/whole_dataset/mouse_638850_registered.csv")
metadata$section <- as.character(metadata$section)

# filter out only relevant sections
metadata <- metadata %>% 
  select(production_cell_id,
         section,
         CCF_level1,
         CCF_level2,
         final_qc_passed,
         hrc_mmc_cluster_alias,
         hrc_mmc_cluster_name,
         hrc_mmc_supertype_name,
         hrc_mmc_subclass_name,
         hrc_mmc_class_name,
         structure_id,
         structure_acronym,
         volume_x,
         volume_y,
         volume_z)

sd_domains <- fread('/scratch/mouse_638850_sd.csv')

metadata <- merge(metadata,
                  sd_domains,
                  by = 'production_cell_id',
                  all.x = T,
                  all.y = F)

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

# add additional cell type information such as nt-type
metadata <- merge(metadata,
                  cl.df,
                  by.x = "hrc_mmc_cluster_name",
                  by.y = "cluster_id_label",
                  all.x = T,
                  all.y = F)

# add region information for cell types
metadata <- merge(metadata,
                  cl.anat.df,
                  by.x = "hrc_mmc_cluster_name",
                  by.y = "cluster_id_label",
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
  filter(!is.na(spatial_domain_level_1)) %>% 
  filter(spatial_domain_level_1 != "LQ") %>% 
  filter(spatial_domain_level_1 != "Borders")

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
  select(hrc_mmc_subclass_name, hrc_mmc_subclass_id, hrc_mmc_class_name) %>% 
  group_by_all() %>% 
  unique %>% 
  arrange(hrc_mmc_class_name, hrc_mmc_subclass_name)


##################### plot example sections #####################
# add more sections
# only plot hemisections

example_sections <- c("1199650941",
                      "1199650950",
                      "1199650965",
                      "1199650975",
                      "1199650999",
                      "1199651012",
                      "1199651024",
                      "1199651039",
                      "1199651057",
                      "1199651072",
                      "1199651084",
                      "1199651103")

plot_data <- metadata_subset %>% 
  filter(section %in% example_sections) %>% 
  filter(volume_x < 5.6)


plot <- ggplot(plot_data,
               aes(x=volume_x,
                   y=volume_y,
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
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd_example_sections.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)


plot <- ggplot(plot_data,
               aes(x=volume_x,
                   y=volume_y,
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
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/sd1_example_sections.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)


plot <- ggplot(plot_data,
               aes(x=volume_x,
                   y=volume_y,
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
  facet_wrap(~fct_rev(section), nrow = 2) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/ccf_broad_example_sections.png", 
       plot = plot, 
       width = 6,
       height = 6, 
       dpi = 160)

##################### create heatmap for subclass dominace #####################

dominance_scores <- metadata_subset %>%
  dplyr::count(spatial_domain_level_2, hrc_mmc_subclass_name) %>%
  group_by(spatial_domain_level_2) %>%
  mutate(Dominance = n/max(n))

dominance_scores$spatial_domain_level_2 <- factor(dominance_scores$spatial_domain_level_2, 
                                      levels = region.mapping$spatial_domain_level_2)

subclass_dominance <- ggplot(dominance_scores, aes(x = hrc_mmc_subclass_name, 
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

# plot class tiles
cols <- setNames(class_colors$class_color, 
                 class_colors$class_id_label)

clplot <- ggplot(data=add.meta) + 
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
        #axis.text.x = element_text(angle = 45, 
        #                           hjust = 1),
        legend.position = "right")


sd2_plot <- ggplot(data=region.mapping) + 
  geom_tile(aes(y=spatial_domain_level_2, 
                x=1, 
                fill=spatial_domain_level_2)) +
  scale_fill_manual(values=spatial_domain_palette, 
                    name = "Spatial domain", 
                    limits = force) +
  theme_void() +
  xlab( "Spatial Domain") +
  theme(axis.title.x = element_text(size=8, 
                                    hjust=1, 
                                    vjust=0.5),
        #axis.text.x = element_text(angle = 45, 
        #                           hjust = 1),
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

############## plot Jaccard overlay with broad landmark clusters ############

# subset data to relevant features
jaccard_df <- metadata %>%
  filter(final_qc_passed == T) %>%
  filter(!is.na(spatial_domain_level_1)) %>% 
  filter(spatial_domain_level_1 != "LQ") %>%
  filter(spatial_domain_level_1 != "Borders") %>%
  filter(!is.na(broad_region)) %>%
  select(production_cell_id,
         spatial_domain_level_1,
         spatial_domain_level_2,
         broad_region,
         registration_landmark)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$broad_region
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$spatial_domain_level_1
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
                   levels = broad_landmarks_color$broad_region) 

tb.df$ref.cl <- factor(tb.df$ref.cl,
                  levels = spatial_domain_1_color$spatial_domain_level_1) 

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

ggsave(filename = "/results/jaccard.pdf", 
       plot = plot, 
       width = 7,
       height = 5, 
       dpi = 160)


############################ plot landmarks ########################


# subset data to relevant features
jaccard_df <- metadata %>%
  filter(final_qc_passed == T) %>%
  filter(!is.na(spatial_domain_level_2)) %>% 
  filter(!str_detect(spatial_domain_level_2, "LQ")) %>%
  filter(!str_detect(spatial_domain_level_2, "Borders")) %>%
  filter(!is.na(registration_landmark)) %>%
  select(production_cell_id,
         spatial_domain_level_2,
         registration_landmark)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$registration_landmark
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$spatial_domain_level_2
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
                       levels = spatial_domain_color$spatial_domain_level_2) 

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

ggsave(filename = "/results/jaccard_sd2.pdf", 
       plot = plot, 
       width = 7,
       height = 5, 
       dpi = 160)


################ compare to registration CCF level1 #########################

# subset data to relevant features
jaccard_df <- metadata %>%
  filter(final_qc_passed == T) %>%
  filter(!is.na(spatial_domain_level_1)) %>% 
  filter(!str_detect(spatial_domain_level_1, "LQ")) %>%
  filter(!str_detect(spatial_domain_level_2, "Borders")) %>%
  filter(!is.na(CCF_level1)) %>%
  select(production_cell_id,
         spatial_domain_level_1,
         CCF_level1)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$CCF_level1
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$spatial_domain_level_1
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
                       levels = spatial_domain_1_color$spatial_domain_level_1) 

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


################ compare to registration CCF level2 #########################

# subset data to relevant features
jaccard_df <- metadata %>%
  filter(final_qc_passed == T) %>%
  filter(!is.na(spatial_domain_level_2)) %>% 
  filter(!str_detect(spatial_domain_level_2, "LQ")) %>%
  filter(!is.na(CCF_level2)) %>%
  select(production_cell_id,
         spatial_domain_level_2,
         CCF_level2)

# convert cell label to rownames
jaccard_df <- jaccard_df %>% 
  tibble::column_to_rownames("production_cell_id")

# convert spatial domain level 1 and parcellation division to factors
cl <- jaccard_df$CCF_level2
names(cl) <- rownames(jaccard_df)
cl <- as.factor(cl)

ref.cl <- jaccard_df$spatial_domain_level_2
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
                       levels = spatial_domain_color$spatial_domain_level_2) 

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
                                   size = 6),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust = 1, 
                                 size = 6),  #remove y axis labels
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
       width = 26,
       height = 10, 
       dpi = 160)
