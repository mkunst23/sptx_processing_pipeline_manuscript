######### script to illustrate mapping strategy ################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(cowplot)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(data.table)
  library(rearrr)
  library(reticulate)
  library(anndata)
  library(ttr)
  library(devtools)
  library(RColorBrewer)
  library(googlesheets4)
  library(plyr)
  library(forcats)
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()

ss <- "https://docs.google.com/spreadsheets/d/1a4URa_t3oc824vjIA8LcyYmJrv0-2Rm7Q5eJq42re2c/edit?usp=sharing"
supertype_colors <- read_sheet(ss, sheet="supertypes")
supertype_color_palette <- setNames(supertype_colors$supertype_color_new, supertype_colors$supertype_id_label)

cluster_colors <- read_sheet(ss, sheet="clusters")
cluster_color_palette <- setNames(cluster_colors$cluster_color, cluster_colors$cluster_label)

# read in metadata file
metadata <- fread("/data/merscope_638850_mouseadult_registered_v2/whole_dataset/mouse_638850_registered.csv")


filter_palette <- c("grey","red")

example_sections = c("1199650984",
                     "1199651018",
                     "1199651094",
                     "1199651109")

example_supertypes <- c("0682 RN Spp1 Glut_1",
                        "0023 L4/5 IT CTX Glut_1",
                        "1180 OPC NN_2",
                        "1076 NLL Gata3 Gly-Gaba_3",
                        "0125 L5 NP CTX Glut_4",
                        "0979 PGRN-PARN-MDRN Hoxb5 Glut_7",
                        "0048 IT AON-TT-DP Glut_3",
                        "1137 SPVI-SPVC Sall3 Nfib Gly-Gaba_2",
                        "0151 OB-in Frmd7 Gaba_2",
                        "0130 OB Eomes Ms4a15 Glut_1",
                        "0396 SI-MPO-LPO Lhx8 Gaba_1",
                        "0295 NDB-SI-ant Prdm12 Gaba_1",
                        "0056 MEA Slc17a7 Glut_2")

# filter out relevant columns
filtered_metadata <- metadata %>%
  filter(basic_qc_passed == T) %>% 
  filter(doublets_qc_passed == T) %>%
  filter(hrc_mmc_supertype_name %in% example_supertypes) %>% 
  select(production_cell_id,
         volume_x,
         volume_y,
         section,
         final_qc_passed,
         hrc_mmc_class_name,
         hrc_mmc_class_avg_correlation,
         hrc_mmc_subclass_name,
         hrc_mmc_subclass_avg_correlation,
         hrc_mmc_supertype_name,
         hrc_mmc_supertype_avg_correlation,
         hrc_mmc_supertype_thr_criteria,
         hrc_mmc_cluster_name,
         hrc_mmc_cluster_avg_correlation,
         hrc_mmc_cluster_thr,
         hrc_mmc_cluster_thr_criteria,
         hrc_mmc_thr)


# extract section threshold
section_threshold <- filtered_metadata %>% 
  select(hrc_mmc_supertype_name,
         hrc_mmc_thr) %>% 
  unique() %>% 
  arrange(hrc_mmc_supertype_name)


# plot distribution of correlation coefficients for selected 
plot <- ggplot(filtered_metadata,
               aes(x = hrc_mmc_supertype_name,
                   y = hrc_mmc_cluster_avg_correlation,
                   fill = hrc_mmc_supertype_name)) +
  geom_violin() +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  scale_fill_manual(values = supertype_color_palette) +
  geom_hline(yintercept = 0.5, 
             linetype = "dashed", 
             color = "black") +
  labs(x = "Supertype Name",
       y = "# of detected spots per cell") +
  guides(color = "none") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = "none")


ggsave(filename = "/results/mapping_filter_old.pdf", 
       plot = plot, 
       width = 10,
       height = 7, 
       dpi = 300)



plot <- ggplot(filtered_metadata,
               aes(x = hrc_mmc_supertype_name,
                   y = hrc_mmc_cluster_avg_correlation,
                   fill = hrc_mmc_supertype_name)) +
  geom_violin() +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  scale_fill_manual(values = supertype_color_palette) +
  geom_hline(yintercept = 0.5, 
             linetype = "dashed", 
             color = "black") +
  labs(x = "Supertype Name",
       y = "# of detected spots per cell") +
  guides(color = "none") +
  geom_segment(aes(x = 0.75,
                   xend = 1.25,
                   y = section_threshold$hrc_mmc_thr[1],
                   yend = section_threshold$hrc_mmc_thr[1],
                   color = "red")) +
  geom_segment(aes(x = 1.75,
                   xend = 2.25,
                   y = section_threshold$hrc_mmc_thr[2],
                   yend = section_threshold$hrc_mmc_thr[2],
                   color = "red")) +
  geom_segment(aes(x = 2.75,
                   xend = 3.25,
                   y = section_threshold$hrc_mmc_thr[3],
                   yend = section_threshold$hrc_mmc_thr[3],
                   color = "red")) +
  geom_segment(aes(x = 3.75,
                   xend = 4.25,
                   y = section_threshold$hrc_mmc_thr[4],
                   yend = section_threshold$hrc_mmc_thr[4],
                   color = "red")) +
  geom_segment(aes(x = 4.75,
                   xend = 5.25,
                   y = section_threshold$hrc_mmc_thr[5],
                   yend = section_threshold$hrc_mmc_thr[5],
                   color = "red")) +
  geom_segment(aes(x = 5.75,
                   xend = 6.25,
                   y = section_threshold$hrc_mmc_thr[6],
                   yend = section_threshold$hrc_mmc_thr[6],
                   color = "red")) +
  geom_segment(aes(x = 6.75,
                   xend = 7.25,
                   y = section_threshold$hrc_mmc_thr[7],
                   yend = section_threshold$hrc_mmc_thr[7],
                   color = "red")) +
  geom_segment(aes(x = 7.75,
                   xend = 8.5,
                   y = section_threshold$hrc_mmc_thr[8],
                   yend = section_threshold$hrc_mmc_thr[8],
                   color = "red")) +
  geom_segment(aes(x =8.75,
                   xend = 9.25,
                   y = section_threshold$hrc_mmc_thr[9],
                   yend = section_threshold$hrc_mmc_thr[9],
                   color = "red")) +
  geom_segment(aes(x =9.75,
                   xend = 10.25,
                   y = section_threshold$hrc_mmc_thr[10],
                   yend = section_threshold$hrc_mmc_thr[10],
                   color = "red")) +
  geom_segment(aes(x =10.75,
                   xend = 11.25,
                   y = section_threshold$hrc_mmc_thr[11],
                   yend = section_threshold$hrc_mmc_thr[11],
                   color = "red")) +
  geom_segment(aes(x =11.75,
                   xend = 12.25,
                   y = section_threshold$hrc_mmc_thr[12],
                   yend = section_threshold$hrc_mmc_thr[12],
                   color = "red")) +
  geom_segment(aes(x =12.75,
                   xend = 13.25,
                   y = section_threshold$hrc_mmc_thr[13],
                   yend = section_threshold$hrc_mmc_thr[13],
                   color = "red")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = "none")


ggsave(filename = "/results/mapping_filter_new.pdf", 
       plot = plot, 
       width = 10,
       height = 7, 
       dpi = 300)


####### plot example sections with old and new filtering threshold #############

# example for initual threshold to high
# supertype: 0682 RN Spp1 Glut_1
# sections: 1199651021 1199651018 1199651015 1199650962

cluster_to_plot <- "0682 RN Spp1 Glut_1"
example_sections <- c("1199651021","1199651018","1199651015","1199650965")

classifier_palette <- c("lightgrey","black","red")

threshold <- 0.5

example_sections_plot <- metadata %>%
  filter(basic_qc_filter == F) %>% 
  filter(doublets_filter == F) %>%
  filter(section %in% example_sections) %>% 
  select(x_coordinate,
         y_coordinate,
         section,
         flat_CDM_supertype_name,
         flat_CDM_cluster_avg_correlation,
         flat_CDM_cluster_thr)

# Create the classifier column
example_sections_plot <- example_sections_plot %>%
  mutate(classifier = case_when(
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation > threshold ~ "hq",
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation <= threshold ~ "lq",
    flat_CDM_supertype_name != cluster_to_plot ~ "bg"
  ))


# plot failed and passed in black and red and the rest in grey

plot <- ggplot(example_sections_plot,
       aes(x=x_coordinate,
           y=y_coordinate,
           color = classifier,
           alpha = classifier,
           size = classifier
       )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  scale_color_manual(values = classifier_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.1, 1, 1)) +
  scale_size_manual(values = c(0.1, 2, 2)) +
  guides(color = "none") +
  guides(alpha = "none") +
  guides(size = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/mapping_0682_RN_Spp1_Glut_1_fixed_filter.png", 
       plot = plot, 
       width = 5,
       height = 15, 
       dpi = 300)


threshold <- section_threshold[9,2]

example_sections_plot <- metadata_sis %>%
  filter(basic_qc_filter == F) %>% 
  filter(doublets_filter == F) %>%
  filter(section %in% example_sections) %>% 
  select(x_coordinate,
         y_coordinate,
         section,
         flat_CDM_supertype_name,
         flat_CDM_cluster_avg_correlation,
         flat_CDM_cluster_thr)

# Create the classifier column
example_sections_plot <- example_sections_plot %>%
  mutate(classifier = case_when(
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation > threshold ~ "hq",
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation <= threshold ~ "lq",
    flat_CDM_supertype_name != cluster_to_plot ~ "bg"
  ))


plot <- ggplot(example_sections_plot,
               aes(x=x_coordinate,
                   y=y_coordinate,
                   color = classifier,
                   alpha = classifier,
                   size = classifier
               )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  scale_color_manual(values = classifier_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.1, 1, 1)) +
  scale_size_manual(values = c(0.1, 2, 2)) +
  guides(color = "none") +
  guides(alpha = "none") +
  guides(size = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/mapping_0682_RN_Spp1_Glut_1_dynamic_filter.png", 
       plot = plot, 
       width = 5,
       height = 15, 
       dpi = 300)

#################### example for initial threshold too low ######################
# supertype: 0048 IT AON-TT-DP Glut_3
# sections: 1199651109 1199651112 1199651115 1199651121

cluster_to_plot <- "0048 IT AON-TT-DP Glut_3"
example_sections <- c("1199651109","1199651112","1199651115","1199651121")

classifier_palette <- c("lightgrey","black","red")

threshold <- 0.5

example_sections_plot <- metadata_sis %>%
  filter(basic_qc_passed == T) %>% 
  filter(doublets_qc_passed == T) %>%
  filter(section %in% example_sections) %>% 
  select(volume_x,
         volume_y,
         section,
         hrc_mmc_supertype_name,
         hrc_mmc_cluster_avg_correlation,
         hrc_mmc_cluster_thr)

# Create the classifier column
example_sections_plot <- example_sections_plot %>%
  mutate(classifier = case_when(
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation > threshold ~ "hq",
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation <= threshold ~ "lq",
    flat_CDM_supertype_name != cluster_to_plot ~ "bg"
  ))

plot <- ggplot(example_sections_plot,
               aes(x=x_coordinate,
                   y=y_coordinate,
                   color = classifier,
                   alpha = classifier,
                   size = classifier
               )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  scale_color_manual(values = classifier_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.1, 1, 1)) +
  scale_size_manual(values = c(0.1, 2, 2)) +
  guides(color = "none") +
  guides(alpha = "none") +
  guides(size = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/mapping_0048_IT_AON-TT-DP_Glut_3_fixed_filter.png", 
       plot = plot, 
       width = 5,
       height = 15, 
       dpi = 300)

threshold <- section_threshold[2,2]

example_sections_plot <- metadata_sis %>%
  filter(basic_qc_filter == F) %>% 
  filter(doublets_filter == F) %>%
  filter(section %in% example_sections) %>% 
  select(x_coordinate,
         y_coordinate,
         section,
         flat_CDM_supertype_name,
         flat_CDM_cluster_avg_correlation,
         flat_CDM_cluster_thr)

# Create the classifier column
example_sections_plot <- example_sections_plot %>%
  mutate(classifier = case_when(
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation > threshold ~ "hq",
    flat_CDM_supertype_name == cluster_to_plot & flat_CDM_cluster_avg_correlation <= threshold ~ "lq",
    flat_CDM_supertype_name != cluster_to_plot ~ "bg"
  ))

plot <- ggplot(example_sections_plot,
               aes(x=x_coordinate,
                   y=y_coordinate,
                   color = classifier,
                   alpha = classifier,
                   size = classifier
               )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  scale_color_manual(values = classifier_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.1, 1, 1)) +
  scale_size_manual(values = c(0.1, 1, 1)) +
  guides(color = "none") +
  guides(alpha = "none") +
  guides(size = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(filename = "/results/mapping_0048_IT_AON-TT-DP_Glut_3_dynamic_filter.png", 
       plot = plot, 
       width = 5,
       height = 15, 
       dpi = 300)


##################### example cases of bimodal distribution ###############

# find cases of bimodal distribution

bimdal_clusters <- metadata_sis %>% 
  filter(is_bimodal_cluster == TRUE) %>% 
  pull(flat_CDM_cluster_name) %>% 
  unique()

example_clusters <- c("5155 SPVI-SPVC Sall3 Lhx1 Gly-Gaba_2",
                      "5174 DCO Il22 Gly-Gaba_3")



# filter out relevant columns
filtered_metadata_sis <- metadata_sis %>%
  filter(basic_qc_filter == F) %>% 
  filter(doublets_filter == F) %>%
  filter(flat_CDM_cluster_name %in% example_clusters) %>% 
  select(cell_label,
         x_coordinate,
         y_coordinate,
         section,
         basic_qc_filter,
         doublets_filter,
         final_filter,
         is_bimodal_cluster,
         flat_CDM_class_name,
         flat_CDM_class_avg_correlation,
         CDM_class_color,
         flat_CDM_subclass_name,
         flat_CDM_subclass_avg_correlation,
         CDM_subclass_color,
         flat_CDM_supertype_name,
         flat_CDM_supertype_avg_correlation,
         flat_CDM_supertype_thr_criteria,
         CDM_supertype_color,
         flat_CDM_cluster_name,
         flat_CDM_cluster_avg_correlation,
         flat_CDM_cluster_thr,
         flat_CDM_cluster_thr_criteria,
         flat_CDM_thr,
         CDM_cluster_color)


# extract section threshold
section_threshold <- filtered_metadata_sis %>% 
  select(flat_CDM_cluster_name,
         flat_CDM_thr) %>% 
  unique() %>% 
  arrange(flat_CDM_cluster_name)


# plot distribution of correlation coefficients for selected 
plot <- ggplot(filtered_metadata_sis,
               aes(x = flat_CDM_cluster_name,
                   y = flat_CDM_cluster_avg_correlation,
                   fill = flat_CDM_cluster_name)) +
  geom_violin() +
  geom_boxplot(width = .2,
               alpha = .6,
               fatten = NULL,
               show.legend = F) +
  stat_summary(fun = "median",
               show.legend = F,
               position = position_dodge(.2)) +
  scale_fill_manual(values = cluster_color_palette) +
  geom_hline(yintercept = 0.5, 
             linetype = "dashed", 
             color = "black") +
  labs(x = "Cluster Name",
       y = "average correlation") + 
  geom_segment(aes(x = 0.75, 
                   xend = 1.25,
                   y = section_threshold[1,2],
                   yend = section_threshold[1,2],
                   color = "red")) +
  geom_segment(aes(x = 1.75, 
                   xend = 2.25,
                   y = section_threshold[2,2],
                   yend = section_threshold[2,2],
                   color = "red")) +theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = "none")


ggsave(filename = "/results/bimodal_filter.pdf", 
       plot = plot, 
       width = 10,
       height = 7, 
       dpi = 300)



cluster_to_plot <- "5174 DCO Il22 Gly-Gaba_3"

example_sections <- c("1199650965",
                      "1199650968",
                      "1199650975",
                      "1199650972")

