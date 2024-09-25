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

# load corrected cell coordinates
load("/scratch/coordinates_sis.rda")
# load data segmented with in-house cellpose model
load("/scratch/metadata_sis.rda")

# add rotated coordinates to metadata
metadata_sis <- merge(metadata_sis,
                      coordinates_sis,
                      by = 0)

metadata_sis <- metadata_sis %>% 
  rownames_to_column( var = "cell_label")

# filter out relevant columns
filtered_metadata_sis <- metadata_sis %>%
  filter(basic_qc_filter == F) %>% 
  filter(doublets_filter == F) %>%
  filter(flat_CDM_supertype_name %in% example_supertypes) %>% 
  select(cell_label,
         x_coordinate,
         y_coordinate,
         section,
         basic_qc_filter,
         doublets_filter,
         final_filter,
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
  select(flat_CDM_supertype_name,
         flat_CDM_thr) %>% 
  unique() %>% 
  arrange(flat_CDM_supertype_name)


# plot distribution of correlation coefficients for selected 
plot <- ggplot(filtered_metadata_sis,
               aes(x = flat_CDM_supertype_name,
                   y = flat_CDM_cluster_avg_correlation,
                   fill = flat_CDM_supertype_name)) +
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
  geom_segment(aes(x = 0.75, 
                   xend = 1.25,
                   y = section_threshold[1,2],
                   yend = section_threshold[1,2],
                   color = "red")) +
  geom_segment(aes(x = 1.75, 
                   xend = 2.25,
                   y = section_threshold[2,2],
                   yend = section_threshold[2,2],
                   color = "red")) +
  geom_segment(aes(x = 2.75, 
                   xend = 3.25,
                   y = section_threshold[3,2],
                   yend = section_threshold[3,2],
                   color = "red")) +
  geom_segment(aes(x = 3.75, 
                   xend = 4.25,
                   y = section_threshold[4,2],
                   yend = section_threshold[4,2],
                   color = "red")) +
  geom_segment(aes(x = 4.75, 
                   xend = 5.25,
                   y = section_threshold[5,2],
                   yend = section_threshold[5,2],
                   color = "red")) +
  geom_segment(aes(x = 5.75, 
                   xend = 6.25,
                   y = section_threshold[6,2],
                   yend = section_threshold[6,2],
                   color = "red")) +
  geom_segment(aes(x = 6.75, 
                   xend = 7.25,
                   y = section_threshold[7,2],
                   yend = section_threshold[7,2],
                   color = "red")) +
  geom_segment(aes(x = 7.75, 
                   xend = 8.5,
                   y = section_threshold[8,2],
                   yend = section_threshold[8,2],
                   color = "red")) +
  geom_segment(aes(x =8.75, 
                   xend = 9.25,
                   y = section_threshold[9,2],
                   yend = section_threshold[9,2],
                   color = "red")) +
  geom_segment(aes(x =9.75, 
                   xend = 10.25,
                   y = section_threshold[10,2],
                   yend = section_threshold[10,2],
                   color = "red")) +
  geom_segment(aes(x =10.75, 
                   xend = 11.25,
                   y = section_threshold[11,2],
                   yend = section_threshold[11,2],
                   color = "red")) +
  geom_segment(aes(x =11.75, 
                   xend = 12.25,
                   y = section_threshold[12,2],
                   yend = section_threshold[12,2],
                   color = "red")) +
  geom_segment(aes(x =12.75, 
                   xend = 13.25,
                   y = section_threshold[13,2],
                   yend = section_threshold[13,2],
                   color = "red")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = "none")


ggsave(filename = "/results/mapping_filter.pdf", 
       plot = plot, 
       width = 10,
       height = 7, 
       dpi = 300)

# supertype:0023 L4/5 IT CTX Glut_1
# sections: 1199651075 1199651106 1199651057

# supertype: 1076 NLL Gata3 Gly-Gaba_3
# sections: 1199650981 1199650978 1199650975



example_sections_plot <- filtered_metadata_sis %>% 
  filter(section %in% example_sections)

ggplot(example_sections_plot,
       aes(x=x_coordinate,
           y=y_coordinate,
           color = flat_CDM_supertype_name,
           alpha = flat_CDM_cluster_thr_criteria
       )) +
  geom_point(size=.5,stroke=0,shape=19,) +
  coord_fixed() +
  ggtitle("Gene Filter") +
  scale_color_manual(values = supertype_color_palette) +
  scale_y_reverse() +
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.5, 1)) +
  guides(color = "none") +
  facet_wrap(~fct_rev(section), ncol = 1) +
  theme_void() +
  theme(strip.text = element_blank())
