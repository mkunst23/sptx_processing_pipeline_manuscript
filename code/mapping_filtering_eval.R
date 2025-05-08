################ visualizing improved mapping filtering ##################

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
  library(devtools)
  library(RColorBrewer)
  library(googlesheets4)
  library(forcats)
})

options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

gs4_deauth()

# read in metadata file
metadata <- fread("/data/merscope_638850_mouseadult_registered_v2/whole_dataset/mouse_638850_registered.csv")

# subset relevant columns
metadata_filtered <- metadata %>% 
  select(production_cell_id,
         basic_qc_passed,
         doublets_qc_passed,
         final_qc_passed,
         hrc_mmc_cluster_avg_correlation,
         hrc_mmc_class_name,
         hrc_mmc_subclass_name,
         hrc_mmc_supertype_name,
         hrc_mmc_cluster_name)

metadata_filtered$hrc_mmc_supertype_name <- as.character(metadata_filtered$hrc_mmc_supertype_name)

# filter out basic qc and doublets
metadata_filtered <- metadata_filtered %>% 
  filter(basic_qc_passed == T) %>% 
  filter(doublets_qc_passed == T) %>% 
  filter(!is.na(hrc_mmc_supertype_name)) %>% 
  filter(hrc_mmc_supertype_name != "NA")

# determine cells that are being filtered out by MAD
bad_cell_mad <- metadata_filtered %>% 
  filter(final_qc_passed == F) %>% 
  pull(production_cell_id) %>% 
  length()

# determine cells that are being filterd out by fixed threshold
bad_cells_old <- metadata_filtered %>% 
  filter(hrc_mmc_cluster_avg_correlation <= 0.5) %>% 
  pull(production_cell_id) %>% 
  length()


failure_rates_by_supertype <- metadata_filtered %>%
  filter(!is.na(hrc_mmc_supertype_name)) %>%
  group_by(hrc_mmc_supertype_name) %>%
  summarise(
    old_failure_rate = mean(hrc_mmc_cluster_avg_correlation <= 0.5, na.rm = TRUE) * 100,
    new_failure_rate = mean(final_qc_passed == FALSE, na.rm = TRUE) * 100
  )


failure_rates_by_supertype <- metadata_filtered |> 
  filter(!is.na(hrc_mmc_supertype_name)) |> 
  group_by(hrc_mmc_supertype_name) |> 
  summarise(
    total_cells = n(),
    cells_old_fail = sum(hrc_mmc_cluster_avg_correlation <= 0.5, na.rm = TRUE),
    cells_new_fail = sum(final_qc_passed == FALSE, na.rm = TRUE),
    old_failure_rate = mean(hrc_mmc_cluster_avg_correlation <= 0.5, na.rm = TRUE) * 100,
    new_failure_rate = mean(final_qc_passed == FALSE, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# convert tot long format
failure_rate_long <- failure_rates_by_supertype %>% 
  pivot_longer(
    cols = c(old_failure_rate,new_failure_rate),
    names_to = "filtering_method",
    values_to = "percentage"
  )


devtools::source_url("https://raw.githubusercontent.com/yjunechoe/geom_paired_raincloud/master/geom_paired_raincloud.R")

plot <- ggplot(failure_rate_long, aes(fct_rev(filtering_method), 
                   percentage, 
                   fill = filtering_method)) +
  geom_paired_raincloud(alpha = .5) +
  geom_point(aes(group = hrc_mmc_supertype_name),
             position = position_nudge(c(.001, -.001)),
             alpha = .5, 
             shape = 16) +
  geom_line(aes(group = hrc_mmc_supertype_name),
            position = position_nudge(c(.13, -.13)),
            linetype = 1,
            linewidth = 0.1) +
  geom_boxplot(position = position_nudge(c(.07, -.07)),
               alpha = .5, width = .04, outlier.shape = " ") +
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save plot
ggsave(filename = "/results/mapping_filter_qc.pdf", 
       plot = plot, 
       width = 12,
       height = 8, 
       dpi = 300)

# add difference to failure rate
failure_rates_by_supertype$diff <- failure_rates_by_supertype$old_failure_rate - failure_rates_by_supertype$new_failure_rate


# save as csv file
fwrite(failure_rates_by_supertype,
       "/results/failure_rate_by_supertype.csv",
       row.names = F)
