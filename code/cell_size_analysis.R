
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
})


options(stringsAsFactors = FALSE)
options(scipen=9999)
options(mc.cores=30)

compute_r2 <- function(data, volume, group) {
  volume <- enquo(volume)
  group <- enquo(group)
  
  # Convert to names for base R formula construction
  volume_str <- rlang::as_name(volume)
  group_str <- rlang::as_name(group)
  
  # Construct the formula
  formula <- as.formula(paste0(volume_str, " ~ ", group_str))
  
  # Fit the model
  model <- aov(formula, data = data)
  
  # Extract vector of the response
  y <- rlang::eval_tidy(volume, data)
  
  # Compute R²
  ss_total <- sum((y - mean(y))^2)
  ss_between <- anova(model)[["Sum Sq"]][1]
  r2 <- ss_between / ss_total
  
  return(r2)
}


# load metadata and extract relevant columns
# read in metadata file
metadata_sis <- fread("/data/merscope_638850_mouseadult_registered_v2/whole_dataset/mouse_638850_registered.csv")


# filter out relevant data
metadata_sis <- metadata_sis %>% 
  filter(final_qc_passed == T) %>% 
  select(hrc_mmc_subclass_name,
         hrc_mmc_supertype_name,
         volume)

load("/scratch/metadata_vpt.rda")

metadata_vpt <- metadata_vpt %>% 
  filter(final_filter == F) %>% 
  select(CDM_subclass_name,
         CDM_supertype_name,
         volume)

r2_method1 <- compute_r2(metadata_vpt,
                         volume,
                         CDM_supertype_name)

r2_method2 <- compute_r2(metadata_sis,
                         volume,
                         hrc_mmc_subclass_name)


observed_diff <- r2_method2 - r2_method1

cat(sprintf("Observed R² difference (Method2 - Method1): %.4f\n", observed_diff))

# -----------------------------
# Step 4: Permutation Test
# -----------------------------
set.seed(123)
n_perm <- 1000
perm_diffs <- numeric(n_perm)

for (i in 1:n_perm) {
  df1_perm <- metadata_vpt
  df2_perm <- metadata_sis
  df1_perm$CDM_supertype_name <- sample(metadata_vpt$CDM_supertype_name)
  df2_perm$hrc_mmc_subclass_name <- sample(metadata_sis$hrc_mmc_subclass_name)
  
  r2_1 <- compute_r2(df1_perm,
                     volume,
                     CDM_supertype_name)
  r2_2 <- compute_r2(df2_perm,
                     volume,
                     hrc_mmc_subclass_name)
  perm_diffs[i] <- r2_2 - r2_1
}

# -----------------------------
# Step 5: Calculate P-Value
# -----------------------------
p_val <- mean(perm_diffs >= observed_diff)
cat(sprintf("Permutation test p-value: %.4f\n", p_val))

# -----------------------------
# Step 6: Plot Results
# -----------------------------
hist(perm_diffs, breaks = 40, col = "gray", main = "Permutation Distribution of R² Differences",
     xlab = "R²(Method 2) - R²(Method 1)")
abline(v = observed_diff, col = "red", lwd = 2)
