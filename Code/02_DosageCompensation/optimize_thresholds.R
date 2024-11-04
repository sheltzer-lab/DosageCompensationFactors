library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(ggplot2)
library(ggpubr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "buffering_ratio.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "DosageCompensation")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_cptac <- read_parquet((here(output_data_dir, 'copy_number_cptac.parquet')))

expr_procan <- read_parquet(here(output_data_dir, "expression_procan.parquet"))
expr_depmap <- read_parquet(here(output_data_dir, "expression_depmap.parquet"))
expr_cptac <- read_parquet(here(output_data_dir, 'expression_cptac.parquet'))
meta_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))

# === Plot Protein Abundance distributions per Copy Number ===
## DepMap
df_depmap <- expr_depmap %>%
  inner_join(y = copy_number, by = c("Model.ID", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one")

df_depmap %>%
  filter(Gene.CopyNumber < 5 & Gene.CopyNumber > 0) %>%
  filter(Protein.Expression.Normalized > -4 & Protein.Expression.Normalized < 4) %>% # For visualization only
  mutate(Gene.CopyNumber = as.factor(Gene.CopyNumber)) %>%
  ggdensity("Protein.Expression.Normalized",
            add = "mean", rug = FALSE,
            color = "Gene.CopyNumber", fill = "Gene.CopyNumber") +
  scale_x_continuous(breaks = c(seq(-0.8, 0.8, 0.2), seq(-4, -1, 0.5), seq(1, 4, 0.5)))

## ProCan
df_procan <- expr_procan %>%
  inner_join(y = copy_number, by = c("Model.ID", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one")

df_procan %>%
  filter(Gene.CopyNumber < 5 & Gene.CopyNumber > 0) %>%
  #filter(Protein.Expression.Normalized > -4 & Protein.Expression.Normalized < 4) %>% # For visualization only
  mutate(Gene.CopyNumber = as.factor(Gene.CopyNumber)) %>%
  ggdensity("Protein.Expression.Normalized",
            add = "mean", rug = FALSE,
            color = "Gene.CopyNumber", fill = "Gene.CopyNumber")
  #+ scale_x_continuous(breaks = c(seq(-0.8, 0.8, 0.2), seq(-4, -1, 0.5), seq(1, 4, 0.5)))

## CPTAC
df_cptac <- expr_cptac %>%
  inner_join(y = copy_number_cptac, by = c("Model.ID", "Gene.ENSEMBL.Id"),
             na_matches = "never", relationship = "many-to-one")

df_cptac %>%
  filter(Gene.CopyNumber < 5 & Gene.CopyNumber > 0) %>%
  filter(Protein.Expression.Normalized > 15 & Protein.Expression.Normalized < 35) %>% # For visualization only
  mutate(Gene.CopyNumber = as.factor(Gene.CopyNumber)) %>%
  ggdensity("Protein.Expression.Normalized",
            add = "mean", rug = FALSE,
            color = "Gene.CopyNumber", fill = "Gene.CopyNumber")
  #+ scale_x_continuous(breaks = c(seq(-0.8, 0.8, 0.2), seq(-4, -1, 0.5), seq(1, 4, 0.5)))

# === Define "optimal" Buffered distributions on DepMap, ProCan & CPTAC ===
# Thresholds defining the "optimal" thresholds for change in protein abundance to disomy per copy number alteration
buffered_log2fc <- list(
  "1" = c(low = log2(0.7), high = log2(1.05)),  # Monosomy, equivalent to protein abundance change between [-30%, +5%]
  "3" = c(low = log2(0.95), high = log2(1.3)),  # Trisomy, equivalent to protein abundance change between [-5%, +30%]
  "4" = c(low = log2(0.95), high = log2(1.6))   # Tetrasomy, equivalent to protein abundance change between [-5%, +60%]
)
# Optimal buffered distributions can alternatively be given by example

datasets <- bind_rows(df_depmap, df_procan, df_cptac) %>%
  drop_na(Protein.Expression.Normalized) %>%
  group_by(Dataset, Model.ID, Gene.Symbol, Gene.CopyNumber) %>%
  summarize(Protein.Expression.Normalized = mean(Protein.Expression.Normalized), .groups = "drop") %>%
  group_by(Dataset, Gene.Symbol) %>%
  mutate(Protein.Expression.Disomic = median(Protein.Expression.Normalized[Gene.CopyNumber == 2])) %>%
  ungroup()

datasets_opt <- datasets %>%
  drop_na(Protein.Expression.Disomic) %>%
  filter(Gene.CopyNumber %in% c(1, 3, 4)) %>%
  mutate(Log2FC = Protein.Expression.Normalized - Protein.Expression.Disomic) %>%
  mutate(Threshold.Optimal = buffered_log2fc[as.character(Gene.CopyNumber)]) %>%
  unnest_wider(Threshold.Optimal) %>%
  mutate(Buffered.Optimal = low < Log2FC & Log2FC < high,
         Buffering.Ratio = buffering_ratio(Protein.Expression.Disomic, Protein.Expression.Normalized,
                                           2, Gene.CopyNumber))

# === Iterate through threshold combinations and determine distance of estimated to optimal Buffered cohorts ===
df_results <- data.frame(
  Dataset = character(),
  Gene.CopyNumber = integer(),
  Threshold.Buffered = numeric(),
  Threshold.AntiScaling = numeric(),
  Buffered.Optimal.Count = integer(),
  Buffered.Estimate.Count = integer(),
  Buffered.Distance.Wasserstein = numeric(),
  Buffered.Distance.Log2FC = numeric()
)

for (thresh_buf in seq(0, 1, 0.1)) {
  for (thresh_as in seq(0, 1, 0.1)) {
    thresholds <- list(Buffered = thresh_buf, AntiScaling = thresh_as)

    df_thresh_res <- datasets_opt %>%
      mutate(Buffering.Class = buffering_class(Buffering.Ratio,
                                                  Protein.Expression.Disomic, Protein.Expression.Normalized,
                                                  2, Gene.CopyNumber, br_cutoffs_ = thresholds),
             Buffered.Estimate = Buffering.Class == "Buffered") %>%
      drop_na()

    buffered_mising <- df_thresh_res %>%
      filter(Buffered.Estimate) %>%
      count(Dataset, Gene.CopyNumber, Buffered.Estimate) %>%
      complete(Dataset, Gene.CopyNumber, Buffered.Estimate) %>%
      pull(n) %>%
      is.na() %>%
      any()

    if (buffered_mising) next

    df_results <- df_thresh_res %>%
      group_by(Dataset, Gene.CopyNumber) %>%
      summarize(Buffered.Optimal.Count = sum(Buffered.Optimal, na.rm = TRUE),
                Buffered.Estimate.Count = sum(Buffered.Estimate, na.rm = TRUE),
                Buffered.Distance.Wasserstein = transport::wasserstein1d(Log2FC[Buffered.Estimate], Log2FC[Buffered.Optimal]),
                Buffered.Distance.Log2FC = mean(Log2FC[Buffered.Estimate]) - mean(Log2FC[Buffered.Optimal]),
                .groups = "drop") %>%
      mutate(Threshold.Buffered = thresh_buf,
             Threshold.AntiScaling = thresh_as) %>%
      bind_rows(df_results)
  }
}

# === Are changes in Wasserstein distance significantly different between datasets? ===

# === Get optimal thresholds ===
df_opt_thresholds <- df_results %>%
  group_by(Threshold.Buffered, Threshold.AntiScaling) %>%
  summarize(Median.Count.Diff = median(Buffered.Estimate.Count - Buffered.Optimal.Count),
            Mean.Wasserstein = mean(Buffered.Distance.Wasserstein),
            Max.Wasserstein = max(Buffered.Distance.Wasserstein),
            Min.Wasserstein = min(Buffered.Distance.Wasserstein),
            Sum.Wasserstein = sum(Buffered.Distance.Wasserstein))

df_opt_thresholds %>%
  filter(Threshold.Buffered == 0.2) %>%
  arrange(Median.Count.Diff, Mean.Wasserstein)

df_opt_thresholds %>%
  filter(Threshold.AntiScaling > 0.2 & Threshold.AntiScaling < 0.4) %>%
  arrange(Median.Count.Diff, Mean.Wasserstein)

# === Plot protein abundance (Buffered vs. Disomic, Buffered vs. Scaling, Buffered vs. AntiScaling)

# === Obtain alternative optimal thresholds ===
# Contraints: Buffered vs. Disomic insignificant, Buffered vs. Scaling significant, Buffered vs. AntiScaling significant
# Optimization: maximize sum of absolute threshold values
