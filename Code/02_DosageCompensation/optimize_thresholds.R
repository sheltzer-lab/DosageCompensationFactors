library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(ggplot2)
library(ggpubr)
library(ggsignif)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "visualization.R"))

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
  mutate(Gene.CopyNumber = if_else(Model.SampleType == "Normal", 2, Gene.CopyNumber)) %>%
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
df_results <- data.frame()

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
                Buffered.Optimal.Share = sum(Buffered.Estimate, na.rm = TRUE) / n(),
                Buffered.Estimate.Count = sum(Buffered.Estimate, na.rm = TRUE),
                Buffered.Estimate.Share = sum(Buffered.Estimate, na.rm = TRUE) / n(),
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
  summarize(Median.Count.Diff = median(abs(Buffered.Estimate.Count - Buffered.Optimal.Count)),
            Median.Buffered.Share = median(Buffered.Estimate.Share),
            Mean.Wasserstein = mean(Buffered.Distance.Wasserstein),
            Max.Wasserstein = max(Buffered.Distance.Wasserstein),
            Min.Wasserstein = min(Buffered.Distance.Wasserstein),
            Sum.Wasserstein = sum(Buffered.Distance.Wasserstein))

# === Plot protein abundance (Buffered vs. Disomic, Buffered vs. Scaling, Buffered vs. AntiScaling)
thresholds <- list(Buffered = 0.2, AntiScaling = 0.3)

datasets_opt_classes <- datasets_opt %>%
  mutate(Buffering.Class = buffering_class(Buffering.Ratio,
                                           Protein.Expression.Disomic, Protein.Expression.Normalized,
                                           2, Gene.CopyNumber, br_cutoffs_ = thresholds)) %>%
  drop_na(Buffering.Class)

datasets_opt_classes %>%
  filter(Dataset == "ProCan") %>%
  group_by(Gene.Symbol, Gene.CopyNumber, Buffering.Class) %>%
  summarize(Protein.Expression = mean(Protein.Expression.Normalized, na.rm = TRUE),
            Disomic = mean(Protein.Expression.Disomic, na.rm = TRUE),
            Log2FC = Protein.Expression - Disomic) %>%
  ggplot() +
  aes(x = Buffering.Class, y = Log2FC) +
  geom_hline(yintercept = 0, color = "red") +
  geom_hline(yintercept = c(log2(1/2), log2(3/2), log2(4/2)), color = default_color) +
  geom_boxplot(outliers = FALSE, alpha = 3/4) +
  facet_grid(~Gene.CopyNumber) +
  geom_signif(comparisons = list(c("Anti-Scaling", "Buffered"),
                                 c("Buffered", "Scaling"),
                                 c("Anti-Scaling", "Scaling")),
              test = wilcox.test, map_signif_level = print_signif, y_position = c(1, 1, 1.5),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  ylim(c(-3,3))

datasets_opt_classes %>%
  filter(Buffering.Class == "Buffered" & Dataset == "ProCan") %>%
  group_by(Gene.Symbol, Gene.CopyNumber) %>%
  summarize(Buffered = mean(Protein.Expression.Normalized),
            Disomic = mean(Protein.Expression.Disomic)) %>%
  pivot_longer(c("Buffered", "Disomic"),
               names_to = "Group", values_to = "Protein.Expression") %>%
  ggplot() +
  aes(x = Group, y = Protein.Expression) +
  geom_boxplot(outliers = FALSE) +
  geom_signif(comparisons = list(c("Buffered", "Disomic")),
              test = wilcox.test, map_signif_level = print_signif, y_position = 8,
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  facet_grid(~Gene.CopyNumber) +
  ylim(c(0,10))

# === Obtain alternative optimal thresholds ===
# Contraints: Buffered vs. Disomic insignificant, Buffered vs. Scaling significant, Buffered vs. AntiScaling significant
# Optimization: maximize sum of absolute threshold values
