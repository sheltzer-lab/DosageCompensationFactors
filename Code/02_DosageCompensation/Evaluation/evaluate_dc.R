library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)
library(purrr)


here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "analysis.R"))
source(here("Code", "evaluation.R"))
source(here("Code", "02_DosageCompensation", "baseline.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "DosageCompensation")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

dc_report <- list()

expr_buf_procan <- read_parquet(here(output_data_dir, 'expression_buffering_procan.parquet'))
expr_buf_depmap <- read_parquet(here(output_data_dir, 'expression_buffering_depmap.parquet'))
expr_depmap <- read_parquet(here(output_data_dir, "expression_depmap.parquet"))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
dc_factors_ids <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet")) %>%
  select(Protein.Uniprot.Accession, Gene.Symbol)


add_debug_data <- function(df_dc) {
    df_dc %>%
      mutate(Debug.ScalingDirection.Gene = sign(2^Protein.Expression.Normalized - 2^Protein.Expression.Baseline) *
        sign(Gene.CopyNumber - Gene.CopyNumber.Baseline),
             Debug.ScalingDirection.Chromosome = sign(2^Protein.Expression.Normalized - 2^Protein.Expression.Baseline) *
               sign(ChromosomeArm.CopyNumber - ChromosomeArm.CopyNumber.Baseline),
             Log2FC.CopyNumber.Gene = log2(Gene.CopyNumber) - log2(Gene.CopyNumber.Baseline),
             Log2FC.CopyNumber.Chromosome = log2(ChromosomeArm.CopyNumber) - log2(ChromosomeArm.CopyNumber.Baseline))
}

expr_buf_procan <- add_debug_data(expr_buf_procan)
expr_buf_depmap <- add_debug_data(expr_buf_depmap)

# Buffering Ratio
(ggplot() +
  geom_density(aes(x = expr_buf_depmap$Buffering.GeneLevel.Ratio, color = "DepMap")) +
  geom_density(aes(x = expr_buf_procan$Buffering.GeneLevel.Ratio, color = "ProCan")) +
  labs(x = "Gene Level Buffering Ratio", color = "Dataset")) %>%
  save_plot("buffering_ratio_distribution.png")

# TODO: Visualize for DepMap

plot_br_distributions <- function(dataset, dataset_name) {
    max_abs_br <- ceiling(max(abs(dataset$Buffering.GeneLevel.Ratio), na.rm = TRUE))

    (dataset %>%
      mutate(AntiScaling = Debug.ScalingDirection.Gene < 0) %>%
      ggplot() +
      aes(x = Buffering.GeneLevel.Ratio, color = AntiScaling) +
      geom_density() +
      geom_vline(xintercept = br_cutoffs$Buffered, color = "orange", linetype = "dashed") +
      scale_x_continuous(limits = c(-max_abs_br, max_abs_br), breaks = seq(-max_abs_br, max_abs_br, 1))) %>%
      save_plot(paste0("buffering_ratio_anti-scaling_", dataset_name, ".png"))

    (dataset %>%
      ggplot() +
      aes(x = Buffering.GeneLevel.Ratio, color = as.factor(Debug.ScalingDirection.Gene)) +
      geom_density() +
      geom_vline(xintercept = br_cutoffs$Buffered, color = "orange", linetype = "dashed") +
      scale_x_continuous(limits = c(-max_abs_br, max_abs_br), breaks = seq(-max_abs_br, max_abs_br, 1))) %>%
      save_plot(paste0("buffering_ratio_scalingdirection_", dataset_name, ".png"))

    (dataset %>%
      ggplot() +
      aes(x = Buffering.GeneLevel.Ratio, color = Buffering.GeneLevel.Class) +
      geom_density() +
      geom_vline(xintercept = br_cutoffs$Buffered, color = "orange", linetype = "dashed") +
      scale_x_continuous(limits = c(-max_abs_br, max_abs_br), breaks = seq(-max_abs_br, max_abs_br, 1))) %>%
      save_plot(paste0("buffering_ratio_classes_", dataset_name, ".png"))

    (dataset %>%
      filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
      mutate(Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
      ggplot() +
      aes(x = Buffering.GeneLevel.Ratio, color = Event) +
      geom_density() +
      geom_vline(xintercept = br_cutoffs$Buffered, color = "orange", linetype = "dashed") +
      scale_x_continuous(limits = c(-max_abs_br, max_abs_br), breaks = seq(-max_abs_br, max_abs_br, 1))) %>%
      save_plot(paste0("buffering_ratio_event_", dataset_name, ".png"))

    (dataset %>%
      filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
      mutate(Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
      ggplot() +
      aes(x = Buffering.GeneLevel.Ratio, color = Buffering.GeneLevel.Class) +
      geom_density() +
      facet_grid(~Event) +
      geom_vline(xintercept = br_cutoffs$Buffered, color = "orange", linetype = "dashed") +
      scale_x_continuous(limits = c(-max_abs_br, max_abs_br), breaks = seq(-max_abs_br, max_abs_br, 1))) %>%
      save_plot(paste0("buffering_ratio_classes_by-event_", dataset_name, ".png"), width = 400)

    return()
}

plot_br_distributions(expr_buf_procan, "procan")
plot_br_distributions(expr_buf_depmap, "depmap")

## Distribution of Confidence Score
# TODO: Facet by Gene vs. ChrArm Level
ggplot() +
  geom_density(aes(x = expr_buf_depmap$Buffering.GeneLevel.Ratio.Confidence, color = "DepMap")) +
  geom_density(aes(x = expr_buf_procan$Buffering.GeneLevel.Ratio.Confidence, color = "ProCan")) +
  labs(x = "Gene Level Buffering Ratio (Confidence Score)", color = "Dataset")

# Buffering Classes
## Similarity between Chromosome-Level Avg Log2FC classes and Gene-Level BR classes
classes_chr_log2fc_avg <- expr_buf_depmap %>%
  filter(abs(ChromosomeArm.CNA) > 0) %>%
  mutate(Event = if_else(ChromosomeArm.CNA > 0, "Gain", "Loss")) %>%
  select(Event, Gene.Symbol, Buffering.ChrArmLevel.Average.Class) %>%
  drop_na() %>%
  group_by(Event, Gene.Symbol) %>%
  summarize(BufferingClass = first(Buffering.ChrArmLevel.Average.Class), .groups = "drop")

classes_gene_br <- expr_buf_depmap %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  select(Event, Gene.Symbol, Buffering.GeneLevel.Class) %>%
  rename(BufferingClass = Buffering.GeneLevel.Class) %>%
  drop_na() %>%
  group_by(Event, Gene.Symbol) %>%
  count(BufferingClass, name = "Count") %>%
  slice_max(Count, n = 1) %>%
  ungroup()

classes_chr_br <- expr_buf_depmap %>%
  filter(abs(ChromosomeArm.CNA) > 0) %>%
  mutate(Event = if_else(ChromosomeArm.CNA > 0, "Gain", "Loss")) %>%
  select(Event, Gene.Symbol, Buffering.ChrArmLevel.Class) %>%
  rename(BufferingClass = Buffering.ChrArmLevel.Class) %>%
  drop_na() %>%
  group_by(Event, Gene.Symbol) %>%
  count(BufferingClass, name = "Count") %>%
  slice_max(Count, n = 1) %>%
  ungroup()

jaccard_index <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return(intersection/union)
}

genes_chr_log2fc_avg <- classes_chr_log2fc_avg %>%
  group_by(Event, BufferingClass) %>%
  summarize(Genes.Chr.Log2FC.Avg = list(Gene.Symbol), .groups = "drop")

genes_chr_br <- classes_chr_br %>%
  group_by(Event, BufferingClass) %>%
  summarize(Genes.Chr.BR = list(Gene.Symbol), .groups = "drop")

genes_gene_br <- classes_gene_br %>%
  group_by(Event, BufferingClass) %>%
  summarize(Genes.BR = list(Gene.Symbol), .groups = "drop")

genes_per_class <- list(genes_chr_log2fc_avg, genes_chr_br, genes_gene_br) %>%
  purrr::reduce(inner_join, by = c("Event", "BufferingClass"), relationship = "one-to-one") %>%
  group_by(Event, BufferingClass) %>%
  mutate(JaccardSimilarity.ChrAvgLog2FC.GeneBR = jaccard_index(unlist(Genes.Chr.Log2FC.Avg), unlist(Genes.BR)),
         JaccardSimilarity.ChrLog2FC.ChrBR = jaccard_index(unlist(Genes.Chr.Log2FC.Avg), unlist(Genes.Chr.BR)),
         JaccardSimilarity.ChrBR.GeneBR = jaccard_index(unlist(Genes.Chr.BR), unlist(Genes.BR))) %>%
  ungroup()

dc_report$ClassSimilarities <- genes_per_class %>%
  select(Event, BufferingClass, starts_with("Jaccard")) %>%
  df2reportlist()

## Alternative: Match Frequency of buffering classes per gene, class, and event
gene_class_match_chr <- expr_buf_depmap %>%
  filter(abs(ChromosomeArm.CNA) > 0) %>%
  mutate(Event = if_else(ChromosomeArm.CNA > 0, "Gain", "Loss")) %>%
  select(Event, Gene.Symbol, Buffering.ChrArmLevel.Average.Class, Buffering.ChrArmLevel.Class) %>%
  mutate(Match = Buffering.ChrArmLevel.Average.Class == Buffering.ChrArmLevel.Class) %>%
  group_by(Event, Buffering.ChrArmLevel.Average.Class) %>%
  add_count(Gene.Symbol, name = "Total") %>%
  ungroup() %>%
  group_by(Event, Gene.Symbol, Buffering.ChrArmLevel.Average.Class) %>%
  summarize(MatchFrequency = sum(Match) / first(Total), .groups = "drop") %>%
  drop_na()

class_summary <- gene_class_match_chr %>%
  group_by(Event, Buffering.ChrArmLevel.Average.Class) %>%
  summarize(AverageMatchFrequency = mean(MatchFrequency), .groups = "drop")

dc_report$ClassMatchFrequency <- df2reportlist(class_summary)

# Check if buffered proteins are significantly downregulated (upon CN gain) in comparison to their scaling counterparts
# Idea: Sets of genes classified as buffered should have lower expression (if CN gained) and should be close to
# baseline expression (-> Anti-Scaling genes classified as buffering may drag expression down artificially)
# Goal: Evaluate thresholds

buffered_expression_chr_gain <- expr_buf_depmap %>%
  select(Gene.Symbol, ChromosomeArm.CNA, Buffering.ChrArmLevel.Class, Protein.Expression.Normalized) %>%
  filter(Buffering.ChrArmLevel.Class != "Anti-Scaling") %>%
  filter(ChromosomeArm.CNA > 0) %>%
  differential_expression(Gene.Symbol, Buffering.ChrArmLevel.Class,
                          Protein.Expression.Normalized, groups = c("Scaling", "Buffered"))

signif_buffered_chr_gain <- buffered_expression_chr_gain %>%
  filter(Log2FC < 0) %>%
  filter(TTest.p.adj < p_threshold) %>%
  nrow() %>%
  (\(x) x/nrow(buffered_expression_chr_gain))

buffered_expression_chr_loss <- expr_buf_depmap %>%
  select(Gene.Symbol, ChromosomeArm.CNA, Buffering.ChrArmLevel.Class, Protein.Expression.Normalized) %>%
  filter(Buffering.ChrArmLevel.Class != "Anti-Scaling") %>%
  filter(ChromosomeArm.CNA < 0) %>%
  differential_expression(Gene.Symbol, Buffering.ChrArmLevel.Class,
                          Protein.Expression.Normalized, groups = c("Scaling", "Buffered"))

signif_buffered_chr_loss <- buffered_expression_chr_loss %>%
  filter(Log2FC > 0) %>%
  filter(TTest.p.adj < p_threshold) %>%
  nrow() %>%
  (\(x) x/nrow(buffered_expression_chr_loss))

dc_report$SignifBufferingShare <- list(ChrGain = signif_buffered_chr_gain,
                                       ChrLoss = signif_buffered_chr_loss)

## Expression Distance to Baseline (set of genes classified as buffered should minimize this distance)
expr_baseline_dist <- expr_buf_depmap %>%
  select(Gene.Symbol, ChromosomeArm.CNA, Buffering.ChrArmLevel.Class,
         Protein.Expression.Normalized, Protein.Expression.Baseline) %>%
  filter(Buffering.ChrArmLevel.Class == "Buffered") %>%
  filter(ChromosomeArm.CNA != 0) %>%
  group_by(Gene.Symbol, ChromosomeArm.CNA) %>%
  summarize(Log2FC = mean(Protein.Expression.Normalized, na.rm = TRUE) -
    mean(Protein.Expression.Baseline, na.rm = TRUE), .groups = "drop")

exp_baseline_dist_avg <- expr_baseline_dist %>%
  group_by(ChromosomeArm.CNA) %>%
  summarize(AvgLog2FC = mean(Log2FC))

antiscaling_as_buffered <- expr_baseline_dist %>%
  filter(sign(ChromosomeArm.CNA) != sign(Log2FC))

dc_report$Buffered_Expression2Baseline_AverageDistance <- df2reportlist(exp_baseline_dist_avg)
dc_report$Buffered_Expression2Baseline_AntiScaling <- df2reportlist(antiscaling_as_buffered)

# Baseline Calculation
num_baseline_samples <- expr_buf_procan %>%
  mutate(RoundedPloidy = round(CellLine.Ploidy)) %>%
  filter(ChromosomeArm.CNA == 0) %>%
#  filter(RoundedPloidy == 2) %>%
  group_by(Gene.Symbol) %>%
  summarize(Observations = sum(is.na(Protein.Expression.Normalized)))

# TODO: Check if min number of observations is ensured by filter_genes()
dc_report$AvailableBaselineObservations <- quantile(num_baseline_samples$Observations)

# Copy Number Baseline
df_cn_eval <- expr_depmap %>%
  select(Gene.Symbol, CellLine.CustomId) %>%
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  select(Gene.Symbol, CellLine.CustomId, Gene.CopyNumber, CellLine.Ploidy,
         ChromosomeArm.CNA, CellLine.AneuploidyScore) %>%
  mutate(PloidyDistance = abs(2 - CellLine.Ploidy),
         CombinedDistance = 0.3 * PloidyDistance / max(PloidyDistance, na.rm = TRUE) +
           0.7 * CellLine.AneuploidyScore / max(CellLine.AneuploidyScore, na.rm = TRUE)) %>%
  calculate_median_baseline(Gene.Symbol, Gene.CopyNumber, target_colname = "MedianAll") %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     distance_col = CellLine.AneuploidyScore, target_colname = "WeightedNeutral.AneuploidyScore",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     distance_col = PloidyDistance, target_colname = "WeightedNeutral.PloidyDistance",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     distance_col = CombinedDistance, target_colname = "WeightedNeutral.CombinedDistance",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     distance_col = CellLine.AneuploidyScore, target_colname = "MeanNeutral", weighted = FALSE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber, summ_func = median,
                     distance_col = CellLine.AneuploidyScore, target_colname = "MedianNeutral", weighted = FALSE) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol, MedianAll, WeightedNeutral.AneuploidyScore, WeightedNeutral.PloidyDistance, WeightedNeutral.CombinedDistance,
         MeanNeutral, MedianNeutral)

dc_report$corr_CN_BaseLine_weightedAS_MedianAll <- cor.test(df_cn_eval$WeightedNeutral.AneuploidyScore,
                                                            df_cn_eval$MedianAll, method = "spearman")

cn_baseline_plot <- df_cn_eval %>%
  pivot_longer(c("MedianAll", "WeightedNeutral.AneuploidyScore", "WeightedNeutral.PloidyDistance", "WeightedNeutral.CombinedDistance",
                 "MeanNeutral", "MedianNeutral"),
               names_to = "Methods",
               values_to = "Gene.CopyNumber.Baseline") %>%
  ggplot() +
  aes(color = Methods, x = Gene.CopyNumber.Baseline) +
  geom_density()

cn_baseline_plot %>%
  save_plot("copynumber_baseline_methods.png")

# Expression Baseline
df_expr_eval <- expr_depmap %>%
  select(Gene.Symbol, CellLine.CustomId, Protein.Expression.Normalized) %>%
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  select(Gene.Symbol, CellLine.CustomId, ChromosomeArm.CNA, CellLine.Ploidy,
         CellLine.AneuploidyScore, Protein.Expression.Normalized) %>%
  mutate(PloidyDistance = abs(2 - CellLine.Ploidy),
         CombinedDistance = 0.3 * PloidyDistance / max(PloidyDistance, na.rm = TRUE) +
           0.7 * CellLine.AneuploidyScore / max(CellLine.AneuploidyScore, na.rm = TRUE)) %>%
  calculate_median_baseline(Gene.Symbol, Protein.Expression.Normalized, target_colname = "MedianAll") %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     distance_col = CellLine.AneuploidyScore, target_colname = "WeightedNeutral.AneuploidyScore",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     distance_col = PloidyDistance, target_colname = "WeightedNeutral.PloidyDistance",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     distance_col = CombinedDistance, target_colname = "WeightedNeutral.CombinedDistance",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     distance_col = CellLine.AneuploidyScore, target_colname = "MeanNeutral", weighted = FALSE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, summ_func = median,
                     distance_col = CellLine.AneuploidyScore, target_colname = "MedianNeutral", weighted = FALSE) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol, MedianAll, WeightedNeutral.AneuploidyScore, WeightedNeutral.PloidyDistance, WeightedNeutral.CombinedDistance,
         MeanNeutral, MedianNeutral)

dc_report$corr_Expr_Baseline_weightedAS_MeanNeutral <- cor.test(df_expr_eval$WeightedNeutral.AneuploidyScore,
                                                                df_expr_eval$MeanNeutral, method = "spearman")

expr_baseline_plot <- df_expr_eval %>%
  pivot_longer(c("MedianAll", "WeightedNeutral.AneuploidyScore", "WeightedNeutral.PloidyDistance", "WeightedNeutral.CombinedDistance",
                 "MeanNeutral", "MedianNeutral"),
               names_to = "Methods",
               values_to = "Protein.Expression.Baseline") %>%
  ggplot() +
  aes(color = Methods, x = Protein.Expression.Baseline) +
  geom_density()

expr_baseline_plot %>%
  save_plot("expression_baseline_methods.png")

# Check ratio of data lost after join with factor data
data_loss_procan <- (expr_buf_procan %>%
  anti_join(y = dc_factors_ids,
            by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never") %>%
  nrow()) / nrow(expr_buf_procan)

dc_report$data_loss_procan <- data_loss_procan

data_loss_depmap <- (expr_buf_depmap %>%
  anti_join(y = dc_factors_ids,
            by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never") %>%
  nrow()) / nrow(expr_buf_depmap)

dc_report$data_loss_depmap <- data_loss_depmap

# Check correlation between DC Scores and Protein Coefficient of Variance
# Buffering Score should be negatively correlated with Protein Coefficient of Variance
dc_cv <- expr_buf_procan %>%
  select(Gene.Symbol, Protein.Expression.Normalized, Buffering.GeneLevel.Ratio, Buffering.GeneLevel.SF) %>%
  group_by(Gene.Symbol) %>%
  mutate(Protein.Expression.CV = sd(2^Protein.Expression.Normalized, na.rm = TRUE) /
    mean(2^Protein.Expression.Normalized, na.rm = TRUE),
         Buffering.Ratio.Average = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
         Buffering.SF.Average = mean(Buffering.GeneLevel.SF, na.rm = TRUE))

dc_report$corr_CV_BufRatio_Average <- cor.test(dc_cv$Protein.Expression.CV, dc_cv$Buffering.Ratio.Average,
                                               method = "spearman")
dc_report$corr_CV_SF_Average <- cor.test(dc_cv$Protein.Expression.CV, dc_cv$Buffering.SF.Average,
                                         method = "spearman")
dc_report$corr_CV_BufRatio <- cor.test(dc_cv$Protein.Expression.CV, dc_cv$Buffering.GeneLevel.Ratio,
                                       method = "spearman")
dc_report$corr_CV_SF <- cor.test(dc_cv$Protein.Expression.CV, dc_cv$Buffering.GeneLevel.SF,
                                 method = "spearman")

## Check correlation between datasets
buf_matched <- match_datasets(expr_buf_procan, expr_buf_depmap)
corr_gene <- dataset_correlation(buf_matched,
                                 Dataset, Buffering.GeneLevel.Ratio,
                                 "Gene")
corr_chr <- dataset_correlation(buf_matched,
                                Dataset, Buffering.ChrArmLevel.Ratio,
                                "Chromosome Arm")
corr_chr_avg <- dataset_correlation(buf_matched,
                                    Dataset, Buffering.ChrArmLevel.Average.Ratio,
                                    "Chromosome Arm (Average)")

dc_dataset_corr_plot <- corr_gene %>%
  bind_rows(corr_chr, corr_chr_avg) %>%
  rename(`Inter-Dataset Correlation` = "Correlation",
         `Buffering Level` = "Comparison") %>%
  jittered_boxplot(`Buffering Level`, `Inter-Dataset Correlation`, alpha = 0.75, jitter_width = 0.2) %>%
  save_plot("dc_dataset_correlation.png", height = 100)

corr_summary <- list(
  Chr = list(mean = mean(corr_chr$Correlation, na.rm = TRUE),
             median = median(corr_chr$Correlation, na.rm = TRUE),
             sd = sd(corr_chr$Correlation, na.rm = TRUE)),
  ChrAvg = list(mean = mean(corr_chr_avg$Correlation, na.rm = TRUE),
                median = median(corr_chr_avg$Correlation, na.rm = TRUE),
                sd = sd(corr_chr_avg$Correlation, na.rm = TRUE)),
  Gene = list(mean = mean(corr_gene$Correlation, na.rm = TRUE),
              median = median(corr_gene$Correlation, na.rm = TRUE),
              sd = sd(corr_gene$Correlation, na.rm = TRUE))
)

dc_report$corr_procan_depmap_summary <- corr_summary

# === Write Report ===
report_file <- here(reports_dir, "dc_report.txt")
write_report(dc_report, report_file)

# Median CN, Weighted Mean Expr:
#   * Chr:    mean = 0.587, median = 0.599, sd = 0.129
#   * ChrAvg: mean = 0.806, median = 0.811, sd = 0.132
#   * Gene:   mean = 0.524, median = 0.534, sd = 0.119
# Weighted CN & Expr (ChrAvg: unweighted expression)
#   * Chr:    mean = 0.575, median = 0.588, sd = 0.127
#   * ChrAvg: mean = 0.824, median = 0.847, sd = 0.122
#   * Gene:   mean = 0.533, median = 0.542, sd = 0.119
# Weighted CN & Expr (Chr+ChrAvg: Constant CN + CNA)
#   * Chr:    mean = 0.586, median = 0.599, sd = 0.120
#   * ChrAvg: mean = 0.926, median = 0.931, sd = 0.032
#   * Gene:   mean = 0.533, median = 0.542, sd = 0.119
# Unweighted Mean CN & Expr (Chr.CN = Ploidy + CNA):
#   * Chr:    mean = 0.486, median = 0.486, sd = 0.120
#   * ChrAvg: mean = 0.811, median = 0.821, sd = 0.122
#   * Gene:   mean = 0.408, median = 0.411, sd = 0.114
# Unweighted Mean CN & Expr (Chr.CN = Baseline + CNA):
#   * Chr:    mean = 0.494, median = 0.501, sd = 0.106
#   * ChrAvg: mean = 0.883, median = 0.886, sd = 0.042
#   * Gene:   mean = 0.408, median = 0.411, sd = 0.114
# === Weighting Methods ===
# Ploidy Distance:
#   * Chr:    mean = 0.575, median = 0.588, sd = 0.127
#   * ChrAvg: mean = 0.926, median = 0.931, sd = 0.032
#   * Gene:   mean = 0.533, median = 0.542, sd = 0.119
# Aneuploidy Score:
#   * Chr:    mean = 0.574, median = 0.592, sd = 0.132
#   * ChrAvg: mean = 0.926, median = 0.932, sd = 0.032
#   * Gene:   mean = 0.537, median = 0.548, sd = 0.122
# Aneuploidy Score (Chr.CN = round(ploidy) + CNA):
#   * Chr:    mean = 0.573, median = 0.590, sd = 0.124
#   * ChrAvg: mean = 0.926, median = 0.932, sd = 0.032
#   * Gene:   mean = 0.537, median = 0.548, sd = 0.122

# === Combine Plots for publishing ===
br_cn_gain <- plot_buffering_ratio_expr(buffering_ratio, cn_diff = 1) + xlab("Expression Difference (CN Gain)")
br_cn_loss <- plot_buffering_ratio_expr(buffering_ratio, cn_diff = -1) +
  xlab("Expression Difference (CN Loss)") + ylab(NULL) + theme(axis.text.y = element_blank())
br_cn_diff <- plot_buffering_ratio_cn(buffering_ratio, tick_distance = 0.2) + ylab(NULL)

cn_baseline_plot_modified <- cn_baseline_plot +
  theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 2))

plot_br <- cowplot::plot_grid(br_cn_gain, br_cn_loss, br_cn_diff,
                              ncol = 3, labels = c("A", "", "B"), align = "h", axis = "tb", rel_widths = c(1.2, 1, 1.1))

plot_eval <- cowplot::plot_grid(cn_baseline_plot_modified, dc_dataset_corr_plot, ncol = 2, labels = c("C", "D"))

plot_publish <- cowplot::plot_grid(plot_br, plot_eval, rel_heights = c(1,0.8), nrow = 2)

cairo_pdf(here(plots_dir, "dc-method_evaluation_publish.pdf"), width = 11)
plot_publish
dev.off()
