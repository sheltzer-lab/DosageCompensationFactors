library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)


here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "analysis.R"))
source(here("Code", "evaluation.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

dc_report <- list()

# === Load Datasets ===
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_wgd <- read_parquet(here(output_data_dir, "copy_number_wgd.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_no_wgd <- read_parquet(here(output_data_dir, "copy_number_no-wgd.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_p0211 <- read_parquet(here(output_data_dir, 'copy_number_p0211.parquet'))

expr_procan <- read_parquet(here(output_data_dir, "expression_procan.parquet"))
expr_depmap <- read_parquet(here(output_data_dir, "expression_depmap.parquet"))
expr_matched_renorm <- read_parquet(here(output_data_dir, "expression_matched_renorm.parquet"))
expr_p0211 <- read_parquet(here(output_data_dir, 'expression_p0211.parquet'))

dc_factors_ids <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet")) %>%
  select(Protein.Uniprot.Accession, Gene.Symbol)

# === Combine Datasets and Calculate Buffering & Dosage Compensation ===

calculate_weights <- function(distances) {
  weights <- (1/(1+distances)) / sum(1/(1+distances), na.rm = TRUE)
  return(weights)
}
calculate_baseline <- function(df, gene_col, chr_arm_cna_col, value_col,
                               target_colname = "Baseline", distance_col = NULL, weighted = FALSE) {
  baseline_expr <- df %>%
    select({ { gene_col } }, { { chr_arm_cna_col } }, { { value_col } }, { { distance_col } }) %>%
    filter({ { chr_arm_cna_col } } == 0) %>%
    group_by({ { gene_col } }) %>%
    mutate(Weights = calculate_weights({ { distance_col } })) %>%
    summarize(!!target_colname := if_else(weighted,
                                          sum({ { value_col } } * Weights, na.rm = TRUE),
                                          mean({ { value_col } }, na.rm = TRUE))) %>%
    ungroup()

  df %>%
    inner_join(y = baseline_expr, by = quo_name(enquo(gene_col)),
               unmatched = "error", na_matches = "never", relationship = "many-to-one")
}

calculate_median_baseline <- function(df, gene_col, value_col, target_colname = "Baseline") {
  df %>%
    group_by({ { gene_col } }) %>%
    mutate(!!target_colname := median({ { value_col } }, na.rm = TRUE)) %>%
    ungroup()
}

calculate_protein_neutral_cv <- function(df, gene_col, chr_arm_cna_col, expr_col) {
  variance <- df %>%
    select({ { gene_col } }, { { chr_arm_cna_col } }, { { expr_col } }) %>%
    filter({ { chr_arm_cna_col } } == 0) %>%
    group_by({ { gene_col } }) %>%
    summarize(`Protein Neutral CV` = sd(2^{ { expr_col } }, na.rm = TRUE) / mean(2^{ { expr_col } }, na.rm = TRUE)) %>%
    ungroup()

  df %>%
    inner_join(y = variance, by = quo_name(enquo(gene_col)),
               unmatched = "error", na_matches = "never", relationship = "many-to-one")
}

filter_genes <- function(df, gene_col, chr_arm_cna_col, expr_col, min_samples = 10) {
  filtered <- df %>%
    group_by({ { gene_col } }, { { chr_arm_cna_col } }) %>%
    mutate(Samples = sum(!is.na({ { expr_col } }))) %>%
    group_by({ { gene_col } }) %>%
    filter(all(Samples >= min_samples)) %>%
    select(-Samples) %>%
    ungroup()

  return(filtered)
}

build_dataset <- function(df, df_copy_number, cellline_col = "CellLine.CustomId", min_samples = 10) {
  test <- df %>%
    inner_join(y = df_copy_number, by = c(cellline_col, "Gene.Symbol"),
               na_matches = "never", relationship = "many-to-one") %>%
    # ToDo: Evaluate if filtering might be unneccessary for gene-level dosage compensation analysis
    filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, min_samples = min_samples) %>%
    # Note: Chromosome arm CNA based on basal (rounded) ploidy of cell line
    mutate(ChromosomeArm.CopyNumber.Baseline = round(CellLine.Ploidy),
           ChromosomeArm.CopyNumber = round(CellLine.Ploidy) + ChromosomeArm.CNA) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                       distance_col = CellLine.AneuploidyScore, target_colname = "Gene.CopyNumber.Baseline",
                       weighted = TRUE) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                       distance_col = CellLine.AneuploidyScore, target_colname = "Protein.Expression.Baseline",
                       weighted = TRUE) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                       distance_col = CellLine.AneuploidyScore, target_colname = "Protein.Expression.Baseline.Unweighted",
                       weighted = FALSE) %>%
    calculate_protein_neutral_cv(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
    group_by(Gene.Symbol, ChromosomeArm.CNA) %>%
    mutate(Protein.Expression.Average = mean(Protein.Expression.Normalized, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Log2FC = Protein.Expression.Normalized - Protein.Expression.Baseline,
           Log2FC.Average = Protein.Expression.Average - Protein.Expression.Baseline.Unweighted) %>%
    # ToDo: Evaluate if building mean positive/negative CNV per gene is interesting
    mutate(Buffering.GeneLevel.Ratio = buffering_ratio(2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                       2^Gene.CopyNumber.Baseline, 2^Gene.CopyNumber),
           Buffering.GeneLevel.SF = scaling_factor(2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                   2^Gene.CopyNumber.Baseline, 2^Gene.CopyNumber),
           Buffering.ChrArmLevel.Ratio = buffering_ratio(2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                         ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber),
           Buffering.ChrArmLevel.Average.Ratio = buffering_ratio(2^Protein.Expression.Baseline.Unweighted, 2^Protein.Expression.Average,
                                                                 2L, 2L + ChromosomeArm.CNA)) %>%
    mutate(Buffering.GeneLevel.Class = buffering_class(Buffering.GeneLevel.Ratio),
           Buffering.GeneLevel.SF.Class = buffering_class_sf(Buffering.GeneLevel.SF),
           Buffering.ChrArmLevel.Class = buffering_class(Buffering.ChrArmLevel.Ratio),
           Buffering.ChrArmLevel.Log2FC.Class = buffering_class_log2fc(Log2FC,
                                                                       cn_base = ChromosomeArm.CopyNumber.Baseline,
                                                                       cn_var = ChromosomeArm.CopyNumber),
           Buffering.ChrArmLevel.Average.Class = buffering_class_log2fc(Log2FC.Average,
                                                                        cn_base = 2L,
                                                                        cn_var = 2L + ChromosomeArm.CNA))
}

# === Process & Write datasets to disk ===

# Note: DepMap copy number data does not cover all cell lines in ProCan (333 cell lines lost here)
expr_buf_procan <- expr_procan %>%
  build_dataset(copy_number) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_procan.parquet'), version = "2.6")

expr_buf_depmap <- expr_depmap %>%
  build_dataset(copy_number) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap.parquet'), version = "2.6")

expr_buf_matched_renorm <- expr_matched_renorm %>%
  build_dataset(copy_number) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'), version = "2.6")

## Whole Genome Doubling
buf_wgd <- expr_depmap %>%
  build_dataset(copy_number_wgd) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap_wgd.parquet'), version = "2.6")

buf_no_wgd <- expr_depmap %>%
  build_dataset(copy_number_no_wgd) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap_no-wgd.parquet'), version = "2.6")

## P0211
expr_buf_p0211 <- expr_p0211 %>%
  build_dataset(copy_number_p0211, cellline_col = "Sample.ID", min_samples = 2) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_p0211.parquet'), version = "2.6")

# === Evaluation ===
## Copy Number
df_cn_eval <- expr_depmap %>%
  select(Gene.Symbol, CellLine.CustomId) %>%
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
               na_matches = "never", relationship = "many-to-one") %>%
  select(Gene.Symbol, CellLine.CustomId, Gene.CopyNumber, CellLine.Ploidy,
         ChromosomeArm.CNA, CellLine.AneuploidyScore) %>%
  calculate_median_baseline(Gene.Symbol, Gene.CopyNumber, target_colname = "MedianAll") %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     distance_col = CellLine.AneuploidyScore, target_colname = "WeightedNeutral.AneuploidyScore",
                     weighted = TRUE) %>%
  mutate(PloidyDistance = abs(2 - CellLine.Ploidy)) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     distance_col = PloidyDistance, target_colname = "WeightedNeutral.PloidyDistance",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     distance_col = CellLine.AneuploidyScore, target_colname = "MeanNeutral", weighted = FALSE) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol, MedianAll, WeightedNeutral.AneuploidyScore, WeightedNeutral.PloidyDistance, MeanNeutral)

dc_report$corr_CNmethod_weightedAS_MedianAll <- cor.test(df_cn_eval$WeightedNeutral.AneuploidyScore,
                                                         df_cn_eval$MedianAll, method = "spearman")

cn_baseline_plot <- df_cn_eval %>%
  pivot_longer(c("MedianAll", "WeightedNeutral.AneuploidyScore", "WeightedNeutral.PloidyDistance", "MeanNeutral"),
               names_to = "Methods",
               values_to = "Gene.CopyNumber.Baseline") %>%
  ggplot() +
  aes(color = Methods, x = Gene.CopyNumber.Baseline) +
  geom_density()

cn_baseline_plot %>%
  save_plot("copynumber_baseline_methods.png")

## Expression
df_expr_eval <- expr_depmap %>%
  select(Gene.Symbol, CellLine.CustomId, Protein.Expression.Normalized) %>%
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
               na_matches = "never", relationship = "many-to-one") %>%
  select(Gene.Symbol, CellLine.CustomId, ChromosomeArm.CNA, CellLine.Ploidy,
         CellLine.AneuploidyScore, Protein.Expression.Normalized) %>%
  calculate_median_baseline(Gene.Symbol, Protein.Expression.Normalized, target_colname = "MedianAll") %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     distance_col = CellLine.AneuploidyScore, target_colname = "WeightedNeutral.AneuploidyScore",
                     weighted = TRUE) %>%
  mutate(PloidyDistance = abs(2 - CellLine.Ploidy)) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     distance_col = PloidyDistance, target_colname = "WeightedNeutral.PloidyDistance",
                     weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     distance_col = CellLine.AneuploidyScore, target_colname = "MeanNeutral", weighted = FALSE) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol, MedianAll, WeightedNeutral.AneuploidyScore, WeightedNeutral.PloidyDistance, MeanNeutral)

dc_report$corr_ExprMethod_weightedAS_MeanNeutral <- cor.test(df_expr_eval$WeightedNeutral.AneuploidyScore,
                                                             df_expr_eval$MeanNeutral, method = "spearman")

expr_baseline_plot <- df_expr_eval %>%
  pivot_longer(c("MedianAll", "WeightedNeutral.AneuploidyScore", "WeightedNeutral.PloidyDistance", "MeanNeutral"),
               names_to = "Methods",
               values_to = "Protein.Expression.Baseline") %>%
  ggplot() +
  aes(color = Methods, x = Protein.Expression.Baseline) +
  geom_density()

expr_baseline_plot %>%
  save_plot("expression_baseline_methods.png")

## Check ratio of data lost after join with factor data
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

## Check correlation between DC Scores and Protein Coefficient of Variance
## Buffering Score should be negatively correlated with Protein Coefficient of Variance
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

corr_gene %>%
  bind_rows(corr_chr, corr_chr_avg) %>%
  jittered_boxplot(Comparison, Correlation) %>%
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
write_list(dc_report, report_file)

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
