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

output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_wgd <- read_parquet(here(output_data_dir, "copy_number_wgd.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_no_wgd <- read_parquet(here(output_data_dir, "copy_number_no-wgd.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)

expr_procan <- read_parquet(here(output_data_dir, "expression_procan.parquet"))
expr_depmap <- read_parquet(here(output_data_dir, "expression_depmap.parquet"))

expr_combined <- read_parquet(here(output_data_dir, "expression_combined.parquet"))
expr_combined_celllines <- read_parquet(here(output_data_dir, "expression_combined_celllines.parquet"))
expr_combined_genes <- read_parquet(here(output_data_dir, "expression_combined_genes.parquet"))

expr_matched <- read_parquet(here(output_data_dir, "expression_matched.parquet"))
expr_matched_renorm <- read_parquet(here(output_data_dir, "expression_matched_renorm.parquet"))

# === Combine Datasets and Calculate Buffering & Dosage Compensation ===

calculate_weights <- function(base, var) {
  distances <- abs(base - var)
  weights <- (1/(1+distances)) / sum(1/(1+distances), na.rm = TRUE)
  return(weights)
}
calculate_baseline <- function(df, gene_col, chr_arm_cna_col, value_col,
                               target_colname = "Baseline", ploidy_col = NULL, weighted = TRUE) {
  baseline_expr <- df %>%
    select({ { gene_col } }, { { chr_arm_cna_col } }, { { value_col } }, { { ploidy_col } }) %>%
    filter({ { chr_arm_cna_col } } == 0) %>%
    group_by({ { gene_col } }) %>%
    mutate(Weights = calculate_weights(2, { { ploidy_col } })) %>%
    summarize(!!target_colname := if_else(weighted,
                                          sum({ { value_col } } * Weights, na.rm = TRUE),
                                          mean({ { value_col } }, na.rm = TRUE))) %>%
    ungroup()

  df %>%
    inner_join(y = baseline_expr, by = quo_name(enquo(gene_col)),
               unmatched = "error", na_matches = "never", relationship = "many-to-one") %>%
    rename()
}

calculate_protein_neutral_cv <- function(df, gene_col, chr_arm_cna_col, expr_col) {
  variance <- df %>%
    select({ { gene_col } }, { { chr_arm_cna_col } }, { { expr_col } }) %>%
    filter({ { chr_arm_cna_col } } == 0) %>%
    group_by({ { gene_col } }) %>%
    summarize(`Protein Neutral CV` = log2(sd({ { expr_col } }, na.rm = TRUE) / abs(mean({ { expr_col } }, na.rm = TRUE)))) %>%
    ungroup()

  df %>%
    inner_join(y = variance, by = quo_name(enquo(gene_col)),
               unmatched = "error", na_matches = "never", relationship = "many-to-one")
}

filter_genes <- function(df, gene_col, chr_arm_cna_col, expr_col) {
  filtered <- df %>%
    group_by({ { gene_col } }, { { chr_arm_cna_col } }) %>%
    mutate(Samples = sum(!is.na({ { expr_col } }))) %>%
    group_by({ { gene_col } }) %>%
    filter(all(Samples > 10)) %>%
    select(-Samples) %>%
    ungroup()

  return(filtered)
}

build_dataset <- function(df, df_copy_number, cellline_col = "CellLine.CustomId") {
  test <- df %>%
    inner_join(y = df_copy_number, by = c(cellline_col, "Gene.Symbol"),
               na_matches = "never", relationship = "many-to-one") %>%
    # ToDo: Evaluate if filtering might be unneccessary for gene-level dosage compensation analysis
    filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, CellLine.Ploidy,
                       ploidy_col = CellLine.Ploidy, target_colname = "ChromosomeArm.CopyNumber.Baseline",
                       weighted = TRUE) %>%
    # Note: Chromosome arm CNA based on ploidy of cell line
    mutate(ChromosomeArm.CopyNumber = CellLine.Ploidy + ChromosomeArm.CNA) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                       ploidy_col = CellLine.Ploidy, target_colname = "Gene.CopyNumber.Baseline",
                       weighted = TRUE) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                       ploidy_col = CellLine.Ploidy, target_colname = "Protein.Expression.Baseline",
                       weighted = TRUE) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                       ploidy_col = CellLine.Ploidy, target_colname = "Protein.Expression.Baseline.Unweighted",
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
           Buffering.ChrArmLevel.Ratio = buffering_ratio(2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                         ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber),
           Buffering.ChrArmLevel.Average.Ratio = buffering_ratio(2^Protein.Expression.Baseline.Unweighted, 2^Protein.Expression.Average,
                                                                 ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber)) %>%
    mutate(Buffering.GeneLevel.Class = buffering_class(Buffering.GeneLevel.Ratio),
           Buffering.ChrArmLevel.Class = buffering_class(Buffering.ChrArmLevel.Ratio),
           Buffering.ChrArmLevel.Log2FC.Class = buffering_class_log2fc(Log2FC,
                                                                       cn_base = 2L,
                                                                       cn_var = 2L + ChromosomeArm.CNA),
           Buffering.ChrArmLevel.Average.Class = buffering_class_log2fc(Log2FC.Average,
                                                                        cn_base = 2L,
                                                                        cn_var = 2L + ChromosomeArm.CNA))
}

# === Process & Write datasets to disk ===

# Note: DepMap copy number data does not cover all cell lines in ProCan (333 cell lines lost here)
expr_procan %>%
  build_dataset(copy_number) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_procan.parquet'), version = "2.6")

expr_depmap %>%
  build_dataset(copy_number) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap.parquet'), version = "2.6")

# expr_combined %>%
#   build_dataset(copy_number) %>%
#   write_parquet(here(output_data_dir, 'expression_buffering_combined.parquet'), version = "2.6")
#
# expr_combined_celllines %>%
#   build_dataset(copy_number) %>%
#   write_parquet(here(output_data_dir, 'expression_buffering_combined_celllines.parquet'), version = "2.6")
#
# expr_combined_genes %>%
#   build_dataset(copy_number) %>%
#   write_parquet(here(output_data_dir, 'expression_buffering_combined_genes.parquet'), version = "2.6")
#
# expr_matched %>%
#   build_dataset(copy_number) %>%
#   write_parquet(here(output_data_dir, 'expression_buffering_matched.parquet'), version = "2.6")

expr_matched_renorm %>%
  build_dataset(copy_number) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'), version = "2.6")

## Whole Genome Doubling
expr_depmap %>%
  build_dataset(copy_number_wgd) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap_wgd.parquet'), version = "2.6")

expr_depmap %>%
  build_dataset(copy_number_no_wgd) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap_no-wgd.parquet'), version = "2.6")

# === Evaluation ===
## Copy Number
df_cn_eval <- expr_depmap %>%
  select(Gene.Symbol, CellLine.CustomId) %>%
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
               na_matches = "never", relationship = "many-to-one") %>%
  select(Gene.Symbol, CellLine.CustomId, Gene.CopyNumber,
         ChromosomeArm.CNA, CellLine.Ploidy) %>%
  group_by(Gene.Symbol) %>%
  mutate(MedianAll = median(Gene.CopyNumber, na.rm = TRUE)) %>%
  ungroup() %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     ploidy_col = CellLine.Ploidy, target_colname = "WeightedNeutral", weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                     ploidy_col = CellLine.Ploidy, target_colname = "MeanNeutral", weighted = FALSE) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol, MedianAll, WeightedNeutral, MeanNeutral)

cor.test(df_cn_eval$WeightedNeutral, df_cn_eval$MedianAll, method = "spearman")

cn_baseline_plot <- df_cn_eval %>%
  pivot_longer(c("MedianAll", "WeightedNeutral", "MeanNeutral"),
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
  select(Gene.Symbol, CellLine.CustomId, ChromosomeArm.CNA,
         CellLine.Ploidy, Protein.Expression.Normalized) %>%
  group_by(Gene.Symbol) %>%
  mutate(MedianAll = median(Protein.Expression.Normalized, na.rm = TRUE)) %>%
  ungroup() %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     ploidy_col = CellLine.Ploidy, target_colname = "WeightedNeutral", weighted = TRUE) %>%
  calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                     ploidy_col = CellLine.Ploidy, target_colname = "MeanNeutral", weighted = FALSE) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol, MedianAll, WeightedNeutral, MeanNeutral)

cor.test(df_expr_eval$WeightedNeutral, df_expr_eval$MeanNeutral, method = "spearman")

expr_baseline_plot <- df_expr_eval %>%
  pivot_longer(c("MedianAll", "WeightedNeutral", "MeanNeutral"),
               names_to = "Methods",
               values_to = "Protein.Expression.Baseline") %>%
  ggplot() +
  aes(color = Methods, x = Protein.Expression.Baseline) +
  geom_density()

expr_baseline_plot %>%
  save_plot("expression_baseline_methods.png")