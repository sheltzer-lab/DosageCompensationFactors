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
source(here("Code", "02_DosageCompensation", "baseline.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "DosageCompensation")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_wgd <- read_parquet(here(output_data_dir, "copy_number_wgd.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_no_wgd <- read_parquet(here(output_data_dir, "copy_number_no-wgd.parquet")) %>%
  select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name)
copy_number_p0211 <- read_parquet(here(output_data_dir, 'copy_number_p0211.parquet'))
copy_number_cptac <- read_parquet((here(output_data_dir, 'copy_number_cptac.parquet')))

expr_procan <- read_parquet(here(output_data_dir, "expression_procan.parquet"))
expr_depmap <- read_parquet(here(output_data_dir, "expression_depmap.parquet"))
expr_matched_renorm <- read_parquet(here(output_data_dir, "expression_matched_renorm.parquet"))
expr_p0211 <- read_parquet(here(output_data_dir, 'expression_p0211.parquet'))
expr_cptac <- read_parquet(here(output_data_dir, 'expression_cptac.parquet'))

# === Combine Datasets and Calculate Buffering & Dosage Compensation ===
# TODO: Rewrite to return vector and allow usage within dplyr::mutate()
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

# ToDo: Evaluate if filtering is neccessary for gene-level dosage compensation analysis
filter_genes <- function(df, gene_col, chr_arm_cna_col, expr_col, min_samples = 10) {
  filtered <- df %>%
    group_by({ { gene_col } }, { { chr_arm_cna_col } }) %>%
    mutate(Samples = sum(!is.na({ { expr_col } }))) %>%
    group_by({ { gene_col } }) %>%
    filter(all(Samples >= min_samples)) %>% # TODO: verify that this works (see diff exp code)
    select(-Samples) %>%
    ungroup()

  return(filtered)
}

baseline_estimation <- function (df) {
  df %>%
    # Note: Chromosome arm CNA based on basal (rounded) ploidy of cell line
    mutate(ChromosomeArm.CopyNumber.Baseline = round(CellLine.Ploidy),
           ChromosomeArm.CopyNumber = round(CellLine.Ploidy) + ChromosomeArm.CNA) %>%
    mutate(PloidyDistance = abs(2 - CellLine.Ploidy),
           CombinedDistance = 0.3 * PloidyDistance / max(PloidyDistance, na.rm = TRUE) +
             0.7 * CellLine.AneuploidyScore / max(CellLine.AneuploidyScore, na.rm = TRUE)) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Gene.CopyNumber,
                       distance_col = CombinedDistance, target_colname = "Gene.CopyNumber.Baseline",
                       weighted = FALSE, summ_func = median) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                       distance_col = CombinedDistance, target_colname = "Protein.Expression.Baseline",
                       weighted = FALSE, summ_func = median) %>%
    calculate_baseline(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                       distance_col = CombinedDistance, target_colname = "Protein.Expression.Baseline.Unweighted",
                       weighted = FALSE, summ_func = mean)
}

calculate_dc <- function(df) {
  df %>%
    group_by(Gene.Symbol, ChromosomeArm.CNA) %>%
    mutate(Protein.Expression.Average = mean(Protein.Expression.Normalized, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Log2FC = Protein.Expression.Normalized - Protein.Expression.Baseline,
           Log2FC.Average = Protein.Expression.Average - Protein.Expression.Baseline.Unweighted) %>%
    # ToDo: Evaluate if building mean positive/negative CNV per gene is interesting
    mutate(Buffering.GeneLevel.Ratio = buffering_ratio(2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                       Gene.CopyNumber.Baseline, Gene.CopyNumber),
           Buffering.GeneLevel.SF = scaling_factor(2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                   Gene.CopyNumber.Baseline, Gene.CopyNumber),
           Buffering.ChrArmLevel.Ratio = buffering_ratio(2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                         ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber),
           Buffering.ChrArmLevel.Average.Ratio = buffering_ratio(2^Protein.Expression.Baseline.Unweighted, 2^Protein.Expression.Average,
                                                                 2L, 2L + ChromosomeArm.CNA),
           Buffering.GeneLevel.Ratio.Confidence = buffering_ratio_confidence(Gene.CopyNumber.Baseline, Gene.CopyNumber,
                                                                             `Protein Neutral CV`),
           Buffering.ChrArmLevel.Ratio.Confidence = buffering_ratio_confidence(ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber,
                                                                               `Protein Neutral CV`)) %>%
    mutate(Buffering.GeneLevel.Class = buffering_class(Buffering.GeneLevel.Ratio,
                                                       2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                       Gene.CopyNumber.Baseline, Gene.CopyNumber),
           Buffering.GeneLevel.SF.Class = buffering_class_sf(Buffering.GeneLevel.SF),
           Buffering.ChrArmLevel.Class = buffering_class(Buffering.ChrArmLevel.Ratio,
                                                         2^Protein.Expression.Baseline, 2^Protein.Expression.Normalized,
                                                         ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber),
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
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, min_samples = 10) %>%
  baseline_estimation() %>%
  calculate_protein_neutral_cv(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
  calculate_dc() %>%
  write_parquet(here(output_data_dir, 'expression_buffering_procan.parquet'), version = "2.6")

expr_buf_depmap <- expr_depmap %>%
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, min_samples = 10) %>%
  baseline_estimation() %>%
  calculate_protein_neutral_cv(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
  calculate_dc() %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap.parquet'), version = "2.6")

expr_buf_matched_renorm <- expr_matched_renorm %>%
  inner_join(y = copy_number, by = c("CellLine.CustomId", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, min_samples = 10) %>%
  baseline_estimation() %>%
  calculate_protein_neutral_cv(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
  calculate_dc() %>%
  write_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'), version = "2.6")

## Whole Genome Doubling
buf_wgd <- expr_depmap %>%
  inner_join(y = copy_number_wgd, by = c("CellLine.CustomId", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, min_samples = 10) %>%
  baseline_estimation() %>%
  calculate_protein_neutral_cv(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
  calculate_dc() %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap_wgd.parquet'), version = "2.6")

buf_no_wgd <- expr_depmap %>%
  inner_join(y = copy_number_no_wgd, by = c("CellLine.CustomId", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, min_samples = 10) %>%
  baseline_estimation() %>%
  calculate_protein_neutral_cv(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
  calculate_dc() %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap_no-wgd.parquet'), version = "2.6")

## P0211
expr_buf_p0211 <- expr_p0211 %>%
  inner_join(y = copy_number_p0211, by = c("Sample.ID", "Gene.Symbol"),
             na_matches = "never", relationship = "many-to-one") %>%
  filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized, min_samples = 2) %>%
  baseline_estimation() %>%
  calculate_protein_neutral_cv(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
  calculate_dc() %>%
  write_parquet(here(output_data_dir, 'expression_buffering_p0211.parquet'), version = "2.6")

## CPTAC
expr_buf_cptac <- expr_cptac %>%
  inner_join(y = copy_number_cptac %>% select(-Model.CancerType, -Protein.Uniprot.Accession, -Gene.Symbol),
             by = c("Model.ID", "Gene.ENSEMBL.Id"), na_matches = "never", relationship = "many-to-one") %>%
  ### Baseline Estimation
  mutate(ChromosomeArm.CNA = NA,
         ChromosomeArm.CopyNumber.Baseline = NA,
         ChromosomeArm.CopyNumber = NA,
         Gene.CopyNumber.Baseline = 2L) %>%
  group_by(Gene.ENSEMBL.Id, Protein.Uniprot.Accession) %>%
  mutate(Protein.Expression.Baseline = median(Protein.Expression.Normalized[Model.SampleType == "Normal" & Gene.CNV == 0],
                                              na.rm = TRUE),
         Protein.Expression.Baseline.Unweighted = NA) %>%
  ungroup() %>%
  # filter(Model.SampleType == "Tumor") %>%
  ### DC Calculation
  calculate_protein_neutral_cv(Gene.ENSEMBL.Id, Gene.CNV, Protein.Expression.Normalized) %>%
  calculate_dc() %>%
  write_parquet(here(output_data_dir, 'expression_buffering_cptac.parquet'), version = "2.6")
