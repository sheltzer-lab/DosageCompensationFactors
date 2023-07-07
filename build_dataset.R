library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)
library(gtsummary)


here::i_am("DosageCompensationFactors.Rproj")

source(here("parameters.R"))
source(here("buffering_ratio.R"))

output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet"))

procan_expr_avg <- read_parquet(here(output_data_dir, "expression_average_goncalves.parquet"))
depmap_expr <- read_parquet(here(output_data_dir, "expression_depmap.parquet"))

# === Combine Datasets and Calculate Buffering & Dosage Compensation ===

calculate_baseline_copynumber <- function(df, gene_col, copynumber_col) {
  df %>%
    group_by({{ gene_col }}) %>%
    mutate(Gene.CopyNumber.Baseline = median({{ copynumber_col }}, na.rm = TRUE)) %>%
    ungroup()
}

calculate_baseline_expression <- function(df, gene_col, chr_arm_cna_col, expr_col, aneuploidy_score_col = NULL) {
  baseline_expr <- df %>%
    select({ { gene_col } }, { { chr_arm_cna_col } }, { { expr_col } }, { { aneuploidy_score_col } }) %>%
    filter({ { chr_arm_cna_col } } == 0) %>%
    group_by({ { gene_col } }) %>%
    mutate(Weights = ifelse(is.null({ { aneuploidy_score_col } }),
                            NA,
                            (1 / (1 + { { aneuploidy_score_col } })) /
                              sum(1 / (1 + { { aneuploidy_score_col } }), na.rm = TRUE))) %>%
    summarize(Protein.Expression.Baseline = ifelse(is.null({ { aneuploidy_score_col } }),
                                                   mean({ { expr_col } }, na.rm = TRUE),
                                                   sum({ { expr_col } } * Weights, na.rm = TRUE)),
              Protein.Expression.Baseline.Unweighted = mean({ { expr_col } }, na.rm = TRUE)) %>%
    ungroup()

  df %>%
    inner_join(y = baseline_expr, by = quo_name(enquo(gene_col)),
               unmatched = "error", na_matches = "never", relationship = "many-to-one")
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

build_dataset <- function(df, cellline_col, df_copy_number, df_dc_factors) {
  df %>%
    inner_join(y = copy_number, by = c(quo_name(enquo(cellline_col)), "Gene.Symbol"),
               na_matches = "never", relationship = "many-to-one") %>%
    # ToDo: Evaluate if filtering might be unneccessary for gene-level dosage compensation analysis
    filter_genes(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
    calculate_baseline_copynumber(Gene.Symbol, Gene.CopyNumber) %>%
    calculate_baseline_expression(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized,
                                  aneuploidy_score_col = CellLine.AneuploidyScore) %>%
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
           Buffering.ChrArmLevel.Class = buffering_class_log2fc(Log2FC,
                                                                cn_base = ChromosomeArm.CopyNumber.Baseline,
                                                                cn_var = ChromosomeArm.CopyNumber),
           Buffering.ChrArmLevel.Average.Class = buffering_class_log2fc(Log2FC.Average,
                                                                        cn_base = ChromosomeArm.CopyNumber.Baseline,
                                                                        cn_var = ChromosomeArm.CopyNumber)) %>%
    left_join(y = dc_factors, by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
              na_matches = "never", relationship = "many-to-one")
}


# === Summarize Distribution of Obersavtions ===

summary_tbl_depmap <- expr_buf_depmap %>%
  tbl_summary(include = c(Gene.CopyNumber, Protein.Expression.Normalized, Protein.Expression.Average,
                          Protein.Expression.Baseline, Protein.Expression.Baseline.Unweighted,
                          Log2FC, Log2FC.Average,
                          Buffering.GeneLevel.Class, Buffering.ChrArmLevel.Class,
                          Buffering.ChrArmLevel.Average.Class),
              by = ChromosomeArm.CNA) %>%
  modify_header(label = "**Chromosome Arm CNA**") %>%
  bold_labels()

summary_tbl_goncalves <- expr_buf_goncalves %>%
  tbl_summary(include = c(Gene.CopyNumber, Protein.Expression.Normalized, Protein.Expression.Average,
                          Protein.Expression.Baseline, Protein.Expression.Baseline.Unweighted,
                          Log2FC, Log2FC.Average,
                          Buffering.GeneLevel.Class, Buffering.ChrArmLevel.Class,
                          Buffering.ChrArmLevel.Average.Class),
              by = ChromosomeArm.CNA) %>%
  modify_header(label = "**Chromosome Arm CNA**") %>%
  bold_labels()

summary_tbl_merged <-   tbl_merge(
    tbls = list(summary_tbl_depmap, summary_tbl_goncalves),
    tab_spanner = c("**DepMap**", "**ProCan (Goncalves)**")
  )

summary_tbl_merged

# === Process & Write datasets to disk ===

procan_expr_avg %>%
  build_dataset(CellLine.SangerModelId, copy_number, dc_factors) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_goncalves.parquet'), version = "2.6")

depmap_expr %>%
  build_dataset(CellLine.DepMapModelId, copy_number, dc_factors) %>%
  write_parquet(here(output_data_dir, 'expression_buffering_depmap.parquet'), version = "2.6")