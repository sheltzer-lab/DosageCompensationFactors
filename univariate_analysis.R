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

source(here("parameters.R"))
source(here("buffering_ratio.R"))

output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_avg <- read_parquet(here(output_data_dir, "expression_average.parquet"))
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet"))

# === Calculate Buffering & Dosage Compensation ===

calculate_baseline_copynumber <- function(df, gene_col, copynumber_col) {
  df %>%
    group_by({{ gene_col }}) %>%
    mutate(Gene.CopyNumber.Baseline = median({{ copynumber_col }}, na.rm = TRUE)) %>%
    ungroup()
}

calculate_baseline_expression <- function(df, gene_col, expr_col) {
  df %>%
    group_by({{ gene_col }}) %>%
    mutate(Protein.Expression.Baseline = median({{ expr_col }}, na.rm = TRUE)) %>%
    ungroup()
}

expr_avg_dc <- expr_avg %>%
  inner_join(y = copy_number, by = c("CellLine.SangerModelId", "Gene.Symbol"),
             na_matches = "never", relationship = "one-to-one") %>%
  calculate_baseline_copynumber(Gene.Symbol, Gene.CopyNumber) %>%
  calculate_baseline_expression(Gene.Symbol, Protein.Expression.Normalized) %>%
  mutate(Buffering.Ratio = buffering_ratio(Protein.Expression.Baseline, Protein.Expression.Normalized,
                                           Gene.CopyNumber.Baseline, Gene.CopyNumber)) %>%
  mutate(Buffering.Class = buffering_class(Buffering.Ratio)) %>%
  left_join(y = dc_factors, by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one")

# === Write Processed dataset to disk ===

write_parquet(expr_avg_dc, here(output_data_dir, 'expression_average_dc.parquet'),
              version = "2.6")