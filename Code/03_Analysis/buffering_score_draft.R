library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)
library(pROC)
library(skimr)
library(cowplot)
library(broom)
library(corrplot)
library(gtsummary)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "03_Analysis", "analysis.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Univariate")
goncalves_plots_dir <- here(plots_base_dir, "Univariate", "Goncalves")
goncalves_chr_plots_dir <- here(goncalves_plots_dir, "ChromosomeArm-Level")
goncalves_gene_plots_dir <- here(goncalves_plots_dir, "Gene-Level")
goncalves_comparison_plots_dir <- here(goncalves_plots_dir, "Comparison")
depmap_plots_dir <- here(plots_base_dir, "Univariate", "DepMap")
depmap_chr_plots_dir <- here(depmap_plots_dir, "ChromosomeArm-Level")
depmap_gene_plots_dir <- here(depmap_plots_dir, "Gene-Level")

dir.create(output_data_dir, recursive = TRUE)
dir.create(goncalves_plots_dir, recursive = TRUE)
dir.create(goncalves_chr_plots_dir, recursive = TRUE)
dir.create(goncalves_gene_plots_dir, recursive = TRUE)
dir.create(goncalves_comparison_plots_dir, recursive = TRUE)
dir.create(depmap_plots_dir, recursive = TRUE)
dir.create(depmap_chr_plots_dir, recursive = TRUE)
dir.create(depmap_gene_plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))


test <- expr_buf %>%
  select(Gene.Symbol, CellLine.Name,
         Gene.CopyNumber.Baseline, Gene.CopyNumber,
         Protein.Expression.Baseline, Protein.Expression.Normalized,
         Buffering.ChrArmLevel.Average.Ratio, Buffering.ChrArmLevel.Average.Class,
         Buffering.GeneLevel.Ratio, Buffering.GeneLevel.Class,
         ChromosomeArm.CNA, Log2FC.Average) %>%
  filter(ChromosomeArm.CNA != 0) %>%
  group_by(Gene.Symbol, ChromosomeArm.CNA) %>%
  mutate(Expression.Log2FC = Protein.Expression.Normalized - Protein.Expression.Baseline,
         CopyNumber.Log2FC = Gene.CopyNumber - Gene.CopyNumber.Baseline) %>%
  mutate(Buffering.Score = buffering_score(Expression.Log2FC, CopyNumber.Log2FC)) %>%
  mutate(Buffering.Score.Class = buffering_score_class(Buffering.Score)) %>%
  ungroup() %>%
  distinct(Gene.Symbol, ChromosomeArm.CNA, Buffering.Score, Buffering.Score.Class,
           Buffering.ChrArmLevel.Average.Class, Log2FC.Average)

summary_tbl <- test %>%
  tbl_summary(include = c(Log2FC.Average, Buffering.Score,
                          Buffering.Score.Class, Buffering.ChrArmLevel.Average.Class),
              by = ChromosomeArm.CNA) %>%
  modify_header(label = "**Chromosome Arm CNA**") %>%
  bold_labels()

calculate_buffering <- function(df, method = "corr", level = "gene", condition = "gain", average = TRUE) {
  if (method == "corr" &
    level == "gene" &
    condition == "gain") {
    df_proc <- df %>%
      filter_cn_gain() %>%
      group_by(Gene.Symbol) %>%
      mutate(Buffering.Score = buffering_score(Protein.Expression.Normalized - Protein.Expression.Baseline,
                                               Gene.CopyNumber - Gene.CopyNumber.Baseline)) %>%
      ungroup() %>%

      mutate(Buffering.Score.Class = buffering_score_class(Buffering.Score),
             UniqueId = Gene.Symbol) %>%
      distinct(Gene.Symbol, .keep_all = TRUE)
    return(df_proc)
  }
}

df_buf_score <- expr_buf %>%
  calculate_buffering()

buf_score_summary <- df_buf_score %>%
  tbl_summary(include = c(Buffering.Score, Buffering.Score.Class)) %>%
  bold_labels()