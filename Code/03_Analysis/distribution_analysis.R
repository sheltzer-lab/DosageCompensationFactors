library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(skimr)
library(openxlsx)
library(gtsummary)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Univariate", "ProCan")

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

# === Summarize Distribution of Obersavtions ===

summary_tbl_depmap <- expr_buf_depmap %>%
  tbl_summary(include = c(Gene.CopyNumber, Protein.Expression.Log2, Protein.Expression.Normalized, Protein.Expression.Average,
                          Protein.Expression.Baseline, Protein.Expression.Baseline.Unweighted,
                          Log2FC, Log2FC.Average,
                          Buffering.GeneLevel.Ratio, Buffering.ChrArmLevel.Ratio, Buffering.ChrArmLevel.Average.Ratio,
                          Buffering.GeneLevel.Class, Buffering.ChrArmLevel.Class, Buffering.ChrArmLevel.Average.Class),
              by = ChromosomeArm.CNA) %>%
  modify_header(label = "**Chromosome Arm CNA**") %>%
  bold_labels()

summary_tbl_procan <- expr_buf_procan %>%
  tbl_summary(include = c(Gene.CopyNumber, Protein.Expression.Log2, Protein.Expression.Normalized, Protein.Expression.Average,
                          Protein.Expression.Baseline, Protein.Expression.Baseline.Unweighted,
                          Log2FC, Log2FC.Average,
                          Buffering.GeneLevel.Ratio, Buffering.ChrArmLevel.Ratio, Buffering.ChrArmLevel.Average.Ratio,
                          Buffering.GeneLevel.Class, Buffering.ChrArmLevel.Class, Buffering.ChrArmLevel.Average.Class),
              by = ChromosomeArm.CNA) %>%
  modify_header(label = "**Chromosome Arm CNA**") %>%
  bold_labels()

summary_tbl_merged <- tbl_merge(tbls = list(summary_tbl_depmap, summary_tbl_procan),
                                tab_spanner = c("**DepMap**", "**ProCan**"))
summary_tbl_merged %>%
  as_gt() %>%
  gt::gtsave(filename = here(tables_dir, "summary_table.html"))

## Plot Copy Number Distribution and Filter Thresholds
cn_diff_quantiles <- quantile(expr_buf_procan$Gene.CopyNumber - expr_buf_procan$Gene.CopyNumber.Baseline,
                              probs = seq(0, 1, 0.01))
cn_dist <- expr_buf_procan %>%
  mutate(Log2FC.CopyNumber =  Gene.CopyNumber - Gene.CopyNumber.Baseline) %>%
  ggplot() +
  aes(x = Log2FC.CopyNumber) +
  geom_density() +
  geom_vline(xintercept = cn_diff_quantiles["5%"], linetype = "dashed", color = "red") +
  geom_vline(xintercept = cn_diff_quantiles["95%"], linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(-2, 6, 0.5))
save_plot(cn_dist, "copy_number_distribution.png", width = 300)