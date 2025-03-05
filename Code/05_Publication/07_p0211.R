# === Plot Buffering in Genes on Rtr13 and RM13 (P0211) ===
library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir
tables_dir <- here(tables_base_dir, "Publication")

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

# === Load Datasets ===

expr_p0211 <- read_parquet(here(output_data_dir, 'expression_p0211.parquet'))
expr_buf_p0211 <- read_parquet(here(output_data_dir, "expression_buffering_p0211.parquet"))

copy_number_p0211 <- read_parquet(here(output_data_dir, 'copy_number_p0211.parquet'))

# === Supplemental Figures ===

## PCA plots
pca_pre_p0211 <- expr_p0211 %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Log2) %>%
  plot_pca()

pca_norm_p0211 <- expr_p0211 %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  plot_pca()

## Protein Log2FC per chromosome
p0211_expr_baseline <- expr_p0211 %>%
  filter(CellLine.Name == "RPE1") %>%
  group_by(Gene.Symbol) %>%
  summarize(Protein.Expression.Baseline = median(Protein.Expression.Normalized, na.rm = TRUE))

chr_heatmap <- expr_p0211 %>%
  left_join(y = p0211_expr_baseline, by = "Gene.Symbol",
            unmatched = "error", na_matches = "never", relationship = "many-to-one") %>%
  left_join(y = copy_number_p0211, by = c("Gene.Symbol", "Sample.ID"),
            unmatched = "error", na_matches = "never", relationship = "many-to-one") %>%
  summarize(Log2FC = Protein.Expression.Normalized - Protein.Expression.Baseline,
            .by = c(Gene.Symbol, Sample.Name, Sample.ID, CellLine.Name, Gene.Chromosome, Gene.StartPosition)) %>%
  bidirectional_heatmap(Log2FC, Sample.Name, Gene.Chromosome,
                        transpose = TRUE, cluster_rows = TRUE)

## Log2FC by replicate
logfc_heatmap <- expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  bidirectional_heatmap(Log2FC, Gene.Symbol, Sample.Name,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE)

## ChrArm Buffering Ratio by replicate
br_heatmap <- expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  mutate(Label = str_replace_all(Buffering.ChrArmLevel.Class, c("Anti-Scaling" = "-", "Buffered" = "*", "Scaling" = ""))) %>%
  bidirectional_heatmap(Buffering.ChrArmLevel.Ratio, Gene.Symbol, Sample.Name,
                        text_col = Label, cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE)

## Average Log2FC by Cell Line
logfc_heatmap_avg <- expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  mutate(Label = str_replace_all(Buffering.ChrArmLevel.Average.Class, c("Anti-Scaling" = "-", "Buffered" = "*", "Scaling" = ""))) %>%
  bidirectional_heatmap(Log2FC.Average, Gene.Symbol, CellLine.Name,
                        text_col = Label, cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE)

# TODO: Export table

