# === Plot Buffering in Genes on Rtr13 and RM13 (P0211) ===
library(here)
library(arrow)
library(dplyr)
library(stringr)
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
  plot_pca()  +
  labs(color = "Cell Line") +
  theme(legend.position = "none")

pca_norm_p0211 <- expr_p0211 %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  plot_pca() +
  labs(color = "Cell Line") +
  theme(legend.position = "none")

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
  mutate(Chromosome = fct_reorder(paste("Chr", Gene.Chromosome), Gene.Chromosome)) %>%
  bidirectional_heatmap(Log2FC, Sample.Name, Chromosome,
                        transpose = TRUE, cluster_rows = TRUE,
                        title = "Average Protein Log2FC per Chromosome", treeheight_row = 25)

## Log2FC by replicate
logfc_heatmap <- expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  bidirectional_heatmap(Log2FC, Gene.Symbol, Sample.Name,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE,
                        title = "Protein Log2FC", treeheight_col = 25, treeheight_row = 25)

## ChrArm Buffering Ratio by replicate
br_heatmap <- expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  mutate(Label = str_replace_all(Buffering.ChrArmLevel.Class, c("Anti-Scaling" = "#", "Buffered" = "=", "Scaling" = ""))) %>%
  bidirectional_heatmap(Buffering.ChrArmLevel.Ratio, Gene.Symbol, Sample.Name,
                        text_col = Label, cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE,
                        title = "Buffering Ratio", treeheight_col = 25, treeheight_row = 25)

## Average Log2FC by Cell Line
logfc_heatmap_avg <- expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  mutate(Label = str_replace_all(Buffering.ChrArmLevel.Average.Class, c("Anti-Scaling" = "#", "Buffered" = "=", "Scaling" = ""))) %>%
  bidirectional_heatmap(Log2FC.Average, Gene.Symbol, CellLine.Name,
                        text_col = Label, cluster_rows = FALSE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE,
                        title = "Average Protein Log2FC", treeheight_col = 25, treeheight_row = 25)

figure_s7_sub1 <- cowplot::plot_grid(chr_heatmap$gtable, pca_norm_p0211,
                                     rel_widths = c(1, 0.5),
                                     nrow = 1, ncol = 2, labels = c("A", "B"))

figure_s7 <- cowplot::plot_grid(figure_s7_sub1, logfc_heatmap$gtable,
                                br_heatmap$gtable, logfc_heatmap_avg$gtable,
                                nrow = 4, ncol = 1,
                                rel_heights = c(1.55, 1, 1, 0.75),
                                labels = c("", "C", "D", "E"))

cairo_pdf(here(plots_dir, "figure_s7.pdf"), width = 12, height = 11)
figure_s7
dev.off()

# === Tables ===
# Remark: Buffering tables exported in 01_buffering_ratio.R
