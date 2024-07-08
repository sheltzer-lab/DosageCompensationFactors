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

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Gene")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))
expr_buf_p0211 <- read_parquet(here(output_data_dir, "expression_buffering_p0211.parquet"))

# === Determine Genes that are significantly Buffered on average ===

# signif_buf_genes <- function(df, buffering_col, gene_col) {
#   df %>%
#   select({ { gene_col } }, { { buffering_col } }) %>%
#   drop_na() %>%
#   group_by({ { gene_col } }) %>%
#   # Avoid having not enough samples for t-test
#   add_count({ { gene_col } }) %>%
#   filter(n > 1) %>%
#   summarize(TTest.p = t.test({ { buffering_col } }, mu = 0)$p.value,
#             Buffering.Ratio.Average = mean({ { buffering_col } })) %>%
#   mutate(TTest.p.adjusted = p.adjust(TTest.p, method = "BY"))
# }
#
# # ToDo: Repeat analysis with Log2FC (Gene Level & Chromosome Level)
#
# test_all <- expr_buf_procan %>%
#   signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)
#
# plot_all <- test_all %>%
#   mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
#   plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
#   save_plot("volcano_buffering_all.png")
#
# test_gain <- expr_buf_procan %>%
#   filter_cn_gain() %>%
#   signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)
#
# plot_gain <- test_gain %>%
#   mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
#   plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
#   save_plot("volcano_buffering_cn-gain.png")
#
# test_loss <- expr_buf_procan %>%
#   filter_cn_loss() %>%
#   signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)
#
# plot_loss <- test_loss %>%
#   mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
#   plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
#   save_plot("volcano_buffering_cn-loss.png")
#
# test_filtered <- expr_buf_procan %>%
#   filter_cn_diff() %>%
#   signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)
#
# plot_filtered <- test_filtered %>%
#   mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
#   plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
#   save_plot("volcano_buffering_cn-filtered.png")

# === Plot Buffering in Genes on Rtr13 and RM13 (P0211) ===
## ChrArm Buffering Ratio by Sample
expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  bidirectional_heatmap(Buffering.ChrArmLevel.Ratio, Gene.Symbol, Sample.Name,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE) %>%
  save_plot("genes_chr13_buffering_p0211.png", width = 500)

## ChrArm Buffering Ratio by Cell Line
expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  bidirectional_heatmap(Buffering.ChrArmLevel.Ratio, Gene.Symbol, CellLine.Name,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE) %>%
  save_plot("genes_chr13_buffering_average_p0211.png", width = 500)

## Absolute Average Log2FC by Cell Line
expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  mutate(Log2FC.Average.Abs = abs(Log2FC.Average)) %>%
  bidirectional_heatmap(Log2FC.Average.Abs, Gene.Symbol, CellLine.Name,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE) %>%
  save_plot("genes_chr13_log2fc_absolute_average_p0211.png", width = 500)

## Average Log2FC by Cell Line
expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  bidirectional_heatmap(Log2FC.Average, Gene.Symbol, CellLine.Name,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE) %>%
  save_plot("genes_chr13_log2fc_average_p0211.png", width = 500)

## Log2FC by Sample
expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  bidirectional_heatmap(Log2FC, Gene.Symbol, Sample.Name,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE) %>%
  save_plot("genes_chr13_log2fc_p0211.png", width = 500)

# === Export Tables ===

expr_buf_p0211 %>%
  select(Sample.Name, Gene.Symbol, Protein.Expression.Normalized) %>%
  pivot_wider(names_from = "Sample.Name", values_from = "Protein.Expression.Normalized", id_cols = "Gene.Symbol") %>%
  write.xlsx(here(tables_base_dir, "p0211_expression_processed_wide.xlsx"),
             colNames = TRUE)

expr_buf_p0211 %>%
  filter(Gene.Chromosome == 13) %>%
  filter(CellLine.Name != "RPE1") %>%
  select(Gene.Symbol, Gene.ChromosomeArm, Sample.Name, CellLine.Name, CellLine.Replicate,
         ChromosomeArm.CopyNumber, ChromosomeArm.CopyNumber.Baseline,
         Protein.Expression.Normalized, Protein.Expression.Baseline,
         Log2FC, Log2FC.Average,
         Buffering.ChrArmLevel.Ratio, Buffering.ChrArmLevel.Class, Buffering.ChrArmLevel.Ratio.Confidence,
         Buffering.ChrArmLevel.Average.Class) %>%
  write.xlsx(here(tables_base_dir, "p0211_dosage_compensation.xlsx"),
             colNames = TRUE)

chr_arms <- expr_buf_p0211 %>%
  distinct(Gene.Symbol, Gene.Chromosome, Gene.ChromosomeArm)

expr_buf_p0211 %>%
  select(CellLine.Name, Sample.Name, Gene.Symbol, Protein.Expression.Normalized) %>%
  filter(CellLine.Name %in% c("RM13", "RPE1")) %>%
  differential_expression(Gene.Symbol, CellLine.Name, Protein.Expression.Normalized,
                          groups = c("RPE1", "RM13")) %>%
  left_join(y = chr_arms, by = "Gene.Symbol") %>%
  write.xlsx(here(tables_base_dir, "p0211_diff_exp_RM13.xlsx"),
             colNames = TRUE)

expr_buf_p0211 %>%
  select(CellLine.Name, Sample.Name, Gene.Symbol, Protein.Expression.Normalized, Gene.Chromosome) %>%
  filter(CellLine.Name %in% c("RM13", "RPE1")) %>%
  filter(Gene.Chromosome == 13) %>%
  differential_expression(Gene.Symbol, CellLine.Name, Protein.Expression.Normalized,
                          groups = c("RPE1", "RM13")) %>%
  left_join(y = chr_arms, by = "Gene.Symbol") %>%
  write.xlsx(here(tables_base_dir, "p0211_diff_exp_RM13_chr13.xlsx"),
             colNames = TRUE)

expr_buf_p0211 %>%
  select(CellLine.Name, Sample.Name, Gene.Symbol, Protein.Expression.Normalized) %>%
  filter(CellLine.Name %in% c("Rtr13", "RPE1")) %>%
  differential_expression(Gene.Symbol, CellLine.Name, Protein.Expression.Normalized,
                          groups = c("RPE1", "Rtr13")) %>%
  left_join(y = chr_arms, by = "Gene.Symbol") %>%
  write.xlsx(here(tables_base_dir, "p0211_diff_exp_Rtr13.xlsx"),
             colNames = TRUE)

expr_buf_p0211 %>%
  select(CellLine.Name, Sample.Name, Gene.Symbol, Protein.Expression.Normalized, Gene.Chromosome) %>%
  filter(CellLine.Name %in% c("Rtr13", "RPE1")) %>%
  filter(Gene.Chromosome == 13) %>%
  differential_expression(Gene.Symbol, CellLine.Name, Protein.Expression.Normalized,
                          groups = c("RPE1", "Rtr13")) %>%
  left_join(y = chr_arms, by = "Gene.Symbol") %>%
  write.xlsx(here(tables_base_dir, "p0211_diff_exp_Rtr13_chr13.xlsx"),
             colNames = TRUE)

# === Low Variance Buffered Genes ===
mean_var_buf <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  select(Dataset, Gene.Symbol, Buffering.GeneLevel.Ratio, Buffering.GeneLevel.Class) %>%
  drop_na(Buffering.GeneLevel.Ratio) %>%
  group_by(Dataset, Gene.Symbol) %>%
  summarize(Gene.BR.Mean = mean(Buffering.GeneLevel.Ratio),
            Gene.BR.SD = sd(Buffering.GeneLevel.Ratio),
            Gene.Buffered.Count = sum(Buffering.GeneLevel.Class == "Buffered"),
            Samples = sum(!is.na(Buffering.GeneLevel.Class)),
            Gene.Buffered.Share = Gene.Buffered.Count / Samples, .groups = "drop")

buf_rank <- mean_var_buf %>%
  filter(Samples > 20) %>%
  filter(Gene.Buffered.Share > 1/3 & Gene.BR.Mean > br_cutoffs$Buffered) %>%
  add_count(Gene.Symbol) %>%
  filter(n >= 3) %>%
  mean_norm_rank(Gene.BR.SD, Dataset, Gene.Symbol) %>%
  rename(Gene.BR.SD.MeanNormRank = MeanNormRank)

mean_var_buf_highlighted <- mean_var_buf %>%
  left_join(y = buf_rank, by = "Gene.Symbol") %>%
  mutate(Top50 = Gene.Symbol %in% slice_min(buf_rank, Gene.BR.SD.MeanNormRank, n = 50)$Gene.Symbol,
         Top10 = Gene.Symbol %in% slice_min(buf_rank, Gene.BR.SD.MeanNormRank, n = 10)$Gene.Symbol) %>%
  arrange(Gene.BR.SD.MeanNormRank)

scatter_signif_buffered <- mean_var_buf_highlighted %>%
  mutate(Label = if_else(Gene.Symbol %in% buf_labels$Gene.Symbol, Gene.Symbol, NA)) %>%
  filter(Dataset == "CPTAC") %>%
  arrange(Top50) %>%
  ggplot() +
  aes(x = Gene.BR.SD, y = Gene.BR.Mean, label = Label, color = Top50, alpha = Top50) +
  geom_point() +
  geom_label_repel(na.rm = TRUE, max.overlaps = 50) +
  scale_color_manual(values = c(default_color, highlight_color)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(x = "Standard Deviation of Buffering Ratio", y = "Mean Buffering Ratio",
       color = "Frequently Buffered", alpha = "Frequently Buffered")

scatter_signif_buffered %>%
  save_plot("frequently_buffered_genes.png", width = 210, height = 160)

mean_var_buf_highlighted %>%
  drop_na() %>%
  write.xlsx(here(tables_base_dir, "frequently_buffered_genes.xlsx"),
             colNames = TRUE)
