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
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))

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

## Buffering Ratio Gain vs. Loss
### P0211
expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  signif_violin_plot(CellLine.Name, Buffering.ChrArmLevel.Ratio) %>%
  save_plot("buffering_cnv_p0211.png", width = 100, height = 150)

### CPTAC
expr_buf_cptac %>%
  filter(Gene.CopyNumber.Baseline == 2) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  signif_violin_plot(CNV, Buffering.GeneLevel.Ratio) %>%
  save_plot("buffering_cnv_cptac.png", width = 100, height = 150)

### DepMap
expr_buf_depmap %>%
  filter(Gene.CopyNumber.Baseline == 2) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  signif_violin_plot(CNV, Buffering.GeneLevel.Ratio) %>%
  save_plot("buffering_cnv_depmap.png", width = 100, height = 150)

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
low_var_buf <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  analyze_low_br_variance()

low_var_buf_gain <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  filter_cn_gain_abs() %>%
  analyze_low_br_variance()

low_var_buf_loss <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  filter_cn_loss_abs() %>%
  analyze_low_br_variance()

scatter_signif_buffered <- low_var_buf %>%
  mutate(Label = if_else(Top10, Gene.Symbol, NA)) %>%
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

scatter_signif_buffered_cna <- bind_rows(low_var_buf_gain %>% mutate(Gene.CNA = "Gain"),
                                         low_var_buf_loss %>% mutate(Gene.CNA = "Loss")) %>%
  mutate(Label = if_else(Top10, Gene.Symbol, NA)) %>%
  filter(Dataset == "CPTAC") %>%
  arrange(Top50) %>%
  ggplot() +
  aes(x = Gene.BR.SD, y = Gene.BR.Mean, label = Label, color = Top50, alpha = Top50) +
  geom_point() +
  geom_label_repel(na.rm = TRUE, max.overlaps = 50) +
  scale_color_manual(values = c(default_color, highlight_color)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(x = "Standard Deviation of Buffering Ratio", y = "Mean Buffering Ratio",
       color = "Frequently Buffered", alpha = "Frequently Buffered") +
  facet_grid(~Gene.CNA)

scatter_signif_buffered_cna %>%
  save_plot("frequently_buffered_genes_by-cna.png", width = 280, height = 160)

low_var_buf %>%
  drop_na() %>%
  write.xlsx(here(tables_base_dir, "frequently_buffered_genes.xlsx"),
             colNames = TRUE)

low_var_buf_gain %>%
  drop_na() %>%
  write.xlsx(here(tables_base_dir, "frequently_buffered_genes_gain.xlsx"),
             colNames = TRUE)

low_var_buf_loss %>%
  drop_na() %>%
  write.xlsx(here(tables_base_dir, "frequently_buffered_genes_loss.xlsx"),
             colNames = TRUE)

## VCP occurs consistently (gain/loss/all)
### Potential target for cancer therapy drugs (Eerl, CB-5083)
### VCP inhibitors synergize with proteasome inhibitors (Bortezomib)

## Overrepresentation Analysis of Top50
ora_buf <- low_var_buf %>%
  filter(Top50) %>%
  filter(Dataset == "CPTAC") %>%
  arrange(Gene.BR.SD.MeanNormRank) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = TRUE)

ora_buf %>%
  plot_terms() %>%
  save_plot("frequently_buffered_genes_ora-terms.png", height = 200)

# === Random Allelic Expression ===
rae_buf <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  add_factors(df_factors = dc_factors) %>%
  select(Model.ID, Protein.Uniprot.Accession, Gene.Symbol, Dataset,
         Gene.CopyNumber, Gene.CopyNumber.Baseline, Protein.Expression.Normalized,
         Buffering.GeneLevel.Ratio, Buffering.GeneLevel.Class, `Random Allelic Expression`) %>%
  drop_na() %>%
  mutate(`Random Allelic Expression` = `Random Allelic Expression` > 0)

rae_buf_plot <- rae_buf %>%
  filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  ggplot() +
  aes(x = `Random Allelic Expression`, y = Buffering.GeneLevel.Ratio, color = `Random Allelic Expression`) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")), test = t.test,
              map_signif_level = print_signif, y_position = 1.5,
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c("FALSE" = default_color, "TRUE" = highlight_color), guide = "none") +
  #scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1)) +
  facet_grid(~Dataset)

rae_buf_plot %>%
  save_plot("random_allelic_expression.png")

rae_buf_plot_gain <- rae_buf %>%
  filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
  filter_cn_gain_abs() %>%
  ggplot() +
  aes(x = `Random Allelic Expression`, y = Buffering.GeneLevel.Ratio, color = `Random Allelic Expression`) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")), test = t.test,
              map_signif_level = print_signif, y_position = 1.5,
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c("FALSE" = default_color, "TRUE" = highlight_color), guide = "none") +
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1)) +
  facet_grid(~Dataset)

rae_buf_plot_gain %>%
  save_plot("random_allelic_expression_cn-gain.png")

rae_buf_plot_loss <- rae_buf %>%
  filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
  filter_cn_loss_abs() %>%
  ggplot() +
  aes(x = `Random Allelic Expression`, y = Buffering.GeneLevel.Ratio, color = `Random Allelic Expression`) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("FALSE", "TRUE")), test = t.test,
              map_signif_level = print_signif, y_position = 1.5,
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c("FALSE" = default_color, "TRUE" = highlight_color), guide = "none") +
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1)) +
  facet_grid(~Dataset)

rae_buf_plot_loss %>%
  save_plot("random_allelic_expression_cn-loss.png")

rae_buf_plot_cnv <- rae_buf %>%
  filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  filter(`Random Allelic Expression`) %>%
  mutate(Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  ggplot() +
  aes(x = Event, y = Buffering.GeneLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("Gain", "Loss")), test = t.test,
              map_signif_level = print_signif, y_position = 1.5,
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  #scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1)) +
  facet_grid(~Dataset)

rae_buf_plot_cnv %>%
  save_plot("random_allelic_expression_gene-cnv.png")

## Check if RAE leads to change in Buffered proteins
rae_buf_counts <- rae_buf %>%
  filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(Buffering = factor(if_else(Buffering.GeneLevel.Class == "Buffered", "Buffering", "Other"),
                            levels = c("Buffering", "Other")),
         RAE = factor(if_else(`Random Allelic Expression`, "RAE", "Other"),
                            levels = c("RAE", "Other"))) %>%
  count(Buffering, RAE) %>%
  drop_na() %>%
  arrange(Buffering, RAE) %>%
  pivot_wider(names_from = "RAE", values_from = "n") %>%
  tibble::column_to_rownames("Buffering") %>%
  as.matrix()

rae_buf_test <- fisher.test(rae_buf_counts)
