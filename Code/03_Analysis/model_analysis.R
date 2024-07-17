library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(skimr)
library(openxlsx)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "CellLine")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
df_celllines <- read_parquet(here(output_data_dir, "celllines.parquet"))
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

# === Analyze Dosage Compensation on Cell Line level ===

# TODO: Calculate confidence score based on SD, Observations & Gene-Level Confidence
analyze_model_buffering <- function(df, buffering_ratio_col, model_col = Model.ID,
                                    aneuploidy_col = CellLine.AneuploidyScore) {
  cellline_buf_avg <- df %>%
    select({ { model_col } }, { { buffering_ratio_col } }, { { aneuploidy_col } }) %>%
    group_by({ { model_col } }) %>%
    summarize(Model.Buffering.Ratio = mean({ { buffering_ratio_col } }, na.rm = TRUE),
              Observations = sum(!is.na({ { buffering_ratio_col } })),
              SD = sd({ { buffering_ratio_col } }, na.rm = TRUE),
              AneuploidyScore = first({ { aneuploidy_col } })) %>%
    filter(Observations > 50) %>%
    mutate(Model.Buffering.Ratio.ZScore = z_score(Model.Buffering.Ratio),
           Rank = as.integer(rank(Model.Buffering.Ratio.ZScore)),
           Model.Buffering.Ratio.Adjusted = lm(Model.Buffering.Ratio ~ AneuploidyScore, data = .)$residuals) %>%
    select(-AneuploidyScore) %>%
    arrange(Rank)
}

# Note: Not required when merging datasets to match genes and cell lines
filter_samples <- function(df_expr_buf) {
  df_expr_buf %>%
    filter(CellLine.AneuploidyScore > 0 | round(CellLine.Ploidy) != 2) %>%  # Remove non-aneuploid cell lines
    filter_cn_diff(remove_between = c(-0.01, 0.02)) %>% # Remove noise from discontinuity points of buffering ratio
    # filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
    filter(Buffering.GeneLevel.Ratio.Confidence > 0.3) %>% # Remove BR estimates with low confidence
    drop_na(Buffering.GeneLevel.Ratio)
}

add_cellline_names <- function (df, df_cl = df_celllines) {
  df_cl <- df_cl %>%
    distinct(Model.ID, CellLine.Name)

  df %>%
    inner_join(y = df_cl, by = "Model.ID", relationship = "many-to-one")
}

## Combine datasets
expr_buf_procan_filtered <- expr_buf_procan %>%
  filter_samples() %>%
  select("Model.ID", "Gene.Symbol", "Protein.Uniprot.Accession",
         "Buffering.GeneLevel.Ratio", "CellLine.AneuploidyScore", "Dataset") %>%
  drop_na()

expr_buf_depmap_filtered <- expr_buf_depmap %>%
  filter_samples() %>%
  select("Model.ID", "Gene.Symbol", "Protein.Uniprot.Accession",
         "Buffering.GeneLevel.Ratio", "CellLine.AneuploidyScore", "Dataset") %>%
  drop_na()

expr_buf_joined <- expr_buf_procan_filtered %>%
  rename(ProCan = Buffering.GeneLevel.Ratio) %>%
  select(-CellLine.AneuploidyScore) %>%
  inner_join(y = expr_buf_depmap_filtered %>% rename(DepMap = Buffering.GeneLevel.Ratio),
             by = c("Model.ID", "Gene.Symbol", "Protein.Uniprot.Accession"),
             relationship = "one-to-one", na_matches = "never")

common_genes <- intersect(unique(expr_buf_procan$Gene.Symbol),
                          unique(expr_buf_depmap$Gene.Symbol))

## Calculate Cell Line level Dosage Compensation
### All genes & cell lines
cellline_buf_procan <- expr_buf_procan_filtered %>%
  analyze_model_buffering(Buffering.GeneLevel.Ratio) %>%
  add_cellline_names() %>%
  mutate(Dataset = "ProCan")

cellline_buf_depmap <- expr_buf_depmap_filtered %>%
  analyze_model_buffering(Buffering.GeneLevel.Ratio) %>%
  add_cellline_names() %>%
  mutate(Dataset = "DepMap")

### Common genes, all cell lines
cellline_buf_gene_filtered_procan <- expr_buf_procan_filtered %>%
  filter(Gene.Symbol %in% common_genes) %>%
  analyze_model_buffering(Buffering.GeneLevel.Ratio) %>%
  add_cellline_names() %>%
  mutate(Dataset = "ProCan")

cellline_buf_gene_filtered_depmap <- expr_buf_depmap_filtered %>%
  filter(Gene.Symbol %in% common_genes) %>%
  analyze_model_buffering(Buffering.GeneLevel.Ratio) %>%
  add_cellline_names() %>%
  mutate(Dataset = "DepMap")

cellline_buf_merged_gene <- cellline_buf_gene_filtered_procan %>%
  bind_rows(cellline_buf_gene_filtered_depmap) %>%
  arrange(CellLine.Name)

### Joined on common genes and cell lines
cellline_buf_joined_procan <- expr_buf_joined %>%
  analyze_model_buffering(ProCan) %>%
  add_cellline_names() %>%
  mutate(Dataset = "ProCan")

cellline_buf_joined_depmap <- expr_buf_joined %>%
  analyze_model_buffering(DepMap) %>%
  add_cellline_names() %>%
  mutate(Dataset = "DepMap")

cellline_buf_joined <- cellline_buf_joined_procan %>%
  bind_rows(cellline_buf_joined_depmap) %>%
  arrange(CellLine.Name)

## Save results
write.xlsx(cellline_buf_procan, here(tables_base_dir, "cellline_buffering_procan.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_depmap, here(tables_base_dir, "cellline_buffering_depmap.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_gene_filtered_procan, here(tables_base_dir, "cellline_buffering_gene_filtered_procan.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_gene_filtered_depmap, here(tables_base_dir, "cellline_buffering_gene_filtered_depmap.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_joined_procan, here(tables_base_dir, "cellline_buffering_filtered_procan.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_joined_depmap, here(tables_base_dir, "cellline_buffering_filtered_depmap.xlsx"),
           colNames = TRUE)

write_parquet(cellline_buf_procan, here(output_data_dir, "cellline_buffering_procan.parquet"),
              version = "2.6")
write_parquet(cellline_buf_depmap, here(output_data_dir, "cellline_buffering_depmap.parquet"),
              version = "2.6")
write_parquet(cellline_buf_gene_filtered_procan, here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"),
              version = "2.6")
write_parquet(cellline_buf_gene_filtered_depmap, here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"),
              version = "2.6")
write_parquet(cellline_buf_joined_procan, here(output_data_dir, "cellline_buffering_filtered_procan.parquet"),
              version = "2.6")
write_parquet(cellline_buf_joined_depmap, here(output_data_dir, "cellline_buffering_filtered_depmap.parquet"),
              version = "2.6")

## Create plots
cellline_buf_waterfall_procan <- cellline_buf_procan %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_procan.png")
cellline_buf_waterfall_depmap <- cellline_buf_depmap %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_depmap.png")

# ToDo: Facetted waterfall plot
cellline_buf_waterfall_filtered_procan <- cellline_buf_joined_procan %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_filtered_procan.png")
cellline_buf_waterfall_filtered_depmap <- cellline_buf_joined_depmap %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_filtered_depmap.png")

cellline_buf_waterfall_gene_filtered_procan <- cellline_buf_gene_filtered_procan %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_gene_filtered_procan.png")
cellline_buf_waterfall_gene_filtered_depmap <- cellline_buf_gene_filtered_depmap %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_gene_filtered_depmap.png")

### Show cell line intersection betwen datasets in Venn diagram
celllines <- list(ProCan = (expr_buf_procan_filtered %>% distinct(Model.ID))$Model.ID,
                  DepMap = (expr_buf_depmap_filtered %>% distinct(Model.ID))$Model.ID)

ggvenn(celllines, columns = c("ProCan", "DepMap"), fill_alpha = 2/3,
       fill_color = c(bidirectional_color_pal[1], bidirectional_color_pal[5]), show_percentage = FALSE) %>%
  save_plot("cellline_venn_filtered.png", height = 120, width = 150)

# === Determine Correlation between Datasets ===
cellline_dist <- cellline_buf_joined %>%
  violin_plot(Dataset, Model.Buffering.Ratio) %>%
  save_plot("cellline_buffering_distribution.png")

# TODO: QQ-Plot (geom_qq) & shapiro.test()

## Full datasets
cellline_buf_full <- bind_rows(
  cellline_buf_procan %>% select(Model.ID, Model.Buffering.Ratio, Dataset),
  cellline_buf_depmap %>% select(Model.ID, Model.Buffering.Ratio, Dataset)
) %>%
  pivot_wider(names_from = "Dataset", values_from = "Model.Buffering.Ratio", id_cols = "Model.ID")

cellline_kendall <- cor.test(x = cellline_buf_full$ProCan,
                             y = cellline_buf_full$DepMap,
                             method = "kendall")

cellline_pearson <- cor.test(x = cellline_buf_full$ProCan,
                             y = cellline_buf_full$DepMap,
                             method = "pearson")

cellline_spearman <- cor.test(x = cellline_buf_full$ProCan,
                             y = cellline_buf_full$DepMap,
                              method = "spearman")

cat(capture.output(cellline_kendall), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cellline_pearson), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = TRUE, sep = "\n")
cat(capture.output(cellline_spearman), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = TRUE, sep = "\n")

## Datasets filtered by common genes and cell lines
cellline_kendall_filtered <- cor.test(x = (cellline_buf_joined %>% filter(Dataset == "ProCan"))$Model.Buffering.Ratio,
                                      y = (cellline_buf_joined %>% filter(Dataset == "DepMap"))$Model.Buffering.Ratio,
                                      method = "kendall")

cellline_pearson_filtered <- cor.test(x = (cellline_buf_joined %>% filter(Dataset == "ProCan"))$Model.Buffering.Ratio,
                                      y = (cellline_buf_joined %>% filter(Dataset == "DepMap"))$Model.Buffering.Ratio,
                                      method = "pearson")

cellline_spearman_filtered <- cor.test(x = (cellline_buf_joined %>% filter(Dataset == "ProCan"))$Model.Buffering.Ratio,
                                       y = (cellline_buf_joined %>% filter(Dataset == "DepMap"))$Model.Buffering.Ratio,
                                       method = "spearman")

cat(capture.output(cellline_kendall_filtered), file = here(reports_dir, "cellline_buffering_correlation_filtered.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cellline_pearson_filtered), file = here(reports_dir, "cellline_buffering_correlation_filtered.txt"),
    append = TRUE, sep = "\n")
cat(capture.output(cellline_spearman_filtered), file = here(reports_dir, "cellline_buffering_correlation_filtered.txt"),
    append = TRUE, sep = "\n")

## Datasets filtered by common genes only
cellline_buf_merged_gene_wide <- cellline_buf_merged_gene %>%
  pivot_wider(names_from = "Dataset", values_from = "Model.Buffering.Ratio", id_cols = "Model.ID")

cellline_kendall_gene <- cor.test(x = cellline_buf_merged_gene_wide$ProCan,
                                  y = cellline_buf_merged_gene_wide$DepMap,
                                  method = "kendall")

cellline_pearson_gene <- cor.test(x = cellline_buf_merged_gene_wide$ProCan,
                                  y = cellline_buf_merged_gene_wide$DepMap,
                                  method = "pearson")

cellline_spearman_gene <- cor.test(x = cellline_buf_merged_gene_wide$ProCan,
                                   y = cellline_buf_merged_gene_wide$DepMap,
                                   method = "spearman")

cat(capture.output(cellline_kendall_gene), file = here(reports_dir, "cellline_buffering_correlation_gene.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cellline_pearson_gene), file = here(reports_dir, "cellline_buffering_correlation_gene.txt"),
    append = TRUE, sep = "\n")
cat(capture.output(cellline_spearman_gene), file = here(reports_dir, "cellline_buffering_correlation_gene.txt"),
    append = TRUE, sep = "\n")

# === Create aggregated Cell Line Ranking ===

cellline_buf_rank <- cellline_buf_merged_gene %>%
  mean_norm_rank(Model.Buffering.Ratio, Dataset, Model.ID)

cellline_buf_mean <- cellline_buf_merged_gene %>%
  standardized_mean(Model.Buffering.Ratio, Dataset, Model.ID)

cellline_buf_agg <- cellline_buf_rank %>%
  inner_join(y = cellline_buf_mean, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never", unmatched = "error") %>%
  rename(Model.Buffering.MeanNormRank = "MeanNormRank",
         Model.Buffering.StandardizedMean = "StandardizedMean") %>%
  add_cellline_names() %>%
  mutate(CellLine.Name = fct_reorder(CellLine.Name, Model.Buffering.MeanNormRank)) %>%
  write_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"), version = "2.6")

write.xlsx(cellline_buf_agg, here(tables_base_dir, "cellline_buffering_agg.xlsx"),
           colNames = TRUE)

## Compare aggregated ranks
agg_corr_kendall <- cor.test(x = cellline_buf_agg$Model.Buffering.MeanNormRank,
                             y = cellline_buf_agg$Model.Buffering.StandardizedMean,
                             method = "kendall")

agg_corr_pearson <- cor.test(x = cellline_buf_agg$Model.Buffering.MeanNormRank,
                             y = cellline_buf_agg$Model.Buffering.StandardizedMean,
                             method = "pearson")

agg_corr_spearman <- cor.test(x = cellline_buf_agg$Model.Buffering.MeanNormRank,
                             y = cellline_buf_agg$Model.Buffering.StandardizedMean,
                              method = "spearman")

cat(capture.output(agg_corr_kendall), file = here(reports_dir, "cellline_buffering_correlation_agg.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(agg_corr_pearson), file = here(reports_dir, "cellline_buffering_correlation_agg.txt"),
    append = TRUE, sep = "\n")
cat(capture.output(agg_corr_spearman), file = here(reports_dir, "cellline_buffering_correlation_agg.txt"),
    append = TRUE, sep = "\n")

## Plot aggregated ranks
plot_agg_top <- cellline_buf_agg %>%
  slice_max(Model.Buffering.MeanNormRank, n = 10) %>%
  vertical_bar_chart(CellLine.Name, Model.Buffering.MeanNormRank,
                     default_fill_color = head(bidirectional_color_pal, n = 1),
                     text_color = "black", bar_label_shift = 0.1, break_steps = 0.25,
                     value_range = c(0,1), line_intercept = 0.5, value_lab = "Mean Normalized Rank") %>%
  save_plot("cellline_buffering_aggregated_top.png")

plot_agg_bot <- cellline_buf_agg %>%
  slice_min(Model.Buffering.MeanNormRank, n = 10) %>%
  vertical_bar_chart(CellLine.Name, Model.Buffering.MeanNormRank,
                     default_fill_color = tail(bidirectional_color_pal, n = 1),
                     text_color = "black", bar_label_shift = 0.1, break_steps = 0.25,
                     value_range = c(0,1), line_intercept = 0.5, value_lab = "Mean Normalized Rank") %>%
  save_plot("cellline_buffering_aggregated_bot.png")

# === Evaluation ===

## Plot protein expression against copy number ratio for high and low buffering cell lines

eval_buf_expr_procan <- cellline_buf_procan %>%
  filter(Observations >= 500) %>%
  split_by_quantiles(Model.Buffering.Ratio, quantile_low = "10%", quantile_high = "90%",
                     target_group_col = "Bufffering.CellLine.Group") %>%
  inner_join(y = expr_buf_procan, by = "CellLine.Name") %>%
  mutate(CopyNumberRatio = Gene.CopyNumber / Gene.CopyNumber.Baseline)

eval_buf_high <- eval_buf_expr_procan %>%
  filter(Bufffering.CellLine.Group == "High") %>%
  scatter_plot_reg_corr(CopyNumberRatio, Protein.Expression.Normalized) %>%
  save_plot("eval_cellline_buffering_high.png")

eval_buf_low <- eval_buf_expr_procan %>%
  filter(Bufffering.CellLine.Group == "Low") %>%
  scatter_plot_reg_corr(CopyNumberRatio, Protein.Expression.Normalized) %>%
  save_plot("eval_cellline_buffering_low.png")

eval_buf_diff_expr <- eval_buf_expr_procan %>%
  signif_violin_plot(Bufffering.CellLine.Group, Protein.Expression.Normalized,
                     facet_col = CopyNumberRatio, test = t.test) %>%
  save_plot("eval_cellline_buffering_diff_expression.png")

eval_buf_diff_cn <- eval_buf_expr_procan %>%
  signif_violin_plot(Bufffering.CellLine.Group, CopyNumberRatio)  %>%
  save_plot("eval_cellline_buffering_diff_cn.png")

### Conclusion: Highly buffered cell lines show protein expression closer to baseline

# === Combine Plots for publishing ===
plot_bracket <- plot_corr_bracket(cellline_pearson_gene)
plot_stack1 <- cowplot::plot_grid(cellline_buf_waterfall_gene_filtered_procan,
                                  cellline_buf_waterfall_gene_filtered_depmap + ylab(NULL),
                                  nrow = 1, ncol = 2, align = "hv", axis = "tblr", labels = c("ProCan", "DepMap"),
                                  label_y = 0.98, label_x = 0.1, rel_widths = c(1, 1))
plot_stack2 <- cowplot::plot_grid(plot_bracket, plot_stack1,
                                  nrow = 2, ncol = 1,
                                  rel_heights = c(0.1, 1))

plot_agg <- cowplot::plot_grid(plot_agg_top + ylab(NULL) + cowplot::theme_minimal_vgrid()
                                 + theme(axis.text.x = element_blank()),
                               plot_agg_bot + cowplot::theme_minimal_vgrid(),
                               nrow = 2, ncol = 1, align = "v", axis = "lr")

plot_publish <- cowplot::plot_grid(plot_stack2, plot_agg,
                                  nrow = 1, ncol = 2, labels = c("A", "B"),
                                  rel_widths = c(1.618, 1))

cairo_pdf(here(plots_dir, "cellline_buffering_gene_filtered_comparison.pdf"), width = 12)
plot_publish
dev.off()

## Poster
cairo_pdf(here(plots_dir, "cellline_buffering_gene_filtered_procan_poster.pdf"))
cellline_buf_waterfall_gene_filtered_procan <- cellline_buf_gene_filtered_procan %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name, font_size = 6)
cellline_buf_waterfall_gene_filtered_procan +
  ylab("Mean Buffering Ratio (z-score)") +
  theme_light(base_size = 20)
dev.off()
