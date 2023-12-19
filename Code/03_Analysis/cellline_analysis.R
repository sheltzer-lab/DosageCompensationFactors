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
## ToDo: Sanitize Cell Line names
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

# === Analyze Dosage Compensation on Cell Line level ===

analyze_cellline_buffering <- function(df, buffering_ratio_col, cellline_col = CellLine.Name) {
  mean_pop <- mean(df[[quo_name(enquo(buffering_ratio_col))]], na.rm = TRUE)
  sd_pop <- sd(df[[quo_name(enquo(buffering_ratio_col))]], na.rm = TRUE)

  cellline_buf_avg <- df %>%
    select({ { cellline_col } }, { { buffering_ratio_col } }) %>%
    group_by({ { cellline_col } }) %>%
    summarize(Buffering.CellLine.Ratio = mean({ { buffering_ratio_col } }, na.rm = TRUE)) %>%
    mutate(Buffering.CellLine.Ratio.ZScore = (Buffering.CellLine.Ratio - mean_pop) / sd_pop,
           Rank = as.integer(rank(Buffering.CellLine.Ratio.ZScore))) %>%
    arrange(Rank)

  return(cellline_buf_avg)
}

# Note: Not required when merging datasets to match genes and cell lines
filter_genes <- function(df, gene_col = Gene.Symbol, value_col = Protein.Expression.Normalized, keep = "90%") {
  df %>%
    group_by({ { gene_col } }) %>%
    mutate(CV = abs(sd({ { value_col } }, na.rm = TRUE) / mean({ { value_col } }, na.rm = TRUE))) %>%
    ungroup() %>%
    filter(CV < quantile(CV, probs = seq(0, 1, 0.01))[keep])
}

# Note: Does not seem to improve correlation
remove_antiscaling <- function(df, buffering_col, gene_col = Gene.Symbol) {
  df %>%
    group_by({ { gene_col } }) %>%
    mutate(AvgBuf = mean({ { buffering_col } }, na.rm = TRUE)) %>%
    filter(buffering_class(AvgBuf) != "Anti-Scaling") %>%
    select(-AvgBuf)
}

## Combine datasets
expr_buf_procan_filtered <- expr_buf_procan %>%
  filter(CellLine.AneuploidyScore > 0 | CellLine.WGD > 0) %>%  # Remove non-aneuploid cell lines
  filter_cn_diff(remove_between = c(-0, 0.0246)) %>% # Remove noise from discontinuity points of buffering ratio
  select("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession", "Buffering.GeneLevel.Ratio") %>%
  drop_na() %>%
  rename(ProCan = "Buffering.GeneLevel.Ratio")
expr_buf_depmap_filtered <- expr_buf_depmap %>%
  filter(CellLine.AneuploidyScore > 0 | CellLine.WGD > 0) %>%  # Remove non-aneuploid cell lines
  filter_cn_diff(remove_between = c(-0, 0.0246)) %>% # Remove noise from discontinuity points of buffering ratio
  select("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession", "Buffering.GeneLevel.Ratio") %>%
  drop_na() %>%
  rename(DepMap = "Buffering.GeneLevel.Ratio")

expr_buf_filtered <- expr_buf_procan_filtered %>%
  inner_join(y = expr_buf_depmap_filtered, by = c("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession"),
             relationship = "one-to-one", na_matches = "never")

common_genes <- intersect(unique(expr_buf_procan$Gene.Symbol),
                          unique(expr_buf_depmap$Gene.Symbol))

## Calculate Cell Line level Dosage Compensation
cellline_buf_procan <- expr_buf_procan %>%
  analyze_cellline_buffering(Buffering.GeneLevel.Ratio)
cellline_buf_depmap <- expr_buf_depmap %>%
  analyze_cellline_buffering(Buffering.GeneLevel.Ratio)

cellline_buf_gene_filtered_procan <- expr_buf_procan_filtered %>%
  filter(Gene.Symbol %in% common_genes) %>%
  analyze_cellline_buffering(ProCan) %>%
  mutate(Dataset = "ProCan")
cellline_buf_gene_filtered_depmap <- expr_buf_depmap_filtered %>%
  filter(Gene.Symbol %in% common_genes) %>%
  analyze_cellline_buffering(DepMap) %>%
  mutate(Dataset = "DepMap")

cellline_buf_filtered_procan <- expr_buf_filtered %>%
  analyze_cellline_buffering(ProCan) %>%
  mutate(Dataset = "ProCan")
cellline_buf_filtered_depmap <- expr_buf_filtered %>%
  analyze_cellline_buffering(DepMap) %>%
  mutate(Dataset = "DepMap")

cellline_buf_merged <- cellline_buf_filtered_procan %>%
  bind_rows(cellline_buf_filtered_depmap) %>%
  arrange(CellLine.Name)

cellline_buf_merged_gene <- cellline_buf_gene_filtered_procan %>%
  bind_rows(cellline_buf_gene_filtered_depmap) %>%
  pivot_wider(id_cols = "CellLine.Name",
              names_from = "Dataset",
              values_from = c("Buffering.CellLine.Ratio", "Buffering.CellLine.Ratio.ZScore", "Rank")) %>%
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
write.xlsx(cellline_buf_filtered_procan, here(tables_base_dir, "cellline_buffering_filtered_procan.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_filtered_depmap, here(tables_base_dir, "cellline_buffering_filtered_depmap.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_merged, here(tables_base_dir, "cellline_buffering_z-scores_merged.xlsx"),
           colNames = TRUE)

write_parquet(cellline_buf_procan, here(output_data_dir, "cellline_buffering_procan.parquet"),
              version = "2.6")
write_parquet(cellline_buf_depmap, here(output_data_dir, "cellline_buffering_depmap.parquet"),
              version = "2.6")
write_parquet(cellline_buf_gene_filtered_procan, here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"),
              version = "2.6")
write_parquet(cellline_buf_gene_filtered_depmap, here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"),
              version = "2.6")
write_parquet(cellline_buf_filtered_procan, here(output_data_dir, "cellline_buffering_filtered_procan.parquet"),
              version = "2.6")
write_parquet(cellline_buf_filtered_depmap, here(output_data_dir, "cellline_buffering_filtered_depmap.parquet"),
              version = "2.6")

## Create plots
cellline_buf_waterfall_procan <- cellline_buf_procan %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_procan.png")
cellline_buf_waterfall_depmap <- cellline_buf_depmap %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_depmap.png")

# ToDo: Facetted waterfall plot
cellline_buf_waterfall_filtered_procan <- cellline_buf_filtered_procan %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_filtered_procan.png")
cellline_buf_waterfall_filtered_depmap <- cellline_buf_filtered_depmap %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_filtered_depmap.png")

### Show cell line intersection betwen datasets in Venn diagram
celllines <- list(ProCan = (expr_buf_procan_filtered %>% distinct(CellLine.Name))$CellLine.Name,
                  DepMap = (expr_buf_depmap_filtered %>% distinct(CellLine.Name))$CellLine.Name)

ggvenn(celllines, columns = c("ProCan", "DepMap"), fill_alpha = 2/3,
       fill_color = c(bidirectional_color_pal[1], bidirectional_color_pal[5]), show_percentage = FALSE) %>%
  save_plot("cellline_venn_filtered.png", height = 120, width = 150)

# === Determine Correlation between Datasets ===
cellline_dist <- cellline_buf_merged %>%
  violin_plot(Dataset, Buffering.CellLine.Ratio.ZScore) %>%
  save_plot("cellline_buffering_distribution.png")


cellline_kendall <- cor.test(x = (cellline_buf_merged %>% filter(Dataset == "ProCan"))$Buffering.CellLine.Ratio,
                             y = (cellline_buf_merged %>% filter(Dataset == "DepMap"))$Buffering.CellLine.Ratio,
                             method = "kendall")

cellline_pearson <- cor.test(x = (cellline_buf_merged %>% filter(Dataset == "ProCan"))$Buffering.CellLine.Ratio,
                             y = (cellline_buf_merged %>% filter(Dataset == "DepMap"))$Buffering.CellLine.Ratio,
                             method = "pearson")

cellline_spearman <- cor.test(x = (cellline_buf_merged %>% filter(Dataset == "ProCan"))$Buffering.CellLine.Ratio,
                              y = (cellline_buf_merged %>% filter(Dataset == "DepMap"))$Buffering.CellLine.Ratio,
                              method = "spearman")

cat(capture.output(cellline_kendall), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cellline_pearson), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = TRUE, sep = "\n")
cat(capture.output(cellline_spearman), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = TRUE, sep = "\n")

## Datasets filtered by common genes only
cellline_kendall_gene <- cor.test(x = cellline_buf_merged_gene$Buffering.CellLine.Ratio_ProCan,
                                  y = cellline_buf_merged_gene$Buffering.CellLine.Ratio_DepMap,
                                  method = "kendall")

cellline_pearson_gene <- cor.test(x = cellline_buf_merged_gene$Buffering.CellLine.Ratio_ProCan,
                                  y = cellline_buf_merged_gene$Buffering.CellLine.Ratio_DepMap,
                                  method = "pearson")

cellline_spearman_gene <- cor.test(x = cellline_buf_merged_gene$Buffering.CellLine.Ratio_ProCan,
                                   y = cellline_buf_merged_gene$Buffering.CellLine.Ratio_DepMap,
                                   method = "spearman")

cat(capture.output(cellline_kendall_gene), file = here(reports_dir, "cellline_buffering_correlation_gene.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cellline_pearson_gene), file = here(reports_dir, "cellline_buffering_correlation_gene.txt"),
    append = TRUE, sep = "\n")
cat(capture.output(cellline_spearman_gene), file = here(reports_dir, "cellline_buffering_correlation_gene.txt"),
    append = TRUE, sep = "\n")

# === Create aggregated Cell Line Ranking ===

cellline_buf_rank <- cellline_buf_filtered_procan %>%
  bind_rows(cellline_buf_filtered_depmap) %>%
  mean_norm_rank(Rank, Dataset, CellLine.Name)

cellline_buf_mean <- cellline_buf_filtered_procan %>%
  bind_rows(cellline_buf_filtered_depmap) %>%
  standardized_mean(Buffering.CellLine.Ratio, Dataset, CellLine.Name)

cellline_buf_agg <- cellline_buf_rank %>%
  inner_join(y = cellline_buf_mean, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never", unmatched = "error") %>%
  rename(Buffering.CellLine.MeanNormRank = "MeanNormRank",
         Buffering.CellLine.StandardizedMean = "StandardizedMean") %>%
  mutate(CellLine.Name = fct_reorder(CellLine.Name, Buffering.CellLine.MeanNormRank)) %>%
  write_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"), version = "2.6")

write.xlsx(cellline_buf_agg, here(tables_base_dir, "cellline_buffering_agg.xlsx"),
           colNames = TRUE)

cellline_buf_agg %>%
  slice_max(Buffering.CellLine.MeanNormRank, n = 20) %>%
  vertical_bar_chart(CellLine.Name, Buffering.CellLine.MeanNormRank,
                     value_range = c(0.5,1), line_intercept = 0.5, value_lab = "Aggregated Rank") %>%
  save_plot("cellline_buffering_aggregated_top.png")

cellline_buf_agg %>%
  slice_min(Buffering.CellLine.MeanNormRank, n = 20) %>%
  vertical_bar_chart(CellLine.Name, Buffering.CellLine.MeanNormRank,
                     value_range = c(0,0.5), line_intercept = 0.5, value_lab = "Aggregated Rank") %>%
  save_plot("cellline_buffering_aggregated_bot.png")

# === Combine Plots for publishing ===
plot_bracket <- plot_corr_bracket(cellline_pearson)
plot_stack1 <- cowplot::plot_grid(cellline_buf_waterfall_filtered_procan, cellline_buf_waterfall_filtered_depmap,
                                  nrow = 1, ncol = 2, align = "h", axis = "lr", labels = c("ProCan", "DepMap"),
                                  label_y = 0.98, label_x = 0.05, rel_widths = c(1, 1))
plot_stack2 <- cowplot::plot_grid(plot_bracket, plot_stack1,
                                  nrow = 2, ncol = 1,
                                  rel_heights = c(0.1, 1))

cairo_pdf(here(plots_dir, "cellline_buffering_filtered_comparison.pdf"), width = 11)
plot_stack2
dev.off()
