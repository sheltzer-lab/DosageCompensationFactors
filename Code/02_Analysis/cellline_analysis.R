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
source(here("Code", "02_Analysis", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "CellLine")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_goncalves <- read_parquet(here(output_data_dir, "expression_buffering_goncalves.parquet"))
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
filter_genes <- function (df, gene_col = Gene.Symbol, value_col = Protein.Expression.Normalized, keep = "90%") {
  df %>%
    group_by({ { gene_col } }) %>%
    mutate(CV = abs(sd({{value_col}}, na.rm = TRUE) / mean({{value_col}}, na.rm = TRUE))) %>%
    ungroup() %>%
    filter(CV < quantile(CV, probs = seq(0, 1, 0.01))[keep])
}

# Note: Does not seem to improve correlation
remove_antiscaling <- function(df, buffering_col, gene_col = Gene.Symbol) {
  df %>%
    group_by({{gene_col}}) %>%
    mutate(AvgBuf = mean({{buffering_col}}, na.rm = TRUE)) %>%
    filter(buffering_class(AvgBuf) != "Anti-Scaling") %>%
    select(-AvgBuf)
}

## Combine datasets
expr_buf_goncalves_filtered <- expr_buf_goncalves %>%
  filter_cn_diff(remove_between = c(-0, 0.0246)) %>% # Remove noise from discontinuity points of buffering ratio
  select("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession", "Buffering.GeneLevel.Ratio") %>%
  drop_na() %>%
  rename(ProCan = "Buffering.GeneLevel.Ratio")
expr_buf_depmap_filtered <- expr_buf_depmap %>%
  filter_cn_diff(remove_between = c(-0, 0.0246)) %>% # Remove noise from discontinuity points of buffering ratio
  select("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession", "Buffering.GeneLevel.Ratio") %>%
  drop_na() %>%
  rename(DepMap = "Buffering.GeneLevel.Ratio")

expr_buf_filtered <- expr_buf_goncalves_filtered %>%
  inner_join(y = expr_buf_depmap_filtered, by = c("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession"),
             relationship = "one-to-one", na_matches = "never")

## Calculate Cell Line level Dosage Compensation
cellline_buf_goncalves <- expr_buf_goncalves %>%
  analyze_cellline_buffering(Buffering.GeneLevel.Ratio)
cellline_buf_depmap <- expr_buf_depmap %>%
  analyze_cellline_buffering(Buffering.GeneLevel.Ratio)

cellline_buf_filtered_goncalves <- expr_buf_filtered %>%
  analyze_cellline_buffering(ProCan)
cellline_buf_filtered_depmap <- expr_buf_filtered %>%
  analyze_cellline_buffering(DepMap)

cellline_buf_merged <- cellline_buf_filtered_goncalves %>%
  select("CellLine.Name", "Buffering.CellLine.Ratio.ZScore") %>%
  inner_join(y = cellline_buf_filtered_depmap %>% select("CellLine.Name", "Buffering.CellLine.Ratio.ZScore"),
             by = "CellLine.Name", relationship = "one-to-one", na_matches = "never") %>%
  arrange(CellLine.Name) %>%
  rename(ProCan = Buffering.CellLine.Ratio.ZScore.x,
         DepMap = Buffering.CellLine.Ratio.ZScore.y)

## Save results
write.xlsx(cellline_buf_goncalves, here(tables_base_dir, "cellline_buffering_goncalves.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_depmap, here(tables_base_dir, "cellline_buffering_depmap.xlsx"),
           colNames = TRUE)

write.xlsx(cellline_buf_filtered_goncalves, here(tables_base_dir, "cellline_buffering_filtered_goncalves.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_filtered_depmap, here(tables_base_dir, "cellline_buffering_filtered_depmap.xlsx"),
           colNames = TRUE)
write.xlsx(cellline_buf_merged, here(tables_base_dir, "cellline_buffering_z-scores_merged.xlsx"),
           colNames = TRUE)

## Create plots
cellline_buf_waterfall_goncalves <- cellline_buf_goncalves %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_goncalves.png")
cellline_buf_waterfall_depmap <- cellline_buf_depmap %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_depmap.png")

cellline_buf_waterfall_filtered_goncalves <- cellline_buf_filtered_goncalves %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_filtered_goncalves.png")
cellline_buf_waterfall_filtered_depmap <- cellline_buf_filtered_depmap %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name) %>%
  save_plot("cellline_buffering_waterfall_filtered_depmap.png")


# === Determine Correlation between Datasets ===
cellline_dist <- cellline_buf_merged %>%
  pivot_longer(c(ProCan, DepMap), names_to = "Dataset", values_to = "Buffering.CellLine.Ratio.ZScore") %>%
  violin_plot(Dataset, Buffering.CellLine.Ratio.ZScore)

ggsave(here(plots_dir, "cellline_buffering_distribution.png"),
       plot = cellline_dist,
       height = 200, width = 200, units = "mm", dpi = 300)

cellline_kendall <- cor.test(x = cellline_buf_merged$ProCan,
                             y = cellline_buf_merged$DepMap,
                             method = "kendall")

cellline_pearson <- cor.test(x = cellline_buf_merged$ProCan,
                             y = cellline_buf_merged$DepMap,
                             method = "pearson")

cat(capture.output(cellline_kendall), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cellline_pearson), file = here(reports_dir, "cellline_buffering_correlation.txt"),
    append = TRUE, sep = "\n")
