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

analyze_cellline_buffering <- function(df, buffering_ratio_col) {
  mean_pop <- mean(df[[quo_name(enquo(buffering_ratio_col))]], na.rm = TRUE)
  sd_pop <- sd(df[[quo_name(enquo(buffering_ratio_col))]], na.rm = TRUE)

  cellline_buf_avg <- df %>%
    select(UniqueId, starts_with("CellLine"), { { buffering_ratio_col } }) %>%
    group_by(CellLine.Name) %>%
    summarize(Buffering.CellLine.Ratio = mean({ { buffering_ratio_col } }, na.rm = TRUE)) %>%
    mutate(Buffering.CellLine.Ratio.ZScore = (Buffering.CellLine.Ratio - mean_pop) / sd_pop,
           Rank = as.integer(rank(Buffering.CellLine.Ratio.ZScore))) %>%
    arrange(Rank)

  return(cellline_buf_avg)
}

waterfall_plot <- function(df, value_col, rank_col, label_col) {
  xlim <- c(0, max(df[[quo_name(enquo(rank_col))]]))
  label_nudge_x <- floor(xlim[2] / 5)

  ylim1 <- c(min(df[[quo_name(enquo(value_col))]]),
             (df %>% filter({ { rank_col } } == label_nudge_x))[[quo_name(enquo(value_col))]])
  ylim2 <- c((df %>% filter({ { rank_col } } == xlim[2] - label_nudge_x))[[quo_name(enquo(value_col))]],
             max(df[[quo_name(enquo(value_col))]]))

  df %>%
    ggplot() +
    aes(x = { { rank_col } }, y = { { value_col } }, label = { { label_col } }) +
    geom_hline(yintercept = 0, color = "red") +
    geom_point(size = 0.3) +
    geom_text_repel(data = df %>% slice_min({ { rank_col } }, n = 5),
                    xlim = xlim, ylim = ylim1, direction = "y", nudge_x = label_nudge_x,
                    seed = 42, color = "darkblue") +
    geom_text_repel(data = df %>% slice_max({ { rank_col } }, n = 5),
                    xlim = xlim, ylim = ylim2, direction = "y", nudge_x = -label_nudge_x,
                    seed = 42, color = "darkred")
}

cellline_buf_goncalves <- expr_buf_goncalves %>%
  analyze_cellline_buffering(Buffering.GeneLevel.Ratio)

write.xlsx(cellline_buf_goncalves, here(tables_base_dir, "cellline_buffering_goncalves.xlsx"),
           colNames = TRUE)

cellline_buf_waterfall_goncalves <- cellline_buf_goncalves %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name)

ggsave(here(plots_dir, "cellline_buffering_waterfall_goncalves.png"),
       plot = cellline_buf_waterfall_goncalves,
       height = 200, width = 200, units = "mm", dpi = 300)

cellline_buf_depmap <- expr_buf_depmap %>%
  analyze_cellline_buffering(Buffering.GeneLevel.Ratio)

write.xlsx(cellline_buf_depmap, here(tables_base_dir, "cellline_buffering_depmap.xlsx"),
           colNames = TRUE)

cellline_buf_waterfall_depmap <- cellline_buf_depmap %>%
  waterfall_plot(Buffering.CellLine.Ratio.ZScore, Rank, CellLine.Name)

ggsave(here(plots_dir, "cellline_buffering_waterfall_depmap.png"),
       plot = cellline_buf_waterfall_depmap,
       height = 200, width = 200, units = "mm", dpi = 300)


# === Rank Correlation ===

cellline_buf_merged <- cellline_buf_goncalves %>%
  select("CellLine.Name", "Buffering.CellLine.Ratio.ZScore") %>%
  inner_join(y = cellline_buf_depmap %>% select("CellLine.Name", "Buffering.CellLine.Ratio.ZScore"),
             by = "CellLine.Name", relationship = "one-to-one", na_matches = "never") %>%
  arrange(CellLine.Name) %>%
  rename(ProCan = Buffering.CellLine.Ratio.ZScore.x,
         DepMap = Buffering.CellLine.Ratio.ZScore.y)

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