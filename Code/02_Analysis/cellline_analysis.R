library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(skimr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "02_Analysis", "analysis.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Univariate", "Goncalves")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_goncalves <- read_parquet(here(output_data_dir, "expression_buffering_goncalves.parquet"))

# === Analyze Dosage Compensation on Cell Line level ===

mean_pop <- mean(expr_buf_goncalves$Buffering.GeneLevel.Ratio, na.rm = TRUE)
sd_pop <- sd(expr_buf_goncalves$Buffering.GeneLevel.Ratio, na.rm = TRUE)

# ToDo: Why are there only 616 cell lines?
cellline_buf_avg <- expr_buf_goncalves %>%
  select(UniqueId, starts_with("CellLine"), Buffering.GeneLevel.Ratio) %>%
  group_by(CellLine.Name) %>%
  summarize(Buffering.CellLine.Ratio = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE)) %>%
  mutate(Buffering.CellLine.Ratio.ZScore = (Buffering.CellLine.Ratio - mean_pop) / sd_pop,
         Rank = as.integer(rank(Buffering.CellLine.Ratio.ZScore)))

## Waterfall plot
xlim <- c(0, max(cellline_buf_avg$Rank))
ylim <- c(min(cellline_buf_avg$Buffering.CellLine.Ratio.ZScore),
          max(cellline_buf_avg$Buffering.CellLine.Ratio.ZScore))

cellline_buf_avg %>%
  ggplot() +
  aes(x = Rank, y = Buffering.CellLine.Ratio.ZScore, label = CellLine.Name) +
  geom_hline(yintercept = 0, color ="red") +
  geom_point(size = 0.3) +
  geom_text_repel(data = cellline_buf_avg %>% slice_max(Rank, n = 5),
                  xlim = xlim, ylim = ylim,
                  force = 10, seed = 42, color = "darkred") +
  geom_text_repel(data = cellline_buf_avg %>% slice_min(Rank, n = 5),
                  xlim = xlim, ylim = ylim,
                  force = 10, seed = 42, color = "darkblue")