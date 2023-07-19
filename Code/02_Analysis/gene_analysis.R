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

test <- expr_buf_goncalves %>%
  select(Gene.Symbol, Buffering.GeneLevel.Ratio) %>%
  drop_na() %>%
  group_by(Gene.Symbol) %>%
  summarize(TTest.p = t.test(Buffering.GeneLevel.Ratio, mu = 0)$p.value,
            Buffering.GeneLevel.Ratio.Average = mean(Buffering.GeneLevel.Ratio)) %>%
  mutate(TTest.p.adjusted = p.adjust(TTest.p, method = "BY"))
