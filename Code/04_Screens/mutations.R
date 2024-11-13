library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Screens", "Mutations")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

df_model_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))
expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac.parquet"))
cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))

cptac_mut <- df_model_cptac %>%
  select(Model.ID, ends_with("mutation")) %>%
  pivot_longer(ends_with("mutation"), names_to = "Gene.Symbol", values_to = "Mutations") %>%
  mutate(Gene.Symbol = str_remove(Gene.Symbol, "_mutation"),
         Mutated = replace_na(Mutations > 0, FALSE))

cptac_buf <- expr_buf_cptac %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNEvent = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  select(Model.ID, Gene.Symbol, Buffering.GeneLevel.Class, CNEvent) %>%
  mutate(Buffered = replace_na(Buffering.GeneLevel.Class == "Buffered", FALSE))

cptac_buf_mut <- cptac_mut %>%
  inner_join(cptac_buf, by = c("Model.ID", "Gene.Symbol")) %>%
  count(Gene.Symbol, Buffered, Mutated, CNEvent) %>%
  complete(Gene.Symbol, Buffered, Mutated, CNEvent, fill = list(n = 0))

## MutEx analysis seperated by gain and loss events
cptac_buf_mut_mutex <- cptac_buf_mut %>%
  group_by(Gene.Symbol, CNEvent) %>%
  group_modify(~ {
    .x %>%
      df2contingency(Buffered, Mutated) %>%
      mutex_score() %>%
      as.data.frame()
  }) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

## MutEx analysis for gain and loss events combined
cptac_buf_mut_mutex_comb <- cptac_mut %>%
  inner_join(cptac_buf, by = c("Model.ID", "Gene.Symbol")) %>%
  count(Gene.Symbol, Buffered, Mutated) %>%
  complete(Gene.Symbol, Buffered, Mutated, fill = list(n = 0)) %>%
  group_by(Gene.Symbol) %>%
  group_modify(~ {
    .x %>%
      df2contingency(Buffered, Mutated) %>%
      mutex_score() %>%
      as.data.frame()
  }) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

# TODO: Use model BR groups against mutational signatures
