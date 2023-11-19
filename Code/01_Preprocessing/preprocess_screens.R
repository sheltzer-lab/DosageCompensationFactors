library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))

copynumber_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

regex_parentheses <- "(.*)\\s+\\((.*)\\)"

df_celllines <- read_parquet(here(output_data_dir, "celllines.parquet"))

df_prism <- read_csv_arrow(here(screens_data_dir, "PRISM_Repurposing_Public_23Q2.csv")) %>%
  rename(CellLine.DepMapModelId = 1) %>%
  pivot_longer(everything() & !CellLine.DepMapModelId,
               names_to = "Drug", values_to = "Drug.MFI.Log2FC") %>%
  mutate(Drug.Name = str_extract(Drug, regex_parentheses, group = 1),
         Drug.ID = str_extract(Drug, regex_parentheses, group = 2)) %>%
  left_join(y = df_celllines, by = "CellLine.DepMapModelId",
            relationship = "many-to-one", na_matches = "never") %>%
  write_parquet(here(output_data_dir, 'drug_screens.parquet'),
                version = "2.6")

df_growth <- read_csv_arrow(here(screens_data_dir, "growth_rate_20220907.csv")) %>%
  select(model_name, model_id, seeding_density, day4_day1_ratio, replicates) %>%
  rename(CellLine.SangerModelId = "model_id",
         CellLine.Name = "model_name",
         CellLine.SeedingDensity = "seeding_density",
         CellLine.GrowthRatio = "day4_day1_ratio",
         CellLine.Replicates = "replicates") %>%
  filter(CellLine.Replicates > 2) %>%
  group_by(CellLine.Name) %>%
  slice_max(CellLine.Replicates, n = 1) %>%
  ungroup() %>%
  left_join(y = df_celllines, by = "CellLine.SangerModelId",
            relationship = "many-to-one", na_matches = "never") %>%
  write_parquet(here(output_data_dir, 'cellline_growth.parquet'),
                version = "2.6")

df_crispr <- read_csv_arrow(here(screens_data_dir, "CRISPRGeneEffect.csv")) %>%
  rename(CellLine.DepMapModelId = ModelID) %>%
  pivot_longer(everything() & !CellLine.DepMapModelId,
               names_to = "Gene", values_to = "CRISPR.EffectScore") %>%
  mutate(Gene.Symbol = str_extract(Gene, regex_parentheses, group = 1),
         Gene.Entrez.Id = str_extract(Gene, regex_parentheses, group = 2)) %>%
  left_join(y = df_celllines, by = "CellLine.DepMapModelId",
            relationship = "many-to-one", na_matches = "never") %>%
  id2uniprot_acc("Gene.Symbol", "hgnc_symbol") %>%
  select(CellLine.CustomId, CellLine.DepMapModelId, CellLine.SangerModelId, CellLine.Name,
         Protein.Uniprot.Accession, Gene.Symbol, CRISPR.EffectScore) %>%
  write_parquet(here(output_data_dir, 'crispr_screens.parquet'),
                version = "2.6")
