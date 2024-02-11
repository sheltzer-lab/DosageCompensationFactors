library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(readxl)
library(ggplot2)


here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))

output_data_dir <- output_data_base_dir
dir.create(output_data_dir, recursive = TRUE)

df_rep <- read_excel(here(external_data_dir, "1-s2.0-S2667237522001709-mmc3.xlsx"), sheet = 3) %>%
  rename(Gene.Symbol = 1,
         ReproducibilityRank.Ovarian = "Ovarian Reproducibility Rank",
         ReproducibilityRank.Colon = "Colon Reproducibility Rank",
         ReproducibilityRank.CCLE = "CCLE Reproducibility Rank",
         ReproducibilityRank.Aggregated = "Aggregated Reproducibility Rank") %>%
  distinct(Gene.Symbol, ReproducibilityRank.Aggregated) %>%
  drop_na() %>%
  updateGeneSymbols() %>%
  arrange(desc(ReproducibilityRank.Aggregated)) %>%
  write_parquet(here(output_data_dir, 'reproducibility_ranks.parquet'), version = "2.6")

rep_quant <- quantile(df_rep$ReproducibilityRank.Aggregated,
                      probs = seq(0, 1, 0.1))
filter_cutoff <- rep_quant["10%"]

df_rep_filtered <- df_rep %>%
  filter(ReproducibilityRank.Aggregated > filter_cutoff) %>%
  write_parquet(here(output_data_dir, 'reproducibility_ranks_filtered.parquet'), version = "2.6")


# === Evaluation ===
df_rep %>%
  ggplot() +
  aes(sample = ReproducibilityRank.Aggregated) +
  geom_hline(yintercept = filter_cutoff, color = "red") +
  geom_qq_line(alpha = 0.5) +
  geom_qq()

df_rep_removed <- df_rep %>%
  filter(ReproducibilityRank.Aggregated <= filter_cutoff) %>%
  arrange(Gene.Symbol)
