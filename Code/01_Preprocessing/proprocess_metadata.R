library(here)
library(tidyr)
library(dplyr)
library(arrow)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))

output_data_dir <- output_data_base_dir
dir.create(output_data_dir, recursive = TRUE)

cancer_genes <- vroom::vroom(here(external_data_dir, "cancerGeneList.tsv")) %>%
  rename(Occurrences = 9,
         Gene.Symbol = 1) %>%
  filter(Occurrences >= 2) %>%
  mutate(CancerDriverMode = case_when(`Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "Yes" ~ "OG/TSG",
                                      `Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "No" ~ "OG",
                                      `Is Oncogene` == "No" & `Is Tumor Suppressor Gene` == "Yes" ~ "TSG",
                                      TRUE ~ "Unknown")) %>%
  select(Gene.Symbol, CancerDriverMode, Occurrences) %>%
  updateGeneSymbols() %>%
  write_parquet(here(output_data_dir, "cancer_genes.parquet"))
