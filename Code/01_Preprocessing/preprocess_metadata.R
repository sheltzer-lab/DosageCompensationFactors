library(here)
library(tidyr)
library(dplyr)
library(vroom)
library(reshape2)
library(arrow)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))

output_data_dir <- output_data_base_dir
dir.create(output_data_dir, recursive = TRUE)

df_celllines <- read_parquet(here(output_data_dir, "celllines.parquet"))

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

# https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2024Q2&filename=OmicsSomaticMutationsMatrixDamaging.csv
mutation_matrix <- vroom(here(external_data_dir, "OmicsSomaticMutationsMatrixDamaging.csv"))

df_mutations <- mutation_matrix %>%
  melt() %>%
  rename(CellLine.DepMapModelId = 1,
         Gene = 2,
         MutationStatus = 3) %>%
  separate_wider_delim(Gene, delim = " ", names = c("Gene.Symbol", "Gene.ID")) %>%
  mutate(MutationStatus = as.logical(MutationStatus)) %>%
  left_join(y = df_celllines, by = "CellLine.DepMapModelId",
            relationship = "many-to-one", na_matches = "never") %>%
  select(Model.ID, Gene.Symbol, MutationStatus) %>%
  write_parquet(here(output_data_dir, "damaging_mutations_depmap.parquet"))
