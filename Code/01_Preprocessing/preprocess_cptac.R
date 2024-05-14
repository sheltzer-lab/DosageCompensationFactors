library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(magrittr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "preprocessing.R"))

expression_data_dir <- here(external_data_dir, "Expression", "CPTAC", "Proteome_BCM_GENCODE_v34_harmonized_v1")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "Expression", "CPTAC")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
file_list <- list.files(expression_data_dir, include.dirs = FALSE, recursive = FALSE)
file_list <- file_list[grep(".+_proteomics_.+.txt", file_list)]

df_list <- list()

for (filename in file_list) {
  metadata <- sub('\\.txt$', '', filename) %>%
    str_split_fixed("_", 9) %>%
    as.data.frame() %>%
    select(V1, V9)

  df_name <- paste(metadata$V1, metadata$V9, sep = "_")

  df_list[[df_name]] <- read.table(here(expression_data_dir, filename),
                                   sep = "\t", dec = ".", header = FALSE, stringsAsFactors = FALSE) %>%
    janitor::row_to_names(1) %>%
    pivot_longer(everything() & !idx,
                 names_to = "Model.ID", values_to = "Protein.Expression.Log2",
                 names_ptypes = character(), values_ptypes = double(),
                 names_transform = as.character, values_transform = as.numeric) %>%
    mutate(Model.CancerType = metadata$V1,
           Model.SampleType = metadata$V9)
}