library(dplyr)
library(here)
library(arrow)

source(here("Code", "parameters.R"))

output_data_dir <- output_data_base_dir

# Source : ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz
uniprot_mapping <- read.table(here(external_data_dir, "HUMAN_9606_idmapping_selected.tab"),
                                sep = "\t", dec = ".",
                                header = FALSE, stringsAsFactors = FALSE) %>%
  select(V1, V2) %>%
  rename(Protein.Uniprot.Accession = V1,
         Protein.Uniprot.Symbol = V2) %>%
  write_parquet(here(output_data_dir, 'uniprot_mapping.parquet'),
                version = "2.6")
