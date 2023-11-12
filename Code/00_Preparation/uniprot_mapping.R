library(dplyr)
library(here)
library(arrow)

source(here("Code", "parameters.R"))

output_data_dir <- output_data_base_dir

# Source : ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
uniprot_mapping <- read.table(here(external_data_dir, "HUMAN_9606_idmapping.dat"),
                                sep = "\t", dec = ".",
                                header = FALSE, stringsAsFactors = FALSE) %>%
  select(V1, V2, V3) %>%
  rename(Protein.Uniprot.Accession = V1,
         Target.ID.Type = V2,
         Target.ID = V3)

id_types <- unique(uniprot_mapping$Target.ID.Type)
selected_types <- c("UniProtKB-ID", "Gene_Name", "Ensembl")

uniprot_mapping_selected <- uniprot_mapping %>%
  filter(Target.ID.Type %in% selected_types) %>%
  group_by(Protein.Uniprot.Accession, Target.ID.Type) %>%
  summarize(Target.ID = first(Target.ID, na_rm = TRUE)) %>%
  pivot_wider(id_cols = Protein.Uniprot.Accession,
              values_from = Target.ID, names_from = Target.ID.Type) %>%
  rename(Protein.Uniprot.Symbol = "UniProtKB-ID",
         Gene.Symbol = "Gene_Name",
         Gene.ENSEMBL.Id = "Ensembl") %>%
  mutate(Gene.ENSEMBL.Id = sub("\\.\\d+$", "", Gene.ENSEMBL.Id)) %>%
  write_parquet(here(output_data_dir, 'uniprot_mapping.parquet'),
                version = "2.6")
