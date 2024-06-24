library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "analysis.R"))
source(here("Code", "visualization.R"))

ubi_dir <- here(external_data_dir, "PTMs", "CPTAC", "Ubiquitylome_CDAP_v1")
output_data_dir <- output_data_base_dir

df_ubi <- read.table(here(ubi_dir, "CPTAC3_Lung_Squamous_Cell_Carcinoma_Ubiquitylome.ubiquitylsite.tmt11.tsv"),
                     sep = "\t", dec = ".", header = FALSE, stringsAsFactors = FALSE) %>%
  janitor::row_to_names(1)

df_ubi_meta <- read.table(here(ubi_dir, "PDC_biospecimen_manifest.tsv"),
                          sep = "\t", dec = ".", header = FALSE, stringsAsFactors = FALSE) %>%
    janitor::row_to_names(1)

aliquot_id_mapping <- df_ubi_meta %>%
  rename(Model.ID = "Case Submitter ID",
         Model.AliquotSubmitterID = "Aliquot Submitter ID") %>%
  select(Model.ID, Model.AliquotSubmitterID)

df_ubi_tidy <- df_ubi %>%
  mutate(across(matches("Log Ratio"), as.numeric)) %>%
  pivot_longer(matches("Log Ratio"),
               names_to = "Model.AliquotSubmitterID",
               values_to = "Protein.Ubiquitination.Log2") %>%
  mutate(Model.AliquotSubmitterID = str_split_i(Model.AliquotSubmitterID, " ", 1)) %>%
  left_join(aliquot_id_mapping, by = "Model.AliquotSubmitterID") %>%
  rename(Gene.Symbol = Gene,
         Protein.Ubiquitination.Site = Ubiquitylsite,
         Protein.Ubiquitination.Peptide = Peptide) %>%
  select(Model.ID, Model.AliquotSubmitterID, Gene.Symbol,
         Protein.Ubiquitination.Peptide, Protein.Ubiquitination.Site, Protein.Ubiquitination.Log2)

df_ubi_tidy %>%
  write_parquet(here(output_data_dir, 'ubiquitination_cptac.parquet'),
              version = "2.6")

# === Experimental ===
library(magrittr)

expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac.parquet"))

df_ubi_processed <- df_ubi_tidy %>%
  mutate(Protein.Ubiquitination.Log2 = as.vector(scale(Protein.Ubiquitination.Log2))) %>%
  group_by(Model.ID, Gene.Symbol) %>%
  summarize(Protein.Ubiquitination.Log2 = mean(Protein.Ubiquitination.Log2, na.rm = TRUE), .groups = "drop")

ubi_buf_cptac <- expr_buf_cptac %>%
  inner_join(y = df_ubi_processed, by = c("Model.ID", "Gene.Symbol")) %>%
  tibble::column_to_rownames("UniqueId")

resid <- residuals(lm(Protein.Ubiquitination.Log2 ~ Protein.Expression.Normalized, ubi_buf_cptac))

ubi_buf_cptac <- ubi_buf_cptac %>%
  mutate(Protein.Ubiquitination.Adj = resid[rownames(ubi_buf_cptac)],
         Buffering.GeneLevel.Ratio = as.vector(scale(Buffering.GeneLevel.Ratio)))

ubi_buf_cptac %>%
  filter_cn_gain_abs() %$%
  cor.test(Protein.Ubiquitination.Adj, Buffering.GeneLevel.Ratio)

ubi_buf_corr <- ubi_buf_cptac %>%
  filter_cn_gain_abs() %>%
  drop_na(Protein.Ubiquitination.Adj, Buffering.GeneLevel.Ratio) %>%
  #filter(Buffering.GeneLevel.Class == "Buffered") %>%
  group_by(Gene.Symbol) %>%
  summarize(Corr = cor(Buffering.GeneLevel.Ratio, Protein.Ubiquitination.Adj),
            Samples = min(sum(!is.na(Buffering.GeneLevel.Ratio)), sum(!is.na(Protein.Ubiquitination.Adj)))) %>%
  filter(Samples > 10)

ubi_buf_cptac %>%
  filter(Gene.Symbol %in% (ubi_buf_corr %>% filter(Corr > 0.4) %>% pull(Gene.Symbol))) %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Ratio, y = Protein.Ubiquitination.Adj, color = Gene.Symbol) +
  geom_point()

ubi_buf_cptac %>%
  filter_cn_gain_abs() %>%
  drop_na(Buffering.GeneLevel.Class) %>%
  filter(Buffering.GeneLevel.Class %in% c("Buffered", "Scaling")) %>%
  signif_violin_plot(Buffering.GeneLevel.Class, Protein.Ubiquitination.Adj,
                     test = t.test)
