library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(readr)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "analysis.R"))
source(here("Code", "visualization.R"))

expression_data_dir <- here(external_data_dir, "Expression", "CPTAC", "Proteome_BCM_GENCODE_v34_harmonized_v1")
meta_data_dir <- here(expression_data_dir, "README")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "Expression", "CPTAC")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
## Expression
expr_file_list <- list.files(expression_data_dir, include.dirs = FALSE, recursive = FALSE)
expr_file_list <- expr_file_list[grep(".+_proteomics_.+.txt", expr_file_list)]

df_list_expr <- list()

for (filename in expr_file_list) {
  metadata <- sub('\\.txt$', '', filename) %>%
    str_split_fixed("_", 9) %>%
    as.data.frame() %>%
    select(V1, V9)

  df_name <- paste(metadata$V1, metadata$V9, sep = "_")

  df_list_expr[[df_name]] <- read.table(here(expression_data_dir, filename),
                                   sep = "\t", dec = ".", header = FALSE, stringsAsFactors = FALSE) %>%
    janitor::row_to_names(1) %>%
    pivot_longer(everything() & !idx,
                 names_to = "Model.ID", values_to = "Protein.Expression.Log2",
                 names_ptypes = character(), values_ptypes = double(),
                 names_transform = as.character, values_transform = as.numeric) %>%
    mutate(Model.CancerType = metadata$V1,
           Model.SampleType = metadata$V9)
}

## Metadata
meta_filenames <- list.files(meta_data_dir, include.dirs = FALSE, recursive = FALSE)
meta_filenames <- meta_filenames[grep(".+_meta.txt", meta_filenames)]

df_list_meta <- list()

for (filename in meta_filenames) {
  cancer_type <- sub('\\.txt$', '', filename) %>%
    str_split_i("_", 1)

  df_list_meta[[cancer_type]] <- read_delim(here(meta_data_dir, filename), delim = "\t") %>%
    filter(idx != "data_type") %>%
    mutate(Model.CancerType = cancer_type)
}

## Copy Number
gdc_cptac_cnv <- read_parquet(here(output_data_dir, 'gdc_cptac-3_cnv_ascat.parquet'))

## Other
uniprot_mapping <- read_parquet(here(output_data_dir, "uniprot_mapping.parquet"))
df_rep_filtered <- read_parquet(here(output_data_dir, "reproducibility_ranks_filtered.parquet"))

# === Combine & Tidy Datasets ===
## Expression
expr_cptac <- bind_rows(df_list_expr) %>%
  unite("Model.SampleID", c("Model.ID", "Model.SampleType"), sep = '_', remove = FALSE) %>%
  rename(Gene.ENSEMBL.Id = idx) %>%
  ## Remove ENSEMBL version
  mutate(Gene.ENSEMBL.Id = sub("\\.\\d+$", "", Gene.ENSEMBL.Id)) %>%
  ## ID Mapping
  left_join(y = uniprot_mapping %>% select("Gene.ENSEMBL.Id", "Protein.Uniprot.Accession", "Gene.Symbol"),
            by = "Gene.ENSEMBL.Id",
            na_matches = "never", relationship = "many-to-one", multiple = "last") %>% # TODO: Check which UniProt ID is obsolete
  updateGeneSymbols() %>%
  unite("UniqueId", c("Model.SampleID", "Gene.ENSEMBL.Id"), sep = '_', remove = FALSE) %>%
  mutate(Dataset = "CPTAC",
         Model.Type = "Tumor Sample")

## Copy Number
copy_number_cptac <- gdc_cptac_cnv %>%
  rename(Gene.ENSEMBL.Id = gene_id) %>%
  mutate(Gene.CopyNumber = if_else(as.numeric(copy_number) < 10, as.numeric(copy_number), NA), # TODO: Temporary fix
         ## Remove ENSEMBL version
         Gene.ENSEMBL.Id = sub("\\.\\d+$", "", Gene.ENSEMBL.Id)) %>%
  ## ID Mapping
  select(Model.ID, Gene.ENSEMBL.Id, Gene.CopyNumber) %>%
  left_join(y = uniprot_mapping %>% select("Gene.ENSEMBL.Id", "Protein.Uniprot.Accession", "Gene.Symbol"),
            by = "Gene.ENSEMBL.Id",
            na_matches = "never", relationship = "many-to-one", multiple = "last") %>% # TODO: Check which UniProt ID is obsolete
  updateGeneSymbols() %>%
  get_chromosome_arms() %>%
  filter(Gene.Chromosome %in% (1:22)) %>% # Only use autosomal genes
  # TODO: Determine source of quality issues (multiple values per gene and sample)
  drop_na(Gene.CopyNumber) %>%
  group_by(Model.ID, Gene.ENSEMBL.Id) %>%
  summarize(Gene.CopyNumber = max(Gene.CopyNumber, na.rm = TRUE), .groups = "drop")

# === Preprocess Datasets ===
## Expression
expr_cptac_processed <- expr_cptac %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  semi_join(df_rep_filtered, by = "Gene.Symbol") %>%  # Remove proteins with low reproducibilty across datasets
  normalize_samples(Model.SampleID, Protein.Expression.Log2, Gene.ENSEMBL.Id,
                    normalized_colname = "Protein.Expression.Normalized")

## Metadata
metadata_cptac_processed <- bind_rows(df_list_meta) %>%
  rename(Model.ID = idx) %>%
  mutate_all(~type.convert(., as.is = TRUE)) %>%
  mutate(Model.TumorPurity = pmin(WES_purity, WGS_purity, na.rm = TRUE))


# === Save Datasets ===
expr_cptac_processed %>%
  write_parquet(here(output_data_dir, 'expression_cptac.parquet'), version = "2.6")

copy_number_cptac %>%
  write_parquet(here(output_data_dir, 'copy_number_cptac.parquet'), version = "2.6")

metadata_cptac_processed %>%
  write_parquet(here(output_data_dir, 'metadata_cptac.parquet'), version = "2.6")

# === Evaluation ===
expr_dist <- plot_expr_dist(expr_cptac_processed) %>%
  save_plot("cptac_expression_distribution.png")

expr_dist_sample <- expr_cptac_processed %>%
  ggplot() +
  aes(x = Protein.Expression.Normalized, color = Model.SampleType) +
  geom_density()

sample_subset <- expr_cptac_processed %>%
  distinct(Model.SampleID, Model.ID, Model.CancerType, Model.SampleType) %>%
  group_by(Model.CancerType, Model.SampleType) %>%
  slice_sample(n = 2)

sample_dist_pre <- expr_cptac_processed %>%
  filter(Model.SampleID %in% sample_subset$Model.SampleID) %>%
  mutate(Model.SampleID = fct_reorder(Model.SampleID, Model.CancerType)) %>%
  vertical_box_plot(Protein.Expression.Log2, Model.SampleID, Model.CancerType) %>%
  save_plot("cptac_sample_distribution_non-norm.png")

sample_dist_post <- expr_cptac_processed %>%
  filter(Model.SampleID %in% sample_subset$Model.SampleID) %>%
  mutate(Model.SampleID = fct_reorder(Model.SampleID, Model.CancerType)) %>%
  vertical_box_plot(Protein.Expression.Normalized, Model.SampleID, Model.CancerType) %>%
  save_plot("cptac_sample_distribution_norm.png")

types_dist_pre <- expr_cptac_processed %>%
  unite("Model.CancerSampleType", c(Model.CancerType, Model.SampleType), sep = '_', remove = FALSE) %>%
  vertical_box_plot(Protein.Expression.Log2, Model.CancerSampleType, Model.CancerType) %>%
  save_plot("cptac_cancer-type_distribution_non-norm.png")

types_dist_post <- expr_cptac_processed %>%
  unite("Model.CancerSampleType", c(Model.CancerType, Model.SampleType), sep = '_', remove = FALSE) %>%
  vertical_box_plot(Protein.Expression.Normalized, Model.CancerSampleType, Model.CancerType) %>%
  save_plot("cptac_cancer-type_distribution_norm.png")

pca_pre <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, Gene.ENSEMBL.Id, Protein.Expression.Log2) %>%
  plot_pca(color_col = Model.CancerType, label_col = NULL) %>%
  save_plot("cptac_pca_non-norm.png")

pca_pre_sample <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, Gene.ENSEMBL.Id, Protein.Expression.Log2) %>%
  plot_pca(color_col = Model.SampleType, label_col = NULL) %>%
  save_plot("cptac_pca_sample_non-norm.png")

pca_norm <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, Gene.ENSEMBL.Id, Protein.Expression.Normalized) %>%
  plot_pca(color_col = Model.CancerType, label_col = NULL) %>%
  save_plot("cptac_pca_norm.png")

pca_norm <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, Gene.ENSEMBL.Id, Protein.Expression.Normalized) %>%
  plot_pca(color_col = Model.SampleType, label_col = NULL) %>%
  save_plot("cptac_pca_sample_norm.png")

cn_dist <- copy_number_cptac %>%
  ggplot() +
  aes(x = Gene.CopyNumber) +
  geom_bar()

ploidy_dist <- copy_number_cptac %>%
  group_by(Model.ID) %>%
  summarize(Ploidy = mean(Gene.CopyNumber, na.rm = TRUE), .groups = "drop") %>%
  ggplot() +
  aes(x = Ploidy) +
  geom_density()

# Tumor Purity
purity_dist <- metadata_cptac_processed %>%
  ggplot() +
  aes(x = Model.TumorPurity) +
  geom_density() +
  scale_x_continuous(breaks = seq(0, 1, 0.1))

fivenum(metadata_cptac_processed$Model.TumorPurity, na.rm = TRUE)

# WGS and WES ploidy inconsistent
ploidy_cptac <- metadata_cptac_processed %>%
  scatter_plot_reg_corr(WES_ploidy, WGS_ploidy)
