library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "preprocessing.R"))

expression_data_dir <- here(external_data_dir, "Expression")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "Expression")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
## Load dataset from Goncalves et al. (DOI: https://doi.org/10.6084/m9.figshare.19345397.v1)
# ToDo: evaluate use of dataset containing 8498 proteins per cell line
procan_expr <- read.table(here(expression_data_dir, "ProCan-DepMapSanger_protein_matrix_6692_averaged.txt"),
                          sep="\t", dec=".", header=FALSE, stringsAsFactors=FALSE)

## Load Proteomics dataset from DepMap
depmap_expr <- read_csv_arrow(here(expression_data_dir, "Broad-DepMap-Proteomics.csv"))

df_celllines <- read_parquet(here(output_data_dir, "celllines.parquet"))
uniprot_mapping <- read_parquet(here(output_data_dir, "uniprot_mapping.parquet"))
df_rep_filtered <- read_parquet(here(output_data_dir, "reproducibility_ranks_filtered.parquet"))

# === Tidy Datasets ===

procan_expr_tidy <- procan_expr %>%
  janitor::row_to_names(1) %>%
  pivot_longer(everything() & !Project_Identifier,
               names_to = "Protein", values_to = "Protein.Expression.Log2") %>%
  mutate(Protein.Expression.Log2 = as.numeric(Protein.Expression.Log2)) %>%
  separate_wider_delim(Protein,
           delim = ";",
           names = c("Protein.Uniprot.Accession", "Protein.Uniprot.Id")) %>%
  separate_wider_delim(Project_Identifier,
                       delim = ";",
                       names = c("CellLine.SangerModelId", "CellLine.Name")) %>%
  # Add cell line identifiers
  select(-CellLine.Name) %>%
  inner_join(y = df_celllines, by = "CellLine.SangerModelId",
             relationship = "many-to-one", na_matches = "never") %>%
  left_join(y = uniprot_mapping %>% select("Protein.Uniprot.Accession", "Gene.Symbol"),
            by = "Protein.Uniprot.Accession",
            na_matches = "never", relationship = "many-to-one") %>%
  updateGeneSymbols() %>%
  unite("UniqueId", c("CellLine.SangerModelId", "Protein.Uniprot.Accession"),
        sep = '_', remove = FALSE) %>%
  mutate(Dataset = "ProCan",
         Model.Type = "Cell Line")

depmap_expr_tidy <- depmap_expr %>%
  rename(CellLine.DepMapModelId = 1) %>%
  pivot_longer(everything() & !CellLine.DepMapModelId,
               names_to = "Protein", values_to = "Protein.Expression.Log2") %>%
  # ToDo: Find more elegant way
  separate_wider_delim(Protein,
           delim = " ",
           names = c("Gene.Symbol", "Protein.Uniprot.Accession")) %>%
  mutate(Protein.Uniprot.Accession = gsub("[()]", "", Protein.Uniprot.Accession)) %>%
  # Add cell line identifiers
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId",
             relationship = "many-to-one", na_matches = "never") %>%
  select(-Gene.Symbol) %>%
  left_join(y = uniprot_mapping %>% select("Protein.Uniprot.Accession", "Gene.Symbol"),
            by = "Protein.Uniprot.Accession",
            na_matches = "never", relationship = "many-to-one") %>%
  updateGeneSymbols() %>%
  unite("UniqueId", c("CellLine.DepMapModelId", "Gene.Symbol", "Protein.Uniprot.Accession"),
        sep = '_', remove = FALSE) %>%
  mutate(Dataset = "DepMap",
         Model.Type = "Cell Line")

# === Preprocess Datasets ===
# ToDo: Evaluate if applying batch effect correction unaveraged dataset and then averaging it is a better approach
procan_expr_processed <- procan_expr_tidy %>%
  filter(str_detect(Protein.Uniprot.Id, "HUMAN")) %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  semi_join(df_rep_filtered, by = "Gene.Symbol") %>%  # Remove proteins with low reproducibilty across datasets
  # ToDo: Reconsider standardization
  # mutate_at(c('Protein.Expression.Log2'), ~(scale(.) %>% as.vector)) %>%
  normalize_samples(CellLine.SangerModelId, Protein.Expression.Log2, Protein.Uniprot.Accession,
                    normalized_colname = "Protein.Expression.Normalized")

depmap_expr_processed <- depmap_expr_tidy %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  semi_join(df_rep_filtered, by = "Gene.Symbol") %>%  # Remove proteins with low reproducibilty across datasets
  group_by(Model.ID, CellLine.SangerModelId, CellLine.DepMapModelId, CellLine.Name, Dataset,
           UniqueId, Gene.Symbol, Protein.Uniprot.Accession) %>%
  summarize(Protein.Expression.Log2 = mean(Protein.Expression.Log2, na.rm = TRUE)) %>%
  ungroup() %>%
  # ToDo: Reconsider standardization
  # mutate_at(c('Protein.Expression.Log2'), ~(scale(.) %>% as.vector)) %>%
  normalize_samples(CellLine.DepMapModelId, Protein.Expression.Log2, Protein.Uniprot.Accession,
                    normalized_colname = "Protein.Expression.Normalized")


# === Create combined datasets ===
## ToDo: Expression baseline & averages are calculated across datasets (problematic, different means)
## ToDo: Consider either renormalizing or merging after DC analysis
## Unmatched
expr_combined <- procan_expr_processed %>%
  bind_rows(depmap_expr_processed)

## Matched Cell Lines
common_celllines <- intersect(unique(procan_expr_processed$Model.ID),
                          unique(depmap_expr_processed$Model.ID))

expr_combined_celllines <- expr_combined %>%
  filter(Model.ID %in% common_celllines)

## Matched Genes
common_genes <- intersect(unique(procan_expr_processed$Gene.Symbol),
                          unique(depmap_expr_processed$Gene.Symbol))

expr_combined_genes <- expr_combined %>%
  filter(Gene.Symbol %in% common_genes)

## Matched Genes per matched Cell Line
expr_matched <- match_datasets(procan_expr_processed, depmap_expr_processed)

## Matched Genes per matched Cell Line (re-normalized)
expr_matched_renorm <- expr_matched %>%
  select(-Protein.Expression.Normalized) %>%
  mutate(CellLine.Name.Unique = paste(Model.ID, Dataset, sep = "_")) %>%
  normalize_samples(CellLine.Name.Unique, Protein.Expression.Log2, Protein.Uniprot.Accession,
                    normalized_colname = "Protein.Expression.Normalized") %>%
  select(-CellLine.Name.Unique)

# === Quality Control ===
## Check correlation between datasets
corr_combined <- dataset_correlation(expr_combined,
                                     Dataset, Protein.Expression.Log2,
                                     "Combined")
corr_matched <- dataset_correlation(expr_matched,
                                    Dataset, Protein.Expression.Log2,
                                    "Matched")
corr_matched_norm <- dataset_correlation(expr_matched,
                                         Dataset, Protein.Expression.Normalized,
                                         "Matched (Normalized)")
corr_matched_renorm <- dataset_correlation(expr_matched_renorm,
                                           Dataset, Protein.Expression.Normalized,
                                           "Matched (Renormalized)")

### Note: Median correlation between ProCan and DepMap lower (0.09) than reported (0.44)
corr_combined %>%
  bind_rows(corr_matched) %>%
  bind_rows(corr_matched_norm) %>%
  bind_rows(corr_matched_renorm) %>%
  jittered_boxplot(Comparison, Correlation) %>%
  save_plot("proteomics_dataset_correlation.png", height = 100)


## Plot distribution
procan_expr_dist <- procan_expr_processed %>%
  plot_expr_dist() %>%
  save_plot("expression_distributions_procan.png", height = 200, width = 300)

depmap_expr_dist <- depmap_expr_processed %>%
  plot_expr_dist() %>%
  save_plot("expression_distributions_depmap.png", height = 200, width = 300)

combined_expr_dist <- expr_combined %>%
  plot_expr_dist() %>%
  save_plot("expression_distributions_combined.png", height = 200, width = 300)

combined_celllines_expr_dist <- expr_combined_celllines %>%
  plot_expr_dist() %>%
  save_plot("expression_distributions_combined_celllines.png", height = 200, width = 300)

combined_genes_expr_dist <- expr_combined_genes %>%
  plot_expr_dist() %>%
  save_plot("expression_distributions_combined_genes.png", height = 200, width = 300)

matched_expr_dist <- expr_matched %>%
  plot_expr_dist() %>%
  save_plot("expression_distributions_matched.png", height = 200, width = 300)

matched_renorm_expr_dist <- expr_matched_renorm %>%
  plot_expr_dist() %>%
  save_plot("expression_distributions_matched_renorm.png", height = 200, width = 300)

## ToDo: Create & Write report

# === Write processed datasets to disk ===

write_parquet(procan_expr_processed, here(output_data_dir, 'expression_procan.parquet'),
              version = "2.6")
write_parquet(depmap_expr_processed, here(output_data_dir, 'expression_depmap.parquet'),
              version = "2.6")

write_parquet(expr_combined, here(output_data_dir, 'expression_combined.parquet'),
              version = "2.6")
write_parquet(expr_combined_celllines, here(output_data_dir, 'expression_combined_celllines.parquet'),
              version = "2.6")
write_parquet(expr_combined_genes, here(output_data_dir, 'expression_combined_genes.parquet'),
              version = "2.6")

write_parquet(expr_matched, here(output_data_dir, 'expression_matched.parquet'),
              version = "2.6")
write_parquet(expr_matched_renorm, here(output_data_dir, 'expression_matched_renorm.parquet'),
              version = "2.6")
