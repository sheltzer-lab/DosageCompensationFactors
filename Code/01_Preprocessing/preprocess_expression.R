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


expression_data_dir <- here(external_data_dir, "Expression")
copynumber_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
## Load dataset from Goncalves et al. (DOI: https://doi.org/10.6084/m9.figshare.19345397.v1)
# ToDo: evaluate use of dataset containing 8498 proteins per cell line
procan_expr <- read_excel(here(expression_data_dir, "ProCan-DepMapSanger_protein_matrix_6692_averaged.xlsx"))
## Load Proteomics dataset from DepMap
depmap_expr <- read_csv_arrow(here(expression_data_dir, "Broad-DepMap-Proteomics.csv"))
### Loaad cell line model file from DepMap
df_celllinenames <- read_csv_arrow(here(copynumber_data_dir, "Model.csv")) %>%
  select(ModelID, CellLineName) %>%
  rename(CellLine.DepMapModelId = "ModelID",
         CellLine.Name = "CellLineName")

# === Tidy Datasets ===

procan_expr_tidy <- procan_expr %>%
  pivot_longer(everything() & !Project_Identifier,
               names_to = "Protein", values_to = "Protein.Expression.Log2") %>%
  separate_wider_delim(Protein,
           delim = ";",
           names = c("Protein.Uniprot.Accession", "Protein.Uniprot.Id")) %>%
  separate_wider_delim(Project_Identifier,
                       delim = ";",
                       names = c("CellLine.SangerModelId", "CellLine.Name")) %>%
  unite("UniqueId", c("CellLine.SangerModelId", "Protein.Uniprot.Accession"),
        sep = '_', remove = FALSE)

depmap_expr_tidy <- depmap_expr %>%
  rename(CellLine.DepMapModelId = 1) %>%
  pivot_longer(everything() & !CellLine.DepMapModelId,
               names_to = "Protein", values_to = "Protein.Expression.Log2") %>%
  # ToDo: Find more elegant way
  separate_wider_delim(Protein,
           delim = " ",
           names = c("Gene.Symbol", "Protein.Uniprot.Accession")) %>%
  mutate(Protein.Uniprot.Accession = gsub("[()]", "", Protein.Uniprot.Accession)) %>%
  unite("UniqueId", c("CellLine.DepMapModelId", "Gene.Symbol", "Protein.Uniprot.Accession"),
        sep = '_', remove = FALSE) %>%
  # This requirement is odd
  unite("UniqueProtId", c("Gene.Symbol", "Protein.Uniprot.Accession"),
        sep = '_', remove = FALSE)


# === Preprocess Datasets ===

normalize_celllines <- function(df, cellline_col, value_col, group_col, normalized_colname = "Normalized") {
  df_pre <- data.frame(df)

  df %>%
    assert_rows(col_concat, is_uniq, {{ cellline_col }}, {{ group_col }}) %>%
    select({{ cellline_col }}, {{ value_col }}, {{ group_col }}) %>%
    pivot_wider(names_from = {{ cellline_col }}, values_from = {{ value_col }}) %>%
    tibble::column_to_rownames(var = quo_name(enquo(group_col))) %>%
    select(where(is.numeric)) %>%
    normalizeBetweenArrays(method = "cyclicloess", cyclic.method = "fast") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = quo_name(enquo(group_col))) %>%
    pivot_longer(everything() & !{{ group_col }},
                 names_to = quo_name(enquo(cellline_col)), values_to = normalized_colname) %>%
    inner_join(y = df_pre,
               by = c(quo_name(enquo(cellline_col)), quo_name(enquo(group_col))),
               relationship = "one-to-one", na_matches = "never")
}

remove_noisefloor <- function(df, value_col, percentile_cutoff = 0.0001) {
  # ToDo: Determine Noise Floor by Signal-To-Noise-Ratio
  noise_floor <- quantile(df[[quo_name(enquo(value_col))]],
                          probs = percentile_cutoff, na.rm = TRUE)[[1]]

  df %>%
    mutate(!!enquo(value_col) := replace({{ value_col }}, {{ value_col }} < noise_floor, NA))
}

# ToDo: Evaluate if applying batch effect correction unaveraged dataset and then averaging it is a better approach
procan_expr_processed <- procan_expr_tidy %>%
  filter(str_detect(Protein.Uniprot.Id, "HUMAN")) %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  # ToDo: Reconsider standardization
  # mutate_at(c('Protein.Expression.Log2'), ~(scale(.) %>% as.vector)) %>%
  normalize_celllines(CellLine.SangerModelId, Protein.Expression.Log2, Protein.Uniprot.Accession,
                      normalized_colname = "Protein.Expression.Normalized") %>%
  mapIds("UNIPROT", "SYMBOL",
         "Protein.Uniprot.Accession", "Gene.Symbol") %>%
  mutate(Dataset = "ProCan")

depmap_expr_processed <- depmap_expr_tidy %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  # ToDo: Reconsider standardization
  # mutate_at(c('Protein.Expression.Log2'), ~(scale(.) %>% as.vector)) %>%
  normalize_celllines(CellLine.DepMapModelId, Protein.Expression.Log2, UniqueProtId,
                      normalized_colname = "Protein.Expression.Normalized") %>%
  left_join(y = df_celllinenames, by = "CellLine.DepMapModelId",
               relationship = "many-to-one", na_matches = "never") %>%
  select(-UniqueProtId) %>%
  mutate(Dataset = "DepMap")


# === Create combined datasets ===
## Unmatched
expr_combined <- procan_expr_processed %>%
  bind_rows(depmap_expr_processed)

## Matched Cell Lines
common_celllines <- intersect(unique(procan_expr_processed$CellLine.Name),
                          unique(depmap_expr_processed$CellLine.Name))

expr_combined_celllines <- expr_combined %>%
  filter(CellLine.Name %in% common_celllines)

## Matched Genes
common_genes <- intersect(unique(procan_expr_processed$Gene.Symbol),
                          unique(depmap_expr_processed$Gene.Symbol))

expr_combined_genes <- expr_combined %>%
  filter(Gene.Symbol %in% common_genes)

## Matched Genes per matched Cell Line
expr_matched_procan <- procan_expr_processed %>%
  semi_join(y = depmap_expr_processed,
            by = c("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession"),
            na_matches = "never")

expr_matched_depmap <- depmap_expr_processed %>%
  semi_join(y = procan_expr_processed,
            by = c("CellLine.Name", "Gene.Symbol", "Protein.Uniprot.Accession"),
            na_matches = "never")

expr_matched <- expr_matched_procan %>%
  bind_rows(expr_matched_depmap)

## Matched Genes per matched Cell Line (re-normalized)
expr_matched_renorm <- expr_matched %>%
  select(-Protein.Expression.Normalized) %>%
  mutate(CellLine.Name.Unique = paste(CellLine.Name, Dataset, sep = "_")) %>%
  normalize_celllines(CellLine.Name.Unique, Protein.Expression.Log2, Protein.Uniprot.Accession,
                      normalized_colname = "Protein.Expression.Normalized") %>%
  select(-CellLine.Name.Unique)

# === Quality Control ===
plot_expr_dist <- function(df) {
  df %>%
    ggplot() +
    geom_density(na.rm = TRUE, aes(Protein.Expression.Log2, color = "Non-Normalized"), show.legend = TRUE) +
    geom_density(na.rm = TRUE, aes(Protein.Expression.Normalized, color = "Normalized"), show.legend = TRUE) +
    xlab("Protein Expression")
}

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