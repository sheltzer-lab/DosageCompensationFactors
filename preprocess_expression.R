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

source(here("parameters.R"))

expression_data_dir <- here(external_data_dir, "Expression")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
## Load PhosphoSitePlus datasets
procan_expr_avg <- read_excel(here(expression_data_dir, "ProCan-DepMapSanger_protein_matrix_6692_averaged.xlsx"))

# === Tidy Dataset ===

procan_expr_avg_tidy <- procan_expr_avg %>%
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

# === Preprocess Dataset ===

normalize_celllines <- function(df, cellline_col, value_col, group_col, normalized_colname = "Normalized") {
  df_pre <- data.frame(df)

  df_test <- df %>%
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
procan_expr_avg_processed <- procan_expr_avg_tidy %>%
  filter(str_detect(Protein.Uniprot.Id, "HUMAN")) %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  normalize_celllines(CellLine.SangerModelId, Protein.Expression.Log2, Protein.Uniprot.Accession,
                      normalized_colname = "Protein.Expression.Normalized")

# === Quality Control ===

expr_dist <- procan_expr_avg_processed %>%
  ggplot() +
  geom_density(na.rm = TRUE, aes(Protein.Expression.Log2, color = "Non-Normalized"), show.legend = TRUE) +
  geom_density(na.rm = TRUE, aes(Protein.Expression.Normalized, color = "Normalized"), show.legend = TRUE) +
  xlab("Protein Expression")

ggsave(here(plots_dir, "expression_distributions.png"), plot = expr_dist,
       height = 200, width = 300, units = "mm", dpi = 300)

## ToDo: Create & Write report

# === Write Processed dataset to disk ===

write_parquet(procan_expr_avg_processed, here(output_data_dir, 'expression_average.parquet'),
              version = "2.6")