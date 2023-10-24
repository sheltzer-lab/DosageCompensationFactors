library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)
library(tibble)
library(pheatmap)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "preprocessing.R"))

expression_data_dir <- here(external_data_dir, "Expression", "P0211")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
p0211_raw <- read.table(here(expression_data_dir, "proteinGroups.txt"), sep="\t", dec=".",
                        header=TRUE, stringsAsFactors=FALSE)
p0211_meta <- read_excel(here(expression_data_dir, "metadata.xlsx"))

# === Tidy Dataset ===

p0211_expr_tidy <- p0211_raw %>%
  select(starts_with("Reporter.intensity.corrected"),
         Potential.contaminant, Only.identified.by.site, Reverse,
         Majority.protein.IDs) %>%
  pivot_longer(starts_with("Reporter.intensity.corrected"),
               names_to = "raw_name", values_to = "Protein.Expression") %>%
  inner_join(y = p0211_meta, by = "raw_name",
             unmatched = "error", relationship = "many-to-one") %>%
  unite("CellLine.CustomId", c("project_id", "cell_line"), remove = FALSE, sep = "_") %>%
  mutate(Sample.ID = as.integer(sample_number),
         Sample.Name = sample_name) %>%
  rename(CellLine.Name = cell_line,
         CellLine.Replicate = replicate,
         ProteinGroup.UniprotIDs = Majority.protein.IDs) %>%
  mutate(Reverse = Reverse == "+" | grepl("REVs", ProteinGroup.UniprotIDs),
         Potential.Contaminant = Potential.contaminant == "+",
         Only.Identified.By.Site = Only.identified.by.site == "+") %>%
  group_by(ProteinGroup.UniprotIDs) %>%
  mutate(Unidentified = all(Protein.Expression == 0) &
           !Potential.Contaminant & !Reverse & !Only.Identified.By.Site,
         Identified.In.All = all(Protein.Expression != 0) &
           !Potential.Contaminant & !Reverse & !Only.Identified.By.Site,
         Identified.In.Some = any(Protein.Expression != 0) &
           !Potential.Contaminant & !Reverse & !Only.Identified.By.Site & !Identified.In.All) %>%
  ungroup() %>%
  select(Sample.ID, Sample.Name, CellLine.CustomId, CellLine.Name, CellLine.Replicate,
         ProteinGroup.UniprotIDs, Protein.Expression,
         Potential.Contaminant, Only.Identified.By.Site, Reverse,
         Unidentified, Identified.In.All, Identified.In.Some)

# === Filter & Normalize Dataset ===

p0211_expr_processed <- p0211_expr_tidy %>%
  filter(Identified.In.All) %>%
  filter(Sample.ID != 10) %>%
  mutate(Protein.Expression.Log2 = if_else(is.finite(log2(Protein.Expression)),
                                           log2(Protein.Expression),
                                           NA)) %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  normalize_samples(Sample.Name, Protein.Expression.Log2, ProteinGroup.UniprotIDs,
                    normalized_colname = "Protein.Expression.Normalized")

# === Evaluation & Quality Control ===

plot_protein_states <- function(df) {
  state_cols <- c("Potential.Contaminant", "Only.Identified.By.Site", "Reverse",
                 "Unidentified", "Identified.In.All", "Identified.In.Some")

  df %>%
    select(Sample.Name, all_of(state_cols)) %>%
    pivot_longer(all_of(state_cols), names_to = "Protein.State", values_to = "Count") %>%
    group_by(Sample.Name, Protein.State) %>%
    summarise(Count = sum(Count)) %>%
    ungroup() %>%
    ggplot() +
    aes(x = Sample.Name, y = Count, fill = Protein.State) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, Inf)) +
    theme(axis.text.x.bottom = element_text(angle = -45, hjust = 0))
}

plot_sample_expr <- function (df, value_col) {
  df %>%
    ggplot() +
    aes(x = {{value_col}}, y = Sample.Name, fill = CellLine.Name) +
    geom_boxplot()
}

plot_protein_states(p0211_expr_tidy)
plot_expr_dist(p0211_expr_processed)
p0211_expr_processed %>%
  plot_sample_expr(Protein.Expression.Log2)
p0211_expr_processed %>%
  plot_sample_expr(Protein.Expression.Normalized)

p0211_expr_processed %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Log2) %>%
  plot_pca()

p0211_expr_processed %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  plot_pca()

mat_log2 <- p0211_expr_processed %>%
  select(Sample.Name, ProteinGroup.UniprotIDs, Protein.Expression.Log2) %>%
  pivot_wider(names_from = Sample.Name, values_from = Protein.Expression.Log2) %>%
  column_to_rownames(var = "ProteinGroup.UniprotIDs")

pheatmap(mat_log2, show_rownames = F,
         scale = "row", na_col = "black", cluster_cols = T, cluster_rows = F,
         color = colorRampPalette(c("blue", "white", "red"))(15))

mat_norm <- p0211_expr_processed %>%
  select(Sample.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  pivot_wider(names_from = Sample.Name, values_from = Protein.Expression.Normalized) %>%
  column_to_rownames(var = "ProteinGroup.UniprotIDs")

pheatmap(mat_norm, show_rownames = F,
         scale = "row", na_col = "black", cluster_cols = T, cluster_rows = F,
         color = colorRampPalette(c("blue", "white", "red"))(15))