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

expression_data_dir <- here(input_data_dir, "P0211")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "P0211")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
p0211_raw <- read.table(here(expression_data_dir, "proteinGroups.txt"), sep="\t", dec=".",
                        header=TRUE, stringsAsFactors=FALSE)
p0211_meta <- read_excel(here(expression_data_dir, "metadata.xlsx"))
uniprot_mapping <- read_parquet(here(output_data_dir, "uniprot_mapping.parquet"))

# === Tidy Dataset ===
state_cols <- c("Identified.In.All", "Identified.In.Some", "Potential.Contaminant", "Reverse",
                "Only.Identified.By.Site", "Unidentified")

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
         ProteinGroup.UniprotIDs, Protein.Expression, all_of(state_cols))

# === Filter & Normalize Dataset ===

p0211_expr_processed <- p0211_expr_tidy %>%
  filter(Identified.In.All) %>%
  filter(Sample.ID != 10) %>%
  mutate(Protein.Expression.Log2 = if_else(is.finite(log2(Protein.Expression)),
                                           log2(Protein.Expression),
                                           NA)) %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  normalize_samples(Sample.Name, Protein.Expression.Log2, ProteinGroup.UniprotIDs,
                    normalized_colname = "Protein.Expression.Normalized") %>%
  # No Batch Effect removal, as Proteomics have been measured in one run
  # and defining batches as cell lines would remove differences between them
  select(-all_of(state_cols))

# === Annotation ===
annotations <- p0211_expr_processed %>%
  select(ProteinGroup.UniprotIDs) %>%
  mutate(Protein.Uniprot.Accession = ProteinGroup.UniprotIDs) %>%
  separate_rows(Protein.Uniprot.Accession, sep = ";") %>%
  mutate(Protein.Uniprot.Accession = str_trim(Protein.Uniprot.Accession)) %>%
  distinct(Protein.Uniprot.Accession, .keep_all = TRUE) %>%
  left_join(y = uniprot_mapping %>% select("Protein.Uniprot.Accession", "Gene.Symbol"),
            by = "Protein.Uniprot.Accession",
            na_matches = "never", relationship = "many-to-one") %>%
  drop_na() %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  get_chromosome_arms() %>%
  filter(Gene.Chromosome %in% (1:22)) %>% # Only use autosomal genes
  mutate(Gene.Chromosome = as.integer(Gene.Chromosome)) %>%
  # Use the first annotated accession listed in protein groups for annotation
  distinct(ProteinGroup.UniprotIDs, .keep_all = TRUE)

p0211_expr_annotated <- p0211_expr_processed %>%
  inner_join(y = annotations, by = "ProteinGroup.UniprotIDs",
             na_matches = "never", relationship = "many-to-one") %>%
  unite("UniqueId", c("Sample.ID", "Protein.Uniprot.Accession"),
        sep = '_', remove = FALSE)

## Add Copy Number Metadata required for DC analysis
id_cols <- c("Sample.ID", "Gene.Symbol")
cn_cols <- c("Gene.CopyNumber", "Gene.Chromosome", "Gene.ChromosomeArm", "Gene.ChromosomeBand",
             "ChromosomeArm.CNA", "Gene.StartPosition", "Gene.EndPosition",
             "CellLine.Ploidy", "CellLine.WGD", "CellLine.AneuploidyScore")

p0211_copy_number <- p0211_expr_annotated %>%
  mutate(Gene.CopyNumber = NA,    # ToDo: Ask for CN data
         CellLine.Ploidy = 2L,
         CellLine.WGD = 0L,
         ChromosomeArm.CNA = case_when(
           CellLine.Name == "RM13" & Gene.Chromosome == 13 ~ -1L,
           CellLine.Name == "Rtr13" & Gene.Chromosome == 13 ~ +1L,
           TRUE ~ 0L
         ),
         CellLine.AneuploidyScore = abs(ChromosomeArm.CNA)*2L + 1) %>% # RPE-1 is pseudo-disomic, trisomy on 10q
  select(all_of(id_cols), all_of(cn_cols))


# === Save Datasets ===
p0211_expr_annotated %>%
  select(-any_of(cn_cols)) %>%
  write_parquet(here(output_data_dir, 'expression_p0211.parquet'), version = "2.6")

p0211_copy_number %>%
  write_parquet(here(output_data_dir, 'copy_number_p0211.parquet'), version = "2.6")

# === Evaluation & Quality Control ===
plot_protein_states <- function(df) {
  df %>%
    select(Sample.Name, all_of(state_cols)) %>%
    pivot_longer(all_of(state_cols), names_to = "Protein.State", values_to = "Count") %>%
    mutate(Protein.State = factor(Protein.State, levels = state_cols)) %>%
    group_by(Sample.Name, Protein.State) %>%
    summarise(Count = sum(Count)) %>%
    ungroup() %>%
    arrange(Protein.State) %>%
    ggplot() +
    aes(x = Sample.Name, y = Count, fill = Protein.State) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_y_continuous(expand = c(0, Inf)) +
    theme(axis.text.x.bottom = element_text(angle = -45, hjust = 0))
}

plot_sample_expr <- function (df, value_col) {
  df %>%
    ggplot() +
    aes(x = {{value_col}}, y = Sample.Name, fill = CellLine.Name) +
    geom_boxplot()
}

plot_protein_states(p0211_expr_tidy) %>%
  save_plot("p0211_protein_states.png")

plot_expr_dist(p0211_expr_processed) %>%
  save_plot("p0211_expression_distribution.png")

p0211_expr_processed %>%
  plot_sample_expr(Protein.Expression.Log2) %>%
  save_plot("p0211_sample_distribution_non-norm.png")

p0211_expr_processed %>%
  plot_sample_expr(Protein.Expression.Normalized) %>%
  save_plot("p0211_sample_distribution_norm.png")

p0211_expr_processed %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Log2) %>%
  plot_pca() %>%
  save_plot("p0211_pca_non-norm.png")

p0211_expr_processed %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  plot_pca() %>%
  save_plot("p0211_pca_norm.png")

mat_log2 <- p0211_expr_processed %>%
  select(Sample.Name, ProteinGroup.UniprotIDs, Protein.Expression.Log2) %>%
  pivot_wider(names_from = Sample.Name, values_from = Protein.Expression.Log2) %>%
  column_to_rownames(var = "ProteinGroup.UniprotIDs")

png(here(plots_dir, "p0211_heatmap_non-norm.png"),
      width = 300, height = 300, units = "mm", res = 300)
pheatmap(mat_log2, show_rownames = FALSE,
         scale = "row", na_col = "black", cluster_cols = TRUE, cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(15))
dev.off()

mat_norm <- p0211_expr_processed %>%
  select(Sample.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  pivot_wider(names_from = Sample.Name, values_from = Protein.Expression.Normalized) %>%
  column_to_rownames(var = "ProteinGroup.UniprotIDs")

png(here(plots_dir, "p0211_heatmap_norm.png"),
      width = 300, height = 300, units = "mm", res = 300)
pheatmap(mat_norm, show_rownames = FALSE,
         scale = "row", na_col = "black", cluster_cols = TRUE, cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(15))
dev.off()


p0211_expr_baseline <- p0211_expr_annotated %>%
  filter(CellLine.Name == "RPE1") %>%
  group_by(Gene.Symbol) %>%
  summarize(Protein.Expression.Baseline = mean(Protein.Expression.Normalized, na.rm = TRUE))

p0211_expr_log2fc <- p0211_expr_annotated %>%
  left_join(y = p0211_expr_baseline, by = "Gene.Symbol",
            unmatched = "error", na_matches = "never", relationship = "many-to-one") %>%
  group_by(Gene.Symbol, Sample.Name, Sample.ID, CellLine.Name, Gene.Chromosome, Gene.StartPosition) %>%
  summarize(Log2FC = Protein.Expression.Normalized - Protein.Expression.Baseline) %>%
  ungroup()

png(here(plots_dir, "p0211_log2fc_chr.png"),
      width = 300, height = 100, units = "mm", res = 300)
p0211_expr_log2fc %>%
  bidirectional_heatmap(Log2FC, Sample.Name, Gene.Chromosome,
                        transpose = TRUE, cluster_rows = TRUE)
dev.off()

for (sample in unique(p0211_expr_log2fc$Sample.ID)) {
  title <- (p0211_expr_log2fc %>%
    filter(Sample.ID == sample) %>%
    distinct(Sample.Name))$Sample.Name

  p0211_expr_log2fc %>%
    filter(Sample.ID == sample) %>%
    bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                          highlight_buckets = 13,
                          threshold_low = log2(1) - log(2),
                          threshold_high = log2(3) - log2(2),
                          x_lab = "Chromosome & Gene Position", title = title) %>%
    save_plot(paste0("p0211_log2fc_chr_sample", sample, ".png"), height = 100)
}