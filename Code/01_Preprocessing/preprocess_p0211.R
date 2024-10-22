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
library(scales)


here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "preprocessing.R"))

expression_data_dir <- here(input_data_dir, "P0211")
chunduri_data_dir <- here(external_data_dir, "Expression", "Chunduri")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "P0211")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
p0211_raw <- read.table(here(expression_data_dir, "proteinGroups.txt"), sep="\t", dec=".",
                        header=TRUE, stringsAsFactors=FALSE)
p0211_meta <- read_excel(here(expression_data_dir, "metadata.xlsx"))
# Chunduri et al. 2021, DOI: https://doi.org/10.1038/s41467-021-25288-x
chunduri_raw <- read.table(here(chunduri_data_dir, "proteinGroups.txt"), sep="\t", dec=".",
                           header=TRUE, stringsAsFactors=FALSE)
uniprot_mapping <- read_parquet(here(output_data_dir, "uniprot_mapping.parquet"))
df_rep_filtered <- read_parquet(here(output_data_dir, "reproducibility_ranks_filtered.parquet"))


# TODO: === Add Chunduri to DatasetReferences.md ===
# TODO: === Rename file to preprocess_engineered.R ===

# === Tidy Dataset ===
state_cols <- c("Identified.In.All", "Identified.In.Some", "Potential.Contaminant", "Reverse",
                "Only.Identified.By.Site", "Unidentified")

## P0211
p0211_expr_tidy <- p0211_raw %>%
  select(starts_with("Reporter.intensity.corrected"),
         Potential.contaminant, Only.identified.by.site, Reverse,
         Majority.protein.IDs) %>%
  pivot_longer(starts_with("Reporter.intensity.corrected"),
               names_to = "raw_name", values_to = "Protein.Expression") %>%
  inner_join(y = p0211_meta, by = "raw_name",
             unmatched = "error", relationship = "many-to-one") %>%
  unite("Model.ID", c("project_id", "cell_line"), remove = FALSE, sep = "_") %>%
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
  select(Sample.ID, Sample.Name, Model.ID, CellLine.Name, CellLine.Replicate,
         ProteinGroup.UniprotIDs, Protein.Expression, all_of(state_cols)) %>%
    mutate(Dataset = "P0211",
           Model.Type = "Cell Line")

## Chunduri et al.
chunduri_meta <- data.frame(
  raw_name = paste(rep("Reporter.intensity.corrected", 6), 1:6, sep = "."),
  Sample.Tag = c("TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"),
  CellLine.Name = c("RPE1 p53 KO", "RM X", "RM 10;18", "RM 13", "RPE1 p53 KD", "RM 19p")
)

chunduri_tidy <- chunduri_raw %>%
  select(starts_with("Reporter.intensity.corrected"),
         Potential.contaminant, Only.identified.by.site, Reverse,
         Majority.protein.IDs) %>%
  select(contains(".REP_"),
         Potential.contaminant, Only.identified.by.site, Reverse,
         Majority.protein.IDs) %>%
  pivot_longer(starts_with("Reporter.intensity.corrected"),
               names_to = "raw_name", values_to = "Protein.Expression") %>%
  separate_wider_delim(raw_name, ".REP_", names = c("raw_name", "CellLine.Replicate")) %>%
  mutate(Dataset = "Chunduri",
         Model.Type = "Cell Line",
         CellLine.Replicate = as.integer(CellLine.Replicate)) %>%
  inner_join(y = chunduri_meta, by = "raw_name",
             unmatched = "error", relationship = "many-to-one") %>%
  unite("Model.ID", c("Dataset", "CellLine.Name", "CellLine.Replicate"), remove = FALSE, sep = "_") %>%
  unite("Sample.Name", c("CellLine.Name", "CellLine.Replicate"), remove = FALSE, sep = "_") %>%
  mutate(Sample.ID = as.integer(as.factor(Sample.Name))) %>%
  rename(ProteinGroup.UniprotIDs = Majority.protein.IDs) %>%
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
  select(Sample.ID, Sample.Name, Model.ID, CellLine.Name, CellLine.Replicate,
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

chunduri_processed <- chunduri_tidy %>%
  filter(Identified.In.All) %>%
  filter(CellLine.Name %in% c("RPE1 p53 KO", "RM 10;18", "RM 13", "RM 19p")) %>%
  mutate(Protein.Expression.Log2 = if_else(is.finite(log2(Protein.Expression)),
                                           log2(Protein.Expression),
                                           NA)) %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  normalize_samples(Sample.Name, Protein.Expression.Log2, ProteinGroup.UniprotIDs,
                    normalized_colname = "Protein.Expression.Normalized") %>%
  select(-all_of(state_cols))

# TODO: Add batch effect correction, as PCA shows clusters per replicate not per cell line

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
  updateGeneSymbols() %>%
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
        sep = '_', remove = FALSE) %>%
  # TODO: Executed before normalization in CCLE and ProCan
  semi_join(df_rep_filtered, by = "Gene.Symbol") # Remove proteins with low reproducibilty across datasets

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

## P0211
prot_states_p0211 <- plot_protein_states(p0211_expr_tidy) %>%
  save_plot("p0211_protein_states.png")

expr_dist_p0211 <- plot_expr_dist(p0211_expr_annotated) %>%
  save_plot("p0211_expression_distribution.png")

sample_dist_pre_p0211 <- p0211_expr_annotated %>%
  plot_sample_expr(Protein.Expression.Log2) %>%
  save_plot("p0211_sample_distribution_non-norm.png")

sample_dist_norm_p0211 <- p0211_expr_annotated %>%
  plot_sample_expr(Protein.Expression.Normalized) %>%
  save_plot("p0211_sample_distribution_norm.png")

pca_pre_p0211 <- p0211_expr_annotated %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Log2) %>%
  plot_pca() %>%
  save_plot("p0211_pca_non-norm.png")

pca_norm_p0211 <- p0211_expr_annotated %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  plot_pca() %>%
  save_plot("p0211_pca_norm.png")

prot_heatmap_pre <- p0211_expr_annotated %>%
  mutate(ScaledExpression = oob_squish(as.vector(scale(Protein.Expression.Log2)), range = c(-2,2))) %>%
  bidirectional_heatmap(ScaledExpression, Sample.Name, ProteinGroup.UniprotIDs,
                        cluster_cols = TRUE, show_rownames = FALSE, palette_length = 11) %>%
  save_plot("p0211_heatmap_non-norm.png", width = 300, height = 300)

prot_heatmap_norm <-p0211_expr_annotated %>%
  mutate(ScaledExpression = oob_squish(as.vector(scale(Protein.Expression.Normalized)), range = c(-2,2))) %>%
  bidirectional_heatmap(ScaledExpression, Sample.Name, Gene.Symbol,
                        cluster_cols = TRUE, show_rownames = FALSE, palette_length = 11) %>%
  save_plot("p0211_heatmap_norm.png", width = 300, height = 300)

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

chr_heatmap <- p0211_expr_log2fc %>%
  bidirectional_heatmap(Log2FC, Sample.Name, Gene.Chromosome,
                        transpose = TRUE, cluster_rows = TRUE) %>%
  save_plot("p0211_log2fc_chr.png", width = 300, height = 100)

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

## Chunduri
prot_states_chunduri <- plot_protein_states(chunduri_tidy) %>%
  save_plot("chunduri_protein_states.png")

expr_dist_chunduri <- plot_expr_dist(chunduri_processed) %>%
  save_plot("chunduri_expression_distribution.png")

sample_dist_pre_chunduri <- chunduri_processed %>%
  plot_sample_expr(Protein.Expression.Log2) %>%
  save_plot("chunduri_sample_distribution_non-norm.png")

sample_dist_norm_chunduri <- chunduri_processed %>%
  plot_sample_expr(Protein.Expression.Normalized) %>%
  save_plot("chunduri_sample_distribution_norm.png")

pca_pre_chunduri <- chunduri_processed %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Log2) %>%
  plot_pca() %>%
  save_plot("chunduri_pca_non-norm.png")

pca_norm_chunduri <- chunduri_processed %>%
  calculate_pca(Sample.Name, CellLine.Name, ProteinGroup.UniprotIDs, Protein.Expression.Normalized) %>%
  plot_pca() %>%
  save_plot("chunduri_pca_norm.png")

# === Publish ===

title <- (p0211_expr_log2fc %>%
    filter(Sample.ID == 7) %>%
    distinct(Sample.Name))$Sample.Name
expr_bucket <- p0211_expr_log2fc %>%
  filter(Sample.ID == 7) %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        highlight_buckets = 13,
                        threshold_low = log2(1) - log(2),
                        threshold_high = log2(3) - log2(2),
                        x_lab = "Chromosome & Gene Position", title = title)

grid1 <- cowplot::plot_grid(pca_norm + theme(legend.position = "none"), expr_bucket,
                            nrow = 1, ncol = 2, rel_widths = c(1, 1.5), labels = c("B", "C"))

grid2 <- cowplot::plot_grid(grid1, chr_heatmap$gtable, nrow = 2, ncol = 1,
                            rel_heights = c(1, 2/3), labels = c("", "D"))

plot_publish <- cowplot::plot_grid(sample_dist_norm + theme(legend.position = "none"), grid2,
                                   nrow = 1, ncol = 2,
                                   rel_widths = c(1, 2), labels = c("A", ""))

cairo_pdf(here(plots_dir, "preprocessing_publish.pdf"), width = 13)
plot_publish
dev.off()
