library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(STRINGdb)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "annotation.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
temp_dir <- temp_base_dir
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)

cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))

diff_exp <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet")) %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & !is.na(CancerDriverMode), Gene.Symbol, NA))

# === Volcano Plot Panel ===
color_mapping <- scale_color_manual(values = color_palettes$DiffExp,
                                    na.value = color_palettes$Missing)

panel_volcano <- diff_exp %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, color_mapping) +
  theme(legend.position = "none") +
  labs(y = "-log10(p.adj)")

# === ORA Panel ===
genes_up <- diff_exp %>%
  filter(Significant == "Up") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

genes_down <- diff_exp %>%
  filter(Significant == "Down") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_up <- genes_up %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("REAC", "GO:MF", "CORUM"), custom_color = color_palettes$DiffExp["Up"])

ora_down <- genes_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("KEGG", "GO:MF", "CORUM"), custom_color = color_palettes$DiffExp["Down"])

# === STRING Panel ===
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700,
                          network_type = "full", input_directory = temp_dir)

string_up <- genes_up %>%
  create_string_network(Gene.Symbol, Log2FC, string_db, min_score = 900)

string_down <- genes_down %>%
  create_string_network(Gene.Symbol, Log2FC, string_db, min_score = 900)

# Use link in string_up & string_down to create plot manually

# === Combine Panels into Figure ===
figure5 <- cowplot::plot_grid(panel_volcano, ora_down, ora_up,
                              labels = c("A", "B", "C"),
                              rel_widths = c(1, 1, 1),
                              nrow = 1, ncol = 3)

cairo_pdf(here(plots_dir, "figure05_01.pdf"), width = 13, height = 5)
figure5
dev.off()
