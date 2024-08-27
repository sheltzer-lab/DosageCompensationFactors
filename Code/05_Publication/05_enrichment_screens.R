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

cancer_genes <- vroom::vroom(here(external_data_dir, "cancerGeneList.tsv")) %>%
  rename(Occurrences = 9,
         Gene.Symbol = 1) %>%
  filter(Occurrences >= 2) %>%
  mutate(CancerDriverMode = case_when(`Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "Yes" ~ "OG/TSG",
                                      `Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "No" ~ "OG",
                                      `Is Oncogene` == "No" & `Is Tumor Suppressor Gene` == "Yes" ~ "TSG",
                                      TRUE ~ "Unknown")) %>%
  select(Gene.Symbol, CancerDriverMode) %>%
  updateGeneSymbols()

diff_exp <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet")) %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Regulation = case_when(Log2FC > log2fc_threshold & TTest.p.adj < p_threshold ~ "Up",
                                Log2FC < -log2fc_threshold & TTest.p.adj < p_threshold ~ "Down",
                                TRUE ~ NA),
         Label = if_else(!is.na(Regulation) & !is.na(CancerDriverMode), Gene.Symbol, NA))


# === Volcano Plot Panel ===
color_mapping <- scale_color_manual(values = c(Down = bidirectional_color_pal[1],
                                               Up = bidirectional_color_pal[5]),
                                    na.value = color_palettes$Missing)

panel_volcano <- diff_exp %>%
  plot_volcano(Log2FC, TTest.p.adj, Label, Regulation, color_mapping) +
  theme(legend.position = "none")

# === ORA Panel ===
genes_up <- diff_exp %>%
  filter(Regulation == "Up") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

genes_down <- diff_exp %>%
  filter(Regulation == "Down") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_up <- genes_up %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms(selected_sources = c("REAC", "GO:MF", "CORUM"), string_trunc = 25)

ora_down <- genes_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms(selected_sources = c("KEGG", "GO:MF", "CORUM"), string_trunc = 25)

# === STRING Panel ===
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700,
                          network_type = "full", input_directory = temp_dir)

string_up <- genes_up %>%
  create_string_network(Gene.Symbol, Log2FC, string_db, min_score = 900)

string_down <- genes_down %>%
  create_string_network(Gene.Symbol, Log2FC, string_db, min_score = 900)

# Use link in string_up & string_down to create plot manually

# === Combine Panels into Figure ===
figure5 <- cowplot::plot_grid(panel_volcano, ora_up, ora_down,
                              labels = c("A", "B", "C"),
                              rel_widths = c(1.2, 1, 1),
                              nrow = 1, ncol = 3)

cairo_pdf(here(plots_dir, "figure05_01.pdf"), width = 15, height = 5)
figure5
dev.off()
