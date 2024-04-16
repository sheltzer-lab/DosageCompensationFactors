library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(skimr)
library(openxlsx)
library(gprofiler2)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Screens", "Enrichment")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))

# === Identify proteins with significant expression differences between cell lines with high and low buffering ===

diff_exp <- cellline_buf_procan %>%
  split_by_quantiles(Buffering.CellLine.Ratio, target_group_col = "CellLine.Buffering.Group") %>%
  inner_join(y = expr_buf_procan, by = "CellLine.Name", relationship = "one-to-many", na_matches = "never") %>%
  select(CellLine.Buffering.Group, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, CellLine.Buffering.Group, Protein.Expression.Normalized,
                          groups = c("Low", "High"))

## Volcano Plots
volcano_plot <- diff_exp %>%
  mutate(Label = if_else(abs(Log2FC) > log2fc_threshold & TTest.p.adj < p_threshold,
                         Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, TTest.p.adj, Label, NULL) %>%
  save_plot("volcano_plot.png")

# === Enrichment Analysis ===
# TODO: Evaluate gprofiler against WebGestalt
enrichment_analysis <- function(genes, ordered = TRUE, p_thresh = p_threshold) {
  gost(query = genes,
       organism = "hsapiens", ordered_query = ordered,
       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
       measure_underrepresentation = FALSE, evcodes = FALSE,
       user_threshold = p_thresh, correction_method = "g_SCS",
       domain_scope = "annotated", custom_bg = NULL,
       numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
}

enrichment_up <- diff_exp %>%
  filter(Log2FC > log2fc_threshold & TTest.p.adj < p_threshold) %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol) %>%
  pull(Gene.Symbol) %>%
  enrichment_analysis()

enrichment_down <- diff_exp %>%
  filter(Log2FC < -log2fc_threshold & TTest.p.adj < p_threshold) %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol) %>%
  pull(Gene.Symbol) %>%
  enrichment_analysis()

#gostplot(enrichment_up, capped = TRUE, interactive = TRUE)
#gostplot(enrichment_down, capped = TRUE, interactive = TRUE)

plot_terms <- function(enrichment, selected_sources = c("CORUM", "KEGG", "REAC", "GO:MF", "GO:BP"), n = 20) {
  enrichment$result %>%
    filter(p_value < p_threshold) %>%
    filter(source %in% selected_sources) %>%
    slice_min(p_value, n = n) %>%
    mutate(`-log10(p)` = -log10(p_value),
           term_name = fct_reorder(term_name, p_value, .desc = TRUE)) %>%
    vertical_bar_chart(term_name, `-log10(p)`,
                       value_range = c(1, max(.$`-log10(p)`)),
                       line_intercept = 0, bar_label_shift = 0.18, break_steps = 2,
                       category_lab = "Enriched Term", value_lab = "Significance (-log10(p))")
}

 enrichment_up %>%
   plot_terms() %>%
   save_plot("enriched_terms_up.png", width = 300)

 enrichment_down %>%
   plot_terms() %>%
   save_plot("enriched_terms_down.png", width = 300)

