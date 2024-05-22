library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(skimr)
library(openxlsx)
library(gprofiler2)
library(WebGestaltR)
library(STRINGdb)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Screens", "Enrichment")
reports_dir <- reports_base_dir
temp_dir <- temp_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)
dir.create(temp_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))

# === Identify proteins with significant expression differences between cell lines with high and low buffering ===

diff_exp <- cellline_buf_procan %>%
  split_by_3_quantiles(Buffering.CellLine.Ratio, target_group_col = "CellLine.Buffering.Group") %>%
  filter(CellLine.Buffering.Group != "Center") %>%
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
# TODO: GSEA
overrepresentation_analysis <- function(genes, ordered = TRUE, p_thresh = p_threshold, ref_background = NULL,
                                        databases = c("GO:BP", "GO:MF", "KEGG", "REAC", "WP", "CORUM")) {
  gost(query = genes,
       organism = "hsapiens", ordered_query = ordered,
       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
       measure_underrepresentation = FALSE, evcodes = FALSE,
       user_threshold = p_thresh, correction_method = "g_SCS",
       domain_scope = "annotated", custom_bg = ref_background,
       numeric_ns = "", sources = databases, as_short_link = FALSE, highlight = TRUE)
}

ora_webgestalt <- function(genes, ref_background, p_thresh = p_threshold,
                           databases = c("pathway_KEGG", "pathway_Reactome", "pathway_Wikipathway_cancer",
                                         "geneontology_Biological_Process", "geneontology_Molecular_Function")) {
  WebGestaltR(enrichMethod = "ORA", organism = "hsapiens", enrichDatabase = databases,
              interestGene = as.vector(genes), interestGeneType = "genesymbol",
              referenceGene = as.vector(ref_background), referenceGeneType = "genesymbol",
              fdrThr = p_thresh, isOutput = FALSE)
}

ref_gene <- unique(diff_exp$Gene.Symbol)

genes_up <- diff_exp %>%
  filter(Log2FC > log2fc_threshold & TTest.p.adj < p_threshold) %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

ora_up <- genes_up %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_up_wg <- genes_up %>%
  pull(Gene.Symbol) %>%
  ora_webgestalt(ref_background = ref_gene)

genes_down <- diff_exp %>%
  filter(Log2FC < -log2fc_threshold & TTest.p.adj < p_threshold) %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_down <- genes_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_down_wg <- genes_down %>%
  pull(Gene.Symbol) %>%
  ora_webgestalt(ref_background = ref_gene)

#gostplot(ora_up, capped = TRUE, interactive = TRUE)
#gostplot(ora_down, capped = TRUE, interactive = TRUE)

# TODO: Use WebGestaltR::affinityPropagation() for reducing redundant terms
# ora_down_reduced <- WebGestaltR::affinityPropagation(ora_down$result$term_name, -log10(ora_down$result$p_value))

plot_terms <- function(ora, selected_sources = c("CORUM", "KEGG", "REAC", "WP", "GO:MF", "GO:BP"),
                       terms_per_source = 5, p_thresh = p_threshold) {
  ora$result %>%
    filter(p_value < p_thresh) %>%
    filter(source %in% selected_sources) %>%
    group_by(source) %>%
    slice_min(p_value, n = terms_per_source) %>%
    ungroup() %>%
    mutate(`-log10(p)` = -log10(p_value),
           term_name = fct_reorder(str_trunc(term_name, 50), p_value, .desc = TRUE)) %>%
    vertical_bar_chart(term_name, `-log10(p)`, color_col = source, color_discrete = TRUE,
                       value_range = c(1, max(.$`-log10(p)`)),
                       line_intercept = 0, bar_label_shift = 0.18, break_steps = 2,
                       category_lab = "Enriched Term", value_lab = "Significance (-log10(p))") +
    facet_wrap(~source, scales = "free", ncol = 1) +
    theme(legend.position = "none")
}

 ora_up %>%
   plot_terms() %>%
   save_plot("overrepresentation_terms_up.png", height = 300)

 ora_down %>%
   plot_terms() %>%
   save_plot("overrepresentation_terms_down.png", height = 300)

# === Network Analysis ===
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700,
                          network_type = "full", input_directory = temp_dir)

create_string_network <- function(df, gene_col, logfc_col, string_db) {
  require(STRINGdb)
  require(dplyr)
  require(magrittr)
  require(scales)

  max_val <- df %>% pull({ { logfc_col } }) %>% abs() %>% max()
  domain <- c(-max_val, max_val)
  color_func <- scales::col_numeric(palette = bidirectional_color_pal, domain = domain)

  df_mapped <- df %>%
    select({ { gene_col } }, { { logfc_col } }) %>%
    as.data.frame() %>%
    string_db$map(quo_name(enquo(gene_col)), removeUnmappedRows = TRUE) %>%
    mutate(Color = color_func({ { logfc_col } }))

  payload <- df_mapped %$%
    string_db$post_payload(STRING_id, colors = Color)

  df_mapped %$%
    return(list(df_mapped = df_mapped,
                payload = payload,
                network = string_db$get_subnetwork(STRING_id),
                link = string_db$get_link(STRING_id, payload_id = payload, required_score = 700)))
}

string_up <- genes_up %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_up_png <- string_db$get_png(string_up$df_mapped$STRING_id, payload_id = string_up$payload,
                                   required_score = 700, file = here(plots_dir, "string_up.png"))

string_down <- genes_down %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_down_png <- string_db$get_png(string_down$df_mapped$STRING_id, payload_id = string_down$payload,
                                     required_score = 700, file = here(plots_dir, "string_down.png"))

# ===  Karyotype & Expression by ChrArm ===
# TODO: Compare Karyotype & Expression per Chromosome between buffered and scaling cell lines
gene_metadata <- expr_buf_procan %>%
  distinct(Gene.Symbol, Gene.Chromosome, Gene.ChromosomeArm, Gene.StartPosition, Gene.EndPosition) %>%
  mutate(Gene.Chromosome = as.integer(Gene.Chromosome))

# TODO: Nested facet with chromosome arm
diff_exp %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position")

# TODO: Check copy number / chromosome CNA between two cohorts

egfr_dc <- expr_buf_procan %>%
  filter(Gene.Symbol == "EGFR") %>%
  mutate(CNDiff = ChromosomeArm.CopyNumber - ChromosomeArm.CopyNumber.Baseline) %>%
  filter(abs(CNDiff) > 0.5) %>%
  mutate(CNV = if_else(CNDiff > 0, "Gain", "Loss")) %>%
  drop_na(Buffering.ChrArmLevel.Ratio) %>%
  signif_beeswarm_plot(CNV, Buffering.ChrArmLevel.Ratio, color_col = CellLine.AneuploidyScore) %>%
  save_plot("EGFR_DC_ChrArm.png", height = 150, width = 150)
