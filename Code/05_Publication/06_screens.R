library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "annotation.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)

cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))

df_crispr_model_buf <- read_parquet(here(output_data_dir, "model_buffering_gene_dependency_depmap.parquet"))
df_crispr_buf <- read_parquet(here(output_data_dir, "buffering_gene_dependency_depmap.parquet"))

# === Dependency-Buffering Correlation Panel
color_mapping_driver <- scale_color_manual(values = dicrete_color_pal1, na.value = color_palettes$Missing)

volcano_dep_corr <- df_crispr_buf %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(abs(CRISPR.DependencyScore.Corr) > 0.2 &
                           CRISPR.DependencyScore.Corr.p < p_threshold &
                           !is.na(CancerDriverMode),
                         Gene.Symbol, NA)) %>%
  arrange(!is.na(CancerDriverMode)) %>%
  plot_volcano(CRISPR.DependencyScore.Corr, CRISPR.DependencyScore.Corr.p, Label, CancerDriverMode,
               color_mapping = color_mapping_driver, value_threshold = 0.2) +
  labs(x = "Correlation (Dependency Score, Buffering Ratio)", color = "Cancer Driver", y = "-log10(p.adj)")

# === Gene Correlation Panels
selected_genes <- c("EGFR", "CDK6", "ELMO2", "BRAT1")
gene_corr_plots <- list()

for (gene in selected_genes) {
  cor_gene <- df_crispr_buf %>%
    filter(Gene.Symbol == gene) %>%
    distinct(CRISPR.DependencyScore.Corr, CRISPR.DependencyScore.Corr.p)

  panel_gene <- df_crispr_buf %>%
    filter(Gene.Symbol == gene) %>%
    ggplot() +
    aes(x = CRISPR.DependencyScore, y = Buffering.GeneLevel.Ratio) +
    geom_point(size = 2, color = default_color) +
    stat_smooth(method = lm, color = highlight_colors[2]) +
    annotate("text", x = 0, y = 4.5, hjust = 0, size = 5,
             label = paste0(print_corr(cor_gene$CRISPR.DependencyScore.Corr), ", ", print_signif(cor_gene$CRISPR.DependencyScore.Corr.p))) +
    annotate("text", x = 1, y = -3.5, hjust = 1, size = 5, label = gene) +
    labs(x = "Dependency Score", y = "Buffering Ratio") +
    theme(legend.position = c("left", "top"),
          legend.position.inside = c(0.01, 0.90),
          legend.justification = c("left", "top"),
          legend.title = element_blank(),
          legend.text = element_text(size = base_size)) +
    lims(x = c(0, 1), y = c(-3.5, 4.5))

  gene_corr_plots[[gene]] <- panel_gene
}

# === Model Buffering Dependency Difference Panel
volcano_dep_diff <- df_crispr_model_buf %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & !is.na(CancerDriverMode), Gene.Symbol, NA)) %>%
  arrange(!is.na(CancerDriverMode)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, CancerDriverMode,
               value_threshold = 0.1, color_mapping = color_mapping_driver) +
  labs(color = "Cancer Driver", x = "Log2FC (Dependency Score)", y = "-log10(p.adj)")

# === Model Buffering Enrichment Panel
ora_up <- df_crispr_model_buf %>%
  filter(Significant == "Up") %>%
  arrange(desc(Log2FC)) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "CORUM"),
                     custom_color = color_palettes$DiffExp["Up"], string_trunc = 60)

ora_down <- df_crispr_model_buf %>%
  filter(Significant == "Down") %>%
  arrange(Log2FC) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "WP"),
                     custom_color = color_palettes$DiffExp["Down"], string_trunc = 60)

# === Combine Panels into Figure ===
figure6_sub1 <- cowplot::plot_grid(volcano_dep_corr, gene_corr_plots$EGFR, gene_corr_plots$CDK6 + ylab(NULL),
                                   gene_corr_plots$ELMO2 + ylab(NULL), gene_corr_plots$BRAT1 + ylab(NULL),
                                   labels = c("A", "B", "C", "D", "E"),
                                   rel_widths = c(1, 0.55, 0.5, 0.5, 0.5),
                                   nrow = 1, ncol = 5)

figure6_sub2 <- cowplot::plot_grid(volcano_dep_diff, ora_down, ora_up,
                                   labels = c("F", "G", "H"),
                                   rel_widths = c(1, 0.8, 0.8),
                                   nrow = 1, ncol = 3)

figure6 <- cowplot::plot_grid(figure6_sub1, figure6_sub2,
                              nrow = 2, ncol = 1, rel_heights = c(0.75, 1))

cairo_pdf(here(plots_dir, "figure06.pdf"), width = 17, height = 9)
figure6
dev.off()
