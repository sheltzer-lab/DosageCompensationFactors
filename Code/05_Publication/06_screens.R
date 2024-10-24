library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)

cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))

df_crispr_model_buf <- read_parquet(here(output_data_dir, "model_buffering_gene_dependency_depmap.parquet"))
df_crispr_buf <- read_parquet(here(output_data_dir, "buffering_gene_dependency_depmap.parquet"))

# === Dependency-Buffering Correlation Panel
color_mapping_driver <- scale_color_manual(values = discrete_color_pal1b, na.value = color_palettes$Missing)

top_corr <- df_crispr_buf %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  filter(CRISPR.DependencyScore.Corr.p < p_threshold) %>%
  slice_max(abs(CRISPR.DependencyScore.Corr), n = 10)

volcano_dep_corr <- df_crispr_buf %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(abs(CRISPR.DependencyScore.Corr) > 0.2 &
                           CRISPR.DependencyScore.Corr.p < p_threshold &
                           (!is.na(CancerDriverMode) | Gene.Symbol %in% top_corr$Gene.Symbol),
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
    annotate("text", x = 0, y = 4.2, hjust = 0, size = 5,
             label = paste0(print_corr(cor_gene$CRISPR.DependencyScore.Corr), ", ", print_signif(cor_gene$CRISPR.DependencyScore.Corr.p))) +
    annotate("text", x = 1, y = -3.2, hjust = 1, size = 5, label = gene) +
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
top_diff <- df_crispr_model_buf %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  filter(Test.p.adj < p_threshold) %>%
  slice_max(abs(Log2FC), n = 10)

essential_buf_only <- df_crispr_model_buf %>%
  filter(Significant == "Up") %>%
  filter(Mean_GroupA < 0.5 & Mean_GroupB > 0.5)

essential_scaling_only <- df_crispr_model_buf %>%
  filter(Significant == "Down") %>%
  filter(Mean_GroupA > 0.5 & Mean_GroupB < 0.5)

max_abs_log2fc <- df_crispr_model_buf %>% pull(Log2FC) %>% abs() %>% max(na.rm = TRUE)

volcano_dep_diff <- df_crispr_model_buf %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(ExclusiveEssentiality = case_when(
    Significant == "Up" & Mean_GroupA < 0.5 & Mean_GroupB > 0.5 ~ "High Buffering",
    Significant == "Down" & Mean_GroupA > 0.5 & Mean_GroupB < 0.5 ~ "Low Buffering",
    TRUE ~ NA
  )) %>%
  mutate(Label = if_else(!is.na(Significant) &
                           (!is.na(CancerDriverMode) |
                             Gene.Symbol %in% top_diff$Gene.Symbol |
                             !is.na(ExclusiveEssentiality)),
                         Gene.Symbol, NA)) %>%
  arrange(!is.na(CancerDriverMode)) %>%
  mutate(`-log10(p)` = -log10(Test.p.adj)) %>%
  ggplot() +
  aes(x = Log2FC, y = `-log10(p)`, label = Label, color = CancerDriverMode) +
  geom_point(aes(alpha = Significant, shape = ExclusiveEssentiality, size = ExclusiveEssentiality)) +
  color_mapping_driver +
  scale_shape_manual(values = c(`High Buffering` = 17, `Low Buffering` = 15), na.value = 16) +
  scale_size_manual(values = c(`High Buffering` = 2, `Low Buffering` = 2), na.value = 1, guide = "none") +
  scale_alpha_manual(values = c(Up = 1, Down = 1), na.value = 0.5, guide = "none") +
  geom_hline(yintercept = -log10(p_threshold),
             linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0.1,
             linetype = "dashed", color = "black") +
  geom_vline(xintercept = -0.1,
             linetype = "dashed", color = "black") +
  geom_label_repel(min.segment.length = 0.01, label.size = 0.15,
                   seed = 42, max.iter = 30000, max.time = 1.5,
                   point.padding = 0.3, label.padding = 0.3, box.padding = 0.3,
                   force = 5, max.overlaps = 20) +
  xlim(c(-max_abs_log2fc, max_abs_log2fc)) +
  labs(color = "Cancer Driver", x = "Log2FC (Dependency Score)",
       y = "-log10(p.adj)", shape = "Exclusively Essential")
#  theme(legend.position = "top", legend.direction = "horizontal",
#        legend.box = "vertical", legend.box.just = "left", legend.margin = margin(0,0,0,0, unit = 'cm')) +
#  guides(colour = guide_legend(title.position = "top", title.hjust = 0),
#         shape = guide_legend(title.position = "top", title.hjust = 0, override.aes = list(alpha = 1, size = 3)))

# === Model Buffering Enrichment Panel
ora_up <- df_crispr_model_buf %>%
  filter(Significant == "Up") %>%
  arrange(desc(Log2FC)) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "CORUM"), string_trunc = 45,
                     custom_color = color_palettes$DiffExpBackground["Up"])

ora_down <- df_crispr_model_buf %>%
  filter(Significant == "Down") %>%
  arrange(Log2FC) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "WP"), string_trunc = 45,
                     custom_color = color_palettes$DiffExpBackground["Down"])

# === Combine Panels into Figure ===
gene_corr_panel <- cowplot::plot_grid(gene_corr_plots$EGFR + xlab(NULL), gene_corr_plots$CDK6 + xlab(NULL) + ylab(NULL),
                                      gene_corr_plots$ELMO2, gene_corr_plots$BRAT1 + ylab(NULL),
                                      nrow = 2, ncol = 2)

#figure6_sub1 <- cowplot::plot_grid(volcano_dep_corr, gene_corr_panel,
#                                   labels = c("A", "B"),
#                                   nrow = 1, ncol = 2)

figure6_sub2 <- cowplot::plot_grid(ora_down, ora_up,
                                   labels = c("B", "C"),
                                   nrow = 1, ncol = 2)

figure6 <- cowplot::plot_grid(volcano_dep_diff, figure6_sub2,
                              labels = c("A", ""),
                              nrow = 2, ncol = 1)

cairo_pdf(here(plots_dir, "figure06.pdf"), width = 8, height = 9)
figure6
dev.off()
