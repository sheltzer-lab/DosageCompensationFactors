library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(openxlsx)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir
tables_dir <- here(tables_base_dir, "Publication")

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))

df_crispr_model_buf <- read_parquet(here(output_data_dir, "model_buffering_gene_dependency_depmap.parquet"))
df_crispr_model_buf_control <- read_parquet(here(output_data_dir, "model_buffering_gene_dependency_depmap_adherent.parquet"))
df_crispr_buf <- read_parquet(here(output_data_dir, "buffering_gene_dependency_depmap.parquet"))
model_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
drug_screens <- read_parquet(here(output_data_dir, "drug_screens.parquet")) %>% select(-CellLine.Name)
drug_targets <- read_parquet(here(output_data_dir, "drug_effect_buffering_target.parquet"))
drug_mechanisms <- read_parquet(here(output_data_dir, "drug_effect_buffering_mechanism.parquet"))
drug_mechanisms_control <- read_parquet(here(output_data_dir, "drug_effect_buffering_mechanism_control.parquet"))
median_response_buf <- read_parquet(here(output_data_dir, "median_drug_effect.parquet"))

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
  labs(color = "Cancer Driver", x = "Dependency Score Log2FC",
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
                     custom_color = color_palettes$DiffExp["Up"])

ora_down <- df_crispr_model_buf %>%
  filter(Significant == "Down") %>%
  arrange(Log2FC) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "WP"), string_trunc = 45,
                     custom_color = color_palettes$DiffExp["Down"])

# === Drug Response Panel
median_response_plot_buf <- median_response_buf %>%
  drop_na(Buffering) %>%
  ggplot() +
  aes(x = Buffering, y = Drug.Effect.Median, color = Buffering) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("High", "Low")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = 0.06,
              size = 0.8, tip_length = 0, color = "black") +
  scale_color_manual(values = c("High" = color_palettes$BufferingClasses[["Buffered"]],
                                "Low" = color_palettes$BufferingClasses[["Scaling"]])) +
  theme(legend.position = "none") +
  ylim(c(-0.35, 0.15)) +
  labs(x = "Sample Buffering", y = "Median Drug Effect")

median_response_plot_control <- median_response_buf %>%
  drop_na(DrugStatus) %>%
  ggplot() +
  aes(x = DrugStatus, y = Drug.Effect.Median) +
  geom_boxplot() +
  theme(legend.position = "none") +
  xlab("Drug Control") +
  ylim(c(-0.35, 0.15)) +
  ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

panel_drug_response <- cowplot::plot_grid(median_response_plot_buf,
                                          median_response_plot_control,
                                          align = "h", axis = "lr", ncol = 2, rel_widths = c(1, 0.8))

panel_drug_response_horiz <- cowplot::plot_grid(median_response_plot_buf +
                                                  labs(x = "Sample\nBuffering", y = NULL) +
                                                  coord_flip() +
                                                  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                                                median_response_plot_control +
                                                  labs(x = "Drug\nControl", y = "Median Drug Effect") +
                                                  coord_flip() +
                                                  theme(axis.text.y = element_text(), axis.ticks.y = element_line()),
                                                align = "v", axis = "tb",
                                                nrow = 2, rel_heights = c(0.8, 1))

# === Drug Mechanism Panel
moa_heatmap_diff <- drug_mechanisms %>%
  filter(CommonEffect & !is.na(EffectiveIn)) %>%
  mutate(Label = map_signif(DrugEffect.Buffering.Corr.p.adj, thresholds = c(0.05, 0.01, 0.001))) %>%
  ggplot() +
  aes(y = Drug.MOA, x = "Drug Effect Log2FC (High - Low Buffering)", label = Label,
      fill = DrugEffect.Buffering.Group.Log2FC, color = DrugEffect.Buffering.Corr) +
  scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab",
                       limits = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6), oob = scales::squish) +
  scale_color_gradientn(colors = bidirectional_color_pal2, space = "Lab",
                        limits = c(-0.2, 0.2), breaks = c(-0.2, -0.1, 0, 0.1, 0.2), oob = scales::squish) +
  scale_x_discrete(position = "top") +
  geom_tile() +
  geom_tile(color = "white") +
  geom_text(color = "black") +
  labs(x = NULL, y = NULL,
       color = "Correlation (Drug Effect, Sample BR)", fill = "Drug Effect Log2FC (High - Low Buffering)") +
  theme_void() +
  guides(fill = guide_colourbar(order = 1),
         colour = guide_colourbar(order = 2)) +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5),
        legend.key.size = unit(16, "points"),
        legend.box = "vertical",
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.margin = margin(0,0,5, -200))

moa_heatmap_corr <- drug_mechanisms %>%
  filter(CommonEffect & !is.na(EffectiveIn)) %>%
  mutate(Label = map_signif(DrugEffect.Buffering.Group.Log2FC.p.adj, thresholds = c(0.05, 0.01, 0.001))) %>%
  ggplot() +
  aes(y = Drug.MOA, x = "Correlation (Drug Effect ~ Sample BR)", fill = DrugEffect.Buffering.Corr, label = Label) +
  scale_fill_gradientn(colors = bidirectional_color_pal2, space = "Lab",
                       limits = c(-0.2, 0.2), breaks = c(-0.2, -0.1, 0, 0.1, 0.2), oob = scales::squish) +
  geom_tile() +
  geom_tile(color = "white") +
  geom_text(color = "black") +
  labs(x = "Drug Mechanism", y = NULL) +
  theme_void() +
  theme(legend.position = "none")

panel_drug_mechanism <- cowplot::plot_grid(moa_heatmap_diff, moa_heatmap_corr,
                                           nrow = 1, ncol = 2, align = "h", axis = "lr",
                                           rel_widths = c(1, 0.17))

# === Drug Target Panel
target_heatmap_diff <- drug_targets %>%
  filter(CommonEffect & !is.na(EffectiveIn)) %>%
  mutate(Label = map_signif(DrugEffect.Buffering.Corr.p.adj, thresholds = c(0.05, 0.01, 0.001))) %>%
  ggplot() +
  aes(x = Drug.Target, y = "Drug Effect Log2FC (High - Low Buffering)", label = Label,
      fill = DrugEffect.Buffering.Group.Log2FC, color = DrugEffect.Buffering.Corr) +
  scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab",
                       limits = c(-0.6, 0.6), breaks = c(-0.6, -0.3, 0, 0.3, 0.6), oob = scales::squish) +
  scale_color_gradientn(colors = bidirectional_color_pal2, space = "Lab",
                       limits = c(-0.2, 0.2), breaks = c(-0.2, -0.1, 0, 0.1, 0.2), oob = scales::squish) +
  scale_x_discrete(position = "top") +
  geom_tile() +
  geom_text(color = "black") +
  labs(x = NULL, y = NULL, color = "Correlation") +
  theme_void() +
  guides(fill = guide_colourbar(order = 1),
         colour = guide_colourbar(order = 2)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        legend.key.size = unit(16, "points"),
        legend.box = "horizontal",
        legend.position = "right",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0))

target_heatmap_corr <- drug_targets %>%
  filter(CommonEffect & !is.na(EffectiveIn)) %>%
  mutate(Label = map_signif(DrugEffect.Buffering.Group.Log2FC.p.adj, thresholds = c(0.05, 0.01, 0.001))) %>%
  ggplot() +
  aes(x = Drug.Target, y = "Correlation (Drug Effect ~ Model BR)", fill = DrugEffect.Buffering.Corr, label = Label) +
  scale_fill_gradientn(colors = bidirectional_color_pal2, space = "Lab",
                       limits = c(-0.2, 0.2), breaks = c(-0.2, -0.1, 0, 0.1, 0.2), oob = scales::squish) +
  geom_tile() +
  geom_text(color = "black") +
  labs(x = "Drug Target", y = NULL) +
  theme_void() +
  theme(legend.position = "none",
        axis.title.x.bottom = element_text(),
        axis.text.y = element_text(hjust = 1, vjust = 0.5))

panel_drug_target <- cowplot::plot_grid(target_heatmap_diff, target_heatmap_corr,
                                        nrow = 2, ncol = 1, align = "v", axis = "tb",
                                        rel_heights = c(1, 0.6))

# === Combine Panels into Figure ===
gene_corr_panel <- cowplot::plot_grid(gene_corr_plots$EGFR + xlab(NULL), gene_corr_plots$CDK6 + xlab(NULL) + ylab(NULL),
                                      gene_corr_plots$ELMO2, gene_corr_plots$BRAT1 + ylab(NULL),
                                      nrow = 2, ncol = 2)

#figure6_sub1 <- cowplot::plot_grid(volcano_dep_corr, gene_corr_panel,
#                                   labels = c("A", "B"),
#                                   nrow = 1, ncol = 2)

figure6_sub2 <- cowplot::plot_grid(volcano_dep_diff, panel_drug_response,
                                   labels = c("A", "D"), rel_widths = c(1, 0.7),
                                   nrow = 1, ncol = 2)

figure6_sub3 <- cowplot::plot_grid(ora_down, ora_up, panel_drug_mechanism,
                                   ncol = 3, labels = c("B", "C", "E"))

figure6 <- cowplot::plot_grid(figure6_sub2, figure6_sub3,
                              nrow = 2, ncol = 1)

cairo_pdf(here(plots_dir, "figure06.pdf"), width = 12, height = 10)
figure6
dev.off()

# === Tables ===
crispr_model_all <- df_crispr_model_buf %>%
  mutate(ExclusiveEssentiality = case_when(
    Significant == "Up" & Mean_GroupA < 0.5 & Mean_GroupB > 0.5 ~ "High Buffering",
    Significant == "Down" & Mean_GroupA > 0.5 & Mean_GroupB < 0.5 ~ "Low Buffering",
    TRUE ~ NA
  )) %>%
  mutate(Dataset = "DepMap") %>%
  bind_rows(df_crispr_model_buf_control %>% mutate(Dataset = "DepMap (adherent control)")) %>%
  mutate(GroupA = "Low Buffering", GroupB = "High Buffering") %>%
  arrange(Significant, ExclusiveEssentiality)

crispr_model_common <- crispr_model_all %>%
  drop_na(Significant) %>%
  filter(!is.na(ExclusiveEssentiality) | Dataset == "DepMap (adherent control)") %>%
  add_count(Gene.Symbol, Significant) %>%
  filter(n == 2) %>%
  select(-n) %>%
  arrange(Significant, Gene.Symbol)

drugs_export <- bind_rows(
  read_parquet(here(output_data_dir, "sensitivity_correlation_depmap_gene_filtered.parquet")) %>% mutate(Dataset = "DepMap"),
  read_parquet(here(output_data_dir, "sensitivity_correlation_procan_gene_filtered.parquet")) %>% mutate(Dataset = "ProCan"),
  read_parquet(here(output_data_dir, "sensitivity_correlation_mean-norm-rank.parquet")) %>% mutate(Dataset = "Cell Lines (aggregated, Mean Normalized Ranks)")
)

drug_mechanisms_control <- drug_mechanisms_control %>%
  mutate(GroupA = "Low Aneuploidy", GroupB = "High Aneuploidy")

wb <- createWorkbook()
sheet_crispr <- addWorksheet(wb, "CRISPR-KO")
sheet_crispr_common <- addWorksheet(wb, "CRISPR-KO (common)")
sheet_drug_effect <- addWorksheet(wb, "Median Drug Effect")
sheet_drugs <- addWorksheet(wb, "Drugs")
sheet_moa <- addWorksheet(wb, "Drug Mechanisms")
sheet_moa_control <- addWorksheet(wb, "Drug Mechanisms (control)")
sheet_target <- addWorksheet(wb, "Drug Targets")
writeDataTable(wb = wb, sheet = sheet_crispr, x = crispr_model_all)
writeDataTable(wb = wb, sheet = sheet_crispr_common, x = crispr_model_common)
writeDataTable(wb = wb, sheet = sheet_drug_effect, x = median_response_buf)
writeDataTable(wb = wb, sheet = sheet_drugs, x = drugs_export)
writeDataTable(wb = wb, sheet = sheet_moa, x = drug_mechanisms %>% mutate(Dataset = "Cell Lines (aggregated, Mean Normalized Ranks)"))
writeDataTable(wb = wb, sheet = sheet_moa_control, x = drug_mechanisms_control)
writeDataTable(wb = wb, sheet = sheet_target, x = drug_targets %>% mutate(Dataset = "Cell Lines (aggregated, Mean Normalized Ranks)"))
saveWorkbook(wb, here(tables_dir, "supplementary_table9.xlsx"), overwrite = TRUE)
