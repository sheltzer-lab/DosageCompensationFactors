library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))
expr_buf_p0211 <- read_parquet(here(output_data_dir, "expression_buffering_p0211.parquet"))
expr_buf_chunduri <- read_parquet(here(output_data_dir, "expression_buffering_chunduri.parquet"))

expr_buf_eng <- bind_rows(expr_buf_p0211, expr_buf_chunduri) %>%
  mutate(Dataset = "Engin.",
         Gene.Chromosome = as.character(Gene.Chromosome))

# === DC Class Illustration Panel ===
dc_class_line <- expr_buf_depmap %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  filter(Gene.CopyNumber.Baseline == 2) %>%
  filter(Gene.CopyNumber %in% c(1, 3, 4)) %>%
  mutate(Gene.CopyNumber = as.factor(Gene.CopyNumber),
         Buffering.GeneLevel.Class = factor(Buffering.GeneLevel.Class,
                                            levels = c("Scaling", "Buffered", "Anti-Scaling"))) %>%
  drop_na(Buffering.GeneLevel.Class, Protein.Expression.Normalized) %>%
  group_by(Buffering.GeneLevel.Class, Gene.CopyNumber) %>%
  summarize(Mean.Protein.Ratio = (mean(2^Protein.Expression.Normalized / 2^Protein.Expression.Baseline) - 1),
            .groups = "drop") %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Class, y = Mean.Protein.Ratio,
      group = Gene.CopyNumber, color = Gene.CopyNumber, label = Gene.CopyNumber) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_path(linewidth = 1) +
  scale_y_continuous(limits = c(-0.6, 2.5),
                     breaks = seq(-0.5, 2.5, 0.5),
                     labels = paste0(seq(-0.5, 2.5, 0.5) * 100, "%")) +
  scale_color_manual(values = color_palettes$CopyNumbers) +
  labs(x = "Buffering Class", y = "Mean Protein Abundance", color = "Gene Copy Number") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit = 'cm'))

# === Categorical Distribution Panel ===
df_share_gene <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gene CN Gain", "Gene CN Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  select(Buffering.GeneLevel.Class, Dataset, CNV) %>%
  drop_na() %>%
  count(Buffering.GeneLevel.Class, Dataset, CNV) %>%
  group_by(Dataset, CNV) %>%
  mutate(Share = (n / sum(n))) %>%
  ungroup()

df_share_chr <- bind_rows(expr_buf_depmap, expr_buf_eng) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "ChrArm Gain", "ChrArm Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  select(Buffering.ChrArmLevel.Class, Dataset, CNV) %>%
  drop_na() %>%
  count(Buffering.ChrArmLevel.Class, Dataset, CNV) %>%
  group_by(Dataset, CNV) %>%
  mutate(Share = (n / sum(n))) %>%
  ungroup()

stacked_buf_cn <- df_share_gene %>%
  ggplot() +
  aes(fill = Buffering.GeneLevel.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_wrap(~CNV) +
  scale_fill_manual(values = color_palettes$BufferingClasses) +
  scale_y_continuous(labels = scales::percent) +
  labs(fill = NULL, y = "Fraction") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit = 'cm'))

stacked_buf_cn_chr <- df_share_chr %>%
  ggplot() +
  aes(fill = Buffering.ChrArmLevel.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_wrap(~CNV) +
  scale_fill_manual(values = color_palettes$BufferingClasses) +
  scale_y_continuous(labels = scales::percent) +
  labs(fill = "Buffering Class") +
  theme(legend.position = "none")

# === DC Class Protein Expression Panel ===
prot_exp_dc_class <- expr_buf_depmap %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  drop_na(Buffering.GeneLevel.Class, Protein.Expression.Normalized) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gene CN Gain", "Gene CN Loss"),
         Buffering.GeneLevel.Class = factor(Buffering.GeneLevel.Class,
                                            levels = c("Scaling", "Buffered", "Anti-Scaling"))) %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Class, y = Protein.Expression.Normalized, color = Buffering.GeneLevel.Class) +
  geom_hline(aes(yintercept = mean(Protein.Expression.Normalized, na.rm = TRUE)), color = default_color) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  geom_signif(comparisons = list(c("Anti-Scaling", "Buffered"),
                                 c("Buffered", "Scaling"),
                                 c("Anti-Scaling", "Scaling")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(1.5, 1.5, 1.8),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  facet_grid(~Gene.CopyNumber.Event) +
  labs(x = "Buffering Class", y = "Protein Abundance") +
  scale_color_manual(values = color_palettes$BufferingClasses) +
  theme(legend.position = "none")

# === Tumor vs Cell Lines Panel ===
dc_dataset_dist <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  filter(abs(Gene.CopyNumber - Gene.CopyNumber.Baseline) > 0.2) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gene CN Gain", "Gene CN Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  select(Dataset, Buffering.GeneLevel.Ratio, Gene.CopyNumber.Event, Gene.Symbol) %>%
  drop_na() %>%
  ggplot() +
  aes(x = Dataset, y = Buffering.GeneLevel.Ratio, color = Dataset) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  geom_signif(comparisons = list(c("DepMap", "ProCan"),
                                 c("CPTAC", "ProCan"),
                                 c("CPTAC", "DepMap")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(1.8, 1.8, 2.1),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  theme(legend.position = "") +
  scale_color_manual(values = color_palettes$Datasets) +
  facet_wrap(~Gene.CopyNumber.Event, ncol = 2) +
  labs(x = "Dataset", y = "Buffering Ratio")

# === Buffering Gain vs. Loss Panel ===
br_by_cna <- bind_rows(expr_buf_cptac, expr_buf_depmap, expr_buf_eng) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CopyNumber.Event = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Gain", "Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  group_by(Gene.Symbol, CopyNumber.Event, Dataset) %>%
  summarize(Buffering.ChrArmLevel.Ratio = mean(Buffering.ChrArmLevel.Ratio)) %>%
  ggplot() +
  aes(x = CopyNumber.Event, y = Buffering.ChrArmLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  geom_signif(comparisons = list(c("Gain", "Loss")),
              map_signif_level = print_signif, y_position = 1.9, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  labs(x = "Chromosome Arm", y = "Mean Buffering Ratio") +
  ylim(c(-0.5, 2.5)) +
  facet_grid(~Dataset)

br_by_cnv <- bind_rows(expr_buf_cptac, expr_buf_depmap) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  group_by(Gene.Symbol, Gene.CopyNumber.Event, Dataset) %>%
  summarize(Buffering.GeneLevel.Ratio = mean(Buffering.GeneLevel.Ratio)) %>%
  ggplot() +
  aes(x = Gene.CopyNumber.Event, y = Buffering.GeneLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  geom_signif(comparisons = list(c("Gain", "Loss")),
              map_signif_level = print_signif, y_position = 1.9, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  labs(x = "Gene Copy Number", y = "Mean Buffering Ratio") +
  ylim(c(-0.5, 2.5)) +
  facet_grid(~Dataset)

# === Low Variance Buffering Panel ===
low_var_buf <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  analyze_low_br_variance()

scatter_signif_buffered <- low_var_buf %>%
  mutate(Label = if_else(Top10, Gene.Symbol, NA)) %>%
  filter(Dataset == "CPTAC") %>%
  arrange(Top50) %>%
  ggplot() +
  aes(x = Gene.BR.SD, y = Gene.BR.Mean, label = Label, color = Top50, alpha = Top50) +
  geom_point() +
  geom_label_repel(na.rm = TRUE, max.overlaps = 50) +
  scale_color_manual(values = c(default_color, highlight_color)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(x = "Standard Deviation of Buffering Ratio", y = "Mean Buffering Ratio",
       color = "Frequently Buffered", alpha = "Frequently Buffered") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit = 'cm')) +
  guides(colour = guide_legend(title.position="top", title.hjust = 0))

ora_buf <- low_var_buf %>%
  filter(Top50) %>%
  filter(Dataset == "CPTAC") %>%
  arrange(Gene.BR.SD.MeanNormRank) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = TRUE)

ora_buf_terms <- ora_buf %>%
  plot_terms_compact(selected_sources = c("GO:MF", "KEGG", "WP"), custom_color = highlight_color)

# === Combine Panels into Figure ===
#dc_illus <- cowplot::draw_image(here(illustrations_dir, "dc_illustration.svg"))
#buf_classes <- cowplot::draw_image(here(illustrations_dir, "buffering_classes.svg"))

#illustrations <- cowplot::plot_grid(cowplot::ggdraw() + dc_illus,
#                                    cowplot::ggdraw() + buf_classes,
#                                    labels = c("A", "B"))

stacked_buf_all <- cowplot::plot_grid(stacked_buf_cn, stacked_buf_cn_chr + ylab(NULL),
                                      ncol = 2, align = "h", axis = "tb",
                                      labels = c("D", "E"), rel_widths = c(1, 0.7))

figure1_sub1 <- cowplot::plot_grid(dc_class_line, stacked_buf_all,
                                   rel_widths = c(0.5, 1), labels = c("C", ""), ncol = 2)

figure1_sub2 <- cowplot::plot_grid(prot_exp_dc_class, dc_dataset_dist,
                                   rel_widths = c(1, 1), labels = c("F", "G"), ncol = 2)

br_by_cnv_all <- cowplot::plot_grid(br_by_cnv, br_by_cna + ylab(NULL),
                                    rel_widths = c(1, 0.9), labels = c("H", "I"), ncol = 2,
                                    align = "h", axis = "tb")

figure1_sub3 <- cowplot::plot_grid(br_by_cnv_all, scatter_signif_buffered, ora_buf_terms,
                                   rel_widths = c(1, 0.8, 0.8), labels = c("", "J", "K"), ncol = 3)

figure1 <- cowplot::plot_grid(figure1_sub1, figure1_sub2, figure1_sub3,
                              nrow = 3, rel_heights = c(0.8, 1, 1))

cairo_pdf(here(plots_dir, "figure01_01.pdf"), width = 13, height = 13)
figure1
dev.off()
