library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(openxlsx)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Screens", "DrugSensitivity")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
drug_screens <- read_parquet(here(output_data_dir, "drug_screens.parquet"))
drug_meta <- read_parquet(here(output_data_dir, "drug_metadata.parquet"))
cellline_buf_filtered_procan <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_procan.parquet"))
cellline_buf_filtered_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_depmap.parquet"))
cellline_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(CellLine.Name, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

# Calculate Correlation between Drug Sensitivity and Cell Line Buffering
analyze_sensitivity <- function(df_celllines, df_drug_screens, cellline_value_col, drug_value_col) {
  test <- df_celllines %>%
    inner_join(y = df_drug_screens, by = "CellLine.Name",
               relationship = "one-to-many", na_matches = "never") %>%
    mutate(CellLine.Value.Median = median({ { cellline_value_col } }, na.rm = TRUE),
           CellLine.Group = if_else({ { cellline_value_col } } > CellLine.Value.Median, "High", "Low")) %>%
    group_by(Drug.ID, CellLine.Group) %>%
    mutate(Group.Mean = mean({{drug_value_col}}, na.rm = TRUE)) %>%
    group_by(Drug.ID, Drug.Name) %>%
    # https://stackoverflow.com/questions/48041504/calculate-pairwise-correlation-in-r-using-dplyrmutate
    summarize(
      # Correlation analysis Buffering vs Sensitivity
      # ToDo: Avoid calculating correlation twice
      Corr.Sensitivity_Buffering = cor.test({ { cellline_value_col } }, { { drug_value_col } },
                                            method = "spearman")$estimate[["rho"]],
      Corr.p = cor.test({ { cellline_value_col } }, { { drug_value_col } },
                        method = "spearman")$p.value,
      # Statistical test sensitivity low vs. high buffering
      BufferingGroup.Sensitivity.Diff = if_else(Group.Mean[CellLine.Group == "High"][1] < Group.Mean[CellLine.Group == "Low"][1],
                                    "Higher", "Lower"),
      BufferingGroup.Diff.Test.p = wilcox.test({{drug_value_col}}[CellLine.Group == "High"],
                                               {{drug_value_col}}[CellLine.Group == "Low"])$p.value,
    ) %>%
    ungroup() %>%
    mutate(Corr.p.adj = p.adjust(Corr.p, method = "BH"),
           BufferingGroup.Diff.Test.p.adj = p.adjust(BufferingGroup.Diff.Test.p, method = "BH")) %>%
    arrange(Drug.Name)
}

drug_dc_corr_procan <- cellline_buf_filtered_procan %>%
  analyze_sensitivity(drug_screens, Buffering.CellLine.Ratio.ZScore, Drug.MFI.Log2FC)
drug_dc_corr_depmap <- cellline_buf_filtered_depmap %>%
  analyze_sensitivity(drug_screens, Buffering.CellLine.Ratio.ZScore, Drug.MFI.Log2FC)
drug_dc_corr_agg_rank <- cellline_buf_agg %>%
  analyze_sensitivity(drug_screens, Buffering.CellLine.MeanNormRank, Drug.MFI.Log2FC)
drug_dc_corr_agg_mean <- cellline_buf_agg %>%
  analyze_sensitivity(drug_screens, Buffering.CellLine.StandardizedMean, Drug.MFI.Log2FC)


# Evaluate Reproducibility of Results
cor_datasets <- cor.test(drug_dc_corr_depmap$Corr.Sensitivity_Buffering,
         drug_dc_corr_procan$Corr.Sensitivity_Buffering,
         method = "pearson")

cor_agg_methods <- cor.test(drug_dc_corr_agg_rank$Corr.Sensitivity_Buffering,
         drug_dc_corr_agg_mean$Corr.Sensitivity_Buffering,
         method = "pearson")

cat(capture.output(cor_datasets), file = here(reports_dir, "drug_correlation_datasets.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cor_agg_methods), file = here(reports_dir, "drug_correlation_aggregation-methods.txt"),
    append = FALSE, sep = "\n")


drug_dc_corr_agg_rank %>% mutate(Method = "StandardizedMean") %>%
  bind_rows(drug_dc_corr_agg_mean %>% mutate(Method = "MeanNormRank")) %>%
  ggplot(aes(x = Corr.Sensitivity_Buffering, color = Method)) +
  geom_density()

# === Write results ===
write.xlsx(drug_dc_corr_procan, here(tables_base_dir, "sensitivity_correlation_procan_gene_filtered.xlsx"),
           colNames = TRUE)
write.xlsx(drug_dc_corr_depmap, here(tables_base_dir, "sensitivity_correlation_depmap_gene_filtered.xlsx"),
           colNames = TRUE)
write.xlsx(drug_dc_corr_agg_rank, here(tables_base_dir, "sensitivity_correlation_mean-norm-rank.xlsx"),
           colNames = TRUE)
write.xlsx(drug_dc_corr_agg_mean, here(tables_base_dir, "sensitivity_correlation_standardized-mean.xlsx"),
           colNames = TRUE)

# === Visualize results for rank aggregated sensitivities ===
top_sensitivity <- drug_dc_corr_agg_rank %>%
  filter(BufferingGroup.Diff.Test.p < p_threshold) %>%
  filter(BufferingGroup.Sensitivity.Diff == "Higher") %>%
  slice_min(Corr.Sensitivity_Buffering, n = 12)
bot_sensitivity <- drug_dc_corr_agg_rank %>%
  filter(BufferingGroup.Diff.Test.p < p_threshold) %>%
  filter(BufferingGroup.Sensitivity.Diff == "Lower") %>%
  slice_max(Corr.Sensitivity_Buffering, n = 12)

df_sensitivity_agg <- cellline_buf_agg %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never") %>%
  mutate(Buffering.CellLine = if_else(Buffering.CellLine.MeanNormRank > 0.5, "High", "Low"))

df_sensitivity_agg %>%
  inner_join(y = top_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  mutate(Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(Corr.Sensitivity_Buffering, 3),
                                                                 nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, desc(Corr.Sensitivity_Buffering))) %>%
  arrange(Buffering.CellLine.MeanNormRank) %>%
  jittered_boxplot(Label, Drug.MFI.Log2FC, Buffering.CellLine.MeanNormRank,
                   alpha = 3/4, jitter_width = 0.25) %>%
  save_plot("correlation_buffering_sensitivity_top.png", width = 300)

df_sensitivity_agg %>%
  inner_join(y = bot_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  mutate(Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(Corr.Sensitivity_Buffering, 3),
                                                                 nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, Corr.Sensitivity_Buffering)) %>%
  arrange(Buffering.CellLine.MeanNormRank) %>%
  jittered_boxplot(Label, Drug.MFI.Log2FC, Buffering.CellLine.MeanNormRank,
                   alpha = 3/4, jitter_width = 0.25) %>%
  save_plot("correlation_buffering_sensitivity_bot.png", width = 300)

df_sensitivity_agg %>%
  semi_join(y = top_sensitivity, by = "Drug.ID") %>%
  signif_violin_plot(Buffering.CellLine, Drug.MFI.Log2FC, Drug.Name) %>%
  save_plot(paste0("buffering_sensitivity_distributions_top.png"), width = 350)

df_sensitivity_agg %>%
  semi_join(y = bot_sensitivity, by = "Drug.ID") %>%
  signif_violin_plot(Buffering.CellLine, Drug.MFI.Log2FC, Drug.Name) %>%
  save_plot(paste0("buffering_sensitivity_distributions_bot.png"), width = 350)

# Scatter- & Violin-plot visualization of selected drugs
selected_drugs <- c("REGORAFENIB", "MPS1-IN-5", "SECLIDEMSTAT", "ATIPRIMOD", "INARIGIVIR", "C-021", "G-749",
                       "NIMORAZOLE", "GELDANAMYCIN", "ZOTAROLIMUS", "LERCANIDIPINE")

for (drug in selected_drugs) {
  df_sensitivity_agg %>%
    filter(Drug.Name == drug) %>%
    scatter_plot_reg_corr(Buffering.CellLine.MeanNormRank, Drug.MFI.Log2FC,
                          point_size = 2, title_prefix = drug) %>%
    save_plot(paste0("buffering_sensitivity_", drug, "_scatterplot.png"))
}

df_sensitivity_agg %>%
  filter(Drug.Name %in% selected_drugs) %>%
  signif_violin_plot(Buffering.CellLine, Drug.MFI.Log2FC, Drug.Name) %>%
  save_plot(paste0("buffering_sensitivity_distributions_selected.png"), width = 350)

# Get mechanism of action and target of top drugs
top_meta <- drug_meta %>%
  semi_join(y = top_sensitivity, by = "Drug.ID")

bot_meta <- drug_meta %>%
  semi_join(y = bot_sensitivity, by = "Drug.ID")

selected_meta <- drug_meta %>%
  filter(Drug.Name %in% selected_drugs)

# Plot correlation by mechanism of action
moa_corr <- cellline_buf_agg %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never") %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID") %>%
  separate_longer_delim(Drug.MOA, delim = ", ") %>%
  mutate(Drug.MOA = trimws(Drug.MOA)) %>%
  group_by(Drug.MOA) %>%
  summarize(
    # ToDo: Avoid calculating correlation twice
    Corr.Sensitivity_Buffering = cor.test(Buffering.CellLine.MeanNormRank, Drug.MFI.Log2FC,
                                          method = "spearman")$estimate[["rho"]],
    Corr.p = cor.test(Buffering.CellLine.MeanNormRank, Drug.MFI.Log2FC,
                      method = "pearson")$p.value
  )

moa_corr %>%
  filter(Corr.p < p_threshold) %>%
  filter(Corr.Sensitivity_Buffering < -0.15) %>%
  slice_min(Corr.Sensitivity_Buffering, n = 20) %>%
  mutate(Drug.MOA = fct_reorder(Drug.MOA, Corr.Sensitivity_Buffering, .desc = TRUE),
         `-log10(p)` = -log10(Corr.p)) %>%
  vertical_bar_chart(Drug.MOA, Corr.Sensitivity_Buffering, Corr.p,
                     value_range = c(-0.25, 0), line_intercept = 0, bar_label_shift = 0.18,
                     category_lab = "Mechanism of Action", value_lab = "Correlation Buffering-Sensitivity",
                     color_lab = "p-value") %>%
  save_plot("mechanism_buffering_sensitivity_top.png")

# Plot correlation by drug target
target_corr <- cellline_buf_agg %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never") %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID") %>%
  separate_longer_delim(Drug.Target, delim = ", ") %>%
  mutate(Drug.Target = trimws(Drug.Target)) %>%
  group_by(Drug.Target) %>%
  summarize(
    # ToDo: Avoid calculating correlation twice
    Corr.Sensitivity_Buffering = cor.test(Buffering.CellLine.MeanNormRank, Drug.MFI.Log2FC,
                                          method = "spearman")$estimate[["rho"]],
    Corr.p = cor.test(Buffering.CellLine.MeanNormRank, Drug.MFI.Log2FC,
                      method = "pearson")$p.value
  )

target_corr %>%
  filter(Corr.p < p_threshold) %>%
  filter(Corr.Sensitivity_Buffering < -0.15) %>%
  slice_min(Corr.Sensitivity_Buffering, n = 20) %>%
  mutate(Drug.Target = fct_reorder(Drug.Target, Corr.Sensitivity_Buffering, .desc = TRUE),
         `-log10(p)` = -log10(Corr.p)) %>%
  vertical_bar_chart(Drug.Target, Corr.Sensitivity_Buffering, Corr.p,
                     value_range = c(-0.25, 0), line_intercept = 0, bar_label_shift = 0.225,
                     category_lab = "Drug Target", value_lab = "Correlation Buffering-Sensitivity",
                     color_lab = "p-value") %>%
  save_plot("target_buffering_sensitivity_top.png")

# Heatmap plotting Drug Sensitivity, Cell Line Buffering and potential confounders
drug_confounder_heatmap <- function(df, desc = FALSE) {
  theme_settings <- theme(axis.line.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.key.size = unit(12, "points"),
                          legend.title = element_text(size = 10),
                          legend.text = element_text(size = 10))

  horizontal_legend <- theme(legend.direction = "horizontal",
                             legend.title = element_blank())

  df <- df %>%
    mutate(Drug.Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(Corr.Sensitivity_Buffering, 3),
                                                                        nsmall = 3, scientific = FALSE), ")")) %>%
    mutate(Drug.Label = fct_reorder(Drug.Label, Corr.Sensitivity_Buffering, .desc = desc),
           CellLine.Name = fct_reorder(CellLine.Name, Buffering.CellLine.MeanNormRank),
           CellLine.WGD = as.factor(CellLine.WGD))

  heatmap_buffering <- df %>%
    distinct(CellLine.Name, Buffering.CellLine.MeanNormRank) %>%
    ggplot() +
    aes(x = CellLine.Name, y = "Buffering Score", fill = Buffering.CellLine.MeanNormRank) +
    geom_raster() +
    scale_fill_viridis_c(option = "mako", labels = c(0, 0.5, 1), breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    cowplot::theme_minimal_hgrid() +
    theme_settings +
    horizontal_legend

  heatmap_drugs <- df %>%
    ggplot() +
    aes(x = CellLine.Name, y = Drug.Label, fill = Drug.MFI.Log2FC) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma", direction = -1) +
    cowplot::theme_minimal_hgrid() +
    theme_settings

  heatmap_wgd <- df %>%
    distinct(CellLine.Name, CellLine.WGD) %>%
    ggplot() +
    aes(x = CellLine.Name, y = "WGD Status", fill = CellLine.WGD) +
    geom_raster() +
    scale_fill_viridis_d(option = "viridis") +
    cowplot::theme_minimal_hgrid() +
    theme_settings +
    horizontal_legend

  heatmap_aneuploidy <- df %>%
    distinct(CellLine.Name, CellLine.AneuploidyScore) %>%
    ggplot() +
    aes(x = CellLine.Name, y = "Aneuploidy Score", fill = CellLine.AneuploidyScore) +
    geom_raster() +
    scale_fill_viridis_c(option = "rocket") +
    cowplot::theme_minimal_hgrid() +
    theme_settings +
    horizontal_legend

  cowplot::plot_grid(
    heatmap_buffering, NULL, heatmap_drugs, NULL, heatmap_wgd, NULL, heatmap_aneuploidy,
    labels = NULL, ncol = 1, align = "v", axis = "lr",
    rel_heights = c(2, -0.5, length(selected_drugs), -0.5, 2, -0.8, 2)
  )
}

cellline_buf_agg %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never") %>%
  inner_join(y = drug_dc_corr_agg_rank, by = c("Drug.ID", "Drug.Name")) %>%
  filter(Drug.Name %in% selected_drugs) %>%
  drug_confounder_heatmap(desc = TRUE) %>%
  save_plot("drug_confounder_heatmap_selected.png")

cellline_buf_agg %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never") %>%
  inner_join(y = top_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  drug_confounder_heatmap(desc = TRUE) %>%
  save_plot("drug_confounder_heatmap_top.png")

cellline_buf_agg %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never") %>%
  inner_join(y = bot_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  drug_confounder_heatmap() %>%
  save_plot("drug_confounder_heatmap_bot.png")


# ToDo: Cleanup Code
