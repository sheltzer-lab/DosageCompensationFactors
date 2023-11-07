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
cellline_buf_filtered_procan <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_procan.parquet"))
cellline_buf_filtered_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_depmap.parquet"))
cellline_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))

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
                        method = "pearson")$p.value,
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
    append = TRUE, sep = "\n")


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
  mutate(Label = paste0(Drug.Name, " (ρ = ", format(round(Corr.Sensitivity_Buffering, 3),
                                                    nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, desc(Corr.Sensitivity_Buffering))) %>%
  arrange(Buffering.CellLine.MeanNormRank) %>%
  jittered_boxplot(Label, Drug.MFI.Log2FC, Buffering.CellLine.MeanNormRank,
                   alpha = 3/4, jitter_width = 0.25) %>%
  save_plot("correlation_buffering_sensitivity_top.png", width = 300)

df_sensitivity_agg %>%
  inner_join(y = bot_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  mutate(Label = paste0(Drug.Name, " (ρ = ", format(round(Corr.Sensitivity_Buffering, 3),
                                                   nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, Corr.Sensitivity_Buffering)) %>%
  arrange(Buffering.CellLine.MeanNormRank) %>%
  jittered_boxplot(Label, Drug.MFI.Log2FC, Buffering.CellLine.MeanNormRank,
                   alpha = 3/4, jitter_width = 0.25) %>%
  save_plot("correlation_buffering_sensitivity_bot.png", width = 300)

# Scatter- & Violin-plot visualization of selected drugs
interesting_drugs <- c("REGORAFENIB", "IDASANUTLIN", "SECLIDEMSTAT", "ATIPRIMOD", "INARIGIVIR", "C-021", "G-749")

for (drug in interesting_drugs) {
  df <- df_sensitivity_agg %>%
    filter(Drug.Name == drug)

  df %>%
    scatter_plot_regression(Buffering.CellLine.MeanNormRank,
                          Drug.MFI.Log2FC,
                          Drug.MFI.Log2FC ~ Buffering.CellLine.MeanNormRank,
                          point_size = 2, label_coords = c(0.5, -2)) %>%
    save_plot(paste0("buffering_sensitivity_", drug, "_scatterplot.png"))
}

df_sensitivity_agg %>%
  filter(Drug.Name %in% interesting_drugs) %>%
  signif_violin_plot(Buffering.CellLine, Drug.MFI.Log2FC, Drug.Name) %>%
  save_plot(paste0("selected_buffering_sensitivity_distributions.png"))
