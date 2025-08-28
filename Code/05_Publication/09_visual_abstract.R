library(here)
library(arrow)
library(dplyr)
library(stringr)
library(ggplot2)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "buffering_ratio.R"))
source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

theme_set(theme_classic(base_size = base_size))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

# === Load Data ===
df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))
model_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(Model.ID, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

df_model_procan <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20240110.csv")) %>%
  rename(CellLine.SangerModelId = "model_id") %>%
  inner_join(y = df_celllines, by = "CellLine.SangerModelId") %>%
  select(-CellLine.Name)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.DepMapModelId = "ModelID") %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId") %>%
  select(-CellLine.Name)

df_model_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))

df_procan <- model_buf_procan %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, OncotreeCode), by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = df_model_procan %>% select(Model.ID, growth_properties), by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  rename(GrowthPattern = growth_properties)

df_cptac <- model_buf_cptac %>%
  inner_join(y = df_model_cptac, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  rename(OncotreeCode = Model.CancerType) %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

diff_exp_procan <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet"))

gsea_all <- read_parquet(here(output_data_dir, "gsea_hallmark_pan-cancer.parquet")) %>% mutate(Comparison = "Buffering")
gsea_all_as <- read_parquet(here(output_data_dir, "gsea_hallmark_pan-cancer_aneuploidy.parquet")) %>% mutate(Comparison = "Aneuploidy")

ssgsea_cptac <- read_parquet(here(output_data_dir, "ssgsea_unfolded_cptac.parquet"))

drug_mechanisms <- read_parquet(here(output_data_dir, "drug_effect_buffering_mechanism.parquet")) %>% mutate(Comparison = "Buffering")
drug_mechanisms_control <- read_parquet(here(output_data_dir, "drug_effect_buffering_mechanism_control.parquet")) %>% mutate(Comparison = "Aneuploidy")

# === Plot Panels ===
## Aneuploidy
panel_as_wgd <- df_procan %>%
  distinct(Model.ID, Model.Buffering.Ratio, CellLine.WGD, CellLine.AneuploidyScore) %>%
  mutate(WGD = factor(if_else(CellLine.WGD > 0, "WGD", "Non-WGD"), levels = c("Non-WGD", "WGD"))) %>%
  arrange(WGD) %>%
  ggplot() +
  aes(x = CellLine.AneuploidyScore, y = Model.Buffering.Ratio, color = WGD) +
  geom_point(size = 2) +
  stat_smooth(method = lm, color = highlight_colors[2]) +
  scale_color_manual(values = c("Non-WGD" = color_palettes$Missing, "WGD" = color_palettes$WGD[["WGD"]])) +
  labs(x = "Aneuploidy", y = "Buffering") +
  theme(legend.position = "none")

## Growth pattern
panel_growth_pattern <- df_procan %>%
  filter(GrowthPattern == "Adherent" | GrowthPattern == "Suspension") %>%
  ggplot() +
  aes(x = GrowthPattern, y = Model.Buffering.Ratio, color = GrowthPattern) +
  geom_boxplot(outliers = FALSE, size = 1, alpha = 0) +
  scale_color_manual(values = c(Adherent = discrete_color_pal2_bright[3],
                                Suspension = discrete_color_pal2_bright[2]), guide = NULL) +
  labs(x = "Growth Pattern", y = "Buffering")

## DiffExp
selected_genes <- c("CTSA", "MCM2", "UBE2N", "ITGA3", "ITGAV", "EGFR", "CTNNB1", "CCT2")

panel_diffexp <- diff_exp_procan %>%
  mutate(Label = if_else(Gene.Symbol %in% selected_genes, Gene.Symbol, NA)) %>%
  ggplot() +
  aes(x = Log2FC, y = -log10(Test.p.adj)) +
  geom_density_2d_filled(contour_var = "ndensity",
                         breaks = c(0, 1e-5, seq(0.1, 1, length.out = 7))) +
  geom_label(aes(label = Label, color = Significant)) +
  scale_color_manual(values = color_palettes$DiffExp, na.value = color_palettes$Missing) +
  scale_fill_manual(values = colorRampPalette(c("#FFFFFF", default_color))(10)) +
  theme(legend.position = "none") +
  labs(y = "-log10(p)", x = "Protein Log2FC")

## Gene Sets
selected_pathways <- c("HALLMARK_E2F_TARGETS" = "E2F Targets",
                       "HALLMARK_MTORC1_SIGNALING" = "mTORC1 Signaling",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Epithelial Mesenchymal\nTransition")

panel_gsea <- bind_rows(gsea_all, gsea_all_as) %>%
  filter(Dataset == "CPTAC") %>%
  filter(pathway %in% names(selected_pathways)) %>%
  mutate(pathway = str_replace_all(pathway, selected_pathways)) %>%
  mutate(pathway = factor(pathway, levels = rev(selected_pathways))) %>%
  arrange(pathway) %>%
  ggplot() +
  aes(x = Comparison, y = pathway, fill = NES) +
  geom_raster() +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = bidirectional_color_pal,
                       space = "Lab", limits = c(-2, 2), oob = scales::squish) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = base_size) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0))

## UPR
ssgsea_cptac_filtered <- df_cptac %>%
  filter(Observations > 500) %>%
  select(starts_with("Model"), Observations) %>%
  inner_join(y = ssgsea_cptac, by = "Model.ID")

panel_upr <- ssgsea_cptac_filtered %>%
  ggplot() +
  aes(x = Model.Buffering.Ratio, y = HALLMARK_UNFOLDED_PROTEIN_RESPONSE) +
  geom_density_2d_filled(bins = 9) +
  stat_smooth(method = lm, color = highlight_colors[2], fill = highlight_colors[2]) +
  scale_fill_manual(values = colorRampPalette(c("#FFFFFF", default_color))(9)) +
  labs(x = "Buffering", y = "Unfolded Protein Response") +
  theme(legend.position = "none")

## Drug Sensitivity Analysis
selected_drugs <- drug_mechanisms %>%
  filter(CommonEffect & !is.na(EffectiveIn)) %>%
  pull(Drug.MOA)

drugs <- bind_rows(
  drug_mechanisms %>% rename(Log2FC = "DrugEffect.Buffering.Group.Log2FC") %>% select(Drug.MOA, Log2FC, Comparison),
  drug_mechanisms_control %>% select(Drug.MOA, Log2FC, Comparison, Test.p.adj)
) %>%
  filter(Drug.MOA %in% selected_drugs) %>%
  mutate(SignifAneup = any(na.omit(Test.p.adj < p_threshold)),
         Drug.MOA = str_replace(Drug.MOA, "INHIBITOR", "Inhibitor"),
         .by = Drug.MOA) %>%
  filter(!SignifAneup)

panel_drugs <- drugs %>%
  ggplot() +
  aes(x = Comparison, y = Drug.MOA, fill = Log2FC) +
  geom_raster() +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = bidirectional_color_pal,
                       space = "Lab", limits = c(-0.3, 0.3), oob = scales::squish) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = base_size) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0))

# === Combine Figures ===
vis_abstract <- cowplot::plot_grid(panel_as_wgd, panel_growth_pattern,
                                   panel_diffexp, panel_gsea, panel_upr, panel_drugs, plot_buffering_ratio_classes(),
                                   nrow = 1, rel_widths = c(1, 1, 1, 1.1, 1, 1, 1.4))


cairo_pdf(here(plots_dir, "visual_abstract.pdf"), width = 23, height = 3.1)
vis_abstract
dev.off()