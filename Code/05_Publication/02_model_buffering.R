library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(viridisLite)
library(mskcc.oncotree)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

tumor_types <- get_tumor_types()

df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))

cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
cellline_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))
cellline_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(Model.ID, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

df_model_procan <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20240110.csv")) %>%
  rename(CellLine.SangerModelId = "model_id") %>%
  inner_join(y = df_celllines, by = "CellLine.SangerModelId") %>%
  mutate(OncotreeCode = nci_to_oncotree(cancer_type_ncit_id, expand = TRUE)$oncotree_code) %>%
  select(-CellLine.Name)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.DepMapModelId = "ModelID") %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId") %>%
  select(-CellLine.Name)

df_model_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))

intersect(unique(df_model_depmap$OncotreeCode), unique(df_cptac$Model.CancerType))
intersect(unique(tumor_types$oncotree_code), unique(df_cptac$Model.CancerType))

df_depmap <- cellline_buf_depmap %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never")

df_procan <- cellline_buf_procan %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never")

df_cptac <- cellline_buf_cptac %>%
  inner_join(y = df_model_cptac, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  rename(OncotreeCode = Model.CancerType)

get_oncotree_parent <- function(df_tumor_types = NULL, start_code = NULL, target_level = 3) {
  require(mskcc.oncotree)

  if (is.null(start_code) || is.na(start_code)) return(NA)

  if (is.null(df_tumor_types))
    df_tumor_types <- mskcc.oncotree::get_tumor_types()

  current <- df_tumor_types %>%
    filter(oncotree_code == start_code)

  if (nrow(current) == 0) return(NA)

  if (current$level <= target_level) return(current$oncotree_code)
  else get_oncotree_parent(df_tumor_types, current$parent, target_level)
}

# TODO: Control for aneuploidy score
df_cancer_heatmap <- bind_rows(df_depmap, df_procan, df_cptac) %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 1))) %>%
  group_by(Dataset, OncotreeCode) %>%
  summarize(Mean.BR = mean(Model.Buffering.Ratio, na.rm = TRUE),
            Mean.AS = mean(CellLine.AneuploidyScore, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(OncotreeCode = fct_reorder(OncotreeCode, Mean.BR)) %>%
  drop_na(OncotreeCode) %>%
  complete(Dataset, OncotreeCode)

cancer_heatmap_br <- df_cancer_heatmap %>%
  ggplot() +
  aes(y = Dataset, x = OncotreeCode, fill = Mean.BR) +
  geom_tile(aes(color = Mean.AS), alpha = 0) +
  geom_tile() +
  scale_color_viridis_c(na.value = default_color, option = "rocket") +
  scale_fill_viridis_c(na.value = default_color) +
  scale_x_discrete(position = "top") +
  theme_void() +
  labs(x = NULL, fill = "Mean Buffering Ratio", color = "Mean Aneuploidy") +
  guides(fill = guide_colourbar(order = 1),
         colour = guide_colourbar(order = 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        legend.key.size = unit(16, "points"),
        legend.box = "horizontal")

cancer_heatmap_as <- df_cancer_heatmap %>%
  group_by(OncotreeCode) %>%
  summarize(Mean.AS = mean(Mean.AS, na.rm = TRUE),
            Dataset = "Mean Aneuploidy") %>%
  ggplot() +
  aes(y = Dataset, x = OncotreeCode, fill = Mean.AS) +
  geom_tile() +
  scale_fill_viridis_c(na.value = default_color, option = "rocket") +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

cowplot::plot_grid(cancer_heatmap_br, cancer_heatmap_as,
                   nrow = 2, ncol = 1, align = "v", axis = "tb",
                   rel_heights = c(1, 0.2))
