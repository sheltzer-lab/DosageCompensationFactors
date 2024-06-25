library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(magrittr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "analysis.R"))
source(here("Code", "visualization.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Screens", "Ubiquitination")

expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac.parquet"))
ubi_cptac <- read_parquet(here(output_data_dir, 'ubiquitination_cptac.parquet'))

ubi_cptac_processed <- ubi_cptac %>%
  mutate(Protein.Ubiquitination.Log2 = as.vector(scale(Protein.Ubiquitination.Log2))) %>%
  group_by(Model.ID, Gene.Symbol) %>%
  summarize(Protein.Ubiquitination.Log2 = mean(Protein.Ubiquitination.Log2, na.rm = TRUE), .groups = "drop")

ubi_buf_cptac <- expr_buf_cptac %>%
  inner_join(y = ubi_cptac_processed, by = c("Model.ID", "Gene.Symbol")) %>%
  tibble::column_to_rownames("UniqueId")

resid <- residuals(lm(Protein.Ubiquitination.Log2 ~ Protein.Expression.Normalized, ubi_buf_cptac))

ubi_buf_cptac <- ubi_buf_cptac %>%
  mutate(Protein.Ubiquitination.Adj = resid[rownames(ubi_buf_cptac)],
         Buffering.GeneLevel.Ratio = as.vector(scale(Buffering.GeneLevel.Ratio)))

ubi_buf_cptac %>%
  filter_cn_gain_abs() %$%
  cor.test(Protein.Ubiquitination.Adj, Buffering.GeneLevel.Ratio)

ubi_buf_corr <- ubi_buf_cptac %>%
  filter_cn_gain_abs() %>%
  drop_na(Protein.Ubiquitination.Adj, Buffering.GeneLevel.Ratio) %>%
  #filter(Buffering.GeneLevel.Class == "Buffered") %>%
  group_by(Gene.Symbol) %>%
  summarize(Corr = cor(Buffering.GeneLevel.Ratio, Protein.Ubiquitination.Adj, use = "na.or.complete"),
            Samples = min(sum(!is.na(Buffering.GeneLevel.Ratio)), sum(!is.na(Protein.Ubiquitination.Adj)))) %>%
  filter(Samples > 10)

ubi_buf_cptac %>%
  filter(Gene.Symbol %in% (ubi_buf_corr %>% filter(Corr > 0.4) %>% pull(Gene.Symbol))) %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Ratio, y = Protein.Ubiquitination.Adj, color = Gene.Symbol) +
  geom_point()

ubi_buf_cptac %>%
  filter(abs(Gene.CopyNumber - Gene.CopyNumber.Baseline) > 0.2) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.GeneLevel.Class) %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Class, y = Protein.Ubiquitination.Log2, color = Buffering.GeneLevel.Class) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("Buffered", "Scaling")),
              map_signif_level = print_signif, y_position = 0, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  facet_wrap(~Gene.CopyNumber.Event) +
  theme(legend.position = "") +
  scale_color_manual(values = c(Buffered = bidirectional_color_pal[1],
                                Scaling = bidirectional_color_pal[5],
                                `Anti-Scaling` = "black"))

# Model-level Analysis
# TODO: Move to cell line level analysis
model_buf_cptac <- expr_buf_cptac %>%
  filter_cn_diff(remove_between = c(-0.01, 0.02)) %>% # Remove noise from discontinuity points of buffering ratio
  filter(Buffering.GeneLevel.Ratio.Confidence > 0.3) %>%
  #  filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
  select(Model.ID, Gene.Symbol, Buffering.GeneLevel.Ratio) %>%
  drop_na() %>%
  group_by(Model.ID) %>%
  summarize(Buffering.Model.Ratio = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Observations = sum(!is.na(Buffering.GeneLevel.Ratio)),
            SD = sd(Buffering.GeneLevel.Ratio, na.rm = TRUE)) %>%
  filter(Observations > 10)

model_buf_cptac %>%
  split_by_quantiles(Buffering.Model.Ratio, target_group_col = "Buffering.Model.Group") %>%
  inner_join(y = ubi_cptac_processed, "Model.ID") %>%
  drop_na(Protein.Ubiquitination.Log2) %>%
  ggplot() +
  aes(x = Buffering.Model.Group, y = Protein.Ubiquitination.Log2,
      color = Buffering.Model.Group) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("High", "Low")),
              map_signif_level = print_signif, y_position = 0, size = 0.8,
              tip_length = 0, extend_line = -0.05, color = "black") +
  theme(legend.position = "") +
  scale_color_manual(values = c(High = bidirectional_color_pal[1],
                                Low = bidirectional_color_pal[5]))

# Univariate Analyis
library(pROC)

# TODO: Reference function from univariate analysis
establish_binary_classification <- function(df, buffering_class_col) {
  require(dplyr)

  df %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = ifelse({ { buffering_class_col } } == "Buffered", 1, 0)) %>%
    mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
    drop_na(Buffered)
}

roc_ubi_gain <- ubi_buf_cptac %>%
  filter_cn_gain_abs() %>%
  establish_binary_classification(Buffering.GeneLevel.Class) %>%
  select(Protein.Ubiquitination.Log2, Buffered) %$%
  roc(Buffered, Protein.Ubiquitination.Log2, na.rm = TRUE)

roc_ubi_loss <- ubi_buf_cptac %>%
  filter_cn_loss_abs() %>%
  establish_binary_classification(Buffering.GeneLevel.Class) %>%
  select(Protein.Ubiquitination.Log2, Buffered) %$%
  roc(Buffered, Protein.Ubiquitination.Log2, na.rm = TRUE)

auc(roc_ubi_gain)
auc(roc_ubi_loss)

roc_ubi_adj_gain <- ubi_buf_cptac %>%
  filter_cn_gain_abs() %>%
  establish_binary_classification(Buffering.GeneLevel.Class) %>%
  select(Protein.Ubiquitination.Adj, Buffered) %$%
  roc(Buffered, Protein.Ubiquitination.Adj, na.rm = TRUE)

roc_ubi_adj_loss <- ubi_buf_cptac %>%
  filter_cn_loss_abs() %>%
  establish_binary_classification(Buffering.GeneLevel.Class) %>%
  select(Protein.Ubiquitination.Adj, Buffered) %$%
  roc(Buffered, Protein.Ubiquitination.Adj, na.rm = TRUE)

auc(roc_ubi_adj_gain)
auc(roc_ubi_adj_loss)

df_rocs <- rocs_to_df(list(Gain_Ubi = roc_ubi_gain,
                           Gain_UbiAdj = roc_ubi_adj_gain,
                           Loss_Ubi = roc_ubi_loss,
                           Loss_UbiAdj = roc_ubi_adj_loss))

plot_rocs(df_rocs) +
  scale_color_manual(values = c(Gain_Ubi = categorical_color_pal[4],
                                Gain_UbiAdj = categorical_color_pal[3],
                                Loss_Ubi = categorical_color_pal[8],
                                Loss_UbiAdj = categorical_color_pal[7]))

# TODO: Use XGBTree model for prediction
