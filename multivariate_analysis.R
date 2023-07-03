library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)
library(pROC)
library(randomForest)


here::i_am("DosageCompensationFactors.Rproj")

source(here("parameters.R"))
source(here("buffering_ratio.R"))

output_data_dir <- output_data_base_dir
goncalves_plots_dir <- here(plots_base_dir, "Univariate", "Goncalves")
depmap_plots_dir <- here(plots_base_dir, "Univariate", "DepMap")

dir.create(output_data_dir, recursive = TRUE)
dir.create(goncalves_plots_dir, recursive = TRUE)
dir.create(depmap_plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_goncalves <- read_parquet(here(output_data_dir, "expression_buffering_goncalves.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

dc_factor_cols <- c(
  "Protein-Protein Interactions", "Protein Half-Life", "Protein Complexes (CORUM)",
  "Mean 3'-UTR Length", "Mean 5'-UTR Length",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites",
  "Methylation Sites", "Acetylation Sites", "Regulatory Sites",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate",
  "Translation Rate", "Protein Length", "mRNA Length",
  "Intrinsic Protein Disorder", "Low Complexity Score", "Homology Score",
  "Loops In Protein Score", "Protein Polyampholyte Score", "Protein Polarity",
  "Non-Exponential Decay Delta", "Mean mRNA Decay Rate"
)


# === Separate Training & Test Dataset ===

training_set_ratio <- 0.8
set.seed(42)

## Chromosome Gain

shuffle_rows <- function (df) {
  df[sample(nrow(df), replace = FALSE),]
}

impute_na <- function (df) {
  df %>%
    mutate_if(is.numeric, \(x) replace_na(x, median(x, na.rm = TRUE)))
}

expr_buf_goncalves_gain <- expr_buf_goncalves %>%
  filter(ChromosomeArm.CNA > 0) %>%
  filter(Buffering.ChrArmLevel.Class == "Buffered" | Buffering.ChrArmLevel.Class == "Scaling") %>%
  mutate(Buffered = if_else(Buffering.ChrArmLevel.Class == "Buffered", 1, 0)) %>%
  mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
  drop_na(Buffered) %>%
  select("UniqueId", Buffered, all_of(dc_factor_cols)) %>%
  impute_na() %>%
  slice_sample(n = 10000) %>%
  shuffle_rows()

train_size <- ceiling(nrow(expr_buf_goncalves_gain) * training_set_ratio)
test_size <- nrow(expr_buf_goncalves_gain) - train_size

train_set <- expr_buf_goncalves_gain[1:train_size,]
test_set <- expr_buf_goncalves_gain[(train_size+1):(train_size+test_size),]

# === Build Random Forest Model ===

predictors <- train_set %>% select(-Buffered, -UniqueId)
responses <- train_set$Buffered

rf_model <- randomForest(
  x = predictors,
  y = responses,
  na.action = na.roughfix
)

rf_model
which.min(rf_model$err.rate)
rf_model$err.rate[which.min(rf_model$err.rate)]
varImpPlot(rf_model)

test_predictors <- test_set %>% select(-Buffered, -UniqueId)
test_responses <- test_set$Buffered
test_predicted <- predict(rf_model, newdata=test_predictors)
test_predicted_prob <- predict(rf_model, newdata=test_predictors, type = "prob")

success <- test_responses == test_predicted
success_rate <- length(success[success]) / length(success)

rf_roc <- roc(response = test_responses, predictor = as.numeric(test_predicted_prob[,1]))
plot(rf_roc)
auc(rf_roc)

# ToDo: Read https://stackoverflow.com/questions/30366143/how-to-compute-roc-and-auc-under-roc-after-training-using-caret-in-r