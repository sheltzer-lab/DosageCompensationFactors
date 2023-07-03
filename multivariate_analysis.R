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
library(caret)
library(randomForest)
library(MLeval)


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

training_set_ratio <- 0.7
set.seed(42)

## Chromosome Gain

shuffle_rows <- function(df) {
  df[sample(nrow(df), replace = FALSE),]
}

impute_na <- function(df) {
  df %>%
    mutate_if(is.numeric, \(x) replace_na(x, median(x, na.rm = TRUE)))
}

balanced_sample <- function(df, class_col, n_per_class) {
  df %>%
    group_by({ { class_col } }) %>%
    slice_sample(n = n_per_class) %>%
    ungroup()
}

rebalance_binary <- function(df, class_col, target_balance = 0.5) {
  total_size <- nrow(df)
  new_samples <- df %>%
    count({ { class_col } }) %>%
    mutate(Ratio = n/total_size) %>%
    mutate(Samples = c(ceiling((target_balance * min(n)) / (1 - target_balance)), min(n))) %>%
    select(-n, -Ratio)

  df %>%
    left_join(y = new_samples, by = quo_name(enquo(class_col))) %>%
    group_by({ { class_col } }) %>%
    group_map(~slice_sample(.x, n = unique(.x$Samples)), .keep = TRUE) %>%
    bind_rows() %>%
    select(-Samples)
}

expr_buf_goncalves_gain <- expr_buf_goncalves %>%
  filter(ChromosomeArm.CNA > 0) %>%
  filter(Buffering.ChrArmLevel.Class == "Buffered" | Buffering.ChrArmLevel.Class == "Scaling") %>%
  mutate(Buffered = factor(Buffering.ChrArmLevel.Class, levels = c("Scaling", "Buffered"))) %>%
  drop_na(Buffered) %>%
  select(Buffered, all_of(dc_factor_cols)) %>%
  impute_na() %>%
  rebalance_binary(Buffered, target_balance = 0.7) %>%
  shuffle_rows() %>%
  janitor::clean_names()


train_size <- ceiling(nrow(expr_buf_goncalves_gain) * training_set_ratio)
test_size <- nrow(expr_buf_goncalves_gain) - train_size

train_set <- expr_buf_goncalves_gain[1:train_size,]
test_set <- expr_buf_goncalves_gain[(train_size+1):(train_size+test_size),]

# === Build Random Forest Model ===
tc <- trainControl(method = "cv",
                   number = 2,
                   savePredictions = TRUE,
                   classProbs = TRUE,
                   verboseIter = TRUE)

forest.model <- caret::train(buffered ~., data = train_set, method = "rf", trControl = tc)

forest.model$finalModel
forest.model.eval <- evalm(forest.model)
forest.model.importance <- varImp(forest.model)
varImpPlot(forest.model$finalModel)

test_predicted_prob <- predict(forest.model, test_set, type = "prob")
rf_roc <- roc(response = test_set$buffered, predictor = as.numeric(test_predicted_prob[,"Buffered"]))
plot(rf_roc, print.thres="best", print.thres.best.method="closest.topleft", print.auc = TRUE)
rf_coords <- coords(rf_roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
print(rf_coords)
auc(rf_roc)


# === Build Neural Network Model ===

tc <- trainControl(method = "cv",
                   number = 3,
                   savePredictions = TRUE,
                   classProbs = TRUE,
                   verboseIter = TRUE)

neural.model <- caret::train(buffered ~., data = train_set, method = "pcaNNet", trControl = tc)

neural.model$finalModel
test_predicted_prob_neural <- predict(neural.model, test_set, type = "prob")
rf_roc <- roc(response = test_set$buffered, predictor = as.numeric(test_predicted_prob_neural[,"Buffered"]))
plot(rf_roc, print.thres="best", print.thres.best.method="closest.topleft", print.auc = TRUE)