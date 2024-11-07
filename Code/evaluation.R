report_list2lines <- function(report) {
  lines <- NULL
  for (item in names(report)) {
    if (is.list(report[[item]])) {
      lines <- c(lines,
                 sprintf("\n--- %s ---\n", item),
                 report_list2lines(report[[item]]),
                 "------\n")
    }
    else {
      lines <- c(lines, sprintf("%s: %s\n", item, report[[item]]))
    }
  }
  return(lines)
}

write_report <- function(report, report_file) {
  file.remove(report_file)
  lines <- c("===== Report =====\n",
             format(Sys.time(), "%a %b %d %X %Y"),
             "\n\n",
             report_list2lines(report))
  writeLines(lines, con = report_file, sep = "")
}

df2reportlist <- function(df, collapse_cols = "\t\t") {
  result <- df %>%
    as.list() %>%
    transpose() %>%
    lapply(\(lst) paste(lst, collapse = collapse_cols))
  result <- append(result, paste(colnames(df), collapse = collapse_cols), after = 0)
  names(result) <- c(0, rownames(df))   # Using "Columns" as label of first row may cause row to shift
  return(result)
}

# Multivariate
evaluate_model <- function(model, test_set, dir = NULL, filename = NULL, cv_eval = FALSE) {
  require(dplyr)
  require(pROC)
  require(caret)

  test_predicted_reponse <- NULL
  test_performance_metrics <- NULL

  if (cv_eval == TRUE) {
    # Use results from cross validation generated during training for evaluation
    # Get ROC curve across all folds of k-fold CV for model with best parameters
    cv_best_preds <- model$pred %>%
      semi_join(y = model$bestTune)

    model_roc <- roc(response = cv_best_preds$obs, predictor = cv_best_preds$Buffered, na.rm = TRUE)
  } else {
    # Use separate test set for evaluation
    test_predicted_prob <- predict(model, test_set, type = "prob")
    model_roc <- roc(response = test_set$buffered, predictor = as.numeric(test_predicted_prob[, "Buffered"]), na.rm = TRUE)

    # Add model perfomance metrics
    lvl <- c("Buffered", "Scaling")
    test_predicted_reponse <- test_predicted_prob %>%
      mutate(
        Prediction = factor(if_else(Buffered > 0.5, "Buffered", "Scaling"), levels = lvl),
        Response = factor(test_set$buffered, levels = lvl)
      )

    test_performance_metrics <- list(
      precision = precision(data = test_predicted_reponse$Prediction, reference = test_predicted_reponse$Response),
      recall = recall(data = test_predicted_reponse$Prediction, reference = test_predicted_reponse$Response),
      F1 = F_meas(data = test_predicted_reponse$Prediction, reference = test_predicted_reponse$Response)
    )
  }

  eval_results <- list(roc = model_roc,
                       predictedResponse = test_predicted_reponse,
                       performanceMetrics = test_performance_metrics)

  if (is.null(filename) || is.null(dir))
    return(eval_results)

  png(here(dir, filename),
      width = 200, height = 200, units = "mm", res = 200)
  plot(eval_results$roc, print.thres = "best", print.thres.best.method = "closest.topleft",
       print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.1)
  dev.off()

  return(eval_results)
}

explain_model <- function(model, dir, filename = NULL) {
  require(dplyr)
  require(caret)

  # See randomForest::importance() and xgboost::xgb.importance() (rescaled, 0-100)
  ## rf: Mean Decrease in Gini Impurity
  ## xgbLinear: Weight of linear coefficients of feature
  ## xgbTree: fractional contribution of each feature
  if (model$method %in% c("rf", "xgbLinear", "xgbTree")) {
    plot <- caret::varImp(model)$importance %>%
      rename(Importance = "Overall") %>%
      tibble::rownames_to_column(var = "Factor") %>%
      mutate(Factor = factor(Factor, levels = Factor[order(Importance)])) %>%
      arrange(Factor) %>%
      vertical_bar_chart(Factor, Importance,
                         value_range = c(0, 100), bar_label_shift = 1,
                         line_intercept = 0, break_steps = 10,
                         category_lab = "Factor", value_lab = "Relative Importance",
                         title = paste0("Model Importance (", model$method, ")"))

    if (is.null(filename))
      return(plot)

    save_plot(plot, filename, dir = dir)
  }
}
