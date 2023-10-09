library(ggplot2)
library(dplyr)
library(rlang)
library(EnhancedVolcano)
library(psych)
library(corrplot)

here::i_am("DosageCompensationFactors.Rproj")

vertical_bar_chart <- function(df, category_col, value_col,
                               error_low_col = NULL, error_high_col = NULL,
                               value_range = c(0.45, 0.6), break_steps = 0.05,
                               line_intercept = 0.5, bar_label_shift = 0.005,
                               title = NULL, category_lab = NULL, value_lab = NULL) {
  df %>%
    ggplot() +
    aes(x = { { category_col } }, y = { { value_col } },
        label = format(round({ { value_col } }, 3), nsmall = 3)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = line_intercept) +
  { if (!quo_is_null(enquo(error_low_col)) & !quo_is_null(enquo(error_high_col)))
    geom_pointrange(aes(x = { { category_col } }, y = { { value_col } },
                        ymin = { { error_low_col } }, ymax = { { error_high_col } }),
                    colour = "orange", fatten = 1) } +
    geom_text(color = "white", y = value_range[1] + bar_label_shift) +
    scale_y_continuous(breaks = seq(value_range[1], value_range[2], break_steps)) +
    ggtitle(title) +
    xlab(category_lab) +
    ylab(value_lab) +
    coord_flip(ylim = c(value_range[1], value_range[2])) +
    theme_light()
}

plot_text_col <- function(df, x_col, label_col, align = "center") {
  align_param <- case_when(
    align == "left" ~ 0,
    align == "center" ~ 0.5,
    align == "right" ~ 1,
    TRUE ~ 0.5
  )

  df %>%
    ggplot() +
    aes(x = { { x_col } }, y = align_param, label = { { label_col } }) +
    geom_text(color = "black", hjust = align_param) +
    xlab("") +
    ylab("") +
    coord_flip(ylim = c(0, 1)) +
    theme_void() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

violin_plot <- function(df, x, y) {
  df %>%
    ggplot() +
    aes(x = { { x } }, y = { { y } }) +
    geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75),
                color = "#4080DB") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_correlation <- function(df, method = "spearman") {
  cor_matrix <- psych::corr.test(df, method = method,
                                 adjust = "none")

  corrplot(cor_matrix$r, p.mat = cor_matrix$p,
         type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45,
         pch.cex = 1, pch.col = "darkgrey")
}


plot_volcano_buffered <- function(df, ratio_col, signif_col, label_col, class_col,
                                  color_mapping = NULL,
                                  value_threshold = log2fc_threshold, signif_threshold = p_threshold,
                                  title = NULL, subtitle = NULL) {
  buffering_levels <- c("N.S.", "Scaling", "Buffered", "Anti-Scaling")
  color_palette <- c("darkgrey", "#0E992C", "#E6C45A", "#E64356")
  names(color_palette) <- buffering_levels
  color_mapping <- scale_colour_manual(name = "Buffering Class", values = color_palette)

  df %>%
    mutate(Buffering.Class = if_else(TTest.p.adjusted > p_threshold, "N.S.", Buffering.Class)) %>%
    mutate(Buffering.Class = factor({ { class_col } }, levels = buffering_levels)) %>%
    mutate(Label = if_else(Buffering.Class %in% c("Buffered", "Anti-Scaling"),
                           { { label_col } }, NA)) %>%
    plot_volcano({ { ratio_col } }, { { signif_col } }, Label, Buffering.Class,
                 color_mapping = color_mapping,
                 value_threshold = value_threshold, signif_threshold = signif_threshold,
                 title = title, subtitle = subtitle)
}

plot_volcano <- function(df, value_col, signif_col, label_col, color_col,
                         color_mapping = NULL,
                         value_threshold = log2fc_threshold, signif_threshold = p_threshold,
                         title = NULL, subtitle = NULL) {
  df %>%
    mutate(`-Log10(p)` = -log10({ { signif_col } })) %>%
    ggplot() +
    aes(x = { { value_col } }, y = `-Log10(p)`,
        label = { { label_col } }, color = { { color_col } }) +
    geom_point(alpha = 0.5, size = 1) +
    color_mapping +
    geom_hline(yintercept = -log10(signif_threshold),
               linetype = "dashed", color = "black") +
    geom_label_repel(min.segment.length = 0.01, label.size = 0.15,
                     seed = 42, max.iter = 30000, max.time = 1.5,
                     point.padding = 0.3, label.padding = 0.3, box.padding = 0.3,
                     force = 2, max.overlaps = 20)
}

scatter_plot_regression <- function(df, x_col, y_col, formula,
                                    label_coords = c(0, 0), title = NULL) {
  df <- df %>%
    select({ { x_col } }, { { y_col } }) %>%
    drop_na()

  regression <- lm(formula, df)
  pred.int <- predict(regression, interval = "prediction")
  regression_summary <- summary(regression)
  df <- cbind(df, pred.int)
  slope <- regression$coefficients[[quo_name(enquo(x_col))]]
  intercept <- regression$coefficients[["(Intercept)"]]

  regression_plot <- df %>%
    ggplot() +
    aes(x = { { x_col } }, y = { { y_col } }) +
    geom_point(alpha = 0.3, size = 0.3) +
    geom_density_2d(color = "white", alpha = 0.6, linewidth = 0.4) +
    stat_smooth(method = lm, color = "blue") +
    geom_line(aes(y = lwr), color = "red", linetype = "dashed") +
    geom_line(aes(y = upr), color = "red", linetype = "dashed") +
    ggplot2::annotate("text", x = label_coords[1], y = label_coords[2], color = "blue",
                      label = paste("y =", format(round(slope, 5), nsmall = 5),
                                    "* x +", format(round(intercept, 5), nsmall = 5),
                                    ", R² = ", format(round(regression_summary$r.squared, 5), nsmall = 5)
                      )) +
    xlab(quo_name(enquo(x_col))) +
    ylab(quo_name(enquo(y_col))) +
    ggtitle(title)

  return(regression_plot)
}

waterfall_plot <- function(df, value_col, rank_col, label_col) {
  xlim <- c(0, max(df[[quo_name(enquo(rank_col))]]))
  label_nudge_x <- floor(xlim[2] / 5)

  ylim1 <- c(min(df[[quo_name(enquo(value_col))]]),
             (df %>% filter({ { rank_col } } == label_nudge_x))[[quo_name(enquo(value_col))]])
  ylim2 <- c((df %>% filter({ { rank_col } } == xlim[2] - label_nudge_x))[[quo_name(enquo(value_col))]],
             max(df[[quo_name(enquo(value_col))]]))

  df %>%
    ggplot() +
    aes(x = { { rank_col } }, y = { { value_col } }, label = { { label_col } }) +
    geom_hline(yintercept = 0, color = "red") +
    geom_point(size = 0.3) +
    geom_text_repel(data = df %>% slice_min({ { rank_col } }, n = 5),
                    xlim = xlim, ylim = ylim1, direction = "y", nudge_x = label_nudge_x,
                    seed = 42, color = "darkblue") +
    geom_text_repel(data = df %>% slice_max({ { rank_col } }, n = 5),
                    xlim = xlim, ylim = ylim2, direction = "y", nudge_x = -label_nudge_x,
                    seed = 42, color = "darkred")
}

save_plot <- function(plot, filename, dir = plots_dir,
                      height = 200, width = 200, dpi = 300) {
  ggsave(here(dir, filename), plot = plot,
       height = height, width = width, units = "mm", dpi = dpi)
  return(plot)
}


plot_rocs <- function(rocs) {
  roc_data <- data.frame()

  # Extract relevant information from ROC objects into data frame
  for (name in names(rocs)) {
    current_roc <- rocs[[name]]
    roc_df <- data.frame(
      Name = rep(name, length(current_roc$specificities)),
      Specificity = current_roc$specificities,
      Sensitivity = current_roc$sensitivities,
      AUC = rep(auc(current_roc), length(current_roc$specificities))
    )
    roc_data <- rbind(roc_data, roc_df)
  }

  df_label <- roc_data %>%
    distinct(Name, AUC) %>%
    arrange(desc(Name)) %>%
    mutate(x = 0.1,
           y = seq(0.05, 0.75, 0.05)[seq_along(unique(roc_data$Name))],
           AUC = paste0("AUC = ", format(round(AUC, 3), nsmall = 3)))

  plot <- roc_data %>%
    arrange(Sensitivity) %>%
    ggplot() +
    aes(x = Specificity, y = Sensitivity, color = Name) +
    geom_abline(slope = 1, intercept = 1, color = "grey") +
    geom_line() +
    geom_label(data = df_label, mapping = aes(color = Name, label = AUC, x = x, y = y)) +
    scale_x_reverse(limits = c(1, 0)) +
    labs(x = "Specificity", y = "Sensitivity", color = "Model")

  return(plot)
}

print_signif <- function(p, digits = 4) {
  paste0("p ", if_else(p < 10^(-digits),
                      paste0("< ", format(10^(-digits), nsmall = digits, scientific = FALSE)),
                      paste0("= ", format(round(p, digits), nsmall = digits, scientific = FALSE))))
}

map_signif <- function (p) {
  case_when(
      p < 0.0001 ~ "***",
      p < 0.001 ~ "**",
      p < 0.01 ~ "*",
      TRUE ~ "N.S."
    )
}