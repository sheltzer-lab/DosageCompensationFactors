library(ggplot2)
library(dplyr)
library(rlang)
library(EnhancedVolcano)

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
  df %>%
    cor(method = method) %>%
    corrplot(type = "upper", order = "hclust",
             tl.col = "black", tl.srt = 45)
}

plot_volcano <- function(df, value_col, signif_col, label_col,
                         value_threshold = 1.2, signif_threshold = 0.05,
                         title = NULL, subtitle = NULL,
                         alpha_low = 0.3, alpha_high = 0.6) {
  # Scale points according to absolute distance from value_threshold
  min_pointsize <- 1
  max_pointsize <- 6
  max_value <- max(df[[quo_name(enquo(value_col))]], na.rm = TRUE)
  min_value <- min(df[[quo_name(enquo(value_col))]], na.rm = TRUE)
  max_abs_value <- max(abs(min_value), abs(max_value))
  value_cutoff_dist <- (abs(df[[quo_name(enquo(value_col))]]) - value_threshold) / max_abs_value
  point_scaling <- ifelse(abs(df[[quo_name(enquo(value_col))]]) > value_threshold,
                          value_cutoff_dist * (max_pointsize - min_pointsize) + min_pointsize,
                          min_pointsize)

  # Define plot
  volc <- df %>%
    EnhancedVolcano(
      lab = df[[quo_name(enquo(label_col))]],
      x = quo_name(enquo(value_col)),
      y = quo_name(enquo(signif_col)),
      title = title,
      subtitle = subtitle,
      caption = bquote("Value cutoff: " ~ .(value_threshold) ~ "; Significance cutoff: " ~ .(signif_threshold)),
      pCutoff = signif_threshold,
      FCcutoff = value_threshold,
      labSize = 3,
      xlim = c(-max_abs_value - 0.1, max_abs_value + 0.1),
      ylim = c(0, max(-log10(df[[quo_name(enquo(signif_col))]]) + 0.1, na.rm = TRUE)),
      pointSize = point_scaling,
      col = c('darkgrey', 'darkgrey', 'pink', 'red'),
      colAlpha = ifelse(abs(df[[quo_name(enquo(value_col))]]) > value_threshold &
                          df[[quo_name(enquo(signif_col))]] < signif_threshold,
                        alpha_high,
                        alpha_low),
      drawConnectors = TRUE,
      widthConnectors = 0.2,
      arrowheads = FALSE
    )
}

scatter_plot_regression <- function(df, x_col, y_col, formula, label_coords = c(0, 0)) {
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
    ylab(quo_name(enquo(y_col)))

  return(regression_plot)
}

save_plot <- function(plot, filename, dir = plots_dir,
                      height = 200, width = 200, dpi = 300) {
  ggsave(here(dir, filename), plot = plot,
       height = height, width = width, units = "mm", dpi = dpi)
  return(plot)
}