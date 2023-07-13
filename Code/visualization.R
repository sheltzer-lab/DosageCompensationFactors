library(ggplot2)
library(dplyr)
library(rlang)

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