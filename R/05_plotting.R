# R/05_plotting.R -------------------------------------------------------
#
# Plot builders shared across the three thermal-tolerance experiments.
# Each function takes a tidy data frame plus axis labels and returns a
# ggplot object. Saving is handled by save_fig() (in R/00_setup.R).
# ------------------------------------------------------------------------


#' Two-panel survivorship curve (first vs. fourth instar)
#'
#' @param boot_df Tibble from `bootstrap_by_pop()` for one experiment,
#'   with an additional `instar` column (values "First" / "Fourth").
#' @param predictor Column in boot_df holding the x-axis variable.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label (default: standard survival label).
#' @param x_breaks Optional numeric vector of tick positions.
#'
#' @return A faceted ggplot.
survivorship_curve <- function(boot_df,
                               predictor,
                               x_lab,
                               y_lab = "Proportion Alive \u00B195% CI",
                               x_breaks = ggplot2::waiver()) {

  ggplot2::ggplot(boot_df) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = .data[[predictor]],
                   ymin = lower, ymax = upper,
                   fill = population),
      alpha = 0.4, linetype = "dotted", color = "black"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data[[predictor]], y = median,
                   color = population)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[[predictor]], y = median,
                   color = population)
    ) +
    ggplot2::facet_wrap(~ instar) +
    ggplot2::scale_color_manual(values = pop_colors, name = "Population") +
    ggplot2::scale_fill_manual( values = pop_colors, name = "Population") +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    theme_tt(base_size = 14)
}


#' Two-panel summary-point plot (e.g., LT50 or final-survival per pop)
#'
#' @param summary_df Tibble with columns: population, instar, and the y
#'   variable specified by `y`, plus lower/upper bounds.
#' @param y Column name for the point value (e.g., "lt50", "median").
#' @param y_lab Y-axis label.
#'
#' @return ggplot.
summary_point_plot <- function(summary_df,
                               y,
                               y_lab) {

  ggplot2::ggplot(summary_df) +
    ggplot2::geom_errorbar(
      ggplot2::aes(x = population, ymin = lower, ymax = upper),
      width = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = population, y = .data[[y]],
                   fill = population, shape = instar),
      size = 5, stroke = 0.3
    ) +
    ggplot2::facet_wrap(~ instar) +
    ggplot2::scale_fill_manual(values = pop_colors, guide = "none") +
    ggplot2::scale_shape_manual(values = c(First = 21, Fourth = 22),
                                guide = "none") +
    ggplot2::labs(x = "Population", y = y_lab) +
    theme_tt(base_size = 14)
}


# ----------------------------------------------------------------------
# lt50_vs_bioclim()
#
# Scatter of LT50 (with 95% CI error bars) vs. a single bioclimatic
# variable, one point per population, colored by pop_id, with a dashed
# linear regression line clipped to the range of the data.
#
# Inputs:
#   reg_df       — data frame with one row per population, containing at
#                  minimum the columns: pop_id, lt50, lower, upper, plus
#                  the bioclim column named in `predictor`.
#   predictor    — name of the bioclim column to use as the x-axis,
#                  e.g. "bio05".
#   x_lab        — x-axis label (the bioclim variable's pretty name).
#
# Returns a ggplot object.
# ----------------------------------------------------------------------
lt50_vs_bioclim <- function(reg_df,
                            predictor,
                            x_lab) {

  # Fit the regression so we can draw the line ourselves and clip it to
  # the data range — using stat_smooth would extend the line past the
  # observed predictor values, which we want to avoid given n=6.
  fit  <- lm(reg_df$lt50 ~ reg_df[[predictor]])
  xseq <- seq(min(reg_df[[predictor]]), max(reg_df[[predictor]]),
              length.out = 100)
  line_df <- data.frame(
    x = xseq,
    y = coef(fit)[1] + coef(fit)[2] * xseq
  )

  ggplot2::ggplot(reg_df) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = y),
      linetype = "dashed", color = "#888888", linewidth = 0.7
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(x = .data[[predictor]],
                   ymin = lower, ymax = upper),
      width = 0
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[[predictor]], y = lt50,
                   fill = pop_id),
      shape = 21, size = 5, stroke = 0.4
    ) +
    ggplot2::scale_fill_manual(values = pop_colors, name = "Population") +
    ggplot2::labs(x = x_lab,
                  y = "Median Lethal Temperature (\u00B0C)") +
    theme_tt(base_size = 14)
}