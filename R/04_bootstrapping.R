# R/04_bootstrapping.R --------------------------------------------------
#
# Bootstrap 95% CIs around the proportion alive for each combination of
# population (or species / region) × predictor value. Used to build
# survivorship curves and final-survival summary plots.
# ------------------------------------------------------------------------


#' Bootstrap survival CIs for every (group, predictor) combination
#'
#' For each unique combination of the grouping variable `by` and the
#' predictor (time or temp), draws `BOOT_N` bootstrap resamples of
#' `alive`, takes the mean of each, and returns 2.5 / 50 / 97.5
#' percentiles.
#'
#' @param data Output of `load_survival_data()`.
#' @param predictor "time" or "temp".
#' @param by Grouping column name (default "population").
#'
#' @return tibble with one row per (group × predictor value):
#'   <by>, <predictor>, lower, median, upper.
bootstrap_by_pop <- function(data,
                             predictor = c("time", "temp"),
                             by        = "population") {

  predictor <- match.arg(predictor)

  data %>%
    dplyr::group_by(.data[[by]], .data[[predictor]]) %>%
    dplyr::summarise(
      ci     = list(boot_ci(alive)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      lower  = vapply(ci, `[[`, numeric(1), 1),
      median = vapply(ci, `[[`, numeric(1), 2),
      upper  = vapply(ci, `[[`, numeric(1), 3)
    ) %>%
    dplyr::select(-ci)
}


#' Extract final-timepoint survival per group (convenience wrapper)
#'
#' @param boot_df Output of `bootstrap_by_pop()`.
#' @param predictor "time" or "temp".
#' @param final_value The predictor value to keep (e.g., 96 for sustained
#'   hot, 50 for cold, 46 for ramping hot).
final_survival <- function(boot_df,
                           predictor   = c("time", "temp"),
                           final_value) {
  predictor <- match.arg(predictor)
  boot_df[boot_df[[predictor]] == final_value, , drop = FALSE]
}
