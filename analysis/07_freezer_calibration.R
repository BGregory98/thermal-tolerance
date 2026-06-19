# analysis/07_freezer_calibration.R -------------------------------------
#
# Freezer cool-down calibration — supplementary in the revised MS
# (per Reviewer 1's suggestion to move Fig. 6 to supplementary).
#
# Purpose: characterizes how water temperature in the treatment cups
# changes over 110 minutes in a –20 °C freezer, establishing that
# minutes 60–110 correspond to ~0 °C exposure (i.e., time at 0 °C =
# total freezer time – 60 min). This is the calibration used by
# analysis/04_cold.R to convert freezer time to "time at 0 °C".
#
# Replaces the "COLD RAMP IN FREEZER" section (lines ~295–426) of
# lightandtemp.R. Several exploratory plots from that script (volume
# comparison, cup-to-cup variation in different test runs) have been
# omitted here because they didn't make it into the manuscript.
# ------------------------------------------------------------------------

# --- Load the final 6-replicate run that's plotted in Figure 6 ---------
# File is six replicates × 110 minutes, one column per replicate.
path <- file.path(DATA_DIR, "ColdRamp_400.csv")

if (!file.exists(path)) {
  message("Freezer calibration file not found (", path,
          "). Skipping analysis/07.")
} else {

  cold <- readr::read_csv(path, show_col_types = FALSE)
  names(cold) <- c("mins", paste0("rep", 1:6))

  cold_long <- cold %>%
    tidyr::pivot_longer(
      cols      = dplyr::starts_with("rep"),
      names_to  = "rep",
      values_to = "temp_c"
    )

  cold_sum <- cold_long %>%
    dplyr::group_by(mins) %>%
    dplyr::summarise(
      mean_c = mean(temp_c, na.rm = TRUE),
      sd_c   = stats::sd(temp_c, na.rm = TRUE),
      .groups = "drop"
    )

  # --- Figure: mean ± SD over time, with cup-removal markers -----------
  # Dashed verticals at the six cup-removal times used in the actual
  # low-temperature trials (60, 70, 80, 90, 100, 110 minutes).
  removal_times <- seq(FREEZER_COOLDOWN_MIN,
                       FREEZER_COOLDOWN_MIN + COLD_TIME_MAX,
                       by = 10)

  p_freezer <- ggplot2::ggplot(cold_sum, ggplot2::aes(x = mins)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_c - sd_c, ymax = mean_c + sd_c),
      fill = "#6FB4E7", alpha = 0.5
    ) +
    ggplot2::geom_line(ggplot2::aes(y = mean_c),
                       linewidth = 1, color = "black") +
    ggplot2::geom_vline(xintercept = removal_times,
                        linetype = "dashed", color = "#888888") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::scale_x_continuous(breaks = seq(0, 110, by = 10)) +
    ggplot2::labs(
      x = "Time in Freezer (minutes)",
      y = "Water Temperature (\u00B0C) \u00B1SD"
    ) +
    theme_tt(base_size = 16)

  save_fig(p_freezer, "FigS_freezer_calibration",
           width = 9, height = 5)
}
