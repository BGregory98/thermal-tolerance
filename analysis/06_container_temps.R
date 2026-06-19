# analysis/06_container_temps.R -----------------------------------------
#
# Environmental container-habitat temperatures measured in Baltimore,
# July–August 2023. Sun-exposed vs. shaded CDC gravid-trap buckets,
# logged every 5 minutes over 22 hours per replicate.
#
# Purpose in the manuscript (as revised): method-validation — confirms
# that the 38 °C sustained-heat treatment is ecologically plausible in
# real container habitats. Per Reviewer 1, this should NOT be framed as
# characterizing the thermal regimes of distant collection sites.
#
# Replaces the "TEMPERATURE" section (lines ~14–260) of lightandtemp.R.
# ------------------------------------------------------------------------

# --- Load ----------------------------------------------------------------
dat <- readr::read_csv(
  file.path(DATA_DIR, "balt713_trimmed.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(date = lubridate::mdy_hm(date)) %>%
  tidyr::pivot_longer(cols = c(sun, shade),
                      names_to  = "exposure",
                      values_to = "temp_c") %>%
  dplyr::mutate(over38 = temp_c > 38)


# --- Summary statistics --------------------------------------------------
summary_by_exposure <- dat %>%
  dplyr::group_by(exposure) %>%
  dplyr::summarise(
    n            = dplyr::n(),
    mean_c       = mean(temp_c, na.rm = TRUE),
    min_c        = min(temp_c,  na.rm = TRUE),
    max_c        = max(temp_c,  na.rm = TRUE),
    sd_c         = stats::sd(temp_c, na.rm = TRUE),
    range_c      = max_c - min_c,
    obs_over_38  = sum(over38, na.rm = TRUE),
    prop_over_38 = mean(over38, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n======== CONTAINER WATER TEMPERATURE (Baltimore, 2023) ========\n")
print(summary_by_exposure)

readr::write_csv(summary_by_exposure,
                 file.path(TAB_DIR, "container_temp_summary.csv"))


# --- Days reaching 38 °C (matches manuscript text) ----------------------
days_over_38 <- dat %>%
  dplyr::filter(over38) %>%
  dplyr::group_by(exposure, rep) %>%
  dplyr::summarise(
    minutes_over_38 = dplyr::n() * 5,  # 5-min logging interval
    .groups = "drop"
  )

cat("\n======== REPLICATES REACHING 38 \u00B0C ========\n")
print(days_over_38)


# --- Figure: ribbon of mean ± SD across replicates, sun vs shade -------
dat_sum <- dat %>%
  dplyr::group_by(measurement, exposure) %>%
  dplyr::summarise(
    mean_c = mean(temp_c, na.rm = TRUE),
    sd_c   = stats::sd(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# Hour labels along the 22-hour (4 PM -> 2 PM) logging window.
tick_labels <- c("4PM", "8PM", "12AM", "4AM", "8AM", "12PM")
tick_breaks <- seq(1, max(dat_sum$measurement), length.out = length(tick_labels))

p_container <- ggplot2::ggplot(dat_sum) +
  ggplot2::geom_ribbon(
    ggplot2::aes(x = measurement,
                 ymin = mean_c - sd_c, ymax = mean_c + sd_c,
                 fill = exposure),
    alpha = 0.3
  ) +
  ggplot2::geom_line(
    ggplot2::aes(x = measurement, y = mean_c, color = exposure),
    linewidth = 1
  ) +
  ggplot2::geom_hline(yintercept = 38, linetype = "dotted") +
  ggplot2::scale_x_continuous(breaks = tick_breaks, labels = tick_labels) +
  ggplot2::scale_color_manual(
    values = c(sun = "#E41A1C", shade = "#555555"),
    labels = c(sun = "Sun-exposed", shade = "Shaded"),
    name   = NULL
  ) +
  ggplot2::scale_fill_manual(
    values = c(sun = "#E41A1C", shade = "#555555"),
    labels = c(sun = "Sun-exposed", shade = "Shaded"),
    name   = NULL
  ) +
  ggplot2::labs(x = "Time", y = "Water Temperature (\u00B0C)") +
  theme_tt(base_size = 16) +
  ggplot2::theme(legend.position = c(0.12, 0.88))

save_fig(p_container, "FigS_container_temps",
         width = 10, height = 5)
