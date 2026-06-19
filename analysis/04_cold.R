# analysis/04_cold.R -----------------------------------------------------
#
# Figure 4: Low-temperature (freezing) exposure (first + fourth instar).
# Primary summary metric is proportion alive at 50 minutes at 0 °C
# (= 110 minutes total freezer time minus the 60-min cool-down).
#
# Analyzes the five main populations (POPS_MAIN from config.R). Uses
# population-level binomial GLMs; no region / species pooling, per
# decisions from the 2026-04-24 advisor meeting.
# ------------------------------------------------------------------------

# --- Load both instars ------------------------------------------------
c1 <- load_survival_data(
  path             = file.path(DATA_DIR, "cold_first.csv"),
  time_var         = "time",
  exclude_controls = FALSE,   # cold CSVs don't have tcycler
  drop_edge_time   = NULL,    # and don't need the 120-min trim
  pops             = POPS_MAIN
)
c4 <- load_survival_data(
  path             = file.path(DATA_DIR, "cold_fourth.csv"),
  time_var         = "time",
  exclude_controls = FALSE,
  drop_edge_time   = NULL,
  pops             = POPS_MAIN
)

# Convert total freezer time → time at ~0 °C.
c1$time <- c1$time - FREEZER_COOLDOWN_MIN
c4$time <- c4$time - FREEZER_COOLDOWN_MIN

c1$instar <- "First"
c4$instar <- "Fourth"


# --- GLM comparisons --------------------------------------------------
cat("\n======== COLD — FIRST INSTAR ========\n")
m1 <- compare_full_reduced(c1, predictor = "time", reduce_further = TRUE)
print(m1$lrt)
cat("\n(additive vs predictor-only):\n")
print(m1$lrt_further)
print(summary(m1$reduced))

cat("\n======== COLD — FOURTH INSTAR ========\n")
m4 <- compare_full_reduced(c4, predictor = "time")
print(m4$lrt)
print(summary(m4$full))


# --- Bootstrap survivorship curves ------------------------------------
boot_1 <- bootstrap_by_pop(c1, predictor = "time")
boot_4 <- bootstrap_by_pop(c4, predictor = "time")
boot_1$instar <- "First"
boot_4$instar <- "Fourth"
boot_c <- dplyr::bind_rows(boot_1, boot_4)
boot_c$instar <- factor(boot_c$instar, levels = c("First", "Fourth"))


# --- Final-survival (50-min at 0 °C) summary -------------------------
fs_df <- final_survival(boot_c, predictor = "time",
                        final_value = COLD_TIME_MAX)

cat("\n======== COLD — 50-min survival at 0 °C ========\n")
print(fs_df)

readr::write_csv(fs_df, file.path(TAB_DIR, "cold_final_survival.csv"))


# --- Figure 4A/B: survivorship curves ---------------------------------
p_curve <- survivorship_curve(
  boot_c,
  predictor = "time",
  x_lab     = "Time at 0\u00B0C (minutes)",
  x_breaks  = seq(0, 50, by = 10)
)
save_fig(p_curve, "Fig4AB_cold_curves", width = 10, height = 5)


# --- Figure 4C/D: 50-min final survival ------------------------------
p_final <- summary_point_plot(
  fs_df,
  y     = "median",
  y_lab = "50-minute Survival\n\u00B195% CI"
) + ggplot2::scale_y_continuous(limits = c(0, 1))
save_fig(p_final, "Fig4CD_cold_final", width = 8.82, height = 4)


# --- Export -----------------------------------------------------------
saveRDS(
  list(final = fs_df, boot = boot_c,
       model_first = m1, model_fourth = m4),
  file.path(CACHE_DIR, "cold_results.rds")
)
