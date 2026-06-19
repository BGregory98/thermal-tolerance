# analysis/03_sustained_hot.R -------------------------------------------
#
# Figure 3: Sustained high-temperature trial at 38 °C (first + fourth
# instar). Primary summary metric is the proportion alive at 96 min.
#
# Analyzes the five main populations (POPS_MAIN from config.R). Uses
# population-level binomial GLMs; no region / species pooling, per
# decisions from the 2026-04-24 advisor meeting.
# ------------------------------------------------------------------------

# --- Load both instars -------------------------------------------------
sh1 <- load_survival_data(
  path     = file.path(DATA_DIR, "chronic_first.csv"),
  time_var = "time",
  pops     = POPS_MAIN
)
sh4 <- load_survival_data(
  path     = file.path(DATA_DIR, "chronic_fourth.csv"),
  time_var = "time",
  pops     = POPS_MAIN
)

sh1$instar <- "First"
sh4$instar <- "Fourth"


# --- GLM comparisons --------------------------------------------------
# For first instar the original analysis reduced twice: full vs additive,
# and additive vs predictor-only.
cat("\n======== SUSTAINED HOT — FIRST INSTAR ========\n")
m1 <- compare_full_reduced(sh1, predictor = "time", reduce_further = TRUE)
print(m1$lrt)
cat("\n(additive vs predictor-only):\n")
print(m1$lrt_further)
print(summary(m1$reduced))

cat("\n======== SUSTAINED HOT — FOURTH INSTAR ========\n")
m4 <- compare_full_reduced(sh4, predictor = "time")
print(m4$lrt)
print(summary(m4$full))


# --- Bootstrap survivorship curves ------------------------------------
boot_1 <- bootstrap_by_pop(sh1, predictor = "time")
boot_4 <- bootstrap_by_pop(sh4, predictor = "time")
boot_1$instar <- "First"
boot_4$instar <- "Fourth"
boot_sh <- dplyr::bind_rows(boot_1, boot_4)
boot_sh$instar <- factor(boot_sh$instar, levels = c("First", "Fourth"))


# --- Final-survival (96-min) summary ---------------------------------
fs_df <- final_survival(boot_sh, predictor = "time",
                        final_value = SUSTAINED_TIME_MAX)

cat("\n======== SUSTAINED HOT — 96-min survival ========\n")
print(fs_df)

readr::write_csv(fs_df, file.path(TAB_DIR, "sustained_hot_final_survival.csv"))


# --- Figure 3A/B: survivorship curves --------------------------------
p_curve <- survivorship_curve(
  boot_sh,
  predictor = "time",
  x_lab     = "Time at 38\u00B0C (minutes)",
  x_breaks  = c(0, 24, 48, 72, 96)
)
save_fig(p_curve, "Fig3AB_sustained_hot_curves", width = 10, height = 5)


# --- Figure 3C/D: 96-min final survival ------------------------------
p_final <- summary_point_plot(
  fs_df,
  y     = "median",
  y_lab = "96-minute Survival\n\u00B195% CI"
) + ggplot2::scale_y_continuous(limits = c(0, 1))
save_fig(p_final, "Fig3CD_sustained_hot_final", width = 8.82, height = 4)


# --- Export -----------------------------------------------------------
saveRDS(
  list(final = fs_df, boot = boot_sh,
       model_first = m1, model_fourth = m4),
  file.path(CACHE_DIR, "sustained_hot_results.rds")
)
