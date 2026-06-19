# analysis/02_ramping_hot.R ---------------------------------------------
#
# Figure 2: High-temperature ramping trial (first + fourth instar).
#
# Analyzes the five main populations (POPS_MAIN from config.R). Uses
# population-level binomial GLMs; no region / species pooling, per
# decisions from the 2026-04-24 advisor meeting.
#
# A separate hybrid-zone extension (analysis/08_hybrid_zone.R) adds
# Baltimore and SERC to the first-instar ramping data for a
# supplementary BIOCLIM regression.
# ------------------------------------------------------------------------

# --- Load both instars -------------------------------------------------
rh1 <- load_survival_data(
  path     = file.path(DATA_DIR, "tramp_cx.csv"),
  time_var = "temp",
  pops     = POPS_MAIN
)
rh4 <- load_survival_data(
  path     = file.path(DATA_DIR, "tramp_cx_4.csv"),
  time_var = "temp",
  pops     = POPS_MAIN
)

rh1$instar <- "First"
rh4$instar <- "Fourth"


# --- GLM comparisons (full vs reduced) --------------------------------
cat("\n======== RAMPING HOT — FIRST INSTAR ========\n")
m1 <- compare_full_reduced(rh1, predictor = "temp")
print(summary(m1$full))
print(m1$lrt)

cat("\n======== RAMPING HOT — FOURTH INSTAR ========\n")
m4 <- compare_full_reduced(rh4, predictor = "temp")
print(summary(m4$full))
print(m4$lrt)


# --- LT50 per population ----------------------------------------------
lt50_1 <- lt50_by_pop(rh1, predictor = "temp")
lt50_4 <- lt50_by_pop(rh4, predictor = "temp")
lt50_1$instar <- "First"
lt50_4$instar <- "Fourth"

lt50_df <- dplyr::bind_rows(lt50_1, lt50_4)
lt50_df$instar <- factor(lt50_df$instar, levels = c("First", "Fourth"))

cat("\n======== RAMPING HOT — LT50 (°C) ========\n")
print(lt50_df)

readr::write_csv(lt50_df, file.path(TAB_DIR, "ramping_hot_LT50.csv"))


# --- Bootstrap survivorship curves ------------------------------------
boot_1 <- bootstrap_by_pop(rh1, predictor = "temp")
boot_4 <- bootstrap_by_pop(rh4, predictor = "temp")
boot_1$instar <- "First"
boot_4$instar <- "Fourth"
boot_rh <- dplyr::bind_rows(boot_1, boot_4)
boot_rh$instar <- factor(boot_rh$instar, levels = c("First", "Fourth"))


# --- Figure 2A/B: survivorship curves ---------------------------------
p_curve <- survivorship_curve(
  boot_rh,
  predictor = "temp",
  x_lab     = "Temperature (\u00B0C)",
  x_breaks  = seq(36, 46, by = 2)
)
save_fig(p_curve, "Fig2AB_ramping_hot_curves", width = 10, height = 5)


# --- Figure 2C/D: LT50 summary ----------------------------------------
p_lt50 <- summary_point_plot(
  lt50_df,
  y     = "lt50",
  y_lab = "Median Lethal Temperature (\u00B0C)\n\u00B195% CI"
)
save_fig(p_lt50, "Fig2CD_ramping_hot_LT50", width = 8.82, height = 4)


# --- Export -----------------------------------------------------------
saveRDS(
  list(lt50 = lt50_df, boot = boot_rh,
       model_first = m1, model_fourth = m4),
  file.path(CACHE_DIR, "ramping_hot_results.rds")
)
