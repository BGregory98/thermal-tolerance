# analysis/08_hybrid_zone.R ---------------------------------------------
#
# SECONDARY ANALYSIS — mid-latitude hybrid-zone extension.
#
# Framing (from 2026-04-24 advisor meeting): the primary first-instar
# ramping experiment revealed the clearest population-level differences
# in LT50, with northern (Cx. pipiens) populations less heat-tolerant
# than southern (Cx. quinquefasciatus) populations. As a secondary
# question, we ask whether mid-latitude populations — from the
# traditional Barr (1957) hybrid zone — show LT50s resembling Cx.
# pipiens, Cx. quinquefasciatus, or intermediate values.
#
# Populations analyzed here (POPS_HZ from config.R):
#   FL, TX        — southern Cx. quinquefasciatus
#   ME, IL-AG1    — northern Cx. pipiens
#   Baltimore     — mid-latitude
#   SERC          — mid-latitude
#
# IL-BG2 is excluded because it is lab-adapted (>150 generations) and
# unlikely to reflect field-relevant thermal biology at this scale.
#
# Data sources:
#   tramp_cx.csv    — first instar ramping for the four main populations
#   tramp_cx_HZ.csv — first instar ramping for Baltimore and SERC
# ------------------------------------------------------------------------

# --- Load and combine the two first-instar ramping datasets -----------
rh_main <- load_survival_data(
  path     = file.path(DATA_DIR, "tramp_cx.csv"),
  time_var = "temp",
  pops     = POPS_HZ   # automatically excludes IL-BG2
)

rh_new <- load_survival_data(
  path     = file.path(DATA_DIR, "tramp_cx_HZ.csv"),
  time_var = "temp",
  pops     = POPS_HZ
)

# Align columns (tramp_cx_HZ should match tramp_cx; rbind-compatible)
common_cols <- intersect(names(rh_main), names(rh_new))
rh_all <- dplyr::bind_rows(
  rh_main[, common_cols],
  rh_new[,  common_cols]
)

# Re-apply the south -> mid -> north ordering using POP_LEVELS_HZ
rh_all$population <- factor(as.character(rh_all$population),
                            levels = POP_LEVELS_HZ)

cat("\n======== HYBRID-ZONE EXTENSION — sample sizes ========\n")
print(table(rh_all$population))


# --- Per-population LT50 ----------------------------------------------
lt50_hz <- lt50_by_pop(rh_all, predictor = "temp")

cat("\n======== HYBRID-ZONE — LT50 per population (°C) ========\n")
print(lt50_hz)

readr::write_csv(lt50_hz, file.path(TAB_DIR, "hybrid_zone_LT50.csv"))


# --- Bootstrap survivorship curves (all 6 populations) ----------------
boot_hz <- bootstrap_by_pop(rh_all, predictor = "temp")

# --- Figure: first-instar ramping survivorship curves ----------------
# Single-panel version (no instar faceting, since only first-instar data
# is available for Baltimore and SERC).
p_curve <- ggplot2::ggplot(boot_hz) +
  ggplot2::geom_ribbon(
    ggplot2::aes(x = temp, ymin = lower, ymax = upper,
                 fill = population),
    alpha = 0.4, linetype = "dotted", color = "black"
  ) +
  ggplot2::geom_line(
    ggplot2::aes(x = temp, y = median, color = population)
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = temp, y = median, color = population)
  ) +
  ggplot2::scale_color_manual(values = pop_colors, name = "Population") +
  ggplot2::scale_fill_manual( values = pop_colors, name = "Population") +
  ggplot2::scale_x_continuous(breaks = seq(36, 46, by = 2)) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(x = "Temperature (\u00B0C)",
                y = "Proportion Alive \u00B195% CI") +
  theme_tt(base_size = 18)

save_fig(p_curve, "FigS_hybrid_zone_ramping_curves",
         width = 10, height = 5)


# --- Figure: LT50 per population ---------------------------------------
p_lt50 <- ggplot2::ggplot(lt50_hz) +
  ggplot2::geom_errorbar(
    ggplot2::aes(x = population, ymin = lower, ymax = upper),
    width = 0.4
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = population, y = lt50, fill = population),
    shape = 21, size = 5, stroke = 0.3
  ) +
  ggplot2::scale_fill_manual(values = pop_colors, guide = "none") +
  ggplot2::labs(x = "Population",
                y = "Median Lethal Temperature (\u00B0C)\n\u00B195% CI") +
  theme_tt(base_size = 14)

save_fig(p_lt50, "FigS_hybrid_zone_LT50",
         width = 8, height = 4)


# --- A priori BIOCLIM regression --------------------------------------
# Per the advisor meeting: use three biologically-motivated BIOCLIM
# variables (bio04 seasonality, bio05 warmest-month max T, bio06
# coldest-month min T) and regress LT50 separately against each.
#
# With n = 6 populations and 3 predictors fit univariately, this is
# exploratory; we report adjusted R² and the direction of effect rather
# than treating the p-values as confirmatory. The framing in the
# Results should make this clear.

# Load BIOCLIM values at the HZ site coordinates.
bioclim_stack <- load_bioclim_us(subset = A_PRIORI_BIOCLIM)

coord_path <- file.path(DATA_DIR, "population_coordinates.csv")
coords_raw <- readr::read_csv(coord_path, show_col_types = FALSE)
# Coord file mixes legacy names ("Maine") with pop_ids ("IL-AG1").
# Translate legacy names; pass pop_ids through unchanged.
lookup <- setNames(pop_info$pop_id, pop_info$csv_name)
coords_raw$pop_id <- ifelse(
  coords_raw$population %in% pop_info$pop_id,
  coords_raw$population,
  unname(lookup[coords_raw$population])
)
coords_hz <- coords_raw[coords_raw$pop_id %in% POPS_HZ, ]

pts <- terra::vect(
  as.matrix(coords_hz[, c("longitude", "latitude")]),
  type = "points",
  crs  = "+proj=longlat +datum=WGS84 +no_defs"
)
pts <- terra::project(pts, terra::crs(bioclim_stack))
env_hz <- terra::extract(bioclim_stack, pts)
env_hz$ID <- NULL
env_hz$pop_id <- coords_hz$pop_id

# Merge LT50 (with CIs) and environmental values
lt50_for_reg <- lt50_hz[, c("population", "lt50", "lower", "upper")]
names(lt50_for_reg)[1] <- "pop_id"
lt50_for_reg$pop_id <- as.character(lt50_for_reg$pop_id)

reg_df <- merge(lt50_for_reg, env_hz, by = "pop_id")


cat("\n======== HYBRID-ZONE — LT50 vs BIOCLIM ========\n")
print(reg_df)

reg_df$pop_id = factor(reg_df$pop_id, levels = c('FL', 'TX', 'SERC', 'Baltimore', 'IL-AG3', 'ME'))

# Fit univariate LMs for each of the three a priori variables.
# Report the Pearson correlation coefficient (r) alongside the slope and
# R^2, and correct the three p-values for multiple testing across the
# 3 a priori variables (Bonferroni is the literal "corrected for 3 tests";
# Benjamini-Hochberg/FDR also reported).
reg_results <- do.call(rbind, lapply(A_PRIORI_BIOCLIM, function(v) {
  if (!v %in% names(reg_df)) {
    return(data.frame(predictor = v, n = 0, beta = NA, r = NA,
                      r2 = NA, adj_r2 = NA, p_value = NA))
  }
  m <- lm(lt50 ~ reg_df[[v]], data = reg_df)
  s <- summary(m)
  data.frame(
    predictor = v,
    n         = nrow(reg_df),
    beta      = coef(m)[2],
    r         = cor(reg_df$lt50, reg_df[[v]],
                    use = "complete.obs", method = "pearson"),
    r2        = s$r.squared,
    adj_r2    = s$adj.r.squared,
    p_value   = coef(s)[2, 4]
  )
}))

# Multiple-comparison correction across the 3 a priori variables
reg_results$p_bonferroni <- p.adjust(reg_results$p_value, method = "bonferroni")
reg_results$p_BH         <- p.adjust(reg_results$p_value, method = "BH")

reg_results$predictor_label <- c(
  bio04 = "Temperature Seasonality",
  bio05 = "Max T Warmest Month",
  bio06 = "Min T Coldest Month"
)[reg_results$predictor]

reg_results <- reg_results[, c("predictor", "predictor_label", "n",
                               "beta", "r", "r2", "adj_r2",
                               "p_value", "p_bonferroni", "p_BH")]
# Present in a priori order (bio04, bio05, bio06), not ranked by fit, so
# the table does not foreground the strongest predictor.
reg_results <- reg_results[match(A_PRIORI_BIOCLIM, reg_results$predictor), ]

cat("\n======== UNIVARIATE LMs (LT50 ~ a priori BIOCLIM), corrected over 3 tests ========\n")
print(reg_results, row.names = FALSE)

readr::write_csv(reg_results,
                 file.path(TAB_DIR, "hybrid_zone_bioclim_LT50.csv"))


# --- Figure: LT50 vs bio05 (per-population scatter w/ regression line)
# This is the figure intended for the supplementary materials. The
# regression line is clipped to the range of observed bio05 values.
if ("bio05" %in% names(reg_df)) {
  p_bio05 <- lt50_vs_bioclim(
    reg_df,
    predictor = "bio05",
    x_lab     = "Max Temperature of Warmest Month, bio05 (\u00B0C)"
  )
  save_fig(p_bio05, "FigS_hybrid_zone_LT50_vs_bio05",
           width = 8, height = 4.5)
}


# --- Export -----------------------------------------------------------
saveRDS(
  list(
    lt50    = lt50_hz,
    boot    = boot_hz,
    env     = env_hz,
    reg_df  = reg_df,
    results = reg_results
  ),
  file.path(CACHE_DIR, "hybrid_zone_results.rds")
)