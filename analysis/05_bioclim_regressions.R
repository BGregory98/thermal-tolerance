# analysis/05_bioclim_regressions.R -------------------------------------
#
# BIOCLIM ~ phenotype regressions. Fills the gap left by the empty
# __Linear_regressions.R in the original codebase.
#
# The manuscript methods describe this as: fit a linear model of each
# thermal-tolerance summary statistic (LT50 for ramping hot;
# final-timepoint survival for sustained hot and cold) on each of the
# 19 BIOCLIM variables, separately, and report the best-correlated
# predictor by adjusted R-squared.
#
# This script:
#   1. Loads (or reconstructs) the cached BIOCLIM stack
#   2. Extracts bioclim values at each population's coordinates
#   3. Joins them to the phenotype summary tables produced by
#      analysis/02, 03, and 04
#   4. Fits all univariate LMs (19 variables × N phenotype summaries)
#   5. Writes ranked tables and saves the per-phenotype winner
# ------------------------------------------------------------------------


# --- Load BIOCLIM and extract at collection points ---------------------
bioclim_stack <- load_bioclim_us()
env <- extract_bioclim_at_pops(bioclim_stack, pops = pop_info_hz)

cat("\n======== BIOCLIM VALUES AT COLLECTION SITES ========\n")
print(env)

readr::write_csv(env, file.path(TAB_DIR, "bioclim_at_pops.csv"))


# --- Build the phenotype summary table --------------------------------
# Pull cached results from the three experiment scripts.
rh <- readRDS(file.path(CACHE_DIR, "ramping_hot_results.rds"))
sh <- readRDS(file.path(CACHE_DIR, "sustained_hot_results.rds"))
co <- readRDS(file.path(CACHE_DIR, "cold_results.rds"))

# Each experiment contributes two phenotypes (first instar + fourth
# instar). Ramping hot uses LT50; sustained hot and cold use the
# bootstrap median at the final timepoint.

phenotypes <- dplyr::bind_rows(
  rh$lt50 %>%
    dplyr::transmute(
      pop_id    = as.character(population),
      instar,
      experiment = "ramping_hot",
      metric     = "LT50",
      value      = lt50
    ),
  sh$final %>%
    dplyr::transmute(
      pop_id    = as.character(population),
      instar    = as.character(instar),
      experiment = "sustained_hot",
      metric     = "final_survival_96min",
      value      = median
    ),
  co$final %>%
    dplyr::transmute(
      pop_id    = as.character(population),
      instar    = as.character(instar),
      experiment = "cold",
      metric     = "final_survival_50min_0C",
      value      = median
    )
)


# --- Merge with environmental data ------------------------------------
bio_cols <- intersect(names(env), names(BIOCLIM_NAMES))

dat <- phenotypes %>%
  dplyr::inner_join(
    env %>% dplyr::select(pop_id, dplyr::all_of(bio_cols)),
    by = "pop_id"
  )


# --- Fit univariate LMs for every phenotype × predictor ---------------
# (Skips phenotypes with fewer than 3 finite observations, which can
# happen if populations are dropped and a cell is empty.)

fit_one <- function(pheno_df, predictor_col) {
  y <- pheno_df$value
  x <- pheno_df[[predictor_col]]
  if (sum(is.finite(x) & is.finite(y)) < 3) return(NULL)

  m <- stats::lm(y ~ x)
  s <- summary(m)
  tibble::tibble(
    predictor   = predictor_col,
    predictor_label = unname(BIOCLIM_NAMES[predictor_col]),
    n           = sum(is.finite(x) & is.finite(y)),
    beta        = stats::coef(m)[2],
    r2          = s$r.squared,
    adj_r2      = s$adj.r.squared,
    p_value     = stats::coef(s)[2, 4]
  )
}

results <- dat %>%
  dplyr::group_by(experiment, metric, instar) %>%
  dplyr::group_modify(function(sub, key) {
    purrr::map_dfr(bio_cols, function(p) fit_one(sub, p))
  }) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(experiment, instar, dplyr::desc(adj_r2))


cat("\n======== TOP BIOCLIM PREDICTORS PER PHENOTYPE ========\n")
top_predictors <- results %>%
  dplyr::group_by(experiment, metric, instar) %>%
  dplyr::slice_max(adj_r2, n = 3, with_ties = FALSE) %>%
  dplyr::ungroup()
print(top_predictors, n = Inf)

readr::write_csv(results,
                 file.path(TAB_DIR, "bioclim_regressions_all.csv"))
readr::write_csv(top_predictors,
                 file.path(TAB_DIR, "bioclim_regressions_top3.csv"))


# --- Optional: also fit species/region-grouped versions ---------------
# Reviewers want the reframed analysis. If you want the same regressions
# but at the species or region level (pooling populations), set
# `RUN_GROUPED <- TRUE`. By default we only do the per-population
# regressions that match the submitted manuscript.

RUN_GROUPED <- FALSE

if (RUN_GROUPED) {
  message("Grouped regressions not yet implemented. With 2 species and ",
          "2 regions, n is too small for linear regression; use t-tests ",
          "or effect-size comparisons instead.")
}
