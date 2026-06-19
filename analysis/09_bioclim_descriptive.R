# analysis/10_bioclim_lt50_hz_correction.R ------------------------------
#
# Full BIOCLIM screen of first-instar ramping LT50 for the six-population
# mid-latitude ("hybrid-zone") extension, with multiple-comparison
# correction. Per advisor request (Fritz, 2026-05):
#
#   * Regress first-instar ramping LT50 on EACH of the 19 primary BIOCLIM
#     variables AND latitude  (= 20 univariate linear models).
#   * Populations: FL, TX, ME, IL-AG1, Baltimore, SERC  (IL-BG2 excluded).
#   * Correct the 20 p-values for multiple testing (Bonferroni for the
#     family-wise bound; Benjamini-Hochberg for the FDR, which is the one
#     we recommend reporting for a screen of this kind).
#   * Report slope, R^2, adj R^2, raw p, and both corrected p-values.
#
# This is the "return to all 19 + latitude" analysis. The a priori
# three-variable version lives in 09_bioclim_descriptive.R; the broader
# per-phenotype loop lives in 05_bioclim_regressions.R. This script is
# deliberately narrow: ONE phenotype (first-instar LT50), the SIX
# extension populations, all 20 predictors, corrected.
#
# Assumes the project config/helpers are already loaded (same as 05 & 09):
#   load_bioclim_us(), extract_bioclim_at_pops(), BIOCLIM_NAMES,
#   pop_info, pop_info_hz, DATA_DIR, TAB_DIR, CACHE_DIR
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(readr)
})

# --- The six populations in this extension (IL-BG2 excluded) ------------
HZ_POPS <- c("FL", "TX", "Baltimore", "SERC", "IL-AG3", "ME")

# Which phenotype + instar to screen. Kept as variables so the same
# script can be re-pointed at another phenotype later if desired.
TARGET_INSTAR <- "First"   # caches code instar as "First"/"Fourth"; matched case-insensitively

# --- 1. BIOCLIM (all 19) at the six sites ------------------------------
bioclim_stack <- load_bioclim_us()
env <- extract_bioclim_at_pops(bioclim_stack, pops = pop_info_hz)
env$pop_id <- as.character(env$pop_id)

# Restrict to the six target populations, in case pop_info_hz is broader.
env <- env[env$pop_id %in% HZ_POPS, , drop = FALSE]
stopifnot("Did not recover all six HZ populations from env" =
            setequal(env$pop_id, HZ_POPS))

# --- 2. Add latitude as the 20th predictor -----------------------------
# Use a latitude column if extract_bioclim_at_pops already returned one;
# otherwise pull from pop_info, falling back to population_coordinates.csv.
if (!"latitude" %in% names(env)) {
  if ("lat" %in% names(env)) {
    env$latitude <- env$lat
  } else if ("lat" %in% names(pop_info)) {
    lat_lookup <- setNames(pop_info$lat, as.character(pop_info$pop_id))
    env$latitude <- unname(lat_lookup[env$pop_id])
  } else {
    coord_path <- file.path(DATA_DIR, "population_coordinates.csv")
    cc <- readr::read_csv(coord_path, show_col_types = FALSE)
    lookup <- setNames(pop_info$pop_id, pop_info$csv_name)
    cc$pop_id <- ifelse(cc$population %in% pop_info$pop_id,
                        cc$population, unname(lookup[cc$population]))
    lat_lookup <- setNames(cc$latitude, cc$pop_id)
    env$latitude <- unname(lat_lookup[env$pop_id])
  }
}
stopifnot("Missing/!finite latitude for one or more HZ populations" =
            all(is.finite(env$latitude)))

# --- 3. First-instar LT50 for the SIX populations ----------------------
# IMPORTANT: ramping_hot_results.rds (`rh`) holds only the five MAIN
# populations and excludes the hybrid-zone sites (Baltimore, SERC). The
# six-population first-instar LT50 used in the mid-latitude extension is
# produced by the hybrid-zone script. We therefore source LT50 from the
# hybrid-zone cache, and fall back to merging the main-pop first instars
# (rh) with the HZ-only pops (Baltimore, SERC) if that cache holds only
# the two HZ sites. Each population's LT50 is fit by per-population GLM,
# so a given population's value is identical whether fit in the 5-pop or
# 6-pop run, making the merge exact.
HZ_LT50_RDS  <- "hybrid_zone_results.rds"  # <-- six-pop extension cache; adjust name if needed
HZ_LT50_ELEM <- "lt50"                     # <-- list element holding the LT50 tibble

first_lt50 <- function(tbl) {
  tbl %>%
    dplyr::mutate(pop_id = as.character(population)) %>%
    dplyr::transmute(pop_id, LT50 = lt50)
}

lt50    <- NULL
hz      <- NULL
hz_path <- file.path(CACHE_DIR, HZ_LT50_RDS)
if (file.exists(hz_path)) {
  hz <- readRDS(hz_path)
  if (!is.null(hz[[HZ_LT50_ELEM]])) {
    cand <- first_lt50(hz[[HZ_LT50_ELEM]]) %>% dplyr::filter(pop_id %in% HZ_POPS)
    if (setequal(cand$pop_id, HZ_POPS)) lt50 <- cand   # cache already holds all six
  }
}

# Fallback: combine the four main pops in the six (from rh) with the two
# HZ-only pops (from the hybrid-zone cache).
if (is.null(lt50)) {
  rh <- readRDS(file.path(CACHE_DIR, "ramping_hot_results.rds"))
  main_first <- first_lt50(rh$lt50) %>% dplyr::filter(pop_id %in% HZ_POPS)
  need       <- setdiff(HZ_POPS, main_first$pop_id)     # expect Baltimore, SERC
  hz_first   <- if (!is.null(hz) && !is.null(hz[[HZ_LT50_ELEM]]))
                  first_lt50(hz[[HZ_LT50_ELEM]]) %>% dplyr::filter(pop_id %in% need)
                else tibble::tibble(pop_id = character(), LT50 = numeric())
  lt50 <- dplyr::bind_rows(main_first, hz_first) %>%
            dplyr::distinct(pop_id, .keep_all = TRUE)
}

# Loud guard: tells you exactly which populations are missing if it fails.
if (!setequal(lt50$pop_id, HZ_POPS)) {
  stop("Could not assemble first-instar LT50 for all six HZ populations.\n",
       "  Have:    ", paste(sort(lt50$pop_id), collapse = ", "), "\n",
       "  Missing: ", paste(sort(setdiff(HZ_POPS, lt50$pop_id)), collapse = ", "), "\n",
       "  Point HZ_LT50_RDS / HZ_LT50_ELEM at the cache holding Baltimore & SERC LT50.")
}

# --- 4. Assemble the modelling table -----------------------------------
bio_cols   <- intersect(names(env), names(BIOCLIM_NAMES))   # the 19
predictors <- c(bio_cols, "latitude")                       # 20 total
stopifnot("Expected 19 BIOCLIM predictors" = length(bio_cols) == 19L)

dat <- dplyr::inner_join(lt50, env[, c("pop_id", predictors)], by = "pop_id")
stopifnot(nrow(dat) == length(HZ_POPS))
# analysis/11_collinearity_filter.R ------------------------------------
#
# Collinearity prune for the bioclim LT50 screen. Drop-in: source or paste
# this AFTER `dat` (pop_id, LT50, + the 20 predictor columns) and
# `predictors` (length-20 character vector of the 19 BIOCLIM names +
# "latitude") are built in 10_bioclim_lt50_hz_correction.R, and INSTEAD OF
# that script's all-20 regression step.
#
# Steps:
#   1. Pearson correlation matrix among the 20 predictors (six sites).
#   2. Report every pair with |r| >= cutoff.
#   3. Retain a non-redundant subset in which every pairwise |r| < cutoff.
#   4. Regress LT50 on each RETAINED predictor and correct the p-values
#      over the retained count (Bonferroni + Benjamini-Hochberg).
#
# The within-pair survivor is chosen by an a priori PRIORITY order, NOT by
# correlation with LT50. Keeping the filter outcome-independent is what
# makes it legitimate to test the survivors and correct for their number.
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(readr)
})

CORR_CUTOFF <- 0.7

stopifnot("`dat` not found"        = exists("dat"),
          "`predictors` not found" = exists("predictors"),
          all(predictors %in% names(dat)))

# --- 1. Correlation matrix among the 20 predictors ---------------------
X  <- as.matrix(dat[, predictors])
cm <- stats::cor(X, method = "pearson")     # 20 x 20, named

readr::write_csv(
  tibble::as_tibble(round(cm, 3), rownames = "predictor"),
  file.path(TAB_DIR, "bioclim_predictor_correlation_matrix.csv"))

# --- 2. Report collinear pairs (|r| >= cutoff) -------------------------
ij <- which(upper.tri(cm) & abs(cm) >= CORR_CUTOFF, arr.ind = TRUE)
pairs_tbl <- if (nrow(ij) == 0) tibble::tibble() else
  tibble::tibble(var1 = rownames(cm)[ij[, "row"]],
                 var2 = colnames(cm)[ij[, "col"]],
                 r    = round(cm[ij], 3)) %>%
  dplyr::arrange(dplyr::desc(abs(r)))

cat("\n======== PREDICTOR PAIRS WITH |r| >=", CORR_CUTOFF, "========\n")
print(pairs_tbl, n = Inf)

# --- 3. Retain a non-redundant subset ----------------------------------
# Priority = the three a priori thermal vars (assay-relevant), then
# latitude, then the remaining BIOCLIM in their natural order. EDIT this
# vector to change which variable survives a collinear cluster (e.g. put
# a precipitation variable first to guarantee it is retained). Do NOT
# order it by correlation with LT50.
priority <- c("latitude", "bio05", "bio06", "bio04",
              setdiff(predictors, c("latitude", "bio05", "bio06", "bio04")))
priority <- priority[priority %in% predictors]

keep <- character(0)
for (v in priority) {
  if (length(keep) == 0L || all(abs(cm[v, keep]) < CORR_CUTOFF))
    keep <- c(keep, v)
}
predictors_kept <- keep
dropped         <- setdiff(predictors, predictors_kept)

cat(sprintf("\nRetained %d of %d predictors (all pairwise |r| < %.2f):\n  %s\n",
            length(predictors_kept), length(predictors), CORR_CUTOFF,
            paste(predictors_kept, collapse = ", ")))
cat("Dropped as collinear with a higher-priority retained variable:\n  ",
    paste(dropped, collapse = ", "), "\n", sep = "")

# --- 4. Regress LT50 on the retained predictors, then correct ----------
# (fit_one is defined in 10_*.R; redefined here so this section is
# self-contained if sourced separately.)
fit_one <- function(df, predictor) {
  y <- df$LT50; x <- df[[predictor]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NULL)
  s  <- summary(stats::lm(y ~ x, data = data.frame(y = y[ok], x = x[ok])))
  ct <- stats::coef(s)
  tibble::tibble(
    predictor = predictor,
    label     = if (predictor %in% names(BIOCLIM_NAMES))
                  unname(BIOCLIM_NAMES[predictor]) else "Latitude",
    n = sum(ok), slope = ct[2, 1],
    r2 = s$r.squared, adj_r2 = s$adj.r.squared, p_raw = ct[2, 4])
}

res_kept <- purrr::map_dfr(predictors_kept, ~ fit_one(dat, .x)) %>%
  dplyr::mutate(
    p_bonferroni = stats::p.adjust(p_raw, method = "bonferroni"),
    p_BH         = stats::p.adjust(p_raw, method = "BH")) %>%
  dplyr::arrange(p_raw)

cat("\n======== LT50 ~ retained (non-collinear) predictors ========\n")
print(res_kept, n = Inf)

readr::write_csv(res_kept,
                 file.path(TAB_DIR, "bioclim_lt50_hz_screen_filtered.csv"))
saveRDS(res_kept,
        file.path(CACHE_DIR, "bioclim_lt50_hz_screen_filtered.rds"))

# --- 5. Copy-paste summary ---------------------------------------------
best   <- res_kept[which.min(res_kept$p_raw), ]
n_raw  <- sum(res_kept$p_raw        < 0.05, na.rm = TRUE)
n_bh   <- sum(res_kept$p_BH         < 0.05, na.rm = TRUE)
n_bonf <- sum(res_kept$p_bonferroni < 0.05, na.rm = TRUE)

cat(sprintf(
"\n--- Prose summary ---
After removing collinear predictors (|r| >= %.1f), %d of the 20 candidate
variables were retained (%s). Across univariate regressions of first-instar
LT50 (n = %d populations) on these retained variables, the strongest
predictor was %s (%s; R2 = %.2f, raw p = %.3f). %d retained variable(s) were
nominally significant at raw p < 0.05; after correcting for the %d retained
tests, %d remained significant under Benjamini-Hochberg (min adj p = %.3f)
and %d under Bonferroni (min adj p = %.3f).\n",
  CORR_CUTOFF, length(predictors_kept),
  paste(predictors_kept, collapse = ", "), best$n,
  best$predictor, best$label, best$r2, best$p_raw,
  n_raw, length(predictors_kept),
  n_bh,   min(res_kept$p_BH, na.rm = TRUE),
  n_bonf, min(res_kept$p_bonferroni, na.rm = TRUE)))

cat("\nNOTE: with n = 6 the inter-predictor correlations are themselves",
    "\nestimated from six points and are noisy, so the 0.7 cutoff is a rough",
    "\nscreen, not a precise threshold. Report which variable was kept from",
    "\neach collinear cluster as an explicit a priori choice.\n")