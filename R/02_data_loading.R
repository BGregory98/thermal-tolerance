# R/02_data_loading.R ---------------------------------------------------
#
# Single entry point for loading thermal-tolerance trial CSVs.
#
# Raw CSVs use legacy `population` labels ("Florida", "Chicago-AG", etc.).
# This loader renames them to the short pop_id codes ("FL", "IL-AG1"),
# attaches species and region metadata, and optionally filters to the
# populations listed in POPS_MAIN.
# ------------------------------------------------------------------------


#' Load a survival-trial CSV and clean it
#'
#' @param path Path to the CSV (cold_first.csv, chronic_fourth.csv, etc.).
#' @param time_var Name of the predictor variable in the CSV. For the
#'   ramping-hot experiment this is `"temp"`; for the sustained-hot and
#'   cold experiments it's `"time"`.
#' @param exclude_pupae Whether to drop rows representing larvae that
#'   pupated during the trial. Older CSVs marked these with `alive == 2`;
#'   newer CSVs use a separate `pupated` column with values "T"/"F"/NA.
#'   Both encodings are handled when this is TRUE (the default).
#' @param exclude_controls Whether to drop rows where `tcycler == 0`
#'   (control cups). Default TRUE. Ignored if the column is absent
#'   (the cold-experiment CSVs don't have a tcycler column).
#' @param drop_edge_time Numeric or NULL. For the sustained-hot data,
#'   rows with `time >= drop_edge_time` were excluded in the original
#'   analysis "due to possible edge effect". Default 120 for consistency
#'   with the submitted manuscript; pass NULL to keep everything.
#' @param pops Character vector of pop_ids to KEEP. Defaults to POPS_MAIN.
#'
#' @return A tibble with columns: population (factor, with levels per
#'   POP_LEVELS), species, region, trial, alive (numeric 0/1), and
#'   the time_var column, plus any other columns that were in the CSV.
load_survival_data <- function(path,
                               time_var          = c("time", "temp"),
                               exclude_pupae     = TRUE,
                               exclude_controls  = TRUE,
                               drop_edge_time    = 120,
                               pops              = POPS_MAIN) {

  time_var <- match.arg(time_var)

  if (!file.exists(path)) {
    stop("Data file not found: ", path,
         "\nCheck DATA_DIR in config.R.")
  }

  d <- utils::read.csv(path, stringsAsFactors = FALSE)

  # --- Repair column names ----------------------------------------------
  # Some CSVs have `tiMaine` where `time` should be — a legacy find/replace
  # bug ("me" -> "Maine" applied globally, mangling "tiMe"). Auto-rename
  # so downstream code that expects `time` just works.
  if ("tiMaine" %in% names(d) && !"time" %in% names(d)) {
    message("  Repairing column name: 'tiMaine' -> 'time' (", basename(path), ")")
    names(d)[names(d) == "tiMaine"] <- "time"
  }

  # NA handling: use complete.cases only on the columns we actually rely on.
  # (The whole-row version dropped rows where an irrelevant optional column
  # like `pupated` was NA, which was too aggressive.)
  required_cols <- intersect(c("population", "alive", "time", "temp"), names(d))
  d <- d[stats::complete.cases(d[, required_cols, drop = FALSE]), ]

  # --- Exclusions --------------------------------------------------------
  # Pupae: older CSVs encoded pupated larvae as `alive == 2`; newer ones
  # add a separate `pupated` column with values "T"/"F"/NA/"". Handle both.
  if (exclude_pupae) {
    if ("alive" %in% names(d)) {
      d <- d[!(d$alive %in% 2), ]
    }
    if ("pupated" %in% names(d)) {
      # Treat "T" as pupated; F, NA, and blank are all non-pupated.
      pup <- toupper(trimws(as.character(d$pupated)))
      d <- d[!(pup %in% "T"), ]
    }
  }

  if (exclude_controls && "tcycler" %in% names(d)) {
    d <- d[d$tcycler != 0, ]
  }

  # Drop the late time point that caused edge effects (sustained-hot only).
  if (!is.null(drop_edge_time) && "time" %in% names(d)) {
    d <- d[d$time < drop_edge_time, ]
  }

  # --- Rename legacy population labels to pop_ids ------------------------
  # e.g. "Florida" -> "FL", "Chicago-AG" -> "IL-AG1"
  lookup <- setNames(pop_info$pop_id, pop_info$csv_name)
  if (!all(d$population %in% names(lookup))) {
    unknown <- unique(d$population[!d$population %in% names(lookup)])
    stop("Unknown population label(s) in ", basename(path), ": ",
         paste(unknown, collapse = ", "),
         "\nAdd them to `pop_info` in config.R or check the CSV.")
  }
  d$population <- unname(lookup[d$population])

  # --- Attach species/region metadata -----------------------------------
  meta <- pop_info[, c("pop_id", "species", "region")]
  names(meta)[1] <- "population"
  d <- merge(d, meta, by = "population", all.x = TRUE)

  # --- Filter to requested populations ----------------------------------
  d <- d[d$population %in% pops, ]
  d$population <- order_pops(d$population)
  d$species    <- factor(d$species, levels = SPECIES_LEVELS)
  d$region     <- factor(d$region,  levels = REGION_LEVELS)

  # --- Ensure numeric alive (some scripts had it as factor) -------------
  d$alive <- as.numeric(as.character(d$alive))

  tibble::as_tibble(d)
}