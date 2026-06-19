# R/06_bioclim.R --------------------------------------------------------
#
# Load CHELSA BIOCLIM rasters, reproject to a common US-wide grid, and
# extract bioclim values at each population's collection point.
#
# Layer names are parsed from filenames (e.g. CHELSA_bio05_...tif ->
# "bio05", CHELSA_gdd5_...tif -> "gdd5"), so the loader works regardless
# of how many CHELSA variables you have in BIOCLIM_DIR — the classic 19
# (bio01..bio19), the full V.2.1 suite (~75 variables), or any subset.
#
# By default we restrict to the classic 19 bio variables for consistency
# with the submitted manuscript; see load_bioclim_us(subset = ...).
#
# Reprojection is slow (often 30+ min on the full 75-variable stack),
# so the result is cached to disk. On subsequent runs the cached stack
# is loaded instead of re-reprojecting.
# ------------------------------------------------------------------------


#' Pretty-print labels for CHELSA variables used in figure captions and
#' regression-output tables. Looked up by the variable code derived from
#' the filename. If a code is not in this table, its own code is used as
#' the label (e.g. "gdd5" stays "gdd5"), so missing labels are not an
#' error — just less pretty.
BIOCLIM_NAMES <- c(
  # Classic BIOCLIM (bio01-19)
  bio01 = "Mean Annual Air Temperature (\u00B0C)",
  bio02 = "Mean Diurnal Air Temperature Range (\u00B0C)",
  bio03 = "Isothermality",
  bio04 = "Temperature Seasonality",
  bio05 = "Max Temperature of Warmest Month (\u00B0C)",
  bio06 = "Min Temperature of Coldest Month (\u00B0C)",
  bio07 = "Annual Temperature Range (\u00B0C)",
  bio08 = "Mean Temperature of Wettest Quarter (\u00B0C)",
  bio09 = "Mean Temperature of Driest Quarter (\u00B0C)",
  bio10 = "Mean Temperature of Warmest Quarter (\u00B0C)",
  bio11 = "Mean Temperature of Coldest Quarter (\u00B0C)",
  bio12 = "Annual Precipitation",
  bio13 = "Precipitation of Wettest Month",
  bio14 = "Precipitation of Driest Month",
  bio15 = "Precipitation Seasonality",
  bio16 = "Precipitation of Wettest Quarter",
  bio17 = "Precipitation of Driest Quarter",
  bio18 = "Precipitation of Warmest Quarter",
  bio19 = "Precipitation of Coldest Quarter"
  # Extended CHELSA variables (gdd*, pet*, vpd*, kg*, hurs*, cmi*, etc.)
  # are supported by the loader but left unlabeled here; add entries
  # above as needed.
)


#' The 19 classic BIOCLIM variable codes, useful as a default subset.
BIOCLIM_CLASSIC_19 <- sprintf("bio%02d", 1:19)


#' Derive a variable code from a CHELSA filename.
#'
#' Parses "CHELSA_<var>_<period>_V.2.1.tif" (or similar) to "<var>".
#' Returns NA for filenames that don't match the CHELSA naming pattern
#' (e.g. a pre-generated stack like `USA_CHELSA_stack.tif`).
chelsa_var_from_filename <- function(filename) {
  base <- basename(filename)
  m <- regmatches(base,
                  regexec("^CHELSA_([A-Za-z0-9]+)_", base))
  vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_,
         character(1))
}


#' Load (and cache) the reprojected BIOCLIM stack for the continental US
#'
#' On first call, reads every `CHELSA_<var>_*.tif` in BIOCLIM_DIR,
#' reprojects each to a common ~1 arcmin grid over the continental US,
#' stacks them, labels each layer with its variable code, and writes
#' the result to `cache/bioclim_us.tif`. Subsequent calls load the
#' cached stack instead.
#'
#' @param subset Character vector of variable codes to keep, or NULL to
#'   keep all. Defaults to `BIOCLIM_CLASSIC_19` (the 19 variables used
#'   in the submitted manuscript). Pass `NULL` to load everything in
#'   BIOCLIM_DIR (useful for exploratory analyses).
#' @param force Logical. If TRUE, rebuild the cached stack from scratch.
#'
#' @return A `terra::SpatRaster` with one layer per retained CHELSA
#'   variable, named with the parsed variable code (e.g. "bio05").
load_bioclim_us <- function(subset = BIOCLIM_CLASSIC_19, force = FALSE) {

  cache_file <- file.path(CACHE_DIR, "bioclim_us.tif")

  # --- Load from cache if available ------------------------------------
  if (file.exists(cache_file) && !force) {
    message("Loading cached BIOCLIM stack: ", cache_file)
    r <- terra::rast(cache_file)

    # terra writes layer names into the file, so they should survive.
    # If the cache was produced by the old version of this function
    # (which left layers unnamed / named NA), rebuild it from source.
    if (any(is.na(names(r))) || all(grepl("^lyr", names(r)))) {
      message("  Cached stack has no/invalid layer names; rebuilding.")
      return(load_bioclim_us(subset = subset, force = TRUE))
    }
    if (!is.null(subset)) {
      keep <- names(r) %in% subset
      if (!any(keep)) {
        stop("None of the requested subset (", paste(subset, collapse = ", "),
             ") is in the cached stack. Cached layers: ",
             paste(names(r), collapse = ", "))
      }
      r <- r[[which(keep)]]
    }
    return(r)
  }

  # --- Build from scratch ----------------------------------------------
  tif_files <- list.files(BIOCLIM_DIR,
                          pattern = "\\.tif$",
                          full.names = TRUE,
                          ignore.case = TRUE)

  # Filter to files that parse as CHELSA variables; drop anything else
  # (e.g. a pre-generated stack called USA_CHELSA_stack.tif).
  var_codes <- chelsa_var_from_filename(tif_files)
  keep_file <- !is.na(var_codes)
  if (any(!keep_file)) {
    message("  Skipping non-CHELSA files in BIOCLIM_DIR: ",
            paste(basename(tif_files[!keep_file]), collapse = ", "))
  }
  tif_files <- tif_files[keep_file]
  var_codes <- var_codes[keep_file]

  # Optionally restrict to a subset of variable codes.
  if (!is.null(subset)) {
    matches       <- var_codes %in% subset
    tif_files     <- tif_files[matches]
    var_codes     <- var_codes[matches]
    missing_codes <- setdiff(subset, var_codes)
    if (length(missing_codes) > 0) {
      warning("Requested variables not found in BIOCLIM_DIR: ",
              paste(missing_codes, collapse = ", "), call. = FALSE)
    }
  }

  if (length(tif_files) == 0) {
    stop("No matching CHELSA .tif files in BIOCLIM_DIR (", BIOCLIM_DIR,
         "). Check your `subset` argument and that the files are named ",
         "CHELSA_<var>_*.tif.")
  }

  # Order by variable code so the output stack is deterministic.
  ord <- order(var_codes)
  tif_files <- tif_files[ord]
  var_codes <- var_codes[ord]

  # ~1 arcmin grid over the continental US (matches the original script).
  template <- terra::rast(
    xmin = -130, xmax = -65, ymin = 23, ymax = 50,
    crs  = "+proj=longlat +datum=WGS84 +no_defs",
    resolution = c(0.008333333 * 2, 0.008333333 * 2)
  )

  message("Reprojecting ", length(tif_files), " CHELSA layers...")
  projected <- lapply(seq_along(tif_files), function(i) {
    message("  [", i, "/", length(tif_files), "] ", basename(tif_files[i]))
    terra::project(terra::rast(tif_files[i]), template)
  })

  stack <- terra::rast(projected)
  names(stack) <- var_codes

  terra::writeRaster(stack, cache_file, overwrite = TRUE)
  message("Cached to ", cache_file)

  stack
}


#' Extract BIOCLIM values at population collection points
#'
#' @param bioclim_stack SpatRaster returned by `load_bioclim_us()`.
#' @param pops Subset of pop_info to extract at (default: `pop_info_main`).
#'
#' @return Tibble with one row per population and one column per bioclim
#'   layer. pop_id / species / region / lat / lon are in the first five
#'   columns.
extract_bioclim_at_pops <- function(bioclim_stack, pops = pop_info_main) {

  pts <- terra::vect(
    as.matrix(pops[, c("lon", "lat")]),
    type = "points",
    crs  = "+proj=longlat +datum=WGS84 +no_defs"
  )
  pts <- terra::project(pts, terra::crs(bioclim_stack))

  vals <- terra::extract(bioclim_stack, pts)
  vals$ID <- NULL  # drop the extract-output ID column

  # Guard against silently-NA layer names (which would break tibble).
  bad <- is.na(names(vals)) | names(vals) == ""
  if (any(bad)) {
    names(vals)[bad] <- paste0("layer_", which(bad))
    warning("Some bioclim layers had missing names; renamed to ",
            paste(names(vals)[bad], collapse = ", "), call. = FALSE)
  }

  out <- cbind(
    pops[, c("pop_id", "species", "region", "lat", "lon")],
    vals
  )
  tibble::as_tibble(out)
}
