# config.R ---------------------------------------------------------------
#
# PROJECT CONFIGURATION — edit this file when:
#   * Moving the project to a new machine
#   * Dropping / adding populations
#   * Changing the grouping scheme
#
# Everything downstream reads from this file. Do not hard-code paths or
# population lists anywhere else.
# ------------------------------------------------------------------------

# ---- 1. Paths ---------------------------------------------------------

PROJECT_ROOT <- normalizePath("C:/Users/benpg/OneDrive/Desktop/UMD/Research/Thermal_Tolerance/scripts/thermal-tolerance", mustWork = FALSE)

DATA_DIR    <- file.path(PROJECT_ROOT, "data")
BIOCLIM_DIR <- file.path(PROJECT_ROOT, "data", "CHELSA")
OUTPUT_DIR  <- file.path(PROJECT_ROOT, "output")
FIG_DIR     <- file.path(OUTPUT_DIR, "figures")
TAB_DIR     <- file.path(OUTPUT_DIR, "tables")
CACHE_DIR   <- file.path(PROJECT_ROOT, "cache")

for (d in c(OUTPUT_DIR, FIG_DIR, TAB_DIR, CACHE_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ---- 2. Populations ---------------------------------------------------
# pop_info is the master population table. ONE row per population.
#
# Columns:
#   pop_id      — short code used throughout the analysis and in figures
#   csv_name    — the string that appears in the `population` column of
#                 the raw CSVs (thermal-trial data files). This is what
#                 gets renamed to pop_id during data loading.
#   species     — taxonomic assignment; "unknown" for mid-latitude
#                 populations that have not been genotyped
#   region      — used descriptively in figures; "north" / "south" for
#                 the five main populations, "mid" for the two
#                 mid-latitude populations added in the hybrid-zone
#                 extension
#   lat, lon    — collection coordinates (decimal degrees, WGS84)
#   color       — hex code for plots (keyed by pop_id)

pop_info <- data.frame(
  pop_id   = c("FL",       "TX",       "ME",       "IL-AG3",   "IL-BG2",
               "Baltimore","SERC"),
  csv_name = c("Florida",  "Texas",    "Maine",    "Chicago-AG","Chicago-BG",
               "Baltimore","SERC"),
  species  = c("quinquefasciatus", "quinquefasciatus",
               "pipiens",  "pipiens",  "pipiens_molestus",
               "unknown",  "unknown"),
  region   = c("south",    "south",    "north",    "north",    "north",
               "mid",      "mid"),
  lat      = c(27.568669,  30.600440,  44.892258,  42.094783,  41.650225,
               39.291699,  38.888256),
  lon      = c(-80.417803, -96.268930, -68.671079, -87.770168, -87.600140,
               -76.627335, -76.552842),
  color    = c("#DC9F1F",  "#ECE24E",  "#3373B0",  "#6FB4E7",  "#C27CA6",
               "#7FC97F",  "#33A02C"),
  stringsAsFactors = FALSE
)

# -- Which populations are in the MAIN thermal-tolerance analysis? -----
# The five populations from the original submitted manuscript. Per the
# meeting with Megan (2026-04-24), we are reverting to this set and NOT
# dropping IL-BG2.

POPS_MAIN <- c("FL", "TX", "ME", "IL-AG3", "IL-BG2")

# -- Which populations are in the HYBRID-ZONE extension? ---------------
# First-instar ramping only. Excludes IL-BG2 (lab-adapted). Adds
# Baltimore and SERC (mid-latitude). Used by
# analysis/08_hybrid_zone.R for a supplementary BIOCLIM regression.

POPS_HZ <- c("FL", "TX", "ME", "IL-AG3", "Baltimore", "SERC")

pop_info_main <- pop_info[pop_info$pop_id %in% POPS_MAIN, ]
pop_info_hz   <- pop_info[pop_info$pop_id %in% POPS_HZ,   ]

# Factor order for plots (south → mid → north by latitude).
POP_LEVELS <- pop_info$pop_id[order(pop_info$lat)]
POP_LEVELS_MAIN <- POP_LEVELS[POP_LEVELS %in% POPS_MAIN]
POP_LEVELS_HZ   <- POP_LEVELS[POP_LEVELS %in% POPS_HZ]

# Named color vector keyed by pop_id.
pop_colors <- setNames(pop_info$color, pop_info$pop_id)


# ---- 3. Grouping metadata (descriptive only) ---------------------------
# These levels are kept for figure legend ordering. Per the meeting,
# populations are analyzed individually — not pooled by region or
# species — so these are not used as predictors in any GLM.

REGION_LEVELS  <- c("south", "mid", "north")
SPECIES_LEVELS <- c("quinquefasciatus", "pipiens", "pipiens_molestus", "unknown")

species_labels <- c(
  quinquefasciatus = "Cx. quinquefasciatus",
  pipiens          = "Cx. pipiens",
  pipiens_molestus = "Cx. pipiens f. molestus",
  unknown          = "unassigned"
)


# ---- 4. Analysis toggles ----------------------------------------------

BOOT_N     <- 10000
BOOT_PROBS <- c(0.025, 0.50, 0.975)

RAMP_TIME_MAX       <- 46
SUSTAINED_TEMP      <- 38
SUSTAINED_TIME_MAX  <- 96
COLD_TIME_MAX       <- 50
FREEZER_COOLDOWN_MIN <- 60


# ---- 5. BIOCLIM variable selection ------------------------------------
# A priori set of non-redundant CHELSA BIOCLIM variables used in the
# hybrid-zone regression (analysis/08) and the descriptive summary at
# collection sites (analysis/09). Chosen on biological grounds rather
# than data-driven correlation filtering.
#
# bio05 — Max Temperature of Warmest Month    (relevant to heat assays)
# bio06 — Min Temperature of Coldest Month    (relevant to cold assays)
# bio04 — Temperature Seasonality             (generalist vs specialist
#                                              thermal phenotype)

A_PRIORI_BIOCLIM <- c("bio04", "bio05", "bio06")


# ---- 6. Figure output settings ----------------------------------------
FIG_WIDTH_IN  <- 9
FIG_HEIGHT_IN <- 5
FIG_DPI       <- 300
FIG_DEVICE    <- "png"


# ---- 7. Sanity check --------------------------------------------------
config_check <- function() {
  stopifnot(
    dir.exists(PROJECT_ROOT),
    all(POPS_MAIN %in% pop_info$pop_id),
    all(POPS_HZ   %in% pop_info$pop_id),
    all(pop_info$region  %in% REGION_LEVELS),
    all(pop_info$species %in% SPECIES_LEVELS)
  )
  message("config.R OK.")
  message("  Main populations: ", paste(POPS_MAIN, collapse = ", "))
  message("  Hybrid-zone populations: ", paste(POPS_HZ, collapse = ", "))
  invisible(TRUE)
}
