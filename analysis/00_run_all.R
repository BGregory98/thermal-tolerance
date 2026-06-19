# analysis/00_run_all.R -------------------------------------------------
#
# Master script. Sourcing this file runs the full pipeline end-to-end.
#
# USAGE:
#   1. setwd() to the project root (the folder that contains config.R).
#   2. source("analysis/00_run_all.R")
# ------------------------------------------------------------------------

setwd("c:/Users/benpg/OneDrive/Desktop/UMD/Research/Thermal_Tolerance/scripts/thermal-tolerance")

# --- Verify we're in the project root ----------------------------------
if (!file.exists("config.R")) {
  stop(
    "config.R not found in the current working directory.\n",
    "  Current working directory: ", getwd(), "\n\n",
    "Before sourcing this file, setwd() to the project root (the\n",
    "folder that contains config.R), then source this file again. e.g.:\n",
    "  setwd(\"/path/to/thermal-tolerance\")\n",
    "  source(\"analysis/00_run_all.R\")",
    call. = FALSE
  )
}

# --- Setup -------------------------------------------------------------
source("config.R")
config_check()

source("R/00_setup.R")
source("R/01_utils.R")
source("R/02_data_loading.R")
source("R/03_survival_models.R")
source("R/04_bootstrapping.R")
source("R/05_plotting.R")
source("R/06_bioclim.R")   # needed again for analyses 08 and 09


# --- Analyses ----------------------------------------------------------
# Each script is self-contained and writes figures / tables to FIG_DIR
# and TAB_DIR. Comment out any you don't need to rebuild.

# Primary analyses — five populations, population-level GLMs
source("analysis/01_population_map.R")       # Figure 1 (all 7 populations)
source("analysis/02_ramping_hot.R")          # Figure 2
source("analysis/03_sustained_hot.R")        # Figure 3
source("analysis/04_cold.R")                 # Figure 4

# Environmental measurements (main + supplementary)
source("analysis/06_container_temps.R")      # Field container temps
source("analysis/07_freezer_calibration.R")  # Supplementary (former Fig 6)

# Descriptive BIOCLIM comparison at the five main collection sites
source("analysis/09_bioclim_descriptive.R")  # Supplementary Table

# Secondary / hybrid-zone extension — first-instar ramping for six
# populations including Baltimore and SERC, with a priori BIOCLIM
# regression against LT50.
source("analysis/08_hybrid_zone.R")          # Supplementary Figures

message("Pipeline complete.")
