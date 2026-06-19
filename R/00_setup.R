# R/00_setup.R -----------------------------------------------------------
#
# Centralized library loading and ggplot theme. Every analysis script
# source()s this file (via analysis/00_run_all.R or directly).
# ------------------------------------------------------------------------

# --- Core tidyverse-ish packages ---------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)  # dplyr, ggplot2, tidyr, readr, etc.
  library(lubridate)
  library(broom)

  # Modeling
  library(MASS)       # dose.p() for LT50 estimation
  library(lmtest)     # lrtest() for likelihood-ratio tests

  # Geospatial (only needed by 06_bioclim.R and 01_population_map.R)
  library(terra)
  library(raster)
  library(sf)
})

# Note: MASS::select() masks dplyr::select(). If you run into this,
# use dplyr::select() explicitly in any script that loads MASS.


# --- Default ggplot theme ----------------------------------------------
# Matches the visual style of the submitted figures (theme_bw, base_size
# around 16–22) but is consistent across plots.

theme_tt <- function(base_size = 16) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      legend.position  = "right"
    )
}


# --- Safe ggsave helper -----------------------------------------------
# Writes to FIG_DIR (from config.R) using the device settings there.
# Avoids every script having to remember to build its own filepath.

save_fig <- function(plot, filename,
                     width  = FIG_WIDTH_IN,
                     height = FIG_HEIGHT_IN,
                     dpi    = FIG_DPI,
                     device = FIG_DEVICE) {
  path <- file.path(FIG_DIR, paste0(filename, ".", device))
  ggplot2::ggsave(
    filename = path, plot = plot,
    width = width, height = height, units = "in",
    dpi = dpi, device = device
  )
  message("Saved: ", path)
  invisible(path)
}
