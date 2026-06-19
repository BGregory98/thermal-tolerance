# Thermal tolerance in *Culex pipiens* Assemblage larvae — analysis code

Analysis code for Gregory et al., "Stage-specific divergence in larval
thermal tolerance among populations of the *Culex pipiens* Assemblage"
(*Journal of Medical Entomology*).

This directory contains the R pipeline for the high-temperature ramping,
sustained high-temperature, and freezing-exposure thermal tolerance
experiments, plus BIOCLIM environmental regressions, the population
reference map, and the two environmental-temperature supporting analyses
(field container temperatures and freezer calibration).

---

## Project layout

```
thermal-tolerance/
├── README.md                      This file
├── LICENSE                        MIT license
├── CITATION.cff                   Machine-readable citation metadata
├── thermal-tolerance.Rproj        Project file (anchors here(), opens in Positron/RStudio)
├── .gitignore                     Excludes data/, output/, cache/, editor cruft
├── config.R                       Paths, populations, grouping, toggles
├── R/                             Reusable helpers (do not edit to
│   ├── 00_setup.R                 change analyses — edit config.R or
│   ├── 01_utils.R                 the analysis/ scripts)
│   ├── 02_data_loading.R
│   ├── 03_survival_models.R
│   ├── 04_bootstrapping.R
│   ├── 05_plotting.R
│   └── 06_bioclim.R
├── analysis/                      One script per figure / analysis
│   ├── 00_run_all.R               Master — sources everything in order
│   ├── 01_population_map.R        Figure 1
│   ├── 02_ramping_hot.R           Figure 2
│   ├── 03_sustained_hot.R         Figure 3
│   ├── 04_cold.R                  Figure 4
│   ├── 05_bioclim_regressions.R   BIOCLIM ~ phenotype tables
│   ├── 06_container_temps.R       Container-habitat temperatures
│   ├── 07_freezer_calibration.R   Supplementary freezer ramp (Fig. 6)
│   ├── 08_hybrid_zone.R           Hybrid-zone ramping (Baltimore, SERC)
│   └── 09_bioclim_descriptive.R   Descriptive BIOCLIM summaries
├── data/                          (you provide — see "Input data" below)
├── output/
│   ├── figures/                   Generated figures (auto-created)
│   └── tables/                    Generated CSV tables (auto-created)
└── cache/                         Cached R objects & reprojected rasters
                                   (auto-created)
```

---

## Input data

Put the following files in `data/` (or edit `DATA_DIR` in `config.R` to
point elsewhere):

| File                      | Contents                                              |
|---------------------------|-------------------------------------------------------|
| `tramp_cx.csv`            | Ramping hot, first instar (per-individual survival)   |
| `tramp_cx_4.csv`          | Ramping hot, fourth instar                            |
| `chronic_first.csv`       | Sustained 38 °C, first instar                         |
| `chronic_fourth.csv`      | Sustained 38 °C, fourth instar                        |
| `cold_first.csv`          | Freezing exposure, first instar                       |
| `cold_fourth.csv`         | Freezing exposure, fourth instar                      |
| `balt713_trimmed.csv`     | Baltimore container-habitat temps (sun/shade, 2023)   |
| `ColdRamp_400.csv`        | Freezer calibration, 6 replicates at 400 mL           |

Expected columns in the survival CSVs:
`population`, `trial`, `alive` (0/1, or 2 if pupated),
`tcycler` (thermal cycler ID; `0` = control), `row`, `column`, and
either `temp` (ramping) or `time` (sustained / cold).

**CHELSA BIOCLIM rasters** go in `data/CHELSA/` (or edit `BIOCLIM_DIR`).
Files must match `CHELSA_*.tif`.

> The `data/` directory is git-ignored, so none of these files are tracked
> in the repository. See **Data availability** below for where to obtain them.

---

## Data availability

The survival-trial CSVs and the environmental-temperature CSVs are archived
on Dryad: <https://doi.org/10.5061/dryad.x69p8d00x>. Download them into
`data/` (and `data/CHELSA/` for the rasters) before running the pipeline.

The CHELSA BIOCLIM rasters are publicly available from the CHELSA project
(<https://chelsa-climate.org>); this repository does not redistribute them.

---

## How to run

### Fresh machine

1. Clone / copy this directory onto your machine.
2. Open `config.R` and confirm or edit:
   - `PROJECT_ROOT`, `DATA_DIR`, `BIOCLIM_DIR`, `OUTPUT_DIR`
   - `POPS_MAIN` — which populations to include in the main analysis
3. Install required packages (see next section).
4. Set R's working directory to the project root and run the master
   script:

```r
setwd("/path/to/thermal-tolerance")
source("analysis/00_run_all.R")
```

This runs every analysis in order and writes figures to `output/figures/`
and tables to `output/tables/`.

### Running one analysis at a time

Each script is self-contained given the shared setup. Source the
master-script header first:

```r
setwd("/path/to/thermal-tolerance")
source("config.R"); config_check()
source("R/00_setup.R")
source("R/01_utils.R")
source("R/02_data_loading.R")
source("R/03_survival_models.R")
source("R/04_bootstrapping.R")
source("R/05_plotting.R")
source("R/06_bioclim.R")

source("analysis/02_ramping_hot.R")   # just Figure 2
```

---

## Required R packages

```r
install.packages(c(
  "tidyverse", "lubridate", "broom",
  "MASS", "lmtest",
  "terra", "raster", "sf",
  "usmap", "here"
))
```

---

## Reproducibility

This pipeline was developed and run under **R 4.5.1**.
To record the exact package versions used for the published results, run
`sessionInfo()` (or `sessioninfo::session_info()`) at the end of a full run
and paste the output here, or commit it as `session_info.txt`.

For a fully pinned environment, consider initializing
[`renv`](https://rstudio.github.io/renv/) (`renv::init()`); the resulting
`renv.lock` makes the package versions exactly reproducible on another
machine.

---

## License

Code in this repository is released under the MIT License. See
[`LICENSE`](LICENSE) for the full text.

---

## How to cite

If you use this code, please cite both the software and the associated
article. Citation metadata is in [`CITATION.cff`](CITATION.cff), and GitHub
renders a "Cite this repository" button from it. The Dryad data DOI is
already recorded there; add the article DOI once it is assigned.
