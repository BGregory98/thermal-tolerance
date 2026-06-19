# R/03_survival_models.R ------------------------------------------------
#
# Shared GLM machinery used by all three thermal-tolerance experiments.
# Every script calls these functions rather than re-inlining the logic.
#
# Design:
#   * `compare_full_reduced()` — runs the full (interaction) and reduced
#     (additive) binomial GLMs, returns both models plus the LRT.
#   * `lt50_by_pop()` — fits a per-population probit GLM and returns
#     LT50 (plus SE, and optionally 95% CI via dose.p).
#   * `final_survival_by_pop()` — returns the bootstrap-CI final-timepoint
#     proportion alive per population (wrapper around bootstrap_by_pop()).
# ------------------------------------------------------------------------


#' Fit full (interaction) and reduced (additive) binomial GLMs, compare via LRT
#'
#' The full model is `alive ~ predictor * population`; the reduced model
#' is `alive ~ predictor + population`. If the interaction is NOT
#' significant, the caller can optionally reduce further to
#' `alive ~ predictor` by re-running with `reduce_further = TRUE`.
#'
#' @param data A data frame from `load_survival_data()`.
#' @param predictor Character, either "time" or "temp".
#' @param by A grouping variable name for the population-like axis.
#'   Usually "population" but could be "species" or "region" for the
#'   reframed reviewer analysis.
#' @param reduce_further If TRUE, also fits `alive ~ predictor` and
#'   compares it against the additive model.
#' @return List with `full`, `reduced`, `lrt`, and (if reduce_further)
#'   `minimal` and `lrt_further` elements.
compare_full_reduced <- function(data,
                                 predictor = c("time", "temp"),
                                 by = "population",
                                 reduce_further = FALSE) {
  predictor <- match.arg(predictor)

  f_full <- as.formula(sprintf("alive ~ %s * %s", predictor, by))
  f_red  <- as.formula(sprintf("alive ~ %s + %s", predictor, by))

  m_full <- glm(f_full, data = data,
                family = stats::binomial, na.action = na.fail)
  m_red  <- glm(f_red,  data = data,
                family = stats::binomial, na.action = na.fail)

  out <- list(
    full    = m_full,
    reduced = m_red,
    lrt     = lmtest::lrtest(m_full, m_red)
  )

  if (reduce_further) {
    f_min <- as.formula(sprintf("alive ~ %s", predictor))
    m_min <- glm(f_min, data = data,
                 family = stats::binomial, na.action = na.fail)
    out$minimal     <- m_min
    out$lrt_further <- lmtest::lrtest(m_red, m_min)
  }

  out
}


#' Compute LT50 per population using a probit GLM and dose.p
#'
#' For each level of `by`, fits `alive ~ predictor` with a probit link
#' and returns the predictor value at which predicted survival = 0.5,
#' along with standard error and a 95% confidence interval.
#'
#' @param data Data frame from `load_survival_data()`.
#' @param predictor "time" (sustained / cold trials) or "temp" (ramping).
#' @param by Grouping column name. Default "population"; pass "species"
#'   or "region" to pool.
#'
#' @return tibble: <by>, lt50, SE, lower, upper (the latter two are
#'   lt50 ± 1.96 * SE).
lt50_by_pop <- function(data,
                        predictor = c("time", "temp"),
                        by        = "population") {

  predictor <- match.arg(predictor)
  f <- as.formula(sprintf("alive ~ %s", predictor))

  groups <- levels(data[[by]])
  groups <- groups[groups %in% unique(as.character(data[[by]]))]

  rows <- lapply(groups, function(g) {
    sub <- data[as.character(data[[by]]) == g, , drop = FALSE]
    if (nrow(sub) < 2) return(NULL)

    m <- glm(f, data = sub,
             family = stats::binomial(link = "probit"),
             na.action = na.fail)
    lt <- MASS::dose.p(m, p = 0.5)

    tibble::tibble(
      !!by := g,
      lt50 = as.numeric(lt),
      SE   = as.numeric(attr(lt, "SE"))
    )
  })

  out <- dplyr::bind_rows(rows)
  out$ci    <- 1.96 * out$SE
  out$lower <- out$lt50 - out$ci
  out$upper <- out$lt50 + out$ci
  out[[by]] <- factor(out[[by]], levels = groups)
  out
}
