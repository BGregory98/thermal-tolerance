# R/01_utils.R -----------------------------------------------------------
#
# General-purpose helpers used across multiple analysis scripts.
# Define each function ONCE here rather than re-defining in each script.
# ------------------------------------------------------------------------


#' Bootstrap 95% confidence interval of a mean
#'
#' Originally written by Megan Fritz. Given a numeric vector `x`, draws
#' `N` bootstrap samples (with replacement, same size as `x`), computes
#' the mean of each, and returns the 2.5%, 50%, and 97.5% percentiles
#' of that distribution.
#'
#' @param x A numeric vector (typically 0/1 survival data for one
#'   population × timepoint/temperature).
#' @param N Number of bootstrap replicates. Defaults to BOOT_N from config.
#' @return Named numeric vector of length 3: lower, median, upper.
boot_ci <- function(x, N = BOOT_N) {
  boot_means <- replicate(N, mean(sample(x, size = length(x), replace = TRUE)))
  stats::quantile(boot_means, probs = BOOT_PROBS, na.rm = TRUE)
}


#' Convert log-odds to probability
#'
#' @param logit Numeric vector of log-odds.
#' @return Probabilities in [0,1].
logit2prob <- function(logit) {
  odds <- exp(logit)
  odds / (1 + odds)
}


#' Apply pop_id factor ordering defined in config.R
#'
#' Ensures consistent population order in every plot without having to
#' re-specify the levels each time.
#'
#' @param x A character or factor vector of population IDs.
#' @return A factor with levels in POP_LEVELS order, unused levels dropped.
order_pops <- function(x) {
  factor(as.character(x), levels = POP_LEVELS)
}
