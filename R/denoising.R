#' Denoise SALT Observatory spectrum via trend filtering
#'
#' @param sci spectroscopic measurements.
#' @param var measurement variances.
#' @param bpm pixel masks.
#' @param min_mask_width (numeric) Parameter that determines the segmentation of
#' the spectrum into smaller pieces divided by gaps (many consecutive masked
#' spectral pixels). This is the minimum number of consecutive pixels that have
#' to be masked to cause a break in the spectrum.
#' @param variability_bands (Logical) If `TRUE` then variability bands are
#' computed via a bootstrap resampling algorithm based on the uncertainties
#' in the spectral pixel measurements.
#' @param B (integer) The number of bootstrap samples that are drawn in order to
#' estimate the trend filtering variability bands. Larger is better, but is more
#' computationally expensive. My trend filtering bootstrap implementation is set
#' up for parallel computing for speedups when multiple cores are available.
#' `B = 100` typically suffices for a most scenarios, but it's good to crank it
#' up a bit for a plot that's going to be published, so you get nice smooth
#' bands.
#' @param max_iter Maximum iterations allowed for the ADMM trend filtering
#' convex optimization (Ramdas and Tibshirani 2016). Increase this if the trend
#' filtering estimate does not appear to have fully converged to a reasonable
#' estimate of the spectral signal. An estimate that has not
#' fully converged will typically increasingly diverge above
#' or below the observed fluxes as you move from left to right.
#' @param obj_tol (integer) The tolerance used in the convex optimization
#' stopping criterion; when the relative change in the
#' objective function is less than this value, the algorithm
#' terminates. Decrease this if the trend
#' filtering estimate does not appear to have fully converged
#' to a reasonable estimate of the signal. An estimate that
#' has not fully converged will typically increasingly
#' diverge above or below the observed fluxes as you move
#' from left to right
#' @return An object of class [`SALT_denoised`][denoise()]. This is a list
#' with the following elements:
#' \item{x_eval}{Input grid used to evaluate the optimized trend filtering
#' estimate on.}
#' \item{tf_estimate}{Optimized trend filtering estimate, evaluated at
#' `x_eval`.}
#'
#' @export denoise
#'
#' @references
#' \enumerate{
#' \item{Politsch et al. (2020a).
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{
#' Trend filtering – I. A modern statistical tool for time-domain astronomy and
#' astronomical spectroscopy}. \emph{MNRAS}, 492(3), p. 4005-4018.} \cr
#' \item{Politsch et al. (2020b).
#' \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{
#' Trend Filtering – II. Denoising astronomical signals with varying degrees of
#' smoothness}. \emph{MNRAS}, 492(3), p. 4019-4032.}}
#'
#' @examples
#' data(SALT_spectrum)
#' @importFrom glmgen trendfilter trendfilter.control.list
#' @importFrom tidyr drop_na tibble
#' @importFrom dplyr %>% arrange filter select n_distinct
#' @importFrom magrittr %$%
denoise <- function(sci,
                    var,
                    bpm,
                    min_mask_width = 20,
                    variability_bands = FALSE,
                    ...) {
  wavelength <- seq(
    from = sci$axDat$crval[1],
    by = sci$axDat$cdelt[1],
    length = sci$axDat$len[1]
  )

  df <- tibble(
    wavelength = wavelength,
    Q = sci$imDat[, , 2],
    Q_var = var$imDat[, , 2],
    Q_mask = bpm$imDat[, , 2],
    U = sci$imDat[, , 3],
    U_var = var$imDat[, , 3],
    U_mask = bpm$imDat[, , 3],
    I = sci$imDat[, , 1],
    I_var = var$imDat[, , 1],
    I_mask = bpm$imDat[, , 1]
  )

  df$mask <- df$Q_mask
  df <- mask.intervals(df, min_mask_width)
  n_segments <- df %>%
    pull(segment) %>%
    n_distinct()

  out <- vector(mode = "list", length = n_segments)
  for (itr in 1:n_segments) {
    df_segment <- df %>% filter(segment == itr)
    tf_out <- trendfilter_interval(
      x = df_segment$wavelength,
      y = df_segment$Q,
      weights = 1 / df_segment$Q_var,
      obj_tol = obj_tol,
      max_iter = max_iter
    )
    out[[itr]] <- tf_out
  }
}

#' @importFrom trendfiltering sure_trendfilter bootstrap_trendfilter
#' @importFrom dplyr %>%
#' @importFrom magrittr %$%
trendfilter_interval <- function(x, y, weights,
                                 B = 100, alpha = 0.05, bootstrap_bands = T,
                                 max_iter = 5000, obj_tol = 1e-12, ...) {
  SURE_tf <- sure_trendfilter(x, y, weights,
    optimization_params = list(
      max_iter = max_iter,
      obj_tol = obj_tol
    )
  )
  if (!bootstrap_bands) {
    return(SURE_tf)
  } else {
    boot_out <- bootstrap_trendfilter(obj = SURE_tf, B = B, alpha = alpha, ...)
    return(boot_out)
  }
}
