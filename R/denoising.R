#' Denoise a SALT Observatory spectrum via quadratic trend filtering and compute
#' uncertainties via a Gaussian bootstrap
#'
#' @param sci spectroscopic measurements.
#' @param var measurement variances.
#' @param bpm pixel masks.
#' @param min_mask_width Parameter that controls the segmentation of a spectrum
#' into ``sub-spectra'', which are then denoised independently. More precisely,
#' `min_mask_width` is the minimum number of consecutively-masked spectral
#' pixels in order to trigger a break in the spectrum. Defaults to
#' `min_mask_width = 20`.
#' @param nx_eval Integer. If nothing is passed to `x_eval`, then it is defined
#' as `x_eval = seq(min(x), max(x), length = nx_eval)`.
#' @param x_eval (Optional) A grid of inputs to evaluate the optimized trend
#' filtering estimate on. May be ignored, in which case the grid is determined
#' by `nx_eval`.
#' @param compute_uncertainties (Boolean) If `TRUE` then variability bands are
#' computed via a parametric bootstrap algorithm based on the uncertainties
#' in the spectral pixel measurements.
#' @param B (Integer) The number of bootstrap samples that are drawn in order to
#' estimate the trend filtering variability bands. Larger is better, but is more
#' computationally expensive. My trend filtering bootstrap implementation is set
#' up for parallel computing for speedups when multiple cores are available.
#' `B = 100` typically suffices for a most scenarios, but it's good to crank it
#' up a bit for a plot that's going to be published, so you get nice smooth
#' bands.
#' @return An object of class [`'SALT_spec'`][denoise_spectrum()]. This is a
#' list with the following elements:
#' \item{x_eval}{Input grid used to evaluate the optimized trend filtering
#' estimate on.}
#' \item{tf_estimate}{Optimized trend filtering estimate, evaluated at
#' `x_eval`.}
#' \item{validation_method}{"SURE"}
#' \item{lambdas}{Vector of hyperparameter values evaluated in the grid search
#' (always returned in descending order).}
#' \item{edfs}{Vector of effective degrees of freedom for all trend filtering
#' estimators fit during validation.}
#' \item{generalization_errors}{Vector of SURE generalization error estimates,
#' corresponding to the descending-ordered `lambdas` vector.}
#' \item{lambda_min}{Hyperparameter value that minimizes the SURE generalization
#' error curve.}
#' \item{i_min}{Index of `lambdas` that minimizes the SURE error curve.}
#' \item{edf_min}{Effective degrees of freedom of the optimized trend
#' filtering estimator.}
#' \item{n_iter}{The number of iterations needed for the ADMM algorithm to
#' converge within the given tolerance, for each hyperparameter value. If many
#' of these are exactly equal to `max_iter`, then their solutions have not
#' converged with the tolerance specified by `obj_tol`. In which case, it is
#' often prudent to increase `max_iter`.}
#' \item{x}{Vector of observed inputs.}
#' \item{y}{Vector of observed outputs.}
#' \item{weights}{Weights for the observed outputs, defined as the reciprocal
#' variance of the additive noise that contaminates the signal.}
#' \item{fitted_values}{Optimized trend filtering estimate, evaluated at the
#' observed inputs `x`.}
#' \item{residuals}{`residuals = y - fitted_values`}
#'
#' @export denoise_spectrum
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
denoise_spectrum <- function(sci,
                             var,
                             bpm,
                             min_mask_width = 20,
                             compute_uncertainties = FALSE,
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
  df <- mask_intervals(df, min_mask_width)
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
