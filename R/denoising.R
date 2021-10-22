#' Denoise a SALT Observatory spectrum via quadratic trend filtering and compute
#' uncertainties via a bootstrap
#'
#' @param stokes Spectropolarimetric measurements in a Stokes parametrization,
#' passed as a 3-column matrix, with the columns corresponding to the I, Q, U
#' Stokes parameters, respectively.
#' @param variances Measurement variances, in a matrix with dimensions matching
#' those of `stokes`.
#' @param masks Pixel masks, in a matrix with dimensions matching those
#' of `stokes`. Nonzero elements flag bad pixels.
#' @param break_at A free parameter that controls the segmentation of a
#' spectrum. More precisely, `break_at` is the minimum number of
#' consecutively-masked spectral pixels that will trigger a break in the
#' spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, it is
#' advisable to examine the resulting segments to ensure that each is
#' sufficiently long for its own denoising analysis. In particular, we discard
#' any segments that have less than `min_pix_segment` unmasked spectral pixels.
#' Defaults to `min_pix_segment = 10`.
#' @param compute_uncertainties (Boolean) If `TRUE` then variability bands are
#' computed for each of the denoised normalized Stokes spectra, via a parametric
#' bootstrap algorithm.
#' @param ... Additional named arguments to be passed to
#' [trendfiltering::sure_trendfilter()] or
#' [trendfiltering::bootstrap_trendfilter()].
#' @return An object of class
#' [`'polarized_spectrum_denoised'`][denoise_polarized_spectrum()].
#' This is a list with the following elements:
#'
#' @export denoise_polarized_spectrum
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
#'
#' wavelength <- seq(
#'   from = sci$axDat$crval[1],
#'   by = sci$axDat$cdelt[1],
#'   length = sci$axDat$len[1]
#' )
#'
#' stokes <- as_tibble(sci$imDat) %>%
#'   rename_with(function(.cols) c("I", "Q", "U"))
#' variances <- as_tibble(var$imDat) %>%
#'   select(1:3) %>%
#'   rename_with(function(.cols) c("I_vars", "Q_vars", "U_vars"))
#' masks <- as_tibble(bpm$imDat) %>%
#'   rename_with(function(.cols) c("I_mask", "Q_mask", "U_mask"))
#' @importFrom trendfiltering sure_trendfilter bootstrap_trendfilter
#' @importFrom glmgen trendfilter trendfilter.control.list
#' @importFrom tidyr drop_na tibble as_tibble
#' @importFrom dplyr %>% arrange filter select n_distinct bind_cols
#' @importFrom magrittr %$%
denoise_polarized_spectrum <- function(wavelength, stokes, variances, masks,
                                       break_at = 10, min_pix_segment = 10,
                                       compute_uncertainties = FALSE,
                                       mc_cores = parallel::detectCores(),
                                       ...) {
  tf_args <- list(...)

  sure_args <- tf_args[
    which(names(tf_args) %in% names(formals(sure_trendfilter)))
  ]

  bootstrap_args <- tf_args[
    which(names(tf_args) %in% names(formals(bootstrap_trendfilter)))
  ]

  wavelength <- as_tibble_col(wavelength, column_name = "wavelength")
  stokes <- as_tibble(stokes) %>%
    rename_with(function(.cols) c("I", "Q", "U"))
  variances <- as_tibble(variances) %>%
    rename_with(function(.cols) c("I_vars", "Q_vars", "U_vars"))
  masks <- as_tibble(masks) %>%
    rename_with(function(.cols) c("I_mask", "Q_mask", "U_mask"))

  df_full <- bind_cols(
    wavelength,
    stokes,
    variances,
    masks
  ) %>%
    arrange(wavelength)

  df_list <- break_spectrum(df_full, break_at, min_pix_segment)

  sure_tf <- mclapply(
    1:(3 * length(df_list)),
    parallel_sure_tf,
    df_list = df_list,
    sure_args = sure_args,
    mc.cores = mc_cores
  )
}


parallel_sure_tf <- function(X, df_list, sure_args) {
  if (X %in% 1:length(df_list)) {
    df_list[[X]] %$% sure_trendfilter(wavelength, I, 1 / I_vars, x_eval = wavelength)
  }

  if (X %in% (length(df_list) + 1):(2 * length(df_list))) {
    df_list[[X - 3]] %$% sure_trendfilter(wavelength, Q, 1 / Q_vars, x_eval = wavelength)
  }

  if (X %in% (2 * length(df_list) + 1):(3 * length(df_list))) {
    df_list[[X - 6]] %$% sure_trendfilter(wavelength, U, 1 / U_vars, x_eval = wavelength)
  }
}
