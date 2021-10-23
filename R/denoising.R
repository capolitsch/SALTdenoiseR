#' Denoise a polarized spectrum via quadratic trend filtering and compute
#' uncertainties via a bootstrap
#'
#' @param wavelength Vector of wavelength measurements.
#' @param stokes Polarized spectrum measurements, passed as a 3-column tibble,
#' data frame, or matrix, with the columns corresponding to the Stokes
#' parameters I, Q, and U, respectively.
#' @param variances Measurement variances, in a tibble, data frame, or matrix
#' with dimensions matching those of `stokes`.
#' @param masks Pixel masks, in a tibble, data frame, or matrix
#' with dimensions matching those of `stokes`. Nonzero elements flag bad pixels.
#' @param break_at A free parameter that controls the segmentation of a
#' spectrum. More precisely, `break_at` is the minimum number of
#' consecutively-masked spectral pixels that will trigger a break in the
#' spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, the
#' resulting segments are examined to ensure that each is sufficiently long for
#' an independent denoising analysis. In particular, any segment that has less
#' than `min_pix_segment` unmasked spectral pixels is discarded. Defaults to
#' `min_pix_segment = 10`.
#' @param compute_uncertainties (Boolean) If `TRUE`, then variability bands are
#' computed for each of the denoised normalized Stokes spectra, via a parametric
#' bootstrap algorithm.
#' @param ... Additional named arguments to be passed to
#' [trendfiltering::sure_trendfilter()] or
#' [trendfiltering::bootstrap_trendfilter()].
#' `bootstrap_algorithm = "parametric"` is set automatically and cannot be
#' overridden.
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
#' library(dplyr, quietly = TRUE)
#'
#' data(polarized_spectrum_WR_star)
#'
#' wavelength <- seq(
#'   from = sci$axDat$crval[1],
#'   by = sci$axDat$cdelt[1],
#'   length = sci$axDat$len[1]
#' )
#'
#' stokes <- as_tibble(sci$imDat)
#' variances <- as_tibble(var$imDat) %>% select(1:3)
#' masks <- as_tibble(bpm$imDat)
#'
#' spec_denoised <- denoise_polarized_spectrum(wavelength, stokes, variances, masks)
#' @importFrom trendfiltering sure_trendfilter bootstrap_trendfilter
#' @importFrom glmgen trendfilter trendfilter.control.list
#' @importFrom tidyr drop_na tibble as_tibble
#' @importFrom dplyr %>% arrange filter select n_distinct bind_cols
#' @importFrom magrittr %$%
#' @importFrom parallel mclapply detectCores
denoise_polarized_spectrum <- function(wavelength,
                                       stokes,
                                       variances,
                                       masks,
                                       break_at = 10,
                                       min_pix_segment = 10,
                                       compute_uncertainties = FALSE,
                                       mc_cores = parallel::detectCores(),
                                       ...) {
  stopifnot(ncol(stokes) == 3)
  if (missing(variances)) {
    stop("Currently, variances must be passed for this denoising analysis.")
  }
  stopifnot(ncol(variances) == 3 & nrow(variances) == nrow(stokes))

  if (missing(masks)) {
    masks <- matrix(rep(0, 3 * nrow(stokes)), ncol = 3)
  } else {
    stopifnot(ncol(masks) == 3 & nrow(masks) == nrow(stokes))
  }

  tf_args <- list(...)
  tf_args$bootstrap_algorithm <- "parametric"

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
  )

  df_list <- break_spectrum(df_full, break_at, min_pix_segment)

  sure_tf <- mclapply(
    1:(3 * length(df_list)),
    parallel_sure_tf,
    df_list = df_list,
    sure_args = sure_args,
    mc.cores = mc_cores
  )

  boot_tf <- mclapply(
    1:length(sure_tf),
    parallel_bootstrap_tf,
    sure_tf = sure_tf,
    bootstrap_args = bootstrap_args,
    mc.cores = mc_cores
  )
}


#' @noRd
parallel_sure_tf <- function(X, df_list, sure_args) {
  if (X %in% 1:length(df_list)) {
    args <- c(
      list(
        x = df_list[[X]]$wavelength,
        y = df_list[[X]]$I,
        weights = 1 / df_list[[X]]$I_vars,
        x_eval = df_list[[X]]$wavelength
      ),
      sure_args
    )
  }

  if (X %in% (length(df_list) + 1):(2 * length(df_list))) {
    args <- c(
      list(
        x = df_list[[X - 3]]$wavelength,
        y = df_list[[X - 3]]$Q,
        weights = 1 / df_list[[X - 3]]$Q_vars,
        x_eval = df_list[[X - 3]]$wavelength
      ),
      sure_args
    )
  }

  if (X %in% (2 * length(df_list) + 1):(3 * length(df_list))) {
    args <- c(
      list(
        x = df_list[[X - 6]]$wavelength,
        y = df_list[[X - 6]]$U,
        weights = 1 / df_list[[X - 6]]$U_vars,
        x_eval = df_list[[X - 6]]$wavelength
      ),
      sure_args
    )
  }

  do.call(sure_trendfilter, args)
}


#' @noRd
parallel_bootstrap_tf <- function(X, sure_tf, bootstrap_args) {
  args <- c(list(obj = sure_tf[[X]]), bootstrap_args)
  do.call(bootstrap_trendfilter, args)
}
