#' Denoise a polarized spectrum
#'
#' \loadmathjax Denoise a polarized spectrum.
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
#' @param sure_args Additional named arguments to be passed to
#' [sure_trendfilter()][trendfiltering::sure_trendfilter()].
#' @param bootstrap_args Additional named arguments to be passed to
#' [bootstrap_trendfilter()][trendfiltering::bootstrap_trendfilter()]. The
#' argument `bootstrap_algorithm = "parametric"` is fixed due to the
#' design of SALT spectra, and we don't allow this argument to be manually
#' overridden. See the
#' [bootstrap_trendfilter()][trendfiltering::bootstrap_trendfilter()]
#' documentation and
#' [Politsch et al. (2020a)](
#' https://academic.oup.com/mnras/article/492/3/4005/5704413)
#' for details on why this particular bootstrap algorithm is appropriate for
#' SALT spectra.
#' @return An object of class `'polarized_spectrum'`. This is a list
#' with the following elements:
#' \item{n_segments}{The number of segments the spectrum was broken into.}
#' \item{denoised_spectra}{A tibble containing denoised spectra for
#' each of the 3 Stokes parameters, as well as Q/I and U/I. The full column set
#' is `c("wavelength","stokes_I","stokes_Q","stokes_U","Q_norm","U_norm")`.}
#' \item{stokes_I_ensemble}{If `compute_uncertainties = TRUE`, the
#' full bootstrap ensemble for the Stokes I parameter,
#' as an \mjseqn{n \times B} matrix, less any columns from post-hoc pruning.
#' If `compute_uncertainties = FALSE`, this will return `NULL`.}
#' \item{stokes_Q_ensemble}{Same as above, but for the Q Stokes parameter.}
#' \item{stokes_U_ensemble}{Same as above, but for the U Stokes parameter.}
#' \item{Q_norm_ensemble}{Same as above, but for Q/I.}
#' \item{U_norm_ensemble}{Same as above, but for U/I.}
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
#' data(polarized_spectrum_WR_star)
#'
#' suppressMessages(library(dplyr))
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
#' @importFrom tibble as_tibble_col
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
                                       sure_args,
                                       bootstrap_args) {
  stopifnot(ncol(stokes) == 3)
  if (missing(variances)) {
    stop("Currently, variances must be passed for this denoising analysis.")
  }
  stopifnot(ncol(variances) == 3 & nrow(variances) == nrow(stokes))
  if (missing(bootstrap_args)) {
    bootstrap_args <- list()
  } else {
    stopifnot(all(names(bootstrap_args) %in% names(formals(bootstrap_trendfilter))))
  }
  if (missing(sure_args)) {
    sure_args <- list()
  } else {
    stopifnot(all(names(sure_args) %in% names(formals(sure_trendfilter))))
  }

  if (missing(masks)) {
    masks <- matrix(rep(0, 3 * nrow(stokes)), ncol = 3)
  } else {
    stopifnot(ncol(masks) == 3 & nrow(masks) == nrow(stokes))
  }

  bootstrap_args$bootstrap_algorithm <- "parametric"
  if (compute_uncertainties) bootstrap_args$return_ensemble <- TRUE
  mc_cores <- parallel::detectCores()

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

  tf_obj <- mclapply(
    X = 1:(3 * length(df_list)),
    parallel_sure_tf,
    df_list = df_list,
    sure_args = c(sure_args, x_eval = df_list[[X]]$wavelength),
    mc.cores = mc_cores
  )

  if (compute_uncertainties) {
    tf_obj <- mclapply(
      1:(3 * length(df_list)),
      parallel_bootstrap_tf,
      sure_tf = tf_obj,
      bootstrap_args = bootstrap_args,
      mc.cores = bootstrap_args$mc_cores
    )
  }

  wavelength_eval <- lapply(
    X = 1:length(df_list),
    FUN = function(X) tf_obj[[X]][["x_eval"]]
  )

  stokes_I_denoised <- lapply(
    X = 1:length(df_list),
    FUN = function(X) tf_obj[[X]][["tf_estimate"]]
  )

  stokes_Q_denoised <- lapply(
    X = (length(df_list) + 1):(2 * length(df_list)),
    FUN = function(X) tf_obj[[X]][["tf_estimate"]]
  )

  stokes_U_denoised <- lapply(
    X = (2 * length(df_list) + 1):(3 * length(df_list)),
    FUN = function(X) tf_obj[[X]][["tf_estimate"]]
  )

  Q_norm_denoised <- lapply(
    X = 1:length(df_list),
    FUN = function(X) stokes_Q_denoised[[X]] / stokes_I_denoised[[X]]
  )

  U_norm_denoised <- lapply(
    X = 1:length(df_list),
    FUN = function(X) stokes_U_denoised[[X]] / stokes_I_denoised[[X]]
  )

  denoised_spectra <- tibble(
    wavelength = wavelength_eval,
    stokes_I = stokes_I_denoised,
    stokes_Q = stokes_Q_denoised,
    stokes_U = stokes_U_denoised,
    Q_norm = Q_norm_denoised,
    U_norm = U_norm_denoised
  )

  if (compute_uncertainties) {
    stokes_I_ensemble <- lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]][["tf_bootstrap_ensemble"]]
    )

    stokes_Q_ensemble <- lapply(
      X = (length(df_list) + 1):(2 * length(df_list)),
      FUN = function(X) tf_obj[[X]][["tf_bootstrap_ensemble"]]
    )

    stokes_U_ensemble <- lapply(
      X = (2 * length(df_list) + 1):(3 * length(df_list)),
      FUN = function(X) tf_obj[[X]][["tf_bootstrap_ensemble"]]
    )

    Q_norm_ensemble <- lapply(
      X = 1:length(df_list),
      FUN = function(X) stokes_Q_ensemble[[X]] / stokes_I_ensemble[[X]]
    )

    U_norm_ensemble <- lapply(
      X = 1:length(df_list),
      FUN = function(X) stokes_U_ensemble[[X]] / stokes_I_ensemble[[X]]
    )
  } else {
    stokes_I_ensemble <- NULL
    stokes_Q_ensemble <- NULL
    stokes_U_ensemble <- NULL
    Q_norm_ensemble <- NULL
    U_norm_ensemble <- NULL
  }

  structure(list(
    n_segments = length(df_list),
    denoised_spectra = denoised_spectra,
    stokes_I_ensemble = stokes_I_ensemble,
    stokes_Q_ensemble = stokes_Q_ensemble,
    stokes_U_ensemble = stokes_U_ensemble,
    Q_norm_ensemble = Q_norm_ensemble,
    U_norm_ensemble = U_norm_ensemble,
    sure_args = sure_args,
    bootstrap_args = bootstrap_args,
    tf_obj
  ),
  class = c("polarized_spectrum", "list")
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
