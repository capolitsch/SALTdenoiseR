#' Denoise a polarized spectrum via quadratic trend filtering
#'
#' \loadmathjax The `denoise_spectrum()` function uses quadratic trend
#' filtering, tuned by Stein's unbiased risk estimate, to optimally denoise a
#' polarized spectrum in each of the \mjseqn{I}, \mjseqn{Q}, \mjseqn{U} Stokes
#' parameters (sometimes alternately denoted \mjseqn{S_0}, \mjseqn{S_1},
#' \mjseqn{S_2}), which then also leads to estimates for the normalized
#' (unit-less) Stokes parameters \mjseqn{Q/I} and \mjseqn{U/I}. Setting
#' `compute_uncertainties = TRUE` generates a bootstrap ensemble of denoised
#' spectra for each Stokes parameter, which allows variability bands to be
#' computed for each denoised Stokes spectrum by then calling
#' [variability_bands()] on the `denoise_spectrum()` output. The default number
#' of bootstrap samples in each ensemble, `B = 100`, can be increased by
#' specifying a new value for `B` in the `bootstrap_args` list.
#'
#' @param wavelength Vector of wavelength measurements.
#' @param flux Spectropolarimetric measurements, passed as a 3-column tibble,
#' data frame, or matrix, with the columns corresponding to the \mjseqn{I},
#' \mjseqn{Q}, \mjseqn{U} Stokes parameters, respectively.
#' @param variances Measurement variances, in a tibble, data frame, or matrix
#' with dimensions matching those of `flux`.
#' @param masks Pixel masks, in a tibble, data frame, or matrix
#' with dimensions matching those of `flux`. Nonzero elements flag bad pixels.
#' @param break_at A free parameter that controls the segmentation of a
#' spectrum. More precisely, `break_at` is the minimum number of
#' consecutively-masked spectral pixels that will trigger a break in the
#' spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, the
#' resulting segments are examined to ensure that each is sufficiently long for
#' an independent denoising analysis. In particular, any segment that has less
#' than `min_pix_segment` unmasked spectral pixels is discarded. Defaults to
#' `min_pix_segment = 10`.
#' @param compute_uncertainties (Boolean) If `TRUE`, then bootstrap ensembles
#' are created for each denoised spectrum via
#' [`bootstrap_trendfilter()`][trendfiltering::bootstrap_trendfilter()], which
#' allows variability bands of any level to quickly be computed by calling
#' [variability_bands()] on the returned object (see examples).
#' @param mc_cores Multi-core computing using the
#' [`parallel`][`parallel::parallel-package`] package: The number of cores to
#' utilize. Defaults to the number of cores detected.
#' @param sure_args (Optional) A named list of arguments to be passed to
#' [`sure_trendfilter()`][trendfiltering::sure_trendfilter()]. The evaluation
#' grid defaults to the observed wavelength grid.
#' @param bootstrap_args (Optional) A named list of arguments to be passed to
#' [`bootstrap_trendfilter()`][trendfiltering::bootstrap_trendfilter()]. The
#' argument `bootstrap_algorithm = "parametric"` is fixed due to the
#' statistical design of SALT spectra, and we don't allow this argument to be
#' manually overridden. See the
#' [`bootstrap_trendfilter()`][trendfiltering::bootstrap_trendfilter()]
#' documentation and
#' [Politsch et al. (2020a)](
#' https://academic.oup.com/mnras/article/492/3/4005/5704413)
#' for details on why this particular bootstrap algorithm is appropriate for
#' SALT spectra.
#'
#' @return An object of class `'polarized_spectrum'`. This is a list with the
#' following elements:
#' \item{n_segments}{The number of segments the spectrum was broken into.}
#' \item{denoised_spectra}{A list of tibbles containing the denoised spectra.
#' The number of tibbles is equal to `n_segments` and each tibble has the
#' column set `c("wavelength","I","Q","U","Q_norm","U_norm")`, where
#' `Q_norm = Q/I` and `U_norm = U/I`.}
#' \item{I_ensemble}{If `compute_uncertainties = TRUE`, then this will
#' be the full bootstrap ensemble of denoised \mjseqn{I} spectra,
#' returned as an \mjseqn{n \times B} matrix, less any columns from post-hoc
#' pruning
#' (see [`bootstrap_trendfilter()`][trendfiltering::bootstrap_trendfilter()]),
#' which is then used by [variability_bands()] to compute
#' variability bands for the denoised \mjseqn{I} spectrum.
#' If `compute_uncertainties = FALSE`, this will return `NULL`.}
#' \item{Q_ensemble}{Same as above, but for the \mjseqn{Q} Stokes
#' parameter.}
#' \item{U_ensemble}{Same as above, but for the \mjseqn{U} Stokes
#' parameter.}
#' \item{Q_norm_ensemble}{Same as above, but for \mjseqn{Q/I}.}
#' \item{U_norm_ensemble}{Same as above, but for \mjseqn{U/I}.}
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
#' # Any SALT Observatory spectrum can be read into R from its FITS file
#' # using the "FITSio" R package, as below. Here, we've stored the `sci`,
#' # `var`, and `bpm` objects in an R data file within the package so we can
#' # simply source them using `data()`, as below.
#' \dontrun{
#' file.name <- "WR006_c1_12345678_stokes.fits"
#' sci <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 1)
#' var <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 2)
#' bpm <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 4)
#' }
#'
#' data(polarized_spectrum_WR_star)
#'
#' library(dplyr)
#'
#' wavelength <- seq(
#'   from = sci$axDat$crval[1],
#'   by = sci$axDat$cdelt[1],
#'   length = sci$axDat$len[1]
#' )
#'
#' flux <- as_tibble(sci$imDat)
#' variances <- as_tibble(var$imDat) %>% select(1:3)
#' masks <- as_tibble(bpm$imDat)
#'
#' \dontrun{
#' spec_denoised <- denoise_spectrum(
#'   wavelength,
#'   flux,
#'   variances,
#'   masks,
#'   compute_uncertainties = TRUE
#' )
#'
#' bands <- variability_bands(spec_denoised, param = "Q_norm", level = 0.95)
#' }
#' @importFrom trendfiltering sure_trendfilter bootstrap_trendfilter
#' @importFrom glmgen trendfilter trendfilter.control.list
#' @importFrom tidyr drop_na tibble as_tibble
#' @importFrom tibble as_tibble_col
#' @importFrom dplyr %>% arrange filter select n_distinct bind_cols
#' @importFrom magrittr %$%
#' @importFrom parallel mclapply detectCores
denoise_spectrum <- function(wavelength,
                             flux,
                             variances,
                             masks,
                             break_at = 10,
                             min_pix_segment = 10,
                             compute_uncertainties = FALSE,
                             mc_cores = parallel::detectCores(),
                             sure_args,
                             bootstrap_args) {
  stopifnot(ncol(flux) == 3)
  if (missing(variances)) {
    stop("Currently, variances must be passed for this denoising analysis.")
  }
  stopifnot(ncol(variances) == 3 & nrow(variances) == nrow(flux))
  if (missing(bootstrap_args)) {
    bootstrap_args <- list()
  } else {
    stopifnot(all(names(bootstrap_args) %in%
      names(formals(bootstrap_trendfilter))))
  }
  if (missing(sure_args)) {
    sure_args <- list()
  } else {
    stopifnot(all(names(sure_args) %in% names(formals(sure_trendfilter))))
  }

  if (missing(masks)) {
    masks <- matrix(rep(0, 3 * nrow(flux)), ncol = 3)
  } else {
    stopifnot(ncol(masks) == 3 & nrow(masks) == nrow(flux))
  }

  bootstrap_args$bootstrap_algorithm <- "parametric"
  if (compute_uncertainties) bootstrap_args$return_ensemble <- TRUE

  wavelength <- as_tibble_col(wavelength, column_name = "wavelength")
  flux <- as_tibble(flux) %>%
    rename_with(function(.cols) c("I", "Q", "U"))
  variances <- as_tibble(variances) %>%
    rename_with(function(.cols) c("I_vars", "Q_vars", "U_vars"))
  masks <- as_tibble(masks) %>%
    rename_with(function(.cols) c("I_mask", "Q_mask", "U_mask"))

  df_full <- bind_cols(
    wavelength,
    flux,
    variances,
    masks
  )

  df_list <- break_spectrum(df_full, break_at, min_pix_segment)

  tf_obj <- mclapply(
    X = 1:(3 * length(df_list)),
    parallel_sure_tf,
    df_list = df_list,
    sure_args = sure_args,
    mc.cores = mc_cores
  )

  if (compute_uncertainties) {
    tf_obj <- mclapply(
      1:(3 * length(df_list)),
      parallel_bootstrap_tf,
      sure_tf = tf_obj,
      bootstrap_args = bootstrap_args,
      mc.cores = mc_cores
    )
  }

  wavelength_eval <- lapply(
    X = 1:length(df_list),
    FUN = function(X) tf_obj[[X]][["x_eval"]]
  )

  I_denoised <- lapply(
    X = 1:length(df_list),
    FUN = function(X) tf_obj[[X]][["tf_estimate"]]
  )

  Q_denoised <- lapply(
    X = (length(df_list) + 1):(2 * length(df_list)),
    FUN = function(X) tf_obj[[X]][["tf_estimate"]]
  )

  U_denoised <- lapply(
    X = (2 * length(df_list) + 1):(3 * length(df_list)),
    FUN = function(X) tf_obj[[X]][["tf_estimate"]]
  )

  Q_norm_denoised <- lapply(
    X = 1:length(df_list),
    FUN = function(X) Q_denoised[[X]] / I_denoised[[X]]
  )

  U_norm_denoised <- lapply(
    X = 1:length(df_list),
    FUN = function(X) U_denoised[[X]] / I_denoised[[X]]
  )

  denoised_spectra <- lapply(
    X = 1:length(df_list),
    FUN = function(X) {
      tibble(
        wavelength = wavelength_eval[[X]],
        I = I_denoised[[X]],
        Q = Q_denoised[[X]],
        U = U_denoised[[X]],
        Q_norm = Q_norm_denoised[[X]],
        U_norm = U_norm_denoised[[X]]
      )
    }
  )

  if (compute_uncertainties) {
    I_ensemble <- lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]][["tf_bootstrap_ensemble"]]
    )

    Q_ensemble <- lapply(
      X = (length(df_list) + 1):(2 * length(df_list)),
      FUN = function(X) tf_obj[[X]][["tf_bootstrap_ensemble"]]
    )

    U_ensemble <- lapply(
      X = (2 * length(df_list) + 1):(3 * length(df_list)),
      FUN = function(X) tf_obj[[X]][["tf_bootstrap_ensemble"]]
    )

    Q_norm_ensemble <- lapply(
      X = 1:length(df_list),
      FUN = function(X) Q_ensemble[[X]] / I_ensemble[[X]]
    )

    U_norm_ensemble <- lapply(
      X = 1:length(df_list),
      FUN = function(X) U_ensemble[[X]] / I_ensemble[[X]]
    )
  } else {
    I_ensemble <- NULL
    Q_ensemble <- NULL
    U_ensemble <- NULL
    Q_norm_ensemble <- NULL
    U_norm_ensemble <- NULL
  }

  structure(list(
    n_segments = length(df_list),
    denoised_spectra = denoised_spectra,
    I_ensemble = I_ensemble,
    Q_ensemble = Q_ensemble,
    U_ensemble = U_ensemble,
    Q_norm_ensemble = Q_norm_ensemble,
    U_norm_ensemble = U_norm_ensemble,
    sure_args = sure_args,
    bootstrap_args = bootstrap_args,
    tf_obj = tf_obj
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

############################################

#' Quantify statistical uncertainty in a denoised spectrum via bootstrap
#' variability bands
#'
#' @param obj An object of class `"polarized_spectrum"` produced by
#' [`denoise_spectrum()`].
#' @param param The denoised spectrum to compute variability bands for. One of
#' `c("I","Q","U","Q_norm","U_norm")`.
#' @param level The level of the pointwise variability bands. Defaults to
#' `level = 0.95`.
#'
#' @return A list of tibbles, with length equal to `obj$n_segments` and each
#' tibble with column set
#' `c("wavelength","bootstrap_lower_band","bootstrap_upper_band")`.
#'
#' @export variability_bands
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
#' flux <- as_tibble(sci$imDat)
#' variances <- as_tibble(var$imDat) %>% select(1:3)
#' masks <- as_tibble(bpm$imDat)
#' \dontrun{
#' spec_denoised <- denoise_spectrum(
#'   wavelength,
#'   flux,
#'   variances,
#'   masks,
#'   compute_uncertainties = TRUE
#' )
#'
#' bands <- variability_bands(spec_denoised, param = "Q_norm", level = 0.95)
#' }
#' @importFrom dplyr case_when tibble
variability_bands <- function(obj, param, level = 0.95) {
  stopifnot(any(class(obj) == "polarized_spectrum"))
  stopifnot(param %in% c("I", "Q", "U", "Q_norm", "U_norm"))
  stopifnot(level > 0 & level < 1)

  ensemble <- case_when(
    param == "I" ~ obj$I_ensemble,
    param == "Q" ~ obj$Q_ensemble,
    param == "U" ~ obj$U_ensemble,
    param == "Q_norm" ~ obj$Q_norm_ensemble,
    param == "U_norm" ~ obj$U_norm_ensemble,
  )

  bootstrap_lower_band <- lapply(
    X = 1:length(ensemble),
    FUN = function(X) {
      apply(
        ensemble[[X]], 1,
        quantile,
        probs = (1 - level) / 2
      )
    }
  )

  bootstrap_upper_band <- lapply(
    X = 1:length(ensemble),
    FUN = function(X) {
      apply(
        ensemble[[X]], 1,
        quantile,
        probs = 1 - (1 - level) / 2
      )
    }
  )

  lapply(
    X = 1:length(ensemble),
    FUN = function(X) {
      tibble(
        wavelength = obj$denoised_spectra[[X]]$wavelength,
        bootstrap_lower_band = bootstrap_lower_band[[X]],
        bootstrap_upper_band = bootstrap_upper_band[[X]]
      )
    }
  )
}
