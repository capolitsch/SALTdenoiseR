#' Denoise a polarized spectrum via quadratic trend filtering
#'
#' \loadmathjax The `denoise_spectrum()` function uses quadratic trend
#' filtering, tuned by Stein's unbiased risk estimate, to optimally denoise a
#' polarized spectrum in each of the \mjseqn{I}, \mjseqn{Q}, \mjseqn{U} Stokes
#' parameters (sometimes denoted \mjseqn{S_0}, \mjseqn{S_1}, \mjseqn{S_2}),
#' which then also leads to estimates for the normalized Stokes parameters
#' \mjseqn{Q/I} and \mjseqn{U/I}. Setting `compute_uncertainties = TRUE`
#' generates a bootstrap ensemble of denoised spectra for each Stokes parameter,
#' which allows variability bands to be computed for each denoised spectrum by
#' then calling [`variability_bands()`] on the `'polarized_spectrum'` object
#' produced by the `denoise_spectrum()` function call.
#'
#' @param wavelength Vector of wavelength measurements.
#' @param flux Spectropolarimetric measurements, passed as a 3-column tibble,
#' data frame, or matrix, with the columns corresponding to the \mjseqn{I},
#' \mjseqn{Q}, \mjseqn{U} Stokes parameters, respectively.
#' @param variances Measurement variances, in a tibble, data frame, or matrix
#' with dimensions matching those of `flux`.
#' @param masks Pixel masks, in a tibble, data frame, or matrix
#' with dimensions matching those of `flux`. Nonzero elements flag bad pixels.
#' @param break_at The minimum number of consecutively-masked spectral pixels
#' that will trigger a break in the spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, the
#' resulting segments are examined to ensure that each is sufficiently long for
#' a non-trivial/well-defined denoising analysis. In particular, any segment
#' that has less than `min_pix_segment` unmasked spectral pixels is discarded.
#' Defaults to `min_pix_segment = 10`.
#' @param compute_uncertainties (Boolean) If `TRUE`, then bootstrap ensembles
#' are created for each denoised spectrum via
#' [`bootstrap_trendfilter()`][trendfiltering::bootstrap_trendfilter()].
#' The ensembles are stored in the function output so that variability bands
#' of any level can quickly be computed by calls to [`variability_bands()`]
#' without redundant overhead (see examples). Defaults to
#' `compute_uncertainties = FALSE`.
#' @param B If `compute_uncertainties = TRUE`, `B` is the number of bootstrap
#' samples in the ensembles used to estimate the uncertainty in the denoiser.
#' Defaults to `B = 100`.
#' @param mc_cores Multi-core computing using the
#' [`parallel`][parallel::parallel-package] package: The number of cores to
#' use. Defaults to the number of cores detected on the machine, minus 4.
#' @param validation_args (Optional) A named list of arguments to be passed to
#' [`sure_trendfilter()`][trendfiltering::sure_trendfilter()]. The evaluation
#' grid argument `x_eval` of
#' [`sure_trendfilter()`][trendfiltering::sure_trendfilter()] defaults to the
#' observed wavelength grid. The arguments `x`, `y`, `weights`, and `k` cannot
#' be overridden.
#'
#' @return An object of class `'polarized_spectrum'`. This is a list with the
#' following elements:
#' \describe{
#' \item{n_segments}{The number of segments the spectrum was broken into by
#' [`break_spectrum()`].}
#' \item{data}{The original data set, as a list of `n_segments` tibbles. Each
#' tibble contains the observed wavelengths (with the superset of masked
#' wavelengths removed), the Stokes \mjseqn{I}, \mjseqn{Q}, \mjseqn{U} flux
#' measurements, and the Stokes \mjseqn{I}, \mjseqn{Q}, \mjseqn{U} flux
#' measurement variances.}
#' \item{denoised_signals}{The set of denoised spectra, as a list of
#' `n_segments` tibbles. Each tibble contains the wavelength evaluation grid for
#' its respective spectrum segment and all of the denoised signal estimates:
#' \mjseqn{I}, \mjseqn{Q}, \mjseqn{U}, \mjseqn{Q/I}, \mjseqn{U/I}.}
#' \item{ensembles}{If `compute_uncertainties = TRUE`, a list of bootstrap
#' ensembles for the denoised \mjseqn{I}, \mjseqn{Q}, and \mjseqn{U} ensembles,
#' respectively. Each ensemble is returned as an \mjseqn{n \times B} matrix,
#' less any columns from post-hoc pruning
#' (see the
#' [`bootstrap_trendfilter()`][trendfiltering::bootstrap_trendfilter()]
#' documentation).
#' If `compute_uncertainties = FALSE`, this will return `NULL`.}
#' \item{I_analysis_summary}{Technical summary of the denoising analysis of the
#' \mjseqn{I} Stokes parameter.}
#' \item{Q_analysis_summary}{Same as above, but for the \mjseqn{Q}
#' Stokes parameter.}
#' \item{U_analysis_summary}{Same as above, but for the \mjseqn{U}
#' Stokes parameter.}
#' }
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
#' @seealso [`variability_bands()`], [`break_spectrum()`]
#'
#' @examples
#' # Any SALT Observatory spectrum can be read into R from its FITS file
#' # using the "FITSio" R package, as below. For convenience, here we've stored
#' # the `sci`, `var`, and `bpm` data structures for a Wolf-Rayet stellar
#' # spectrum in an R data file within this package, so we can source them using
#' # a simple call to `data()`.
#'
#' \dontrun{
#' install.packages("FITSio")
#' path_to_FITS_file <- "<your_local_path_to_the_FITS_file>"
#' file.name <- "WR006_c1_12345678_stokes.fits"
#' sci <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 1)
#' var <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 2)
#' bpm <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 4)
#' }
#'
#' ####
#'
#' data(polarized_spectrum_WR_star)
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
#' spec_denoised <- denoise_spectrum(
#'   wavelength,
#'   flux,
#'   variances,
#'   masks,
#'   compute_uncertainties = TRUE
#' )
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
                             break_at = 10L,
                             min_pix_segment = 10L,
                             compute_uncertainties = FALSE,
                             B = 100L,
                             mc_cores = parallel::detectCores() - 4,
                             validation_args) {
  stopifnot(ncol(flux) == 3)
  stopifnot(length(wavelength) == nrow(flux))

  if (missing(variances)) {
    stop("Currently, variances must be passed for this denoising analysis.")
  }

  stopifnot(ncol(variances) == 3 & nrow(variances) == nrow(flux))

  if (missing(masks)) {
    masks <- matrix(rep(0, 3 * nrow(flux)), ncol = 3)
  } else {
    stopifnot(ncol(masks) == 3 & nrow(masks) == nrow(flux))
  }

  if (missing(validation_args)) {
    validation_args <- list()
  } else {
    stopifnot(
      all(
        names(validation_args) %in%
          names(formals(sure_trendfilter))
        )
      )
  }

  wavelength <- wavelength %>% as_tibble_col(column_name = "wavelength")
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
    validation_args = validation_args,
    mc.cores = mc_cores
  )

  if (compute_uncertainties) {
    bootstrap_args <- list(
      bootstrap_algorithm = "parametric",
      B = B,
      mc_cores = mc_cores,
      return_ensemble = compute_uncertainties
    )

    tf_obj <- mclapply(
      1:(3 * length(df_list)),
      parallel_bootstrap_tf,
      sure_tf = tf_obj,
      bootstrap_args = bootstrap_args,
      mc.cores = 1
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

  data <- lapply(
    X = 1:length(df_list),
    FUN = function(X) {
      tibble(
        wavelength = df_list[[X]]$wavelength,
        I = df_list[[X]]$I,
        Q = df_list[[X]]$Q,
        U = df_list[[X]]$U,
        I_vars = df_list[[X]]$I_vars,
        Q_vars = df_list[[X]]$Q_vars,
        U_vars = df_list[[X]]$U_vars,
        I_fitted_values = tf_obj[[X]]$fitted_values,
        Q_fitted_values = tf_obj[[length(df_list) + X]]$fitted_values,
        U_fitted_values = tf_obj[[2 * length(df_list) + X]]$fitted_values,
        I_residuals = tf_obj[[X]]$residuals,
        Q_residuals = tf_obj[[length(df_list) + X]]$residuals,
        U_residuals = tf_obj[[2 * length(df_list) + X]]$residuals
      )
    }
  )

  denoised_signals <- lapply(
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

    ensembles <- list(
      I = I_ensemble,
      Q = Q_ensemble,
      U = U_ensemble
    )
  } else {
    ensembles <- NULL
  }

  I_summary <- structure(list(
    lambdas = tf_obj[[1]]$lambdas,
    edfs = tf_obj[[1]]$edfs,
    generalization_errors = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$generalization_errors
    ),
    lambda_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$lambda_min
    ) %>% unlist(),
    i_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$i_min
    ) %>% unlist(),
    edf_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$edf_min
    ) %>% unlist(),
    n_iter = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$n_iter
    ),
    training_errors = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$training_errors
    ),
    optimisms = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$optimisms
    ),
    admm_params = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$admm_params
    ),
    edf_boots = if (compute_uncertainties) {
      lapply(
        X = 1:length(df_list),
        FUN = function(X) tf_obj[[X]]$edf_boots
      )
    } else {
      NULL
    },
    n_iter_boots = if (compute_uncertainties) {
      lapply(
        X = 1:length(df_list),
        FUN = function(X) tf_obj[[X]]$n_iter_boots
      )
    } else {
      NULL
    },
    x_scale = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$x_scale
    ) %>% unlist(),
    y_scale = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$y_scale
    ) %>% unlist(),
    data_scaled = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$data_scaled
    )
  ),
  class = c("stokes_spectrum", "sure_tf", "list")
  )

  Q_summary <- structure(list(
    lambdas = tf_obj[[length(df_list) + 1]]$lambdas,
    edfs = tf_obj[[length(df_list) + 1]]$edfs,
    generalization_errors = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$generalization_errors
    ),
    lambda_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$lambda_min
    ) %>% unlist(),
    i_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$i_min
    ) %>% unlist(),
    edf_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$edf_min
    ) %>% unlist(),
    n_iter = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$n_iter
    ),
    training_errors = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$training_errors
    ),
    optimisms = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$optimisms
    ),
    admm_params = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$admm_params
    ),
    edf_boots = if (compute_uncertainties) {
      lapply(
        X = 1:length(df_list),
        FUN = function(X) tf_obj[[length(df_list) + X]]$edf_boots
      )
    } else {
      NULL
    },
    n_iter_boots = if (compute_uncertainties) {
      lapply(
        X = 1:length(df_list),
        FUN = function(X) tf_obj[[length(df_list) + X]]$n_iter_boots
      )
    } else {
      NULL
    },
    x_scale = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$x_scale
    ) %>% unlist(),
    y_scale = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$y_scale
    ) %>% unlist(),
    data_scaled = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$data_scaled
    )
  ),
  class = c("stokes_spectrum", "sure_tf", "list")
  )

  U_summary <- structure(list(
    lambdas = tf_obj[[length(df_list) + 1]]$lambdas,
    edfs = tf_obj[[length(df_list) + 1]]$edfs,
    generalization_errors = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$generalization_errors
    ),
    lambda_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$lambda_min
    ) %>% unlist(),
    i_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$i_min
    ) %>% unlist(),
    edf_min = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$edf_min
    ) %>% unlist(),
    n_iter = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$n_iter
    ),
    training_errors = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$training_errors
    ),
    optimisms = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$optimisms
    ),
    admm_params = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$admm_params
    ),
    edf_boots = if (compute_uncertainties) {
      lapply(
        X = 1:length(df_list),
        FUN = function(X) tf_obj[[2 * length(df_list) + X]]$edf_boots
      )
    } else {
      NULL
    },
    n_iter_boots = if (compute_uncertainties) {
      lapply(
        X = 1:length(df_list),
        FUN = function(X) tf_obj[[2 * length(df_list) + X]]$n_iter_boots
      )
    } else {
      NULL
    },
    x_scale = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$x_scale
    ) %>% unlist(),
    y_scale = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$y_scale
    ) %>% unlist(),
    data_scaled = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$data_scaled
    )
  ),
  class = c("stokes_spectrum", "sure_tf", "list")
  )

  structure(list(
    n_segments = length(df_list),
    data = data,
    denoised_signals = denoised_signals,
    ensembles = ensembles,
    I_analysis_summary = I_summary,
    Q_analysis_summary = Q_summary,
    U_analysis_summary = U_summary
  ),
  class = c("polarized_spectrum", "list")
  )
}


#' @noRd
parallel_sure_tf <- function(X, df_list, validation_args) {
  if (X %in% 1:length(df_list)) {
    args <- c(
      list(
        x = df_list[[X]]$wavelength,
        y = df_list[[X]]$I,
        weights = 1 / df_list[[X]]$I_vars,
        x_eval = df_list[[X]]$wavelength
      ),
      validation_args
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
      validation_args
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
      validation_args
    )
  }

  do.call(sure_trendfilter, args)
}


#' @noRd
parallel_bootstrap_tf <- function(X, sure_tf, bootstrap_args) {
  args <- c(list(obj = sure_tf[[X]]), bootstrap_args)
  do.call(bootstrap_trendfilter, args)
}
