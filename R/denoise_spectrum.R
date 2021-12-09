#' Denoise a polarized spectrum via trend filtering
#'
#' \loadmathjax The `denoise_spectrum()` function uses quadratic trend
#' filtering, tuned by Stein's unbiased risk estimate, to optimally denoise each
#' of the \mjseqn{IQU} Stokes parameters of a polarized spectrum, in turn also
#' leading to denoised estimates for the normalized Stokes parameters
#' \mjseqn{Q/I} and \mjseqn{U/I}. Setting `compute_uncertainties = TRUE`
#' generates a bootstrap ensemble of denoised spectra for each Stokes parameter,
#' which allows variability bands to be computed for each denoised spectrum by
#' then calling [`bands()`] on the `denoise_spectrum()` output.
#'
#' @param wavelength Vector of wavelength measurements.
#' @param flux Spectropolarimetric measurements, passed as a 3-column tibble,
#' data frame, or matrix, with the columns corresponding to the \mjseqn{IQU}
#' Stokes parameters, respectively.
#' @param variances Measurement variances, in a tibble, data frame, or matrix
#' with dimensions matching those of `flux`.
#' @param mask Pixel masks, in a tibble, data frame, or matrix with dimensions
#' matching those of `flux`. Nonzero elements flag bad pixels.
#' @param break_at The minimum number of consecutively-masked spectral pixels
#' that will trigger a break in the spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, the
#' resulting segments are examined to ensure that each is sufficiently long for
#' a non-trivial/well-defined denoising analysis. In particular, any segment
#' that has fewer than `min_pix_segment` unmasked spectral pixels is discarded.
#' Defaults to `min_pix_segment = 10`.
#' @param compute_uncertainties (Boolean) If `TRUE`, then bootstrap ensembles
#' are created for each denoised spectrum via
#' [trendfiltering::bootstrap_trendfilter()].
#' The ensembles are stored in the function output so that variability bands
#' of any level can quickly be computed by calls to [`bands()`] without
#' redundant overhead calculations (see examples). Defaults to
#' `compute_uncertainties = FALSE`.
#' @param mc_cores Multi-core computing using the
#' [`parallel`][parallel::parallel-package] package: The number of cores to
#' use. Defaults to the number of cores detected on the machine, minus 4.
#' @param ... (Optional) Named arguments to be passed to
#' [trendfiltering::sure_trendfilter()] or
#' [trendfiltering::bootstrap_trendfilter()]. The
#' arguments `x`, `y`, `weights`, `k`, and `algorithm` cannot be overridden.
#'
#' @return An object of class `'polarized_spectrum'`. This is a list with the
#' following elements:
#' \describe{
#' \item{n_segments}{Number of segments the spectrum was broken into by
#' [`break_spectrum()`].}
#' \item{data}{The original data set, as a list of `n_segments` tibbles. Each
#' tibble contains the observed wavelengths (with the union set of masked
#' wavelengths removed), the Stokes \mjseqn{IQU} flux measurements, and the
#' measurement variances.}
#' \item{denoised_signals}{The set of denoised spectra, as a list of
#' `n_segments` tibbles. Each tibble contains the wavelength evaluation grid for
#' its respective segment and all of the denoised signal estimates:
#' \mjseqn{I}, \mjseqn{Q}, \mjseqn{U}, \mjseqn{Q/I}, \mjseqn{U/I}.}
#' \item{ensembles}{If `compute_uncertainties = TRUE`, a list of bootstrap
#' ensembles for the denoised \mjseqn{I}, \mjseqn{Q}, \mjseqn{U} parameters,
#' respectively. Each ensemble is returned as an \mjseqn{n \times B} matrix.
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
#' @seealso [`bands()`], [`break_spectrum()`]
#'
#' @examples
#' # Any SALT Observatory spectrum can be read into R from its FITS file
#' # using the "FITSio" R package, as below. In the example below, we will
#' # analyze the polarized spectrum of a SALT Wolf-Rayet star. For convenience,
#' # here we've stored the header and data units (HDU) from the FITS file into
#' # three variables (`sci`, `var`, and `bpm`) that are stored in an R data file
#' # within this package, so we can source them using a simple call to `data()`.
#'
#' \dontrun{
#' install.packages("FITSio")
#' path_to_FITS_file <- "<your_local_path_to_FITS_file>"
#' file.name <- "WR006_c1_12345678_stokes.fits"
#' sci <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 1)
#' var <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 2)
#' bpm <- FITSio::readFITS(paste0(path_to_FITS_file, file.name), hdu = 4)}
#'
#' ####
#'
#' library(dplyr)
#' data(polarized_spectrum_WR_star)
#'
#' wavelength <- seq(
#'   from = sci$axDat$crval[1],
#'   by = sci$axDat$cdelt[1],
#'   length = sci$axDat$len[1]
#' )
#'
#' flux <- as_tibble(sci$imDat)
#' variance <- as_tibble(var$imDat) %>% select(1:3)
#' mask <- as_tibble(bpm$imDat)
#'
#' spec_denoised <- denoise_spectrum(
#'   wavelength,
#'   flux,
#'   variance,
#'   mask,
#'   compute_uncertainties = TRUE
#' )
#' @importFrom trendfiltering sure_trendfilter bootstrap_trendfilter
#' @importFrom dplyr rename_with arrange filter select n_distinct bind_cols
#' @importFrom magrittr %>% %$% %<>%
#' @importFrom tidyr drop_na tibble as_tibble
#' @importFrom tibble as_tibble_col
#' @importFrom parallel mclapply detectCores
denoise_spectrum <- function(wavelength,
                             flux,
                             variance,
                             mask,
                             wavelength_eval = wavelength,
                             compute_uncertainties = FALSE,
                             break_at = 10L,
                             min_pix_segment = 10L,
                             mc_cores = parallel::detectCores() - 4,
                             ...) {
  stopifnot(ncol(flux) == 3 & nrow(flux) == length(wavelength))
  stopifnot(ncol(variance) == 3 & nrow(variance) == nrow(flux))

  if (missing(mask)) {
    mask <- matrix(rep_len(0, 3 * nrow(flux)), ncol = 3)
  } else {
    stopifnot(ncol(mask) == 3 & nrow(mask) == nrow(flux))
  }

  extra_args <- list(...)

  wavelength <- wavelength %>% as_tibble_col(column_name = "wavelength")
  flux <- as_tibble(flux) %>%
    rename_with(function(.cols) c("I", "Q", "U"))
  variances <- as_tibble(variances) %>%
    rename_with(function(.cols) c("I_vars", "Q_vars", "U_vars"))
  mask <- as_tibble(mask) %>%
    rename_with(function(.cols) c("I_mask", "Q_mask", "U_mask"))

  df_list <- bind_cols(wavelength, flux, variance, mask) %>%
    break_spectrum(break_at, min_pix_segment)

  tf_obj <- mclapply(
    X = 1:(3 * length(df_list)),
    parallel_sure_tf,
    df_list = df_list,
    extra_args = extra_args,
    mc.cores = mc_cores
  )

  if (compute_uncertainties) {
    bootstrap_args <- list(
      algorithm = "parametric",
      B = B,
      x_eval = wavelength_eval,
      mc_cores = mc_cores
    )

    tf_obj <- mclapply(
      1:(3 * length(df_list)),
      parallel_bootstrap_tf,
      sure_tf = tf_obj,
      bootstrap_args = bootstrap_args,
      mc.cores = 1
    )
  }

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
    lambda = tf_obj[[1]]$lambda,
    edf = tf_obj[[1]]$edf,
    error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$error
    ),
    se_error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$se_error
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
    training_error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$training_error
    ),
    optimism = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$optimism
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
    lambda = tf_obj[[length(df_list) + 1]]$lambda,
    edf = tf_obj[[length(df_list) + 1]]$edf,
    error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$error
    ),
    se_error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$se_error
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
    training_error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$training_error
    ),
    optimism = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[length(df_list) + X]]$optimism
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
    lambda = tf_obj[[length(df_list) + 1]]$lambda,
    edf = tf_obj[[length(df_list) + 1]]$edf,
    error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$error
    ),
    se_error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[X]]$se_error
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
    training_error = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$training_error
    ),
    optimism = lapply(
      X = 1:length(df_list),
      FUN = function(X) tf_obj[[2 * length(df_list) + X]]$optimism
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
parallel_sure_tf <- function(X, df_list, extra_args) {
  if (X %in% 1:length(df_list)) {
    args <- c(
      list(
        x = df_list[[X]]$wavelength,
        y = df_list[[X]]$I,
        weights = 1 / df_list[[X]]$I_vars
      ),
      extra_args
    )
  }

  if (X %in% (length(df_list) + 1):(2 * length(df_list))) {
    args <- c(
      list(
        x = df_list[[X - 3]]$wavelength,
        y = df_list[[X - 3]]$Q,
        weights = 1 / df_list[[X - 3]]$Q_vars
      ),
      extra_args
    )
  }

  if (X %in% (2 * length(df_list) + 1):(3 * length(df_list))) {
    args <- c(
      list(
        x = df_list[[X - 6]]$wavelength,
        y = df_list[[X - 6]]$U,
        weights = 1 / df_list[[X - 6]]$U_vars
      ),
      extra_args
    )
  }

  do.call(sure_trendfilter, args)
}


#' @noRd
parallel_bootstrap_tf <- function(X, sure_tf, bootstrap_args) {
  args <- c(list(obj = sure_tf[[X]]), bootstrap_args)
  do.call(bootstrap_trendfilter, args)
}
