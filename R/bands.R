#' Quantify the statistical uncertainty in a denoised spectrum by computing
#' bootstrap variability bands
#'
#' See the parametric bootstrap algorithm in [Politsch et al. (2020a)](
#' https://academic.oup.com/mnras/article/492/3/4005/5704413) for details.
#'
#' @param obj An object of class `"polarized_spectrum"`, produced by
#' [`denoise_spectrum()`].
#' @param param A string specifying which spectrum
#' to compute variability bands for. One of `c("I","Q","U","Q_norm","U_norm")`.
#' @param level The level of the pointwise variability bands. Defaults to
#' `level = 0.95`.
#'
#' @return A list of `obj$n_segments` tibbles, each with the column set
#' `c("wavelength","bootstrap_lower_band","bootstrap_upper_band")`.
#'
#' @export bands
#'
#' @seealso [`denoise_spectrum()`]
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
#'
#' spec_bands <- bands(spec_denoised, param = "Q_norm")
#' @importFrom dplyr case_when tibble
bands <- function(obj, param, level = 0.95) {
  stopifnot(any(class(obj) == "polarized_spectrum"))
  stopifnot(param %in% c("I", "Q", "U", "Q_norm", "U_norm"))
  stopifnot(level > 0 & level < 1)

  if (is.null(obj$ensembles)) {
    stop(paste(
      "obj does not have the bootstrap ensembles necessary to compute\n",
      "variability bands. Re-run denoise_spectrum() with\n",
      "`compute_uncertainties = TRUE`."
    ))
  }

  if (param %in% c("I", "Q", "U")) {
    ensemble <- case_when(
      param == "I" ~ obj$ensembles$I,
      param == "Q" ~ obj$ensembles$Q,
      param == "U" ~ obj$ensembles$U
    )
  } else {
    if (param == "Q_norm") {
      ensemble <- lapply(
        X = 1:obj$n_segments,
        FUN = function(X) obj$ensembles$Q[[X]] / obj$ensembles$I[[X]]
      )
    }

    if (param == "U_norm") {
      ensemble <- lapply(
        X = 1:obj$n_segments,
        FUN = function(X) obj$ensembles$U[[X]] / obj$ensembles$I[[X]]
      )
    }
  }

  bootstrap_lower_band <- lapply(
    X = 1:length(ensemble),
    FUN = function(X) {
      apply(
        ensemble[[X]],
        1,
        quantile,
        probs = (1 - level) / 2
      )
    }
  )

  bootstrap_upper_band <- lapply(
    X = 1:length(ensemble),
    FUN = function(X) {
      apply(
        ensemble[[X]],
        1,
        quantile,
        probs = 1 - (1 - level) / 2
      )
    }
  )

  lapply(
    X = 1:length(ensemble),
    FUN = function(X) {
      tibble(
        wavelength = obj$denoised_signals[[X]]$wavelength,
        bootstrap_lower_band = bootstrap_lower_band[[X]],
        bootstrap_upper_band = bootstrap_upper_band[[X]]
      )
    }
  )
}
