#' Denoise SALT Observatory spectra via trend filtering
#'
#' @param min_mask_width (numeric) Parameter that determines the segmentation of
#' the spectrum into smaller pieces divided by gaps (many consecutive masked
#' spectral pixels). This is the minimum number of consecutive pixels that have
#' to be masked to cause a break in the spectrum.
#' @param variability_bands Compute and plot variability bands? `FALSE` saves
#' significant computing time
#' @param save_plot_pdf` save a plot of the spectrum with trend filtering
#' results superposed as a pdf in the path_to_stokes_files directory?
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
#' @importFrom dplyr %>% arrange filter select
#' @importFrom magrittr %$%
denoise <- function(min.mask.width = 20,
                    variability.bands = FALSE,
                    B = 100,
                    obj_tol = 1e-14,
                    max_iter = 1000,
                    mc_cores = parallel::detectCores()) {

  w <- seq(from = sci$axDat$crval[1], by = sci$axDat$cdelt[1], length = sci$axDat$len[1])

  df <- tibble(
    wavelength = w,
    Q = sci$imDat[, , 2],
    Q.var = var$imDat[, , 2],
    Q.mask = bpm$imDat[, , 2],
    U = sci$imDat[, , 3],
    U.var = var$imDat[, , 3],
    U.mask = bpm$imDat[, , 3],
    I = sci$imDat[, , 1],
    I.var = var$imDat[, , 1],
    I.mask = bpm$imDat[, , 1]
  )

  df$mask <- df$Q.mask
  df <- mask.intervals(df, min.mask.width)
  n.segments <- df %>%
    pull(segment) %>%
    n_distinct()

  out <- vector(mode = "list", length = n.segments)
  for (itr in 1:n.segments) {
    df.segment <- df %>% filter(segment == itr)
    tf.out <- trendfilter.interval(
      x = df.segment$wavelength,
      y = df.segment$Q,
      weights = 1 / df.segment$Q.var,
      obj_tol = obj_tol,
      max_iter = max_iter
    )
    out[[itr]] <- tf.out
  }
}


#' @importFrom dplyr %>% filter pull
mask.intervals <- function(df, min.mask.width = 20) {
  masked <- df %>%
    filter(mask == 1) %>%
    pull(wavelength)
  masked.intervals <- list()
  itr <- 1
  diffs <- sapply(X = 2:length(masked), FUN = function(X) masked[X] - masked[X - 1])
  while (max(diffs) >= min.mask.width) {
    ind <- min(which(diffs >= min.mask.width))
    masked.intervals[[itr]] <- masked[1:(ind)]
    masked <- setdiff(masked, masked[1:(ind)])
    itr <- itr + 1
    diffs <- sapply(X = 2:length(masked), FUN = function(X) masked[X] - masked[X - 1])
  }
  if (itr > 1) {
    inds <- unlist(lapply(X = 1:length(masked.intervals), FUN = function(X) length(masked.intervals[[X]]) > min.mask.width))
    if (length(masked.intervals[inds]) > 0) {
      masked.intervals <- lapply(X = 1:length(masked.intervals[inds]), FUN = function(X) masked.intervals[inds][[X]])
    } else {
      masked.intervals <- NULL
    }
  } else {
    masked.intervals <- NULL
  }

  if (length(masked.intervals) > 0) {
    inds <- c(
      df %>% pull(wavelength) %>% min() - 1,
      sapply(X = 1:length(masked.intervals), FUN = function(X) min(masked.intervals[[X]])),
      df %>% pull(wavelength) %>% max()
    )
  } else {
    inds <- c(df %>% pull(wavelength) %>% min() - 1, df %>% pull(wavelength) %>% max())
  }

  good.wavelength.intervals <- lapply(X = 1:(length(inds) - 1), FUN = function(X) (inds[X] + 1):(inds[X + 1] - 1))

  df$segment <- NA
  for (itr in 1:length(good.wavelength.intervals)) {
    i <- match(good.wavelength.intervals[[itr]], df$wavelength)
    df$segment[i] <- itr
  }

  df %>%
    filter(mask == 0) %>%
    select(-mask)
}

#' @importFrom trendfiltering sure_trendfilter bootstrap_trendfilter
#' @importFrom dplyr %>%
#' @importFrom magrittr %$%
trendfilter.interval <- function(x, y, weights,
                                 B = 100, alpha = 0.05, bootstrap.bands = T,
                                 max_iter = 5000, obj_tol = 1e-12, ...) {
  SURE.out <- sure_trendfilter(x, y, weights,
    optimization.params = list(
      max_iter = max_iter,
      obj_tol = obj_tol
    )
  )
  if (!bootstrap.bands) {
    return(SURE.out)
  } else {
    boot.out <- bootstrap_trendfilter(obj = SURE.out, B = B, alpha = alpha, ...)
    return(boot.out)
  }
}
