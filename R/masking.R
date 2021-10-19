#' Function(s) for handling pixel masks
#'
#' @param df Input data frame
#' @param min_mask_width Parameter that controls the segmentation of a spectrum
#' into ``sub-spectra'', which are then denoised independently. More precisely,
#' `min_mask_width` is the minimum number of consecutive spectral pixels that
#' have to be masked to cause a break in the spectrum.
#'
#' @return
#'
#' @noRd

#' @importFrom dplyr %>% filter pull
mask_intervals <- function(df, min_mask_width = 20) {
  masked <- df %>%
    filter(!!mask) %>%
    pull(wavelength)

  masked_intervals <- list()
  itr <- 1
  diffs <- sapply(
    X = 2:length(masked),
    FUN = function(X) masked[X] - masked[X - 1]
  )

  while (max(diffs) >= min_mask_width) {
    ind <- min(which(diffs >= min_mask_width))
    masked_intervals[[itr]] <- masked[1:(ind)]
    masked <- setdiff(masked, masked[1:(ind)])
    itr <- itr + 1
    diffs <- sapply(
      X = 2:length(masked),
      FUN = function(X) masked[X] - masked[X - 1]
    )
  }

  if (itr > 1) {
    inds <- lapply(
      X = 1:length(masked_intervals),
      FUN = function(X) length(masked_intervals[[X]]) > min_mask_width
    ) %>% unlist()

    if (length(masked_intervals[inds]) > 0) {
      masked_intervals <- lapply(
        X = 1:length(masked_intervals[inds]),
        FUN = function(X) masked_intervals[inds][[X]]
      )
    } else {
      masked_intervals <- NULL
    }
  } else {
    masked_intervals <- NULL
  }

  if (length(masked_intervals) > 0) {
    inds <- c(
      df %>% pull(wavelength) %>% min() - 1,
      sapply(
        X = 1:length(masked_intervals),
        FUN = function(X) min(masked_intervals[[X]])
      ),
      df %>% pull(wavelength) %>% max()
    )
  } else {
    inds <- c(
      df %>% pull(wavelength) %>% min() - 1,
      df %>% pull(wavelength) %>% max()
    )
  }

  good_wavelength_intervals <- lapply(
    X = 1:(length(inds) - 1),
    FUN = function(X) (inds[X] + 1):(inds[X + 1] - 1)
  )

  df$segment <- NA
  for (itr in 1:length(good_wavelength_intervals)) {
    i <- match(good_wavelength_intervals[[itr]], df$wavelength)
    df$segment[i] <- itr
  }

  df %>%
    filter(mask == 0) %>%
    select(-mask)
}
