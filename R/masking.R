#' Function(s) for handling pixel masks
#'
#' @param df A tibble with minimal column set:
#' `c("wavelength","I","Q","U")`, and probably also weights and masks for each
#' of the Stokes parameters.
#' @param spectrum_break Parameter that controls the segmentation of a spectrum
#' into ``sub-spectra'', which are then denoised independently. More precisely,
#' `spectrum_break` is the minimum number of consecutive spectral pixels that
#' have to be masked to cause a break in the spectrum.
#'
#' @return Tibble
#' @noRd

#' @importFrom dplyr %>% filter pull
#' @importFrom matrixStats rowMaxs
#' @importFrom dplyr %>% arrange filter select
#' @importFrom magrittr %$%
#' @importFrom data.table rleid
mask_intervals <- function(df_full, stokes_params = c("I", "Q", "U"),
                           spectrum_break = 20, min_segment_length = 10) {
  df <- df_full %>%
    select(starts_with(c("wavelength", stokes_params))) %>%
    mutate(mask = df_full %>% select(ends_with("mask")) %>%
      as.matrix() %>%
      rowMaxs()) %>%
    select(-ends_with("_mask"))

  # Drop leading and trailing masks, if present
  if (df$mask[1] != 0) {
    df <- df %>% filter(rleid(mask) != 1)
  }
  if (df$mask[nrow(df)] != 0) {
    df <- df %>% filter(rev(rleid(rev(mask))) != 1)
  }

  segments <- df %>%
    rleid(mask) %>%
    table() %>%
    as.double()
  breaks <- which(which(segments >= spectrum_break) %% 2 == 0)
}
