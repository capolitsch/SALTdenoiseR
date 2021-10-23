#' Break a polarized spectrum into segments when one or more wide masks exist
#'
#' Break a polarized spectrum into segments when one or more wide masks (i.e.
#' sequences of many consecutively-masked pixels) exist. The set of masked
#' pixels used to define the segmentation is taken to be the set of pixels that
#' are masked for *any* of the I, Q, or U Stokes parameters. Note that
#' [`denoise_polarized_spectrum()`] calls `break_spectrum()` internally. So
#' while `break_spectrum()` is an exported function, it does not need to be
#' called directly in most cases.
#'
#' @param df_full A tibble with minimal column set
#' `c("wavelength", "I_mask", "Q_mask", "U_mask")`. The values of the mask
#' columns are 0 (or FALSE) for unmasked pixels and 1 (or TRUE) for masked
#' pixels.
#' @param break_at The minimum number of consecutively-masked spectral pixels
#' that will trigger a break in the spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, the
#' resulting segments are examined to ensure that each is sufficiently long for
#' a non-trivial/well-defined denoising analysis. In particular, any segment
#' that has less than `min_pix_segment` unmasked spectral pixels is discarded.
#' Defaults to `min_pix_segment = 10`.
#' @return A list of tibbles that define the segmented spectrum. Each
#' tibble in the list has the same format as the input tibble `df_full`, minus
#' the mask columns. The output tibbles only retain the pixels that are unmasked
#' for every I, Q, U Stokes parameter, hence the removal of the mask columns.
#' @export break_spectrum
#' @examples
#' data(polarized_spectrum_WR_star)
#'
#' suppressMessages(library(dplyr))
#'
#' wavelength <- seq(
#'   from = sci$axDat$crval[1],
#'   by = sci$axDat$cdelt[1],
#'   length = sci$axDat$len[1]
#' ) %>%
#'   tibble::as_tibble_col(column_name = "wavelength")
#'
#' stokes <- as_tibble(sci$imDat) %>%
#'   rename_with(function(.cols) c("I", "Q", "U"))
#'
#' variances <- as_tibble(var$imDat) %>%
#'   select(1:3) %>%
#'   rename_with(function(.cols) c("I_vars", "Q_vars", "U_vars"))
#'
#' masks <- as_tibble(bpm$imDat) %>%
#'   rename_with(function(.cols) c("I_mask", "Q_mask", "U_mask"))
#'
#' df_full <- bind_cols(
#'   wavelength,
#'   stokes,
#'   variances,
#'   masks
#' )
#'
#' df_list <- break_spectrum(df_full)
#' @importFrom dplyr %>% filter pull
#' @importFrom matrixStats rowMaxs
#' @importFrom dplyr %>% arrange filter select
#' @importFrom magrittr %$%
#' @importFrom purrr discard
#' @importFrom data.table rleid
break_spectrum <- function(df_full, break_at = 10, min_pix_segment = 10) {
  stopifnot(
    all(c("wavelength", "I_mask", "Q_mask", "U_mask") %in% names(df_full))
  )

  df <- df_full %>%
    select(starts_with(c("wavelength", "I", "Q", "U"))) %>%
    mutate(mask = df_full %>% select(ends_with("mask")) %>%
      as.matrix() %>%
      rowMaxs()) %>%
    select(-ends_with("_mask")) %>%
    arrange(wavelength)

  # Drop leading and trailing masks, if present
  if (df$mask[1] != 0) {
    df <- df %>% filter(rleid(mask) != 1)
  }
  if (df$mask[nrow(df)] != 0) {
    df <- df %>% filter(rev(rleid(rev(mask))) != 1)
  }

  segments <- rleid(df$mask) %>%
    table() %>%
    as.double()

  changepoints <- c(0, cumsum(segments))

  breaks <- c(
    0,
    seq(2, length(segments), 2)[
      which(segments[seq(2, length(segments), 2)] >= break_at)
    ],
    length(changepoints)
  )

  lapply(
    X = 1:(length(breaks) - 1),
    FUN = function(X) {
      df_segment <- df[(changepoints[(breaks[X] + 1)] + 1):changepoints[breaks[X + 1]], ] %>%
        filter(mask == 0) %>%
        select(-mask)
      if (nrow(df_segment) >= min_pix_segment) {
        return(df_segment)
      } else {
        return(NULL)
      }
    }
  ) %>%
    discard(is.null)
}
