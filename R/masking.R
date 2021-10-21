#' Break a spectrum into segments when chunks of many consecutively-masked
#' pixels are present
#'
#' @param df_full A tibble with minimal column set:
#' `c("wavelength","I","Q","U","I_mask","Q_mask","U_mask")`, and probably also
#' measurement variances for each of the Stokes parameters.
#' @param break_at A free parameter that controls the segmentation of a
#' spectrum. More precisely, `break_at` is the minimum number of
#' consecutively-masked spectral pixels that will trigger a break in the
#' spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, it is
#' advisable to examine the resulting segments to ensure that each is
#' sufficiently long for its own denoising analysis. In particular, we discard
#' any segments that have less than `min_pix_segment` unmasked spectral pixels.
#' Defaults to `min_pix_segment = 10`.
#' @return A list of tibbles corresponding to the segmented spectrum, as defined
#' by the above procedure and input argument choices. Each tibble in the list
#' has the same format as the input tibble `df_full`, less the mask columns.
#' The output tibbles only retain unmasked pixels, hence the removal of the
#' mask columns.
#' @export break_spectrum

#' @importFrom dplyr %>% filter pull
#' @importFrom matrixStats rowMaxs
#' @importFrom dplyr %>% arrange filter select
#' @importFrom magrittr %$%
#' @importFrom purrr discard
#' @importFrom data.table rleid
break_spectrum <- function(df_full, break_at = 10, min_pix_segment = 10) {
  df <- df_full %>%
    select(starts_with(c("wavelength", "I", "Q", "U"))) %>%
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

  segments <- rleid(df$mask) %>%
    table() %>%
    as.double()

  changepoints <- c(0, cumsum(segments))

  breaks <- c(
    0,
    seq(2, length(segments), 2)[
      which(segments[seq(2, length(segments), 2)] >= spectrum_break)
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
