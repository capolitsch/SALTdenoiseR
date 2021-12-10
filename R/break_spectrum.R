#' Break a spectrum into segments when one or more wide masks exist
#'
#' Break a spectrum into segments when one or more wide masks (i.e.
#' sequences of many consecutively-masked pixels) exist. The spectrum can be
#' polarized or unpolarized. When the spectrum is polarized, the set of masked
#' spectral pixels used to define the spectrum's segmentation is
#' taken to be the union set of pixels that are masked for *any* of the
#' \mjseqn{IQU} Stokes parameters. Note that [`denoise_spectrum()`] calls
#' `break_spectrum()` internally. So while `break_spectrum()` is an exported
#' function, it typically does not need to be called directly.
#'
#' @param df A tibble that includes a column titled `"wavelength"` and one or
#' more columns with names ending in `"mask"`, e.g.
#' `c("I_mask", "Q_mask", "U_mask")`. The values of the mask columns should be
#' equal to 0 (or FALSE) for unmasked pixels and nonzero (or TRUE) for masked
#' pixels.
#' @param break_at The minimum number of consecutively-masked spectral pixels
#' that will trigger a break in the spectrum. Defaults to `break_at = 10`.
#' @param min_pix_segment After the segmentation procedure is complete, the
#' resulting segments are examined to ensure that each is sufficiently long for
#' a non-trivial/well-defined denoising analysis. In particular, any segment
#' that has fewer than `min_pix_segment` unmasked spectral pixels is discarded.
#' Defaults to `min_pix_segment = 10`.
#' @return A list of tibbles, with one for each segment of the segmented
#' spectrum. Masked rows are now dropped and each tibble in the list has the
#' same format as the input `df`, minus the mask columns.
#' @export break_spectrum
#'
#' @seealso [`denoise_spectrum()`]
#'
#' @examples
#' library(dplyr)
#' data(polarized_spectrum_WR_star)
#'
#' wavelength <- seq(
#'   from = sci$axDat$crval[1],
#'   by = sci$axDat$cdelt[1],
#'   length = sci$axDat$len[1]
#' ) %>%
#'   tibble::as_tibble_col(column_name = "wavelength")
#'
#' flux <- as_tibble(sci$imDat) %>%
#'   rename_with(function(.cols) c("I", "Q", "U"))
#'
#' variance <- as_tibble(var$imDat) %>%
#'   select(1:3) %>%
#'   rename_with(function(.cols) c("I_var", "Q_var", "U_var"))
#'
#' mask <- as_tibble(bpm$imDat) %>%
#'   rename_with(function(.cols) c("I_mask", "Q_mask", "U_mask"))
#'
#' df <- bind_cols(wavelength, flux, variance, mask)
#'
#' df_list <- break_spectrum(df)
#' @importFrom dplyr filter pull arrange select ends_with first last slice
#' @importFrom magrittr %$% %>% %<>%
#' @importFrom purrr discard
#' @importFrom data.table rleid
#' @importFrom matrixStats rowMaxs
break_spectrum <- function(df, break_at = 10L, min_pix_segment = 10L) {
  stopifnot(any(class(df) %in% c("tibble", "tbl", "data.frame", "data.table")))
  stopifnot("wavelength" %in% names(df))
  stopifnot(df %>% select(ends_with("mask")) %>% ncol() >= 1)
  stopifnot(break_at >= 0)
  stopifnot(min_pix_segment >= 5)

  df %<>% mutate(union_mask = df %>% select(ends_with("mask")) %>%
    as.matrix() %>%
    abs() %>%
    rowMaxs()) %>%
    arrange(wavelength)

  if (first(df$union_mask) != 0) df %<>% filter(rleid(union_mask) != 1)
  if (last(df$union_mask) != 0) df %<>% filter(rev(rleid(rev(union_mask))) != 1)

  segments <- rleid(df$union_mask) %>%
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
      df_segment <- df %>%
        slice((changepoints[(breaks[X]+1)]+1):changepoints[breaks[X+1]]) %>%
        filter(union_mask == 0) %>%
        select(-ends_with("mask"))
      if (nrow(df_segment) >= min_pix_segment) {
        return(df_segment)
      } else {
        return(NULL)
      }
    }
  ) %>%
    discard(is.null)
}
