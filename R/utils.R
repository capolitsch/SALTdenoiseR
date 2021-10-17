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
