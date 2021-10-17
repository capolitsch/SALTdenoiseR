min.mask.width <- 20
variability.bands <- FALSE
B <- 100
obj_tol <- 1e-14
max_iter <- 1000
mc_cores <- parallel::detectCores()
data("SALT_spectrum")

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
