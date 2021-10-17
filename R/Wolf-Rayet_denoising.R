remotes::install_github("capolitsch/trendfiltering", auth_token = "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")
library(trendfiltering)
library(tidyverse)
library(magrittr)

min.mask.width <- 20
variability.bands <- FALSE
B <- 100
obj_tol <- 1e-14
max_iter <- 1000
mc_cores <- parallel::detectCores()

w <- seq(from = sci$axDat$crval[1], by = sci$axDat$cdelt[1], length = sci$axDat$len[1])

df <- tibble(wavelength = w,
             Q = sci$imDat[, , 2],
             Q.var = var$imDat[, , 2],
             Q.mask = bpm$imDat[, , 2],
             U = sci$imDat[, , 3],
             U.var = var$imDat[, , 3],
             U.mask = bpm$imDat[, , 3],
             I = sci$imDat[, , 1],
             I.var = var$imDat[, , 1],
             I.mask = bpm$imDat[, , 1])

df$mask <- df$Q.mask
df <- mask.intervals(df, min.mask.width)
n.segments <- df %>% pull(segment) %>% n_distinct

out <- vector(mode = "list", length = n.segments)
for (itr in 1:n.segments) {
  df.segment <- df %>% filter(segment == itr)
  tf.out <- trendfilter.interval(x = df.segment$wavelength,
                                 y = df.segment$Q,
                                 weights = 1 / df.segment$Q.var,
                                 obj_tol = obj_tol,
                                 max_iter = max_iter)
  out[[itr]] <- tf.out
}
