# SALTdenoiseR package

## Installation
``` r
install.packages("remotes")
remotes::install_github("capolitsch/SALTdenoiseR")
```
## Key references:

1. Politsch et al. (2020a). Trend Filtering – I. A modern statistical tool for time-domain astronomy and Astronomical Spectroscopy. 
*Monthly Notices of the Royal Astronomical Society*, 492(3), p. 4005-4018. [[Publisher](https://academic.oup.com/mnras/article/492/3/4005/5704413)] [[arXiv](https://arxiv.org/abs/1908.07151)]

2. Politsch et al. (2020b). Trend Filtering – II. Denoising astronomical signals with varying degrees of smoothness. 
*Monthly Notices of the Royal Astronomical Society*, 492(3), p. 4019-4032. [[Publisher](https://academic.oup.com/mnras/article/492/3/4019/5704414)] [[arXiv](https://arxiv.org/abs/2001.03552)]


## Denoising SALT spectra

### Summary

My preferred approach to denoising these spectra involves breaking the spectra 
into pieces when there are significant gaps due to masking, and then 
independently running a trend filtering analysis on each segment of a broken
spectrum. This has two benefits: 

1. It's faster. In some cases, *much* faster. Which is helpful because
my bootstrap implementation can sometimes be slow, depending on the
setting and your computing specs. I can probably make the bootstrap
more efficient if speed becomes a problem.
2. It allows for even greater local adaptivity of the trend filtering
estimator, because each sub-spectrum gets it's own optimized
hyperparameter. This is helpful when the signal-to-noise ratio may
have non-trivial variation across the sub-spectra.
3. Lastly, and really most importantly, we shouldn't be estimating
the spectrum -- and certainly not providing uncertainties -- in
windows where a non-trivial number of spectral pixels are
consecutively masked.
       
This segmenting procedure is determined by the `min_mask_width` argument
of [`denoise_spectrum()`](https://capolitsch.github.io/SALTdenoiseR/reference/denoise_spectrum.html).
