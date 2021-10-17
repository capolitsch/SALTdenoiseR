# SALTdenoiseR R package

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

# Denoising spectra of Wolf-Rayet stars

## Summary
My preferred approach involves breaking the spectra 
into pieces when there are significant gaps due to masking, and then 
independently running the trend filtering analysis on each sub-spectrum. 
This has two benefits: 

1. It's faster. In some cases, much much faster. Which is helpful because
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
       
This segmenting procedure is determined by the `min.mask.width` parameter 
described below.

This particular script implements a "naive" bootstrap algorithm, where the
measurements of the `I` Stokes parameter are taken as fixed, and thus the
uncertainty in `flux = Q / I` is assumed to purely arise from `Q`. Since `Q`
is assumed to have a Gaussian distribution, so does flux in this setting.

Since `I` is not fixed, the main consequence is that the variability bands
produced in this script are more narrow than they should be. Use 
`ratio_bootstrap.R` for a better implementation that also
accounts for the uncertainty in `I`. As discussed, the next (and prob final)
improvement would come from account for the covariance between `Q` and `I`,
and any autocorrelations within `Q` and `I`.

Some parameters & variables that may need to be altered occasionally

* `min.mask.width`:    (numeric) Parameter that determines the segmentation of the
                       spectrum into smaller pieces divided by gaps (many 
                       consecutive masked spectral pixels). This is the minimum 
                       number of consecutive pixels that have to be masked to 
                       cause a break in the spectrum.
* `variability bands`: (boolean) Compute and plot variability bands? `FALSE` saves
                       significant computing time
* `save.plot.pdf`: (boolean) Save a plot of the spectrum with trend filtering 
                   results superposed as a pdf in the path_to_stokes_files
                   directory?
* `B`:             (integer) The number of bootstrap samples that are drawn 
                   in order to estimate the trend filtering variability bands. 
                   Larger is better, but is more computationally expensive. 
                   My trend filtering bootstrap implementation is set up for 
                   parallel computing for speedups when multiple cores are
                   available. `B = 100` typically suffices for a most 
                   scenarios, but it's good to crank it up a bit for a plot
                   that's going to be published, so you get nice smooth
                   bands.
* `obj_tol`:       (numeric) Maximum iterations allowed for the ADMM trend
                   filtering convex optimization (Ramdas and Tibshirani 2016). 
                   Increase this if the trend filtering estimate 
                   does not appear to have fully converged to a reasonable 
                   estimate of the spectral signal. An estimate that has not 
                   fully converged will typically increasingly diverge above 
                   or below the observed fluxes as you move from left to 
                   right.
* `max_iter`:      (integer) The tolerance used in the convex optimization 
                   stopping criterion; when the relative change in the 
                   objective function is less than this value, the algorithm 
                   terminates. Decrease this if the trend 
                   filtering estimate does not appear to have fully converged 
                   to a reasonable estimate of the signal. An estimate that 
                   has not fully converged will typically increasingly 
                   diverge above or below the observed fluxes as you move 
                   from left to right
