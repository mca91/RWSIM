Install using `devtools::install_github("https://github.com/mca91/RWSIM/", subdir = "RWSIM")`.

You can also `git clone` and then build and install locally by running the shell script `buildpackage.sh`.

## Example for `ARMA_sim()`:

```r
library(RWSIM)

# AR(1): rho = .7
ARMA_sim(
  ar = .7,              # coefficients in AR polynomial (vector)
  ma_coefs = 0,         # coefficients in MA polynomial (vector)
  innovs = rnorm(100)   # innovations
 )
  
# MA(2): theta_1 = .2, theta_2 = .1
ARMA_sim(
  ar = 0,              
  ma_coefs = c(.2, .1),        
  innovs = rnorm(100)   
 )
```

For processes with higher persistence we may want use some burn-in samples:

```r
# ARMA(1, 2): rho = .7, theta_1 = .3, theta_2 = .1
ARMA_sim(
  ar = .7, 
  ma_coefs = c(.3, .1), 
  innovs = rnorm(300)
)[-c(1:150)]
```

Let's check using simulation.

```r
library(tidyverse)

set.seed(123)

est <- map(1:1000,
    ~ tseries::arma(
        x = ARMA_sim(ar = .7, ma_coefs = c(.3, .1), innovs = rnorm(300))[-c(1:150)], 
        order = c(1, 2), 
        include.intercept = F, 
      )$coef
     )

est %>% 
  reduce(rbind) %>% 
  colMeans()
  
#      ar1       ma1       ma2 
# 0.6836856 0.3117469 0.1090300  
```

## Example for `S2_ar`:

In the exaple below we compute the long-run variance using an estimate of the spectral density at frequency zero. `S2_AR()` will include the first `k` lagged differences of `dat` in the auxiliary regression.

```r
S2_AR(
    dat = t(rnorm(100)),
    k = 5,
    model = "nc",
    remove_lags = 0
)
```

We may use a subset of these `k` differences by suppliying an index vector for the lags to be excluded to `remove_lags`:

```r
S2_AR(
    dat = t(rnorm(100)),
    k = 5,
    model = "nc",
    remove_lags = c(1, 3, 4)
)
```

**Note to self:**

I had some trouble getting this to run on Apple Silicon Mac with the homebrew version of R. For now it runs fine using gfortran as shipped with the CRAN binaries for R. Edited Makevars accordingly (`file.edit("~/.R/Makevars")`):

```
FLIBS   = -L/opt/R/arm64/gfortran/lib
F77     = /opt/R/arm64/gfortran/bin/gfortran
FC      = /opt/R/arm64/gfortran/bin/gfortran
```
