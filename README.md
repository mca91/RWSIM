Install using `devtools::install_github("https://github.com/mca91/RWSIM/", subdir = "RWSIM")`.

You can also `git clone` and then build and install locally by running the shell script `buildpackage.sh`.

Example for `ARMA_sim()`:


```r
library(RWSIM)

# AR(1): rho = .7
ARMA_sim(
  ar = .7,              # coefficients in AR polynomial (vector)
  ma_coefs = 0,        # coefficients in MA polynomial (vector)
  innovs = rnorm(100)   # innovations
 )
  
# MA(2): theta_1 = .2, theta_2 = .1
ARMA_sim(
  ar = 0,              
  ma_coefs = (.2, .1),        
  innovs = rnorm(100)   
 )
```

For processes with higher persistence we may want use some burn-in samples:

```r
# ARMA(1, 2): rho = .7, theta_1 = .3, theta_2 = .1
ARMA_sim(
  ar = .7, 
  ma_coefs = c(0.3, .1), 
  innovs = rnorm(300)
)[-c(1:150)]
```

Let's check using simulation.

```r
library(tidyverse)

est <- map(1:1000,
    ~ tseries::arma(
        x = ARMA_sim(ar = .7, ma_coefs = c(0.3, .1), innovs = rnorm(300))[-c(1:150)], 
        order = c(1, 2), 
        include.intercept = F, 
      )$coef
     )

est %>% 
  reduce(rbind) %>% 
  colMeans()
```

**Note to self:**

I had some trouble getting this to run on Apple Silicon Mac with the homebrew version of R. For now it runs fine using gfortran as shipped with the CRAN binaries for R. Edited Makevars accordingly (`file.edit("~/.R/Makevars")`):

```
FLIBS   = -L/opt/R/arm64/gfortran/lib
F77     = /opt/R/arm64/gfortran/bin/gfortran
FC      = /opt/R/arm64/gfortran/bin/gfortran
```
