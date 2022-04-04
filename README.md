Install using `devtools::install_github("https://github.com/mca91/RWSIM/", subdir = "RWSIM")`.

You can also `git clone` and then build and install locally by running the shell script `buildpackage.sh`.

<br>

### Example for `DF_reg()`:

Run (augmented) Dickey-Fuller regression and obtain t-statistic for H_0: \rho = 0 (the coefficient on y_{t-1}). 

Note that I've added the option to drop lags in `1:p` (where `p` is the maximum lag) using the `remove_lags` argument which has no default so you need to set `remove_lags = 0` for all lags to be included.

```r
library(RWSIM)

set.seed(123)
y <- t(rnorm(100))
DF_Reg(y, p = 10, model = "c", remove_lags = 0)
```

```
[1] -2.525997
```

### Example for `DF_Reg_field()`

Like `DF_reg()` and thus also the same call signature but different output:

**3x1 array**

```
     [,1]      
[1,] residuals from ADF regression         n-p-1 x 1 matrix
[2,] estimate for coefficient rho          1 x 1     matrix
[3,] coefficients on the \Delta y_{t-j}    p x 1     matrix
```

```r
set.seed(123)
y <- t(rnorm(100))

DF <- DF_Reg_field(y, p = 10, model = "c", remove_lags = 0) 
DF
```

```
     [,1]      
[1,] numeric,89 
[2,] 0.8720768 
[3,] numeric,10
```

You may use `list` subsetting to extract entries:

```r
DF[[3]]
```

```
             [,1]
 [1,] -0.005815618
 [2,] -0.081509051
 [3,]  0.107343211
 [4,] -0.007601859
 [5,]  0.034466930
 [6,]  0.031101228
 [7,]  0.103760300
 [8,] -0.018097687
 [9,] -0.082248752
[10,] -0.131725638
```

### Example for `ARMA_sim()`:

```r
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
```


#      ar1       ma1       ma2 
# 0.6836856 0.3117469 0.1090300  
```

<br>

### Example for `S2_ar`:

In the exaple below we compute the long-run variance using an estimate of the spectral density at frequency zero. `S2_AR()` will include the first `k` lagged differences of `dat` in the auxiliary regression.

```r
S2_AR(
    dat = t(rnorm(100)),
    k = 5,
    model = "nc",
    remove_lags = 0
)
```

We may use a subset of these `k` differences by supplying an index vector for the lags to be **excluded** to `remove_lags`:

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
FLIBS   = -L/opt/R/arm64/gfortran/lib -mtune=native
F77     = /opt/R/arm64/gfortran/bin/gfortran
FC      = /opt/R/arm64/gfortran/bin/gfortran
```
