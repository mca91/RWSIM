### How to install

Install using `devtools::install_github("https://github.com/mca91/RWSIM/", subdir = "RWSIM")`.

You can also `git clone` and then build locally and install by running the shell script `buildpackage.sh`.

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

This is equivalent to what `urca::ur.df()` yields

```r
urca::ur.df(t(y), type = "drift", lags = 10)
```

```
############################################################### 
# Augmented Dickey-Fuller Test Unit Root / Cointegration Test # 
############################################################### 

The value of the test statistic is: -2.526 3.2035 
```

(the first value is the DF t-ratio for RW with drift, the second one is the test statistic for the joint hypothesis of an RW model _without_ drift) but we are ~ 60 times faster and more memory efficient :-).

```r
bench::mark(
    urca = urca::ur.df(t(y), type = "drift", lags = 10),
    DF_Reg(y, model = "c", p = 10, remove_lags = 0),
    check = F, relative = T
)
```

```
  expression                                        min median `itr/sec` mem_alloc 
  <bch:expr>                                      <dbl>  <dbl>     <dbl>     <dbl>    
1 urca                                             54.5   56.1       1        65.0 
2 DF_Reg(y, model = "c", p = 10, remove_lags = 0)   1      1        56.8       1 
```

<br>

### Example for `DF_Reg_field()`

Same call signature like `DF_reg()` but different output:

**5x1 array**

```
     [,1]      
[1,] residuals from ADF regression                      n-p-1 x 1  matrix
[2,] estimate for coefficient rho                       1 x 1      matrix
[3,] coefficient estimates on the \Delta y_{t-j}        p x 1      matrix
[4,] coefficient estimates on deterministics            d x 1      matrix
[5,] first differences of dependent variable            n-p-1 x 1  matrix
```

NB: for `model = "nc"` `DF_Reg_field(...)[[4]]` will return an empty (0 x 0) matrix due to way the field is initialised by Armadillo / translated to R.

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
[4,] 0.08116683
[5,] numeric,89
```

You may use `list` subsetting to extract entries. Let's get coefficients on lags of the regressand:

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

Now with modified lag structure:

```r
DF_Reg_field(y, p = 10, model = "c", remove_lags = c(2, 4, 6, 8))[[3]]
```

```
            [,1]
[1,]  0.07951090
[2,]  0.15649299
[3,]  0.01993999
[4,]  0.09818481
[5,] -0.06412958
[6,] -0.12236955
```

<br>

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

```
#      ar1       ma1       ma2 
# 0.6836856 0.3117469 0.1090300  
```

#### Update 17.04.22: ADF regression GDP

New call signature:

```c++
    arma::vec ar_coefs,           // vector of AR coefficients
    arma::vec ma_coefs,           // vector of MA coefficients
    const arma::vec& innovs,      // vector of innovations
    const bool& cumsum = false,   // should cumsum of series be returned?
    const double& rho = 1         // for ADF-regression-type GDP
```

`ARMA_sim()` now supports simulating from ADF regression GDPs via the additional arguments `cumsum` and `rho`. Using `rho = <some coefficient>` and `cumsum = T` , the outcome will be _levels_ of a process with ADF regression representation 

<img src="https://render.githubusercontent.com/render/math?math={\displaystyle \Delta y_t = (\rho - 1) \cdot y_{t-1} %2B \sum_{j=1}^p \gamma_j \Delta y_{t-j} %2B e_i}#gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\displaystyle \color{white} \Delta y_t = (\rho - 1) \cdot y_{t-1} %2B \sum_{j=1}^p \gamma_j \Delta y_{t-j} %2B e_i}#gh-dark-mode-only">

where the AR polynomial is specified via the coefficients passed to `ar_coefs`.

Let's check by estimating an ADF regression using `urca::ur.df()`:

```r
set.seed(321)
y <- c(ARMA_sim(ar_coefs = c(.4, .3, .2), ma = 0, rnorm(300), cumsum = T, rho = .8))
urca::ur.df(y = y, type = "none", lags = 3)@testreg$coef[, 1]
```

```
#    z.lag.1 z.diff.lag1 z.diff.lag2 z.diff.lag3 
# -0.2330476   0.4443002   0.3192869   0.2276719 
```

... and because we like simulation:

```r
purrr::map_dfr(
    1:5e3, 
    ~ { 
        y <- c(ARMA_sim(ar_coefs = c(.4, .3, .2), ma = 0, rnorm(150), rho = .8))
        urca::ur.df(y = y, type = "none", lags = 3)@testreg$coef[,1]
    }
) %>% colMeans()
```

```
#   z.lag.1  z.diff.lag1 z.diff.lag2 z.diff.lag3 
# -0.2010864   0.3930918   0.2897235   0.1983707 
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

---

**Note to self:**

I had some trouble getting this to run on Apple Silicon Mac with the homebrew version of R. For now it runs fine using gfortran as shipped with the CRAN binaries for R. Edited Makevars accordingly (`file.edit("~/.R/Makevars")`):

```
FLIBS   = -L/opt/R/arm64/gfortran/lib -mtune=native
F77     = /opt/R/arm64/gfortran/bin/gfortran
FC      = /opt/R/arm64/gfortran/bin/gfortran
```

```r
dir.create("~/.R/")
file.create("~/.R/Makevars")
file.edit("~/.R/Makevars"
```
