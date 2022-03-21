Install using `devtools::install_github("https://github.com/mca91/RWSIM/", subdir = "RWSIM")`.

You can also `git clone` and then build and install locally by running the shell script `buildpackage.sh`.

**Note to self:**

I had some trouble getting this to run on Apple Silicon Mac with the homebrew version of R. For now it runs fine using gfortran as shipped with the CRAN binaries R. Edited Makevars accordingly (`file.edit("~/.R/Makevars")`):

```
# R native
FLIBS   = -L/opt/R/arm64/gfortran/lib
F77     = /opt/R/arm64/gfortran/bin/gfortran
FC      = /opt/R/arm64/gfortran/bin/gfortran
