#!/bin/bash

R -e "Rcpp::compileAttributes('RWSIM')"

R CMD build RWSIM

R -e "devtools::install('RWSIM')"