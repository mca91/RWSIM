// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// *****************************************************************************//


// simple armadillo version of R's base::seq()
// [[Rcpp::export]]
arma::vec seq_cpp(const int& lo,
                  const int& hi) {

  int n = hi - lo + 1;
  arma::vec sequence(n);
  for(int i = 0; i < n; i++) {
    sequence(i) = lo + i;
  }
  return sequence;
}


// function which transforms time series (panels) data to matrix
// of lags of the input series as needed for prewhitening.
// The returned matrix consists of a column of ones followed by
// p lags of the input series.
// [[Rcpp::export]]
arma::mat mlag(const arma::mat& dat,
               const int& p,
               const bool& drop) {

  // initialize output matrix with desiered dimensions and
  // fill with ones
  arma::mat out(dat.n_rows * (dat.n_cols-p), p+1, fill::ones);
  // append lag series up to order p
  if(p != 0) {
    for(int d=1; d <= p; d++) {
      out.col(d) = vectorise(dat.cols(p-d, dat.n_cols-1-d));
    }
  } else {
    out = join_rows(out, vectorise(dat));
  }
  // drop column of ones if desired
  if(drop) {
    out.shed_col(0);
  }
  return out;
}

// function which computes a regressor matrix whose columns
// are lags of differences up to p-th order of the input
// time series / panel.
// A column of ones is prepended if drop = false
// [[Rcpp::export]]
arma::mat mdiff(const arma::mat& dat,
                const int& p,
                const bool& drop) {

  // output matrix dimensions
  int M = dat.n_cols - p, N = dat.n_rows;
  // initialize output matrix
  arma::mat out(M * N, 1, fill::zeros);
  // compute first difference series
  arma::mat temp(dat.n_rows, dat.n_cols - 1, fill::zeros);
  for(size_t j=1; j <= dat.n_cols-1; j++) {
    temp.col(j-1) = dat.col(j) - dat.col(j-1);
  }
  // obtain the regressor matrix with p lags
  out = mlag(temp, p, false);
  // drop column of ones if desired
  if(drop) {
    out.shed_col(0);
  }
  return out;
}


// function which computes a regressor matrix for DF regressions
// [[Rcpp::export]]
arma::mat DF_Reg_Mat(const arma::mat& y,
                     const int& p,
                     const std::string& model = "nc",
                     const bool& omit_y_lag = false) {

  // output matrix
  arma::mat X;
  // vector of time indices
  arma::vec time = seq_cpp(1, y.n_cols);
  // compute first lag of the input series
  arma::mat ylag = y.cols(0, y.n_cols-2).t();

  if(p != 0) {
    // compute matrix of lags of the first differences of the input series
    arma::mat LFD = mdiff(y, p, true);
    // drop observations
    ylag.shed_rows(0, p-1);

    // regression with constant
    if(model == "c") {
      arma::vec j(LFD.n_rows, fill::ones);
      X = join_rows(ylag, LFD);
      X = join_rows(X, j);
      // regression with constant and linear time trend
    } else if(model == "ct") {
      arma::vec j(LFD.n_rows, fill::ones);
      time.shed_rows(0, p);
      X = join_rows(ylag,LFD);
      X = join_rows(X, j);
      X = join_rows(X, time);
    } else {
      X = join_rows(ylag, LFD);
    }
  } else {
    // regression with constant
    if(model == "c") {
      arma::vec j(ylag.n_rows, fill::ones);
      X = join_rows(ylag, j);
      // regression with constant and linear time trend
    } else if(model == "ct") {
      arma::vec j(ylag.n_rows, fill::ones);
      time.shed_rows(0, p);
      X = join_rows(ylag, j);
      X = join_rows(X, time);
    } else {
      X = ylag;
    }
  }
  if(X.n_cols > 1 && omit_y_lag) {
    X.shed_col(0);
  }
  return X;
}


// fast function which computes *time series* (Augmented) Dickey-Fuller regressions
// returns a field object with residuals, estimated coefficients on lagged differences
// and sigmahat^2
// [[Rcpp::export]]
arma::field<arma::mat> DF_Reg_field(
    const arma::mat& Y,
    const int& p,
    const std::string& model,
    const arma::uvec& remove_lags         // cherry-pick your lags
  ) {

  // initialize output field (set #rows conditional on p)
  arma::field<arma::mat> F;
  if(p == 0) {
    F.set_size(2, 1);
  } else {
    F.set_size(3, 1);
  }

  // transform data for regression
  arma::mat X = DF_Reg_Mat(Y, p, model);

  int nz = accu(remove_lags != 0);
  if(nz != 0) X.shed_cols(remove_lags);


  arma::mat y = mdiff(Y, 0, true);
  y = y.rows(p, y.n_rows - 1);

  // OLS and computation of residuals and sigmahat^2
  int n = X.n_rows, k = X.n_cols;
  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X * coef;
  arma::vec sig2 = arma::trans(resid)*resid / (n - k);

  F(0, 0) = resid;
  F(1, 0) = sig2;

  // p vector of estimated coefficients on _lagged differences_ of input series
  if(p > 0) {
    arma::colvec betas = coef.rows(1, p - nz);
    F(2, 0) = betas;
  }

  return F;
}

// function which computes the AR estimate of the spectral density at frequency zero
// for a time series of AR(k) residuals
// [[Rcpp::export]]
double S2_AR(
    const arma::mat& dat,
    const int& k,
    const std::string& model,
    const arma::uvec& remove_lags         // cherry-pick your lags
  ) {

  arma::field<arma::mat> F = DF_Reg_field(dat, k, model, remove_lags);

  // assign outputs from (A)DF regression
  arma::colvec res = F(0, 0);
  double sigmahat2 = arma::as_scalar(F(1, 0));
  // initialise vector for estimated coefficients on lagged differences
  arma::colvec betas;

  double sum_betasq = 1.0;

  if(k > 0) {
    betas = F(2, 0);
    double s = accu(betas);
    sum_betasq = (1 - s) * (1 - s);
  }

  return sigmahat2 / sum_betasq;

}

// fast function which computes *time series* (Augmented) Dickey-Fuller regressions and
// returns t-ratio on the coefficient of y_t-1
// [[Rcpp::export]]
double DF_Reg(
    const arma::mat& Y,                   // time series
    const int& p,                         // maximum lag order
    const std::string& model,             // deterministic component
    const arma::uvec& remove_lags         // cherry-pick your lags
  ) {

  arma::mat X = DF_Reg_Mat(Y, p, model);
  if(accu(remove_lags) != 0) X.shed_cols(remove_lags);

  int n = X.n_rows, k = X.n_cols;

  arma::mat y = mdiff(Y, 0, true);
  y = y.rows(p, y.n_rows-1);

  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X * coef;

  double sig2 = arma::as_scalar(arma::trans(resid) * resid / (n - k));
  arma::colvec stderrest = arma::sqrt(
    sig2 * arma::diagvec(arma::inv(arma::trans(X) * X))
  );

  double tstats = coef(0) / stderrest(0) ;

  return tstats;
}


// a fast OLS function which returns an arma::field object.
//
// First entry is a column vector of k(+1) estimated coefficients
// Second entry is a N x 1 residual matrix
//
// The dependent variable may be a columm vector. Regressors must be supplied as
// N X K(+1) matrix
// [[Rcpp::export]]
arma::field<arma::mat> OLSRes(const arma::vec& y,
                              const arma::mat& X) {

  arma::field<arma::mat> F(1, 2);

  arma::colvec coef = arma::solve(X, y);
  arma::mat resid = y - X * coef;

  F(0,0) = coef;
  F(0,1) = resid;
  return F;
}

// simple but fast function which runs an ARMA recursion on a vector of innovations
// [[Rcpp::export]]
arma::mat ARMA_sim(arma::vec ar_coefs, // vector of AR coefficients
                   arma::vec ma_coefs, // vector of MA coefficients
                   arma::vec innovs // vector of innovations
) {

  size_t n = innovs.n_elem; //number of observations
  size_t p, q;  //AR and MA orders

  arma::vec one(1, fill::eye); // one in MA polynomial

   if((ar_coefs.n_elem == 1) & (as_scalar(ar_coefs.row(0)) == 0)) {
     p = 0;
   } else {
     p = ar_coefs.n_elem;
   }

   if((ma_coefs.n_elem == 1) & (as_scalar(ma_coefs.row(0)) == 0)) {
     q = 0;
     ma_coefs = one;
   } else {
     q = ma_coefs.n_elem;
     ma_coefs = join_cols(one, ma_coefs); // add one to ma coefficients
   }

  arma::vec e(n+q, fill::zeros); // for n innovations + q starting values in MA recursion
  e(span(q, n+q-1)) = innovs; // combine n innovations and q zero starting value for MA recursion

  arma::vec ma(n, fill::zeros); // for n realizations of MA process
  arma::vec u(n+p, fill::zeros); // for n realizations of ARMA recursion + p starting values in AR recursion

  arma::vec theta = reverse(ma_coefs); // assign coefficients in MA polynomial in reversed order
  arma::vec phi = reverse(ar_coefs); // assign coefficients in AR polynomial in reversed order

    // run MA recursion
    for(size_t t=0; t<n; t++) {
      ma.row(t) = theta.t() * e(span(t, t+q));
    }

    // run ARMA recursion
    if(p != 0) {
      for(size_t t=p; t<n+p; t++) {
        u.row(t) = phi.t() * u(span(t-p, t-1)) + ma.row(t-p);
      }
      u.shed_rows(0, p-1); // drop starting values
      return std::move(u);
    } else {
      return std::move(ma);
    }

}

