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

// function which computes a special sum needed
// in the detrending procedure proposed by Chang (2004)
// [[Rcpp::export]]
arma::vec t_col_sum(const arma::mat& dat, const int& t) {
  arma::vec out(dat.n_rows, fill::zeros);
  for(int i = 1; i<=t; i++) {
    out += double(i) * dat.col(i);
  }
  return(out);
}

// function which implements the detrending procedure suggested by
// Chang (2002)
// [[Rcpp::export]]
arma::mat detrend_recursive(const arma::mat& dat) {
  arma::mat detr(dat.n_rows, dat.n_cols-1);
  for(int i = 1; i<=dat.n_cols-1; i++) {
    detr.col(i-1) = dat.col(i) +  2.0/i * sum(dat.cols(1,i), 1) - 6.0/(i * (i + 1)) * t_col_sum(dat, i);
  }
  return(detr);
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

// function which OLS-detrends a time series supplied in a 1xT matrix
// and returns a matrix of the same dimension
// [[Rcpp::export]]
arma::mat Detrend(arma::mat Y, std::string model = "c") {
  if(model == "nc") {
    return Y;
  }
  // setup regressor matrix
  arma::vec j(Y.n_cols, fill::ones);
  arma::vec time = seq_cpp(1, Y.n_cols);
  arma::mat X = join_rows(j, time);
  if(model == "c") {
    X.shed_col(1);
  }
  arma::vec y = vectorise(Y);
  // run regression and obtain the detrended series (the residuals)
  arma::colvec beta = arma::solve(X, y);
  arma::mat Ydetr = y - X * beta;
  return Ydetr.t();
}

// function which performs local GLS detrending of a time series
// [[Rcpp::export]]
arma::mat GLS_Detrend(arma::mat dat, std::string model) {

  int T = dat.n_cols;

  arma::vec X_cbar = vectorise(dat);
  arma::mat z_cbar(T, 1, fill::zeros);
  arma::mat res(1, T, fill::zeros);

  if(model == "c") {
    double cbar = 7;

    for(int i = 1; i<T; i++) {
      X_cbar(i, 0) = dat(0, i) - (1-cbar/T) * dat(0, i-1);
    }

    z_cbar.fill(1-(1-cbar/T));
    z_cbar(0, 0) = 1;

    // run regression and obtain residuals
    arma::colvec beta = arma::solve(z_cbar, X_cbar);
    res = dat - as_scalar(beta);

  } else if(model == "ct") {
    double cbar = 13.5;
    arma::vec trd = seq_cpp(1, T);

    z_cbar.set_size(T, 2);
    z_cbar.col(0).fill(1-(1-cbar/T));
    z_cbar(0, 0) = 1;
    z_cbar(0, 1) = 1;

    for(int i = 1; i<T; i++) {
      X_cbar(i, 0) = dat(0, i) - (1-cbar/T) * dat(0, i-1);
      z_cbar(i, 1) = i+1 - (1-cbar/T) * i;
    }

    // run regression and obtain residuals
    arma::colvec beta = arma::solve(z_cbar, X_cbar);
    res = dat - as_scalar(beta(0, 0)) - beta(1, 0) * trd.t();

  }

  return res;

}

// fast function which computes *time series* (Augmented) Dickey-Fuller regressions
// returns a field object with residuals, estimated coefficients on lagged differences
// and sigmahat^2
// [[Rcpp::export]]
arma::field<arma::mat> DF_Reg_field(
    const arma::mat& Y,
    const int& p,
    const std::string& model,
    const arma::uvec& remove_lags // cherry-pick your lags
  ) {

  // initialize output field
  arma::field<arma::mat> F;
  F.set_size(5, 1);

  //initialise pmax, nz
  int pmax = p, nz = 0;

  // initialise uvec of lags to be removed
  arma::uvec tbr;

  if(pmax != 0) {
    // _how many_ lags need to be excluded?
    nz = accu(remove_lags != 0);
    // initialise the "full" lags vector
    arma::vec lags = seq_cpp(1, pmax);
    // modify the full lags vector accordingly
    if(nz != 0) lags.shed_rows(remove_lags - 1);
    // modify "pmax" if lags are removed
    pmax = lags.max();
    // lags to be removed
    tbr = find(remove_lags < pmax);
  }

  // obtain ADF regression matrix
  arma::mat X = DF_Reg_Mat(Y, pmax, model);

  // remove lags from the ADF regression matrix
  if(nz != 0) X.shed_cols(remove_lags.rows(tbr));

  // obtain dependent variable for ADF regression
  arma::mat y = mdiff(Y, 0, true);
  y = y.rows(p, y.n_rows - 1);

  // run OLS, compute residuals
  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X * coef;

  // assign outcomes to output field object
  arma::colvec betas;
  if(p != 0) betas = coef.rows(1, p - nz);

  F(0, 0) = resid;
  F(1, 0) = coef.row(0);
  F(2, 0) = betas;
  if(model != "nc") {
    F(3, 0) = coef.rows(pmax + 1 - nz, coef.n_rows - 1);
  }
  F(4, 0) = y;

  return F;
}

// function which computes the AR estimate of the spectral density at frequency zero
// for a time series of AR(k) residuals
// [[Rcpp::export]]
double S2_AR(
    const arma::mat& dat,
    const int& k,                         // maximum lag order (after removing lags!)
    const std::string& model,
    const arma::uvec& remove_lags         // cherry-pick your lags
  ) {

  arma::field<arma::mat> F = DF_Reg_field(dat, k, model, remove_lags);

  // assign outputs from (A)DF regression
  arma::colvec res = F(0, 0);

  int n = dat.n_cols;
  double sigmahat2 =  arma::as_scalar(
    // (w)
    arma::trans(res) * res / (n - F(2,0).n_rows)
    );

  //double sigmahat2 = arma::as_scalar(F(1, 0));
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
    const int& p,                         // maximum lag order (after removing lags!)
    const std::string& model,             // deterministic component
    const arma::uvec& remove_lags         // cherry-pick your lags
  ) {

  //initialise pmax, nz
  int pmax = p, nz = 0;

  // initialise uvec of lags to be removed
  arma::uvec tbr;

  if(pmax != 0) {
    // _how many_ lags need to be excluded?
    nz = accu(remove_lags != 0);
    // initialise the "full" lags index vector
    arma::vec lags = seq_cpp(1, pmax);
    // modify the full lags vector accordingly
    if(nz != 0) lags.shed_rows(remove_lags - 1);
    // modify "pmax" if lags are removed
    pmax = lags.max();
    // lags to be removed
    tbr = find(remove_lags < pmax);
  }

  arma::mat X = DF_Reg_Mat(Y, pmax, model);

  // remove lags from the ADF regression matrix
  if(nz != 0) X.shed_cols(remove_lags.rows(tbr));

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

// compute DFGLS t-test
// [[Rcpp::export]]
double ERS(const arma::mat& dat,
           const int& p,
           const std::string& model,
           const arma::uvec& remove_lags
           ) {

  arma::mat dat_detr = GLS_Detrend(dat, model);
  double t = 0;
  if(model == "c") {
    t = DF_Reg(dat_detr, p, "nc", remove_lags);
  } else if(model == "ct") {
    t = DF_Reg(dat_detr, p, "c", remove_lags);
  }
  return t;
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
arma::mat ARMA_sim(
    arma::vec ar_coefs,           // vector of AR coefficients
    arma::vec ma_coefs,           // vector of MA coefficients
    const arma::vec& innovs,      // vector of innovations
    const bool& cumsum = false,   // should cumsum of series be returned?
    const double& rho = 1,        // for ADF-regression-type GDP
    const double& delta = 0       // trend coefficient in ADF-regression-type GDP
) {

  size_t n = innovs.n_elem; //number of observations
  size_t p, q;  //AR and MA orders

  arma::vec one(1, fill::eye); // one in MA polynomial

  // initialise linear trend vector
  arma::vec tr = seq_cpp(1, n);

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
  e(span(q, n+q-1)) = innovs; // join n innovations and q zero starting value for MA recursion

  arma::vec ma(n, fill::zeros); // for n realizations of MA process
  arma::vec u(n+p, fill::zeros); // for n realizations of ARMA recursion + p starting values in AR recursion

  arma::vec theta = reverse(ma_coefs); // assign coefficients in MA polynomial in reversed order
  arma::vec phi = reverse(ar_coefs); // assign coefficients in AR polynomial in reversed order

    // run MA recursion
    for(size_t t=0; t<n; t++) {
      ma.row(t) = theta.t() * e(span(t, t+q));
    }

    // run ARMA recursion, return (cumsum of) series
    if(p != 0) {
      for(size_t t=p; t<n+p; t++) {
        u.row(t) = (rho - 1.0) * accu(u) + phi.t() * u(span(t-p, t-1)) + delta * tr.row(t-p) + ma.row(t-p);
      }
      if(cumsum) {
        u = arma::cumsum(u);
      }
      u.shed_rows(0, p-1); // drop the p starting values in AR(MA) output vector
      return std::move(u);
    } else {
      if(cumsum) {
        return arma::cumsum(ma);
      } else {
        return std::move(ma);
      }
    }

}


// function which computes (M)AIC/(M)BIC for an ADF model
// [[Rcpp::export]]
double IC(const arma::mat& Y,
           const int& p,
           const int& pmax,
           const std::string& model,
           const arma::uvec& remove_lags,
           const std::string& penalty = "BIC",
           const bool& modified = false
           ) {

  // time series length
  double T = Y.n_elem;
  // do the regression stuff
  arma::field<arma::mat> reg_results = DF_Reg_field(Y, p, model, remove_lags);
  // obtain residuals, trimmed levels, rho estimate
  arma::mat e = reg_results(0, 0);
  arma::mat yLag = Y.cols(pmax + 1, Y.n_cols - 2);
  arma::mat hat_rho = reg_results(1, 0);

  if(pmax > p) e.shed_rows(0, pmax-p-1);

  // number of estimated coefficients on lagged differences in ADF regression
  double k = reg_results(2, 0).n_elem;

  // estimate sigma
  double hat_sigma_sq = arma::as_scalar(arma::trans(e) * e / (T - pmax));

  // compute tau
  double tau = 0;
  if(modified) tau = 1 / (hat_sigma_sq) * dot(hat_rho, hat_rho) * dot(yLag, yLag);

  // compute penalty term
  double C;
  if(penalty == "BIC") {
    C = log(T);
  } else {
    C = 2; // AIC
  }

  // compute (M)IC
  double IC = log(hat_sigma_sq) + k * (C + tau) / T;

  return IC;
}


// Function for computing Stock's M-tests
// [[Rcpp::export]]
arma::mat Mtests(
    const arma::mat dat,
    const int p,
    const std::string model,
    const arma::uvec& remove_lags         // cherry-pick your lags
  ) {

  int L = dat.n_cols;
  double T = (double) L;
  arma::mat ylag = dat.cols(0, L-2);

  // S2_AR
  double S2 = S2_AR(dat, p, model, remove_lags);
  S2 = S2 * (T - 2 - p) / (T-1);
  T = T-1;

  double MZ_alpha = (1/T * dat(0, L-1) * dat(0, L-1) - S2) / (2 * 1/(T*T) * accu(ylag * ylag.t()));
  double MSB = sqrt( 1/(T*T) * accu(ylag * ylag.t()) / S2);
  double MZ_t = MZ_alpha * MSB;

  arma::mat out(1, 3);
  out(0, 0) = MZ_alpha;
  out(0, 1) = MSB;
  out(0, 2) = MZ_t;

  return out;

}

// fast function which computes a time series DF-regression and returns the normalized bias statistic
// [[Rcpp::export]]
double DF_Reg_NB(arma::mat y, int p, std::string model) {

  size_t T = y.n_cols;

  arma::mat X = DF_Reg_Mat(y, p, model);

  y = mdiff(y, 0, true);
  y = y.rows(p, y.n_rows-1);

  arma::colvec coef = arma::solve(X, y);

  double nbs = T * coef(0);

  return nbs;
}


