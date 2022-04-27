// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// seq_cpp
arma::vec seq_cpp(const int& lo, const int& hi);
RcppExport SEXP _RWSIM_seq_cpp(SEXP loSEXP, SEXP hiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type lo(loSEXP);
    Rcpp::traits::input_parameter< const int& >::type hi(hiSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_cpp(lo, hi));
    return rcpp_result_gen;
END_RCPP
}
// mlag
arma::mat mlag(const arma::mat& dat, const int& p, const bool& drop);
RcppExport SEXP _RWSIM_mlag(SEXP datSEXP, SEXP pSEXP, SEXP dropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const bool& >::type drop(dropSEXP);
    rcpp_result_gen = Rcpp::wrap(mlag(dat, p, drop));
    return rcpp_result_gen;
END_RCPP
}
// mdiff
arma::mat mdiff(const arma::mat& dat, const int& p, const bool& drop);
RcppExport SEXP _RWSIM_mdiff(SEXP datSEXP, SEXP pSEXP, SEXP dropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const bool& >::type drop(dropSEXP);
    rcpp_result_gen = Rcpp::wrap(mdiff(dat, p, drop));
    return rcpp_result_gen;
END_RCPP
}
// DF_Reg_Mat
arma::mat DF_Reg_Mat(const arma::mat& y, const int& p, const std::string& model, const bool& omit_y_lag);
RcppExport SEXP _RWSIM_DF_Reg_Mat(SEXP ySEXP, SEXP pSEXP, SEXP modelSEXP, SEXP omit_y_lagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const bool& >::type omit_y_lag(omit_y_lagSEXP);
    rcpp_result_gen = Rcpp::wrap(DF_Reg_Mat(y, p, model, omit_y_lag));
    return rcpp_result_gen;
END_RCPP
}
// DF_Reg_field
arma::field<arma::mat> DF_Reg_field(const arma::mat& Y, const int& p, const std::string& model, const arma::uvec& remove_lags);
RcppExport SEXP _RWSIM_DF_Reg_field(SEXP YSEXP, SEXP pSEXP, SEXP modelSEXP, SEXP remove_lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type remove_lags(remove_lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(DF_Reg_field(Y, p, model, remove_lags));
    return rcpp_result_gen;
END_RCPP
}
// S2_AR
double S2_AR(const arma::mat& dat, const int& k, const std::string& model, const arma::uvec& remove_lags);
RcppExport SEXP _RWSIM_S2_AR(SEXP datSEXP, SEXP kSEXP, SEXP modelSEXP, SEXP remove_lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type remove_lags(remove_lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_AR(dat, k, model, remove_lags));
    return rcpp_result_gen;
END_RCPP
}
// DF_Reg
double DF_Reg(const arma::mat& Y, const int& p, const std::string& model, const arma::uvec& remove_lags);
RcppExport SEXP _RWSIM_DF_Reg(SEXP YSEXP, SEXP pSEXP, SEXP modelSEXP, SEXP remove_lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type remove_lags(remove_lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(DF_Reg(Y, p, model, remove_lags));
    return rcpp_result_gen;
END_RCPP
}
// OLSRes
arma::field<arma::mat> OLSRes(const arma::vec& y, const arma::mat& X);
RcppExport SEXP _RWSIM_OLSRes(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(OLSRes(y, X));
    return rcpp_result_gen;
END_RCPP
}
// ARMA_sim
arma::mat ARMA_sim(arma::vec ar_coefs, arma::vec ma_coefs, const arma::vec& innovs, const bool& cumsum, const double& rho);
RcppExport SEXP _RWSIM_ARMA_sim(SEXP ar_coefsSEXP, SEXP ma_coefsSEXP, SEXP innovsSEXP, SEXP cumsumSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ar_coefs(ar_coefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ma_coefs(ma_coefsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type innovs(innovsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type cumsum(cumsumSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(ARMA_sim(ar_coefs, ma_coefs, innovs, cumsum, rho));
    return rcpp_result_gen;
END_RCPP
}
// test
arma::uvec test(int p, const arma::uvec& remove_lags);
RcppExport SEXP _RWSIM_test(SEXP pSEXP, SEXP remove_lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type remove_lags(remove_lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(test(p, remove_lags));
    return rcpp_result_gen;
END_RCPP
}
// BIC
double BIC(const arma::mat& Y, const int& p, const std::string& model, const arma::uvec& remove_lags);
RcppExport SEXP _RWSIM_BIC(SEXP YSEXP, SEXP pSEXP, SEXP modelSEXP, SEXP remove_lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type remove_lags(remove_lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(BIC(Y, p, model, remove_lags));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RWSIM_seq_cpp", (DL_FUNC) &_RWSIM_seq_cpp, 2},
    {"_RWSIM_mlag", (DL_FUNC) &_RWSIM_mlag, 3},
    {"_RWSIM_mdiff", (DL_FUNC) &_RWSIM_mdiff, 3},
    {"_RWSIM_DF_Reg_Mat", (DL_FUNC) &_RWSIM_DF_Reg_Mat, 4},
    {"_RWSIM_DF_Reg_field", (DL_FUNC) &_RWSIM_DF_Reg_field, 4},
    {"_RWSIM_S2_AR", (DL_FUNC) &_RWSIM_S2_AR, 4},
    {"_RWSIM_DF_Reg", (DL_FUNC) &_RWSIM_DF_Reg, 4},
    {"_RWSIM_OLSRes", (DL_FUNC) &_RWSIM_OLSRes, 2},
    {"_RWSIM_ARMA_sim", (DL_FUNC) &_RWSIM_ARMA_sim, 5},
    {"_RWSIM_test", (DL_FUNC) &_RWSIM_test, 2},
    {"_RWSIM_BIC", (DL_FUNC) &_RWSIM_BIC, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_RWSIM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
