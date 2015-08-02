// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvrnorm
arma::mat mvrnorm(int n, arma::vec mu, arma::mat Sigma);
RcppExport SEXP fmscr_mvrnorm(SEXP nSEXP, SEXP muSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    __result = Rcpp::wrap(mvrnorm(n, mu, Sigma));
    return __result;
END_RCPP
}
// sample_quantile
double sample_quantile(arma::colvec x, double p);
RcppExport SEXP fmscr_sample_quantile(SEXP xSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    __result = Rcpp::wrap(sample_quantile(x, p));
    return __result;
END_RCPP
}
// MSE_trim
double MSE_trim(arma::colvec x, double truth, double trim);
RcppExport SEXP fmscr_MSE_trim(SEXP xSEXP, SEXP truthSEXP, SEXP trimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< double >::type trim(trimSEXP);
    __result = Rcpp::wrap(MSE_trim(x, truth, trim));
    return __result;
END_RCPP
}
// MAD
double MAD(arma::colvec x, double truth);
RcppExport SEXP fmscr_MAD(SEXP xSEXP, SEXP truthSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type truth(truthSEXP);
    __result = Rcpp::wrap(MAD(x, truth));
    return __result;
END_RCPP
}
// coverage_prob
double coverage_prob(arma::mat conf_intervals, double truth);
RcppExport SEXP fmscr_coverage_prob(SEXP conf_intervalsSEXP, SEXP truthSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type conf_intervals(conf_intervalsSEXP);
    Rcpp::traits::input_parameter< double >::type truth(truthSEXP);
    __result = Rcpp::wrap(coverage_prob(conf_intervals, truth));
    return __result;
END_RCPP
}
// median_width
double median_width(arma::mat conf_intervals);
RcppExport SEXP fmscr_median_width(SEXP conf_intervalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type conf_intervals(conf_intervalsSEXP);
    __result = Rcpp::wrap(median_width(conf_intervals));
    return __result;
END_RCPP
}
// shortest_CI
arma::rowvec shortest_CI(arma::vec x, double size, double inc);
RcppExport SEXP fmscr_shortest_CI(SEXP xSEXP, SEXP sizeSEXP, SEXP incSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type inc(incSEXP);
    __result = Rcpp::wrap(shortest_CI(x, size, inc));
    return __result;
END_RCPP
}
// clip
NumericVector clip(NumericVector x, double a, double b);
RcppExport SEXP fmscr_clip(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    __result = Rcpp::wrap(clip(x, a, b));
    return __result;
END_RCPP
}
// OLSvsIV_limit_sim
List OLSvsIV_limit_sim(double tau, double pi_sq, int n_sim);
RcppExport SEXP fmscr_OLSvsIV_limit_sim(SEXP tauSEXP, SEXP pi_sqSEXP, SEXP n_simSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type pi_sq(pi_sqSEXP);
    Rcpp::traits::input_parameter< int >::type n_sim(n_simSEXP);
    __result = Rcpp::wrap(OLSvsIV_limit_sim(tau, pi_sq, n_sim));
    return __result;
END_RCPP
}
// OLSvsIV_nonsimCI
List OLSvsIV_nonsimCI(double tau, double pi_sq, double size, int n_sim);
RcppExport SEXP fmscr_OLSvsIV_nonsimCI(SEXP tauSEXP, SEXP pi_sqSEXP, SEXP sizeSEXP, SEXP n_simSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type pi_sq(pi_sqSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type n_sim(n_simSEXP);
    __result = Rcpp::wrap(OLSvsIV_nonsimCI(tau, pi_sq, size, n_sim));
    return __result;
END_RCPP
}
// OLSvsIV_second_step
List OLSvsIV_second_step(double tau, double pi_sq, double size, double inc, int n_sim);
RcppExport SEXP fmscr_OLSvsIV_second_step(SEXP tauSEXP, SEXP pi_sqSEXP, SEXP sizeSEXP, SEXP incSEXP, SEXP n_simSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type pi_sq(pi_sqSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type inc(incSEXP);
    Rcpp::traits::input_parameter< int >::type n_sim(n_simSEXP);
    __result = Rcpp::wrap(OLSvsIV_second_step(tau, pi_sq, size, inc, n_sim));
    return __result;
END_RCPP
}
// OLSvsIV_onestepCI
List OLSvsIV_onestepCI(double tau, double pi_sq, double size, double inc, int n_sim_outer, int n_sim_inner);
RcppExport SEXP fmscr_OLSvsIV_onestepCI(SEXP tauSEXP, SEXP pi_sqSEXP, SEXP sizeSEXP, SEXP incSEXP, SEXP n_sim_outerSEXP, SEXP n_sim_innerSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type pi_sq(pi_sqSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type inc(incSEXP);
    Rcpp::traits::input_parameter< int >::type n_sim_outer(n_sim_outerSEXP);
    Rcpp::traits::input_parameter< int >::type n_sim_inner(n_sim_innerSEXP);
    __result = Rcpp::wrap(OLSvsIV_onestepCI(tau, pi_sq, size, inc, n_sim_outer, n_sim_inner));
    return __result;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP fmscr_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpp_hello());
    return __result;
END_RCPP
}
