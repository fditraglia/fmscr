#include <RcppArmadillo.h>
#include "helper_functions.h"
using namespace Rcpp;


class limit_sim_OLS_IV {
  public:
    arma::colvec ols, tsls, tauhat, mu;
    arma::mat S;
    limit_sim_OLS_IV(double, double, int); //Class constructor
};
// Class constructor
limit_sim_OLS_IV::limit_sim_OLS_IV(double tau, double pi_sq, int n_sim = 10000){
  double tau_var = (1 - pi_sq) / pi_sq;
  S << 1 << 1 << 0 << arma::endr
    << 1 << 1 / pi_sq << -tau_var << arma::endr
    << 0 << -tau_var << tau_var << arma::endr;
  mu = arma::vec(3);
  mu(0) = tau; mu(1) = 0; mu(2) = tau;
  arma::mat sims = mvrnorm(n_sim, mu, S);
  ols = sims.col(0); tsls = sims.col(1); tauhat = sims.col(2);
}


//' Testing interface to limit_sim_OLS_IV class.
//'
//' @param tau Controls degree of endogeneity of OLS.
//' @param pi_sq First-stage R-squared (strengh of instruments).
//' @param n_sim Number of simulation draws.
//' @return List of draws from the limit experiment for OLS, TSLS and the
//' sample estimator of tau, as well as the mean vector and covariance matrix
//' used in the simulation.
//' @details This is a testing interface to the C++ class that is used to draw
//' from the limit experiment in the OLS versus TSLS simulation experiment from
//' section 5.1 of the paper.
//' @examples
//' foo <- OLSvsIV_limit_sim(5, 0.1)
//' cov(data.frame(ols = foo$ols, tsls = foo$tsls, tauhat = foo$tauhat))
// [[Rcpp::export]]
List OLSvsIV_limit_sim(double tau, double pi_sq, int n_sim = 10000){
  limit_sim_OLS_IV out(tau, pi_sq, n_sim);
  return(List::create(Named("ols") = out.ols,
                      Named("tsls") = out.tsls,
                      Named("tauhat") = out.tauhat,
                      Named("S") = out.S,
                      Named("mu") = out.mu));
}


// [[Rcpp::export]]
List OLSvsIV_nonsimCI(double tau, double pi_sq, double size = 0.05,
                      int n_sim = 10000){
  limit_sim_OLS_IV sims(tau, pi_sq, n_sim);
  double tau_var = (1 - pi_sq) / pi_sq;
  arma::vec w = tau_var / arma::pow(sims.tauhat, 2);
  arma::uvec use_ols = find(w >= 0.5);
  arma::uvec use_tsls = find(w < 0.5);
  arma::vec fmsc(n_sim);
  fmsc.elem(use_ols) = sims.ols.elem(use_ols);
  fmsc.elem(use_tsls) = sims.tsls.elem(use_tsls);

  double cval = R::qnorm(1 - size / 2, 0.0, 1.0, 1, 0);
  arma::vec olsL = sims.ols - cval;
  arma::vec olsU = sims.ols + cval;
  arma::vec tslsL = sims.tsls - (1 / sqrt(pi_sq)) * cval;
  arma::vec tslsU = sims.tsls + (1 / sqrt(pi_sq)) * cval;

  arma::vec naiveL(n_sim);
  naiveL.elem(use_ols) = olsL.elem(use_ols);
  naiveL.elem(use_tsls) = tslsL.elem(use_tsls);
  arma::vec naiveU(n_sim);
  naiveU.elem(use_ols) = olsU.elem(use_ols);
  naiveU.elem(use_tsls) = tslsU.elem(use_tsls);

  double naiveCover = coverage_prob(arma::join_rows(naiveL, naiveU), 0);
  double olsCover = coverage_prob(arma::join_rows(olsL, olsU), 0);
  double tslsCover = coverage_prob(arma::join_rows(tslsL, tslsU), 0);

  double naiveWidth = median_width(arma::join_rows(naiveL, naiveU));
  double olsWidth = median_width(arma::join_rows(olsL, olsU));
  double tslsWidth = median_width(arma::join_rows(tslsL, tslsU));

  return(List::create(Named("tslsCover") = tslsCover,
                      Named("olsCover") = olsCover,
                      Named("naiveCover") = naiveCover,
                      Named("tslsWidth") = tslsWidth,
                      Named("olsWidth") = olsWidth,
                    Named("naiveWidth") = naiveWidth));
}
