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
//' @return Dataframe of draws from the limit experiment for OLS, TSLS and the
//' sample estimator of tau.
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
