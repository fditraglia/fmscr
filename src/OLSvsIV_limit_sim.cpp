#include <RcppArmadillo.h>
#include "helper_functions.h"
using namespace Rcpp;

// Class to simulate from the limit experiment for the example from Section 5.1
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


//' Coverage and Width of OLS, TSLS and Naive CIs
//'
//' @param tau Controls degree of endogeneity of OLS.
//' @param pi_sq First-stage R-squared (strengh of instruments).
//' @param size One minus the norminal coverage probability of the intervals.
//' @param n_sim Number of simulation draws from the limit experiment.
//' @return List containing empirical coverage probabilities and median width
//' of confidence intervals for the OLS and TSLS estimators and the same for a
//' "naive" confidence interval for the FMSC-selected estimator.
//' @details This function gives results based on simulations from the limit
//' experiment for the OLS versus TSLS example in Section 5.1 of the paper.
//' The confidence intervals computed here are for the non-simulation based
//' procedures: the OLS estimator, the TSLS estimator and a naive interval for
//' the post-FMSC estimator. This naive interval is constructed from the
//' textbook interval for whichever estimator the FMSC selects: if OLS is
//' selected it uses the standard OLS interval, and if TSLS is selected it uses
//' the TSLS interval. This procedure can perform very badly depending on
//' parameter values. Note that the median widths in this example are not
//' particularly interesting: the width of each OLS and TSLS interval is fixed
//' across all simulations and the median width of the naive interval equals
//' that of OLS when OLS is chosen more than 50 percent of the time and equals
//' that of TSLS otherwise. The median widths are provided merely for
//' consistency with other functions for which this quantity is more
//' interesting, namely the simulation-based intervals that try to correct some
//' of the deficiencies of the naive interval.
//' @examples
//' foo <- OLSvsIV_nonsimCI(tau = 3, pi_sq = 0.1)
//' as.data.frame(foo)
// [[Rcpp::export]]
List OLSvsIV_nonsimCI(double tau, double pi_sq, double size = 0.05,
                      int n_sim = 50000){
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

// Class to implement "second-step" of two-stage CI procedure from Section 4.3
class second_step_OLS_IV {
  public:
    arma::rowvec fmscCI, avgCI;
    second_step_OLS_IV(double, double, double, double, int); //Class constructor
  private:
    arma::vec w, fmsc, wavg, avg;
    arma::uvec use_ols, use_tsls;
};
// Class constructor
second_step_OLS_IV::second_step_OLS_IV(double tau, double pi_sq,
                                       double size = 0.05, double inc = 0.001,
                                       int n_sim = 1000){
  limit_sim_OLS_IV sims(tau, pi_sq, n_sim);
  double tau_var = (1 - pi_sq) / pi_sq;
  w = tau_var / arma::pow(sims.tauhat, 2);
  use_ols = find(w >= 0.5);
  use_tsls = find(w < 0.5);
  fmsc = arma::vec(n_sim);
  fmsc.elem(use_ols) = sims.ols.elem(use_ols);
  fmsc.elem(use_tsls) = sims.tsls.elem(use_tsls);
  wavg = arma::clamp(w, 0.0, 1.0);
  avg = wavg % sims.ols + (1 - wavg) % sims.tsls; // element-wise mult.
  fmscCI = shortest_CI(fmsc, size, inc);
  avgCI = shortest_CI(avg, size, inc);
}


//' Testing interface to second_step_OLS_IV class.
//'
//' @param tau Controls degree of endogeneity of OLS.
//' @param pi_sq First-stage R-squared (strengh of instruments).
//' @param size One minus the norminal coverage probability of the intervals.
//' @param inc Step size for grid over which width is minimized.
//' @param n_sim Number of simulation draws.
//' @return List with two elements: simulation-based two-sided confidence
//' interval for the post-FMSC estimator evaluated at tau and the same for the
//' moment average estimator.
//' @details This is a testing interface to the C++ class that is use, among
//' other things, as the second-step in the two-step simulation based confidence
//' interval construction from Section 4.3 of the paper, based on drawing from
//' the limit experiment. In this implementation it treats all population
//' parameters that enter the limit as know with the exception of tau. (That is,
//' it abstracts from sampling uncertainty.)
//' @examples
//' as.data.frame(OLSvsIV_second_step(tau = 3, pi_sq = 0.1))
// [[Rcpp::export]]
List OLSvsIV_second_step(double tau, double pi_sq, double size = 0.05,
                         double inc = 0.001, int n_sim = 1000){
  second_step_OLS_IV out(tau, pi_sq, size, inc, n_sim);
  return(List::create(Named("fmsc") = out.fmscCI,
                      Named("avg") = out.avgCI));
}

