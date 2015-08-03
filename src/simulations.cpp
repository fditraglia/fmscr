#include <RcppArmadillo.h>
#include <stdexcept>
#include "helper_functions.h"
using namespace Rcpp;

//[[Rcpp::export]]
List sim_OLSvsIV(double rho, double pi_sq, int N, int n_reps = 1000){
  RNGScope scope;
  if(rho >= sqrt(1 - pi_sq)){
    throw std::invalid_argument("p.d. cov. matrix requires rho^2 < 1 - pi_sq");
  }
  arma::mat L(2,2);
  L << 1 << 0 << arma::endr
    << rho << sqrt(1 - pi_sq - pow(rho, 2)) << arma::endr;
  arma::vec s_x_sq(n_reps);
  arma::vec s_v_sq(n_reps);
  arma::vec g_sq(n_reps);
  arma::vec s_e_sq_ols(n_reps);
  arma::vec s_e_sq_tsls(n_reps);
  arma::vec SE_ols(n_reps);
  arma::vec SE_tsls(n_reps);
  arma::vec ols_estimate(n_reps);
  arma::vec tsls_estimate(n_reps);
  arma::vec tau(n_reps);
  arma::vec tau_var(n_reps);

  for(int i = 0; i < n_reps; i++){
    arma::mat errors = trans(L * reshape(arma::vec(rnorm(2 * N)), 2, N));
    arma::mat Z = reshape(1/sqrt(3) * arma::vec(rnorm(3 * N)), N, 3);
    arma::vec x = sqrt(pi_sq) * arma::sum(Z,1) + errors.col(1);
    arma::vec y = 0.5 * x + errors.col(0);
    double xx = arma::dot(x, x);
    ols_estimate(i) = arma::dot(x, y) / xx;
    arma::vec ols_resid = y - x * ols_estimate(i);
    arma::vec first_stage_coefs = arma::solve(Z, x);
    tsls_estimate(i) = as_scalar(solve(Z * first_stage_coefs, y));
    arma::vec tsls_resid = y - x * tsls_estimate(i);
    arma::mat Qz, Rz, Rzinv;
    qr_econ(Qz, Rz, Z);
    Rzinv = arma::inv(trimatu(Rz));
    arma::mat zz_inv = Rzinv * trans(Rzinv);
    arma::mat zz = trans(trimatu(Rz)) * trimatu(Rz);
    arma::mat zx = trans(Z) * x;
    g_sq(i) = as_scalar(trans(zx) * zz_inv * zx) / N;
    s_e_sq_ols(i) = dot(ols_resid, ols_resid) / N;
    s_e_sq_tsls(i) = dot(tsls_resid, tsls_resid) / N;
    s_x_sq(i) = xx / N;
    s_v_sq(i) = s_x_sq(i) - g_sq(i);
    tau(i) = dot(x, tsls_resid) / sqrt(N);
    tau_var(i) = s_e_sq_tsls(i) * s_x_sq(i) * (s_x_sq(i) / g_sq(i) - 1);
    SE_ols(i) = sqrt(s_e_sq_ols(i) / (N * s_x_sq(i)));
    SE_tsls(i) = sqrt(s_e_sq_tsls(i) / (N * g_sq(i)));
  }
  return(List::create(Named("s_x_sq") = s_x_sq,
                      Named("s_v_sq") = s_v_sq,
                      Named("g_sq") = g_sq,
                      Named("s_e_sq_ols") = s_e_sq_ols,
                      Named("s_e_sq_tsls") = s_e_sq_tsls,
                      Named("SE_ols") = SE_ols,
                      Named("SE_tsls") = SE_tsls,
                      Named("b_ols") = ols_estimate,
                      Named("b_tsls") = tsls_estimate,
                      Named("tau_hat") = tau,
                      Named("tau_var") = tau_var));
}
