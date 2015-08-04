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
  arma::vec w_avg(n_reps);

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
    double tau_squared_est = pow(tau(i), 2) - tau_var(i);
    double sq_bias_est;
    if((tau_squared_est / pow(s_x_sq(i), 2)) >= 0){
    sq_bias_est = tau_squared_est / pow(s_x_sq(i), 2);
    } else {
    sq_bias_est = 0;
    }
    double var_diff = s_e_sq_tsls(i) * (1 / g_sq(i) - 1 / s_x_sq(i));
    w_avg(i) = 1 / (1 + (sq_bias_est / var_diff));
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
                      Named("tau_var") = tau_var,
                      Named("w_avg") = w_avg));
}

//[[Rcpp::export]]
List sim_chooseIVs(double rho, double gamma, int N, int n_reps = 1000){
  RNGScope scope;
  arma::mat V1(3, 3);
  V1 << 1.0 << 0.5 - gamma * rho << rho << arma::endr
     << 0.5 - gamma * rho << 8.0/9.0 - pow(gamma, 2) << 0.0 << arma::endr
     << rho << 0.0 << 1.0 << arma::endr;
  arma::mat L = chol(V1, "lower");

  arma::vec s_w_sq(n_reps);
  arma::vec gamma_hat(n_reps);
  arma::vec q_sq_valid(n_reps);
  arma::vec q_sq_full(n_reps);
  arma::vec s_e_sq_valid(n_reps);
  arma::vec s_e_sq_full(n_reps);
  arma::vec SE_valid(n_reps);
  arma::vec SE_full(n_reps);
  arma::vec b_valid(n_reps);
  arma::vec b_full(n_reps);
  arma::vec tau_hat(n_reps);
  arma::vec tau_var(n_reps);

  for(int i = 0; i < n_reps; i++){
    arma::mat evw = trans(L * reshape(arma::vec(rnorm(3 * N)), 3, N));
    arma::vec w = evw.col(2);
    arma::mat Z = reshape(1/sqrt(3) * arma::vec(rnorm(3 * N)), N, 3);
    arma::vec x = (1.0 / 3.0) * arma::sum(Z,1) + gamma * w + evw.col(1);
    arma::vec y = 0.5 * x + evw.col(0);
    //orthogonalize w with respect to Z since we "pretend" we don't know
    //that w and Z are independent
    w = w - Z * arma::solve(Z, w);
    arma::vec pi_hat = arma::solve(Z, x);
    double ww = arma::dot(w, w);
    gamma_hat(i) = arma::dot(x, w) / ww;
    arma::vec x_valid = Z * pi_hat;
    b_valid(i) = dot(x_valid, y) / dot(x_valid, x_valid);
    arma::vec resid_valid = y - b_valid(i) * x;
    s_e_sq_valid(i) = dot(resid_valid, resid_valid) / N;
    arma::vec x_full = x_valid + gamma_hat(i) * w;
    b_full(i) = dot(x_full, y) / dot(x_full, x_full);
    arma::vec resid_full = y - b_full(i) * x;
    s_e_sq_full(i) = dot(resid_full, resid_full) / N;
    s_w_sq(i) = ww / N;
    arma::mat Szz = Z.t() * Z / N;
    q_sq_valid(i) = as_scalar(pi_hat.t() * Szz * pi_hat);
    q_sq_full(i) = q_sq_valid(i) + pow(gamma_hat(i), 2) * s_w_sq(i);
    SE_valid(i) = sqrt(s_e_sq_valid(i) / (N * q_sq_valid(i)));
    SE_full(i) = sqrt(s_e_sq_full(i) / (N * q_sq_full(i)));
    tau_hat(i) = dot(w, resid_valid) / sqrt(N);
    tau_var(i) = s_w_sq(i) * s_e_sq_valid(i) * q_sq_full(i) / q_sq_valid(i);
  }
  return List::create(Named("s_w_sq") = s_w_sq,
                      Named("gamma") = gamma_hat,
                      Named("q_sq_valid") = q_sq_valid,
                      Named("q_sq_full") = q_sq_full,
                      Named("s_e_sq_valid") = s_e_sq_valid,
                      Named("s_e_sq_full") = s_e_sq_full,
                      Named("SE_valid") = SE_valid,
                      Named("SE_full") = SE_full,
                      Named("b_valid") = b_valid,
                      Named("b_full") = b_full,
                      Named("tau_hat") = tau_hat,
                      Named("tau_var") = tau_var);
}
