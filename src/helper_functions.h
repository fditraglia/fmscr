#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

arma::mat mvrnorm(int n, arma::vec mu, arma::mat Sigma);
double sample_quantile(arma::colvec x, double p);
double MSE_trim(arma::colvec x, double truth, double trim = 0.0);
double MAD(arma::colvec x, double truth);
double coverage_prob(arma::mat conf_intervals, double truth);
double median_width(arma::mat conf_intervals);
arma::rowvec shortest_CI(arma::vec x, double size = 0.05, double inc = 0.001);

#endif
