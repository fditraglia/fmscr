#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

arma::mat mvrnorm(int n, arma::vec mu, arma::mat Sigma);
double sample_quantile(arma::colvec x, double p);

#endif
