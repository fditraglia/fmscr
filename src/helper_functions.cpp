#include <RcppArmadillo.h>
using namespace Rcpp;

//' Generate draws from a multivariate normal distribution
//'
//' @param n Number of samples.
//' @param mu Mean vector.
//' @param Sigma Covariance matrix.
//' @return Draws from the normal distribution.
//' @details This is essentially a stripped-down version of the mvrnorm function
//' from the MASS library in R. Through the magic of Rcpp we're transforming the
//' same standard normal draws as the R version. However, since Armadillo
//' follows a different convention from R in its definition of the
//' eign-decomposition, the output of this function will *not* be the same as
//' that of its R counterpart. Since we access R's function for generating
//' normal draws, we can set the seed from R.
//' @examples
//' mvrnorm(10, c(0,0), diag(1, 2, 2))
// [[Rcpp::export]]
arma::mat mvrnorm(int n, arma::vec mu, arma::mat Sigma){
  RNGScope scope;
  int p = Sigma.n_cols;
  arma::mat X = reshape(arma::vec(rnorm(p * n)), p, n);
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, Sigma);
  X = eigvec * diagmat(sqrt(eigval)) * X;
  X.each_col() += mu;
  return(X.t());
}

//' Calculate a sample quantile
//'
//' @param x Vector of data.
//' @param p Probability for desired quantile (e.g. 0.5 for median)
//' @return Sample quantile
//' @details There are many competing definitions of sample quantiles
//' (see Hyndman & Fan, 1996). Here we simply use the R default definition,
//' which corresponds to Definition 7 in Hyndman & Fan. See ?quantile in R for
//' more details.
//' @examples
//' foo <- rnorm(1000)
//' sample_quantile(foo, 0.16)
// [[Rcpp::export]]
double sample_quantile(arma::colvec x, double p){
  int n = x.n_elem;
  double m = 1 - p;
  int j = floor(n * p + m);
  double g = n * p + m - j;
  arma::colvec x_sort = sort(x);
  return((1 - g) * x_sort(j - 1) + g * x_sort(j));
}
