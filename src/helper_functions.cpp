#include <RcppArmadillo.h>
using namespace Rcpp;

//' Generate draws from a multivariate normal distribution
//'
//' @param n Number of samples.
//' @param mu Mean vector.
//' @param Sigma Covariance matrix.
//' @return Matrix of draws from the normal distribution: each row is a draw.
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

//' Calculate (trimmed) mean-squared error.
//'
//' @param x Vector of estimates.
//' @param true True parameter value.
//' @param trim Fraction of estimates to discard (half from each tail) before
//' calculating MSE (defaults to zero)
//' @return (trimmed) mean-squared error
//' @examples
//' x <- rnorm(1000) + 0.5
//' MSE_trim(x, 0)
//' MSE_trim(x, 0, 0.1)
// [[Rcpp::export]]
double MSE_trim(arma::colvec x, double truth, double trim = 0.0){
  int k = x.n_elem;
  int tail_drop = ceil(k * trim / 2);
  arma::colvec x_trimmed = sort(x);
  x_trimmed = x_trimmed(arma::span(tail_drop, k - tail_drop - 1));
  arma::colvec truth_vec = truth * arma::ones(x_trimmed.n_elem);
  arma::colvec errors = x_trimmed - truth_vec;
  double MSE = dot(errors, errors) / errors.n_elem;
  return(MSE);
}

//' Calculate median absolute deviation
//'
//' @param x Vector of estimates.
//' @param truth True value of the parameter.
//' @return Median absolute deviation.
//' @examples
//' x <- rnorm(1000) + 0.5
//' MAD(x, 0)
// [[Rcpp::export]]
double MAD(arma::colvec x, double truth){
  arma::colvec truth_vec = truth * arma::ones(x.n_rows);
  arma::colvec abs_dev = abs(x - truth_vec);
  return(median(abs_dev));
}


//' Calculate the empirical coverage probability of a matrix of confidence
//' intervals.
//'
//' @param conf_intervals Matrix of confidence intervals in which each row is a
//' CI, the first column is the lower confidence limit, and the second column is
//' the upper confidence limit.
//' @param truth True value of the parameter for which the CIs were constructed.
//' @return Empirical coverage probability.
//' @examples
//' xbar <- replicate(1000, mean(rnorm(25)))
//' ME <- qnorm(0.975) / 5
//' CIs <- cbind(xbar - ME, xbar + ME)
//' coverage_prob(CIs, 0)
// [[Rcpp::export]]
double coverage_prob(arma::mat conf_intervals, double truth){
  arma::colvec truth_vec = truth * arma::ones(conf_intervals.n_rows);
  arma::colvec cover_lower = arma::conv_to<arma::colvec>
                    ::from(conf_intervals.col(0) < truth_vec);
  arma::colvec cover_upper = arma::conv_to<arma::colvec>
                    ::from(conf_intervals.col(1) > truth_vec);
  arma::colvec cover = cover_lower % cover_upper;
  return(sum(cover) / cover.n_elem);
}


//' Calculate empirical median width of a matrix of confidence intervals.
//'
//' @param conf_intervals Matrix of confidence intervals in which each row is a
//' CI, the first column is the lower confidence limit, and the second column is
//' the upper confidence limit.
//' @return Empirical median width of the confidence intervals.
//' @examples
//' xbar <- replicate(1000, mean(rnorm(25)))
//' ME <- qnorm(0.975) / 5
//' CIs <- cbind(xbar - ME, xbar + ME)
//' median_width(CIs)
// [[Rcpp::export]]
double median_width(arma::mat conf_intervals){
  arma::colvec width = conf_intervals.col(1) - conf_intervals.col(0);
  return(median(width));
}


//' Calculate shortest two-sided confidence interval
//'
//' @param x Vector of simulations from sampling distribution of an estimator.
//' @param size One minus the desired coverage probability.
//' @param inc Step size for grid over which width is minimized.
//' @return Shortest (1 - size) * 100 percent interval based on the data provided.
//' @examples
//' x <- rnorm(1000)
//' shortest_CI(x)
// [[Rcpp::export]]
arma::rowvec shortest_CI(arma::vec x, double size = 0.05, double inc = 0.001){
  int n_inc = size / inc - 1;
  arma::vec aL = arma::linspace(inc, size - inc, n_inc);
  arma::vec aU = 1 - (size - aL);
  arma::vec L(n_inc);
  arma::vec U(n_inc);
  for(int i = 0; i < n_inc; i++){
    L(i) = sample_quantile(x, aL(i));
    U(i) = sample_quantile(x, aU(i));
  }
  arma::vec widths = U - L;
  arma::uword min_index;
  widths.min(min_index);
  arma::rowvec out(2);
  out(0) = L(min_index);
  out(1) = U(min_index);
  return(out);
}



