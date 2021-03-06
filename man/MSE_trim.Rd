% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{MSE_trim}
\alias{MSE_trim}
\title{Calculate (trimmed) mean-squared error.}
\usage{
MSE_trim(x, truth, trim = 0)
}
\arguments{
\item{x}{Vector of estimates.}

\item{trim}{Fraction of estimates to discard (half from each tail) before
calculating MSE (defaults to zero)}

\item{true}{True parameter value.}
}
\value{
(trimmed) mean-squared error
}
\description{
Calculate (trimmed) mean-squared error.
}
\examples{
x <- rnorm(1000) + 0.5
MSE_trim(x, 0)
MSE_trim(x, 0, 0.1)
}

