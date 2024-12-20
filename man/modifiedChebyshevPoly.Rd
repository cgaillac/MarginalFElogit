% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modifiedChebyshevPoly.R
\name{modifiedChebyshevPoly}
\alias{modifiedChebyshevPoly}
\title{Computes the coefficients of T_0(2u - 1), ..., T_n(2u-1), where T_k is the
usual k-th Chebyshev polynomial of the 1st kind}
\usage{
modifiedChebyshevPoly(n)
}
\arguments{
\item{n}{maximal k for which we compute the coefficients of T_k(2u - 1)}
}
\value{
the result as a (n + 1) x (n + 1) matrix where each row is one
of the polynomials. The first column corresponds to the constant coefficient,
and rows are filled with 0s at the end for T_0, ..., T_(n-1)
}
\description{
Computes the coefficients of T_0(2u - 1), ..., T_n(2u-1), where T_k is the
usual k-th Chebyshev polynomial of the 1st kind
}
