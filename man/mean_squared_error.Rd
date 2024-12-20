% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean_squared_error.R
\name{mean_squared_error}
\alias{mean_squared_error}
\title{Evaluates the mean square error loss function of a linear regression.}
\usage{
mean_squared_error(beta, X, Y)
}
\arguments{
\item{beta}{a vector of size dimX, the value at which the loss is computed
(i.e. we predict Y using X'beta)}

\item{X}{an array of size n x Tmax x dimX containing the values of the covariates X}

\item{Y}{a matrix of size n x Tmax containing the values of the dependent variable Y}
}
\value{
the value of the mean squared error loss function of the linear regression at beta
}
\description{
Evaluates the mean square error loss function of a linear regression.
}
