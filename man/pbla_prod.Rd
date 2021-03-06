% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbla_prod.R
\name{pbla_prod}
\alias{pbla_prod}
\title{Product PBLA}
\usage{
pbla_prod(r, beta, gamma, N, A = 1)
}
\arguments{
\item{r}{numeric vector of increasing removal times}

\item{beta}{numeric rate}

\item{gamma}{numeric rate}

\item{N}{integer population size}

\item{A}{integer patient zeros}
}
\value{
negative log likelihood
}
\description{
Based on product independence, compute pair-based likelihood approximation. Assumes exponential infectious periods.
}
