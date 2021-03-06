% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbla_f_gsem.R
\name{pbla_f_gsem}
\alias{pbla_f_gsem}
\title{f-based PBLA (General SEM)}
\usage{
pbla_f_gsem(r, beta, gamma, N, A = 1, lag = 0)
}
\arguments{
\item{r}{numeric vector of removal times}

\item{beta}{numeric rate}

\item{gamma}{numeric rate}

\item{N}{integer population size}

\item{A}{integer patient zeros}

\item{lag}{numeric fixed lag}
}
\value{
negative log likelihood
}
\description{
Via 3.3.1., compute pair-based likelihood approximation. Supports exponential infectious periods.
}
