\name{pc_basis}
\alias{pc_basis}
\title{Orthonormal polynomial basis over ordered categories}
\description{Orthonormal polynomial basis over ordered categories.}
\usage{pc_basis(a, degree = 2)}
\arguments{
\item{a}{increasing representative values of the \eqn{K} categories.}
\item{degree}{number of leading basis columns on which covariates act (the model uses 2: slope and curvature).}
}
\value{A list with \code{a}, \code{B} (\eqn{K \times} \code{degree}, the columns on which covariates and the unit effects
\eqn{u_{i1}, u_{i2}} act) and \code{Bh} (the remaining \eqn{K - 1 - } \code{degree} columns, which carry the higher-order unit
effects \eqn{z_{ij}}). The constant column is omitted: it is absorbed by the unit-period intercept.}
\examples{pc_basis(c(0.5, 2.5, 9.5, 19.5, 29.5, 39.5))$B}
