\name{pc_tolerance}
\alias{pc_tolerance}
\title{Tolerance (near-null) summaries of posterior draws}
\description{Tolerance (near-null) summaries of posterior draws.}
\usage{pc_tolerance(x, delta = c(1, 2, 5), distance = FALSE)}
\arguments{
\item{x}{vector or matrix of posterior draws (one column per quantity).}
\item{delta}{tolerances, in the units of \code{x}, fixed in advance.}
\item{distance}{\code{TRUE} for a non-negative distance such as \eqn{W_1}: reports \eqn{P(x < \delta)}; otherwise
\eqn{P(|x| < \delta)}.}
}
\details{Exact zero is on the boundary of a distance and is not tested; the tolerance probability, the equal-tailed interval and the
97.5\% quantile of \eqn{|x|} describe how large the quantity can plausibly be.}
\value{A data frame with the mean, sd, median, 95\% interval, \eqn{P(x > 0)}, the 97.5\% quantile of \eqn{|x|} and one tolerance
probability per \code{delta}.}
\examples{pc_tolerance(cbind(E = rnorm(4000, 0.8, 0.4)), delta = c(1, 2))}
