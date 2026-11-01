\name{pc_W1}
\alias{pc_W1}
\alias{pc_M}
\alias{pc_softmax}
\title{Distances between distributions on ordered categories}
\description{Distances between distributions on ordered categories.}
\usage{
pc_W1(P, Q, a)
pc_M(P, Q, a)
pc_softmax(eta)
}
\arguments{
\item{P}{matrix of distributions, one per row (rows sum to 1).}
\item{Q}{reference distributions: a matrix with the same dimensions as \code{P}, or one vector used for every row.}
\item{a}{increasing representative values of the \eqn{K} categories.}
\item{eta}{matrix of linear predictors, one composition per row.}
}
\details{Probability inputs must be finite, with non-negative normalised compositions, matching dimensions and strictly increasing support. These requirements are checked by \code{pc_W1}; callers of \code{pc_M} and \code{pc_softmax} must ensure their input contracts. Softmax also permits negative-infinite logits for structural-zero cells, provided every row has a finite maximum and no missing or positive-infinite entries. Softmax expects a matrix; a vector is interpreted as a one-column matrix.}
\value{
\code{pc_W1}: row-wise Wasserstein-1 distance \eqn{\sum_{k<K} |P_k - Q_k| (a_{k+1} - a_k)} with cumulative sums \eqn{P_k, Q_k}.
\code{pc_M}: row-wise signed mean displacement \eqn{\sum_k a_k (p_k - q_k)}. Always \eqn{W_1 \ge |M|}, with equality when the
cumulative distributions do not cross. \code{pc_softmax}: row-wise softmax.
}
\examples{
a <- c(0.5, 2.5, 9.5, 19.5)
p <- c(0.1, 0.2, 0.3, 0.4); q <- c(0.25, 0.25, 0.25, 0.25)
pc_W1(rbind(p), q, a); pc_M(rbind(p), q, a)
}
