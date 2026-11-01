\name{pc_contrasts}
\alias{pc_contrasts}
\alias{pc_regionyear}
\title{Posterior contrasts between predicted compositions}
\description{Posterior contrasts between predicted compositions.}
\usage{
pc_contrasts(draws, obj, points, contrasts, reference = NULL,
             version = c("marginal", "typical"), R = 200, seed = NULL)
pc_regionyear(draws, obj)
}
\arguments{
\item{draws}{output of \code{\link{pc_draws}}.}
\item{obj}{the \code{\link{pc_data}} object.}
\item{points}{named list of covariate vectors at which compositions are predicted (named by covariate; omitted covariates are 0).}
\item{contrasts}{named list. Each element is \code{list(type = "M", weights = c(point = w, ...))}, the linear combination
\eqn{\sum w_j \bar a(p_j)} of the means of the predicted compositions, or \code{list(type = "W1", pair = c("p1", "p2"))}, the
Wasserstein-1 distance between two predicted compositions.}
\item{reference}{standard reference distribution \eqn{q^\ast} on which contrasts are evaluated (default: the average of unit references over observed unit-period rows, with each row equally weighted; this does not weight by case totals. Usually supply an explicitly chosen reference). For a pooled population reference with unit population totals \eqn{N_i}, supply \eqn{q_k^\ast = \sum_i N_i q_{ik}/\sum_i N_i}. This differs from \code{colMeans(obj$q)} when populations differ; it is separate from equal unit/period weighting in marginal standardisation.}
\item{version}{\code{"marginal"}: for each posterior draw, average over \code{R} Monte Carlo replicates, each drawing a period with its
effects, the effects of a randomly chosen unit and a fresh cell-level effect; the same replicates are used at every point, so
contrasts are paired. \code{"typical"}: all period effects, unit and cell effects set to zero. The two differ because compositions are non-linear in the linear predictor. Marginal W1 is the mean of paired replicate distances, not the distance between averaged compositions. Typical W1 is averaged over periods.}
\item{R}{Monte Carlo replicates per draw (marginal version).}
\item{seed}{seed set before the Monte Carlo replicates.}
}
\value{\code{pc_contrasts}: an \eqn{S \times} (number of contrasts) matrix of posterior draws. \code{pc_regionyear}: a list of
\eqn{N \times S} matrices of \eqn{W_1} and \eqn{M} between each fitted unit-period composition and its unit's reference.}
\seealso{\code{\link{pc_tolerance}}}
