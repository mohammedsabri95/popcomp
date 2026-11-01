\name{pc_fit}
\alias{pc_fit}
\alias{pc_formula}
\title{Specify and fit the population-referenced count model with R-INLA}
\description{Specify and fit the population-referenced count model with R-INLA.}
\usage{
pc_formula(obj, graph = NULL, spatial = c("bym2", "iid"), higher = TRUE, eps = TRUE,
           prior_prec = c(1, 0.01), prior_phi = c(0.5, 2/3))
pc_fit(obj, formula = pc_formula(obj, ...), ...,
       control.fixed = list(mean = 0, prec = 0.1),
       control.compute = list(config = TRUE), num.threads = "1:1")
}
\arguments{
\item{obj}{output of \code{\link{pc_data}}.}
\item{graph}{\eqn{I \times I} symmetric 0/1 adjacency matrix of the units (needed for \code{spatial = "bym2"}).}
\item{spatial}{BYM2 or iid unit effect on the slope term.}
\item{higher}{include the iid unit effects \eqn{z_{ij}} on the higher-order basis columns (one shared precision).}
\item{eps}{include the cell-level overdispersion effect \eqn{\varepsilon_{itk}}. Whether to include this term is a modelling choice; omission does not represent cell-level extra-Poisson variation.}
\item{prior_prec}{penalised-complexity prior \eqn{P(\sigma > u) = \alpha} for all random-effect precisions, as \code{c(u, alpha)}.}
\item{prior_phi}{PC prior \eqn{P(\phi < u) = \alpha} for the BYM2 mixing parameter.}
\item{...}{Arguments passed to \code{pc_formula} when constructing the default formula.}
\item{formula}{a formula from \code{pc_formula}; arguments in \code{...} are passed to it.}
\item{control.fixed, control.compute, num.threads}{passed to \code{INLA::inla}. \code{config = TRUE} is required by
\code{\link{pc_draws}}.}
}
\details{The model is
\deqn{C_{itk} \sim \mathrm{Poisson}(\lambda_{itk}), \quad \log\lambda_{itk} = \alpha_{it} + \log q_{ik} + \log n_{it} +
\sum_{j=1}^{2} B_j(a_k)\{\gamma_{jt} + x_{it}'\beta_j + u_{ij}\} + \sum_{j\ge 3} B_j(a_k) z_{ij} + \varepsilon_{itk},}
with \eqn{\alpha_{it}} an iid Gaussian effect with fixed precision \eqn{\exp(-10)}, fixed effects
with default \eqn{N(0, 10)} priors and PC priors on the remaining random-effect precisions. The finite Gaussian intercept prior gives a diffuse-prior approximation to the conditional multinomial posterior, rather than exact Bayesian equivalence. Conditional count probabilities are multinomial, but integrating this finite prior leaves a composition-dependent factor. An improper flat log-intercept prior gives exact equivalence for positive totals when the remaining posterior is proper. Approximation accuracy requires assessment for the intended application. An age-invariant covariate effect would cancel from the composition,
which is why covariates act only through the basis. If the default fit fails, \code{pc_fit} retries once with empirical-Bayes
integration.}
\value{\code{pc_formula} returns a formula with its model environment and attribute \code{"pc"}. \code{pc_fit} returns an \code{inla} object with attribute \code{"pc"}.}
