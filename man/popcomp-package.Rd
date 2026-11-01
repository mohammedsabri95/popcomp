\name{popcomp-package}
\alias{popcomp-package}
\alias{popcomp}
\docType{package}
\title{Population-referenced count models for case compositions}
\description{
Models the composition of cases over \eqn{K} ordered categories (for example age groups with representative values
\eqn{a_1 < \dots < a_K}) in units \eqn{i} and periods \eqn{t}, relative to a reference distribution \eqn{q_i} (for example the
population age structure), as a function of exposures \eqn{x_{it}}. Counts are fitted through the multinomial--Poisson
transformation in R-INLA, with a diffuse Gaussian intercept per unit-period (a diffuse-prior approximation to conditional multinomial inference, not exact Bayesian equivalence; see \link{pc_fit}), \eqn{\log q_{ik}} as
an offset, covariates acting on an orthonormal basis over the categories, a BYM2 (or iid) unit effect on the slope term, iid unit
effects on the curvature and higher-order terms, and a cell-level effect for overdispersion.
Contrasts in the mean and Wasserstein-1 distances between predicted compositions are computed as posterior functionals
(\code{\link{pc_contrasts}}) and summarised with tolerance probabilities (\code{\link{pc_tolerance}}).

Workflow: \code{\link{pc_basis}} -> \code{\link{pc_data}} -> \code{\link{pc_formula}} / \code{\link{pc_fit}} -> \code{\link{pc_draws}}
-> \code{\link{pc_contrasts}}, \code{\link{pc_regionyear}} -> \code{\link{pc_tolerance}}. A runnable example is in
\code{system.file("examples", "minimal_example.R", package = "popcomp")}.
}
