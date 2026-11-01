\name{pc_draws}
\alias{pc_draws}
\title{Joint posterior draws of the composition parameters}
\description{Joint posterior draws of the composition parameters.}
\usage{pc_draws(fit, obj, S = 4000, seed = 1L, predictor = TRUE)}
\arguments{
\item{fit}{output of \code{\link{pc_fit}} (fitted with \code{config = TRUE}).}
\item{obj}{the \code{\link{pc_data}} object used for the fit.}
\item{S}{number of joint posterior draws (\code{INLA::inla.posterior.sample}).}
\item{seed}{seed passed to \code{inla.posterior.sample}; R's seed is also set to it, because \code{inla.posterior.sample} picks the
hyperparameter configurations with R's \code{sample()}. Reproducibility has been checked within the recorded R/INLA environment with one-thread sampling; equality across platforms or versions is not guaranteed.}
\item{predictor}{also keep the linear predictor (needed by \code{\link{pc_regionyear}}).}
}
\value{A list with the draws of \eqn{\beta_1, \beta_2}, the period effects \eqn{\gamma_{1t}, \gamma_{2t}}, the unit effects
\eqn{u_{i1}, u_{i2}, z_{ij}}, the cell-level standard deviation and (optionally) the linear predictor, The coefficient, period, unit and predictor arrays have one column per draw. The list also contains scalar \code{S} and vector \code{sigma_eps} (length \code{S}), or \code{NULL} when cell effects are excluded; \code{z} and \code{predictor} can likewise be \code{NULL}.}
