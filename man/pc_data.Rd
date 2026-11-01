\name{pc_data}
\alias{pc_data}
\title{Long-format data for the population-referenced count model}
\description{Long-format data for the population-referenced count model.}
\usage{
pc_data(Y, unit, time, q, X, a, basis = pc_basis(a),
        time_labels = NULL, offset_log_n = TRUE)
}
\arguments{
\item{Y}{\eqn{N \times K} matrix of finite non-negative integer counts, one row per unit-period with a positive total (rows with \eqn{n_{it} = 0} carry no
information on composition and must be removed).}
\item{unit, time}{integer indices (\code{1..I}, \code{1..T}) of each row of \code{Y}.}
\item{q}{\eqn{I \times K} matrix of non-negative reference distributions (rows sum to 1), one per unit.}
\item{X}{\eqn{N \times p} matrix of covariates with column names, one row per row of \code{Y}. Interactions must be supplied as
columns.}
\item{a}{representative values of the categories.}
\item{basis}{output of \code{\link{pc_basis}} with \code{degree = 2}.}
\item{time_labels}{labels of the periods, used in the names of the period effects (default \code{1..T}).}
\item{offset_log_n}{add \eqn{\log n_{it}} to the offset. It cancels in every composition but re-centres the free intercepts near zero,
which can help the mode search. The intercept has finite Gaussian precision, so invariance of the fitted posterior to this option is not exact.}
}
\details{A zero reference probability is a structural zero: the corresponding count must be zero. Such cells produce a negative-infinite offset and zero fitted probability. Positive counts at zero reference mass are incompatible with the model. Do not add pseudocounts automatically. Caller requirements include finite covariates, positive integer time indices, unique syntactically valid covariate and time labels, and finite strictly increasing category values. Counts, covariates, reference probabilities and positive integer time indices are checked. The caller must additionally ensure valid unique covariate and period labels. Zero-total rows must be removed explicitly.}
\value{An object of class \code{"pc_data"}: the long data frame (\code{data}) with one row per cell and the design information
used by \code{\link{pc_formula}}, \code{\link{pc_draws}} and \code{\link{pc_contrasts}}.}
\seealso{\code{\link{pc_fit}}}
