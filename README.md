![popcomp research banner](assets/popcomp-banner.svg)

# popcomp

Bayesian count models and distributional contrasts for population-referenced compositions.

**Version 0.1.0** is available from the [popcomp GitHub repository](https://github.com/mohammedsabri95/popcomp).
The repository contains package source and simulated examples; application data, Stan archives
and internal validation reports are not included.

## Overview

`popcomp` studies how the composition of disease cases across ordered categories, such as age groups,
varies with exposures relative to a population reference composition. The reference helps distinguish
case composition from the age structure of the population in which cases arise. This is a model of
composition; it does not by itself estimate overall disease incidence or identify causal effects.

Small numbers of cases make observed proportions noisy, and zero cells are common. The model works
with counts and uses Bayesian posterior inference to represent uncertainty in covariate and random
effects. Rows with no cases carry no information about composition and must be removed explicitly.
Posterior draws propagate parameter uncertainty into signed mean-age differences, Wasserstein-1
(W1) distances, and probabilities that a contrast lies within a prespecified tolerance.

The motivating application concerns malaria and environmental exposures in the Brazilian Amazon.
The package includes simulated examples only; no empirical malaria findings are asserted here.

## Framework

![Conceptual popcomp statistical workflow](assets/popcomp-workflow.svg)

The banner and workflow are conceptual illustrations, not empirical results. The conditional
composition target motivates the implementation; the finite-prior Poisson approximation is
specified below. These illustrations do not establish causal effects or numerical accuracy.

## Main methodological components

- **Conditional count modelling:** a population offset and unit-period intercepts implement a
  Poisson representation motivated by the multinomial-Poisson transformation. The implemented
  intercept has fixed Gaussian precision `exp(-10)`. This finite prior yields a diffuse-prior
  approximation to the conditional multinomial posterior, not exact Bayesian equivalence.
  Conditional count probabilities are multinomial; integrating a finite Gaussian intercept prior
  leaves an additional composition-dependent factor. See `help("pc_fit")`.
- **Age-basis exposure effects:** covariates act on two orthonormal polynomial basis columns,
  with additional unit effects on higher-order columns. Representative ages must be strictly increasing.
- **Spatial dependence:** BYM2 effects represent spatial variation in the first basis coefficient;
  iid alternatives and additional iid effects are available.
- **Posterior inference:** R-INLA fits the count model and supplies joint posterior draws.
  A failed default fit triggers one retry using empirical-Bayes integration; inspect its message.
- **Distributional contrasts:** W1 measures displacement on the supplied ordered support in the same
  units as that support. A signed mean-age contrast indicates direction; W1 is non-negative and
  can detect shape differences with no change in mean.
- **Standardisation:** typical-unit contrasts set unit and cell effects to zero and average over
  periods. Marginal contrasts average paired replicate functionals over sampled periods, units
  and cell effects. Marginal W1 is an average of paired distances, not W1 between averaged compositions.
- **Tolerance inference:** `pc_tolerance()` reports posterior means, standard deviations, equal-tailed
  95% intervals and `P(abs(contrast) < delta)` or `P(W1 < delta)`. Choose meaningful tolerances in advance.

## Installation

Install into a **new, suitable R library**. Validation used local source installation; the
GitHub installation command below has not been independently checked.

```r
# CRAN dependency and optional GitHub installer (stats ships with R).
install.packages(c("Matrix", "remotes"), repos = "https://cloud.r-project.org")

# INLA is needed for model fitting and posterior sampling.
install.packages("INLA", repos = c(
  CRAN = "https://cloud.r-project.org",
  INLA = "https://inla.r-inla-download.org/R/stable"
))

# Install popcomp from GitHub:
remotes::install_github("mohammedsabri95/popcomp", upgrade = "never")
```

These INLA repository instructions follow the [official installation guidance](https://www.r-inla.org/download/).
The stable repository may supply a different version over time. Reproducing the validation estimates
requires the recorded INLA 25.10.19 environment, not merely whichever stable release is newest.
INLA remains in `Suggests` because the transport and data utilities work without it and both
INLA-dependent functions guard their use with `requireNamespace()`.

For the existing workstation environment, run all R commands through conda:

```sh
conda run -n malaria_r45 R CMD INSTALL /path/to/popcomp
conda run -n malaria_r45 Rscript -e 'library(popcomp); packageVersion("popcomp")'
```

## Getting started

Run the small, fully simulated workflow:

```r
library(popcomp)
source(system.file("examples", "getting_started.R", package = "popcomp"))
```

The script uses nine units on a 3-by-3 adjacency grid, two periods and six age groups. It sets
simulation seed 1, posterior seed 2 and marginal standardisation seed 3. It generates reference
compositions and counts, then executes:

```r
# Y, rows, q, x, a and W are generated in the script above.
obj <- pc_data(Y, unit = rows$unit, time = rows$time, q = q,
               X = cbind(x = x[rows$unit]), a = a)
fit <- pc_fit(obj, graph = W, spatial = "bym2", num.threads = "1:1")
dr <- pc_draws(fit, obj, S = 100, seed = 2)
pts <- list(hi = c(x = quantile(x, .75, names = FALSE)),
            lo = c(x = quantile(x, .25, names = FALSE)))
con <- list(dM = list(type = "M", weights = c(hi = 1, lo = -1)),
            W1 = list(type = "W1", pair = c("hi", "lo")))
post <- pc_contrasts(dr, obj, pts, con, reference = colMeans(q),
                     version = "marginal", R = 30, seed = 3)
pc_tolerance(post[, "dM", drop = FALSE], delta = c(1, 2))
pc_tolerance(post[, "W1", drop = FALSE], delta = c(1, 2), distance = TRUE)
```

Here `dM` is the mean-age difference between the high and low exposure predictions in years.
`W1` is their distributional displacement in years. For example, `p_within_2` is a posterior
probability of a displacement smaller than two years, not a frequentist p-value. The small
posterior sample and Monte Carlo sizes are for illustration; assess simulation error for research.
The original, larger example remains available as `examples/minimal_example.R` in the installed package.
Both simulated examples use `colMeans(q)`, an equally weighted illustrative reference. This is not the
application's pooled population reference. For eligible regions with population totals `population_total`,
supply `reference = colSums(q * population_total) / sum(population_total)`. The population reference
weights and the equal unit/period weights in marginal standardisation serve different purposes.

Supply finite non-negative integer counts, non-negative normalised references with zero cases wherever reference mass is zero, finite covariates,
valid positive integer indices and strictly increasing ages. These count, covariate, time and W1
contracts are checked; contrast types must be `"M"` or `"W1"`. Supply unique syntactically valid
covariate and period labels. The default contrast reference weights
each observed unit-period row equally, not each case. Pass `reference` explicitly for your estimand.

## Main functions

| Export | Purpose | Main input | Main output |
|---|---|---|---|
| `pc_basis` | Construct an orthonormal age basis | Ordered representative ages | `a`, `B`, `Bh` list |
| `pc_data` | Prepare count model data | Counts, indices, population reference, covariates | `pc_data` object |
| `pc_formula` | Specify random effects and priors | Data object, adjacency and model options | Formula with model environment |
| `pc_fit` | Fit with INLA | Data object and formula/options | INLA fit with model attributes |
| `pc_draws` | Sample composition parameters | Fit, data object, sample size, seed | Parameter and optional predictor arrays |
| `pc_contrasts` | Standardised exposure contrasts | Draws, points, contrast definitions, reference | Draw-by-contrast matrix |
| `pc_regionyear` | Compare fitted and reference compositions | Predictor draws and data object | Unit-period-by-draw W1 and mean matrices |
| `pc_tolerance` | Summarise posterior contrasts | Draw matrix, tolerances | Summary data frame |
| `pc_W1` | Ordered-support W1 distance | Probability rows, reference, ages | Row-wise distances |
| `pc_M` | Signed mean displacement | Probability rows, reference, ages | Row-wise mean differences |
| `pc_softmax` | Convert logits to compositions (negative infinity permits structural zeros) | Matrix with one composition per row | Probability matrix |

## Citation

Use `citation("popcomp")` for the software citation. The draft citation is Mohammed Sabri (2026),
*popcomp: Bayesian Count Models for Population-Referenced Compositions*, R package version 0.1.0,
available from the [GitHub repository](https://github.com/mohammedsabri95/popcomp). A version DOI will be added only after a citable archive record exists.

The associated manuscript, *Bayesian count models and distributional contrasts for population-referenced
age compositions, with an application to malaria and deforestation in the Brazilian Amazon*,
is in preparation for **Environmetrics**. It is not described as accepted or published.

## Reproducibility

Validation uses Ubuntu 22.04.5 LTS, R 4.5.2 and INLA 25.10.19 in conda environment `malaria_r45`.
INLA fits and posterior sampling use `num.threads = "1:1"`; BLAS/OpenMP thread counts are set to one
in the validation commands. `pc_draws()` sets both R's random seed and INLA's sampler seed.
INLA internal numerical auto-tuning can also vary the fitted approximation despite one-thread
settings. Disabling it improves repeatability in the diagnostic example but does not by itself
establish accuracy. Numerical configuration, sampling error and prior sensitivity must be assessed
separately. Identical seeds do not promise identical results across configurations, platforms or versions.

Raw validation artifacts, application outputs and local review reports are excluded from the
public repository and R source archive. The package examples generate simulated data only.
Manual `.Rd` help is maintained without roxygen2. No vignette dependencies are required.

## Authorship and licence

Mohammed Sabri, King Abdullah University of Science and Technology (KAUST), Saudi Arabia.
Author and maintainer: Mohammed Sabri; email: mohammed.sabri@kaust.edu.sa. No ORCID is asserted.
popcomp v0.1.0 is distributed under the MIT licence. Copyright (c) 2026 Mohammed Sabri.
See [LICENSE](LICENSE) for the R package attribution and [LICENSE.md](LICENSE.md) for the complete MIT licence.
The package declaration is `License: MIT + file LICENSE`.
