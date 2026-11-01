# Minimal example for popcomp: simulate case counts over 6 ordered age groups in 16 units (a 4 x 4 grid) and 3 periods, with a
# unit-level exposure x that shifts cases towards older ages, then fit the count model and summarise the contrast in mean age and
# the W1 distance between the compositions predicted at the 75th and 25th percentiles of x. Runs in about a minute.
library(popcomp)
set.seed(1)
a <- c(0.5, 7, 20, 35, 50, 70)                         # representative ages of 6 groups
K <- length(a); I <- 16; TT <- 3
bs <- pc_basis(a)
# reference (population) distributions: younger in some units, older in others
q <- t(vapply(seq_len(I), function(i) { w <- exp(-a / runif(1, 20, 35)); w / sum(w) }, numeric(K)))
x <- as.vector(scale(rnorm(I)))                         # standardised unit-level exposure
# 4 x 4 grid adjacency (rook)
g <- expand.grid(r = 1:4, c = 1:4); W <- 1 * (as.matrix(dist(g)) == 1)
# data-generating model: slope on the basis increases with x; cell-level overdispersion
beta1 <- 0.6; rows <- expand.grid(unit = seq_len(I), time = seq_len(TT))
n <- rpois(nrow(rows), 400) + 1
Y <- t(vapply(seq_len(nrow(rows)), function(r) {
  i <- rows$unit[r]
  eta <- log(q[i, ]) + bs$B[, 1] * (beta1 * x[i] + rnorm(1, 0, 0.2)) + rnorm(K, 0, 0.15)
  as.vector(rmultinom(1, n[r], exp(eta) / sum(exp(eta))))
}, numeric(K)))

obj <- pc_data(Y, unit = rows$unit, time = rows$time, q = q, X = cbind(x = x[rows$unit]), a = a)
fit <- pc_fit(obj, graph = W, spatial = "bym2")
dr  <- pc_draws(fit, obj, S = 500, seed = 2)

qstar <- colMeans(q)
pts <- list(hi = c(x = quantile(x, 0.75, names = FALSE)), lo = c(x = quantile(x, 0.25, names = FALSE)))
con <- list(dM = list(type = "M", weights = c(hi = 1, lo = -1)), W1 = list(type = "W1", pair = c("hi", "lo")))
draws_marg <- pc_contrasts(dr, obj, pts, con, reference = qstar, version = "marginal", R = 100, seed = 3)
draws_typ  <- pc_contrasts(dr, obj, pts, con, reference = qstar, version = "typical")

print(pc_tolerance(draws_marg[, "dM", drop = FALSE], delta = c(1, 2)))            # signed contrast in mean age (years)
print(pc_tolerance(draws_marg[, "W1", drop = FALSE], delta = c(1, 2), distance = TRUE))
print(pc_tolerance(draws_typ, delta = 2))
ry <- pc_regionyear(dr, obj)                                                      # unit-period W1 and M vs own reference
print(summary(rowMeans(ry$W1 - abs(ry$M))))                                       # W1 - |M| >= 0: the part a mean shift misses
