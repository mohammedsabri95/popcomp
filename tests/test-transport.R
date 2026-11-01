# Base-R tests (no testthat dependency): W1 and M against direct definitions; W1 >= |M|; equality without crossing.
library(popcomp)
set.seed(11)
for (rep in 1:200) {
  K <- sample(3:12, 1); a <- sort(runif(K, 0, 90))
  p <- prop.table(runif(K)); q <- prop.table(runif(K))
  w1 <- sum(abs(cumsum(p - q))[-K] * diff(a)); m <- sum(a * (p - q))
  stopifnot(abs(pc_W1(rbind(p), q, a) - w1) < 1e-12, abs(pc_M(rbind(p), q, a) - m) < 1e-12, pc_W1(rbind(p), q, a) >= abs(m) - 1e-12)
}
a <- c(1, 2, 4); p <- c(0.2, 0.3, 0.5); q <- c(0.3, 0.3, 0.4)                  # p stochastically larger: no crossing
stopifnot(abs(pc_W1(rbind(p), q, a) - pc_M(rbind(p), q, a)) < 1e-12)
stopifnot(abs(rowSums(pc_softmax(matrix(rnorm(20), 4))) - 1) < 1e-12)
b <- pc_basis(c(0.5, 2.5, 9.5, 19.5, 29.5)); stopifnot(max(abs(crossprod(cbind(b$B, b$Bh)) - diag(4))) < 1e-12)
tol <- pc_tolerance(cbind(E = c(-3, -1, 0.5, 1.5, 2.5)), delta = 2); stopifnot(tol$p_within_2 == 0.6)
cat("popcomp transport tests passed\n")
