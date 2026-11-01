# Posterior functionals: joint posterior draws, contrasts in the mean and in W1 between predicted compositions (marginal over the
# random effects or for a typical unit), unit-period W1 and M, and tolerance summaries.

.rows_by <- function(lat, tag) { r <- grep(paste0("^", tag, ":"), rownames(lat)); r[order(as.integer(sub(".*:", "", rownames(lat)[r])))] }

#' Joint posterior draws of the parameters that define the compositions
pc_draws <- function(fit, obj, S = 4000, seed = 1L, predictor = TRUE) {
  if (!requireNamespace("INLA", quietly = TRUE)) stop("pc_draws() needs the INLA package")
  info <- attr(fit, "pc"); if (is.null(info)) info <- list(higher = TRUE, eps = TRUE)
  fn <- obj$fixed; nc <- obj$N * obj$K; J <- ncol(obj$Bh)
  sel <- c(setNames(as.list(rep(1, length(fn))), fn), list(r1 = seq_len(obj$I), r2 = 0))
  if (info$higher) sel$rz <- nc + seq_len(obj$I * J)                 # the z part of the (Z z, z) block
  if (predictor) sel$Predictor <- 0
  # inla.posterior.sample() picks the hyperparameter configurations with R's sample(), so R's seed is set too: without it, two
  # calls with the same 'seed' give different draws
  set.seed(seed)
  smp <- INLA::inla.posterior.sample(S, fit, selection = sel, num.threads = "1:1", seed = seed)
  lat <- do.call(cbind, lapply(smp, `[[`, "latent")); hp <- do.call(rbind, lapply(smp, `[[`, "hyperpar"))
  fx <- lat[vapply(fn, function(nm) .rows_by(lat, nm), 1L), , drop = FALSE]; rownames(fx) <- fn
  out <- list(S = S,
              beta1 = fx[sprintf("B1_%s", obj$covs), , drop = FALSE], beta2 = fx[sprintf("B2_%s", obj$covs), , drop = FALSE],
              gamma1 = fx[sprintf("B1_y%s", obj$time_labels), , drop = FALSE], gamma2 = fx[sprintf("B2_y%s", obj$time_labels), , drop = FALSE],
              u1 = lat[.rows_by(lat, "r1")[seq_len(obj$I)], , drop = FALSE], u2 = lat[.rows_by(lat, "r2"), , drop = FALSE],
              z = if (info$higher) lat[.rows_by(lat, "rz"), , drop = FALSE] else NULL,
              sigma_eps = if (info$eps) 1 / sqrt(hp[, "Precision for cell"]) else NULL,
              predictor = if (predictor) lat[.rows_by(lat, "Predictor"), , drop = FALSE] else NULL)
  stopifnot(nrow(out$u1) == obj$I, nrow(out$u2) == obj$I, is.null(out$z) || nrow(out$z) == obj$I * J)
  out
}

#' Posterior draws of contrasts between predicted compositions
#' points:    named list of covariate vectors (named by obj$covs; missing names are 0)
#' contrasts: named list; each element is either list(type = "M", weights = c(point = weight, ...)), a linear combination of the
#'            mean of the predicted compositions, or list(type = "W1", pair = c(point1, point2)), the W1 distance between two of them
#' version:   "marginal" averages over R Monte Carlo replicates of (period, unit, cell effect); "typical" uses all period effects with
#'            unit and cell effects at zero
pc_contrasts <- function(draws, obj, points, contrasts, reference = NULL, version = c("marginal", "typical"), R = 200, seed = NULL) {
  version <- match.arg(version)
  valid_type <- vapply(contrasts, function(cn)
    is.list(cn) && is.character(cn$type) && length(cn$type) == 1L &&
      !is.na(cn$type) && cn$type %in% c("M", "W1"), logical(1))
  if (!all(valid_type))
    stop("each contrast type must be exactly 'M' or 'W1'", call. = FALSE)
  if (is.null(reference)) reference <- colSums(obj$q * tabulate(obj$unit, obj$I)) / obj$N   # default only; usually pass q*
  K <- obj$K; a <- obj$a; B <- obj$B; J <- ncol(obj$Bh); I <- obj$I; TT <- obj$TT
  xv <- lapply(points, function(x) { v <- setNames(numeric(length(obj$covs)), obj$covs); v[names(x)] <- x; v })
  if (!is.null(seed)) set.seed(seed)
  one <- function(s) {
    g1 <- draws$gamma1[, s]; g2 <- draws$gamma2[, s]
    if (version == "typical") {
      base <- matrix(log(reference), TT, K, byrow = TRUE) + outer(g1, B[, 1]) + outer(g2, B[, 2])
    } else {
      tr <- sample.int(TT, R, TRUE); ir <- sample.int(I, R, TRUE)
      eps <- if (is.null(draws$sigma_eps)) 0 else matrix(rnorm(R * K, 0, draws$sigma_eps[s]), R, K)
      base <- matrix(log(reference), R, K, byrow = TRUE) + outer(g1[tr] + draws$u1[ir, s], B[, 1]) + outer(g2[tr] + draws$u2[ir, s], B[, 2])
      if (!is.null(draws$z)) base <- base + (matrix(draws$z[, s], I, J, byrow = TRUE) %*% t(obj$Bh))[ir, , drop = FALSE]
      base <- base + eps
    }
    b1 <- draws$beta1[, s]; b2 <- draws$beta2[, s]
    P <- lapply(xv, function(x) pc_softmax(sweep(base, 2, B[, 1] * sum(b1 * x) + B[, 2] * sum(b2 * x), "+")))
    vapply(contrasts, function(cn) {
      if (cn$type == "M") sum(vapply(names(cn$weights), function(p) cn$weights[[p]] * mean(as.vector(P[[p]] %*% a)), 0))
      else mean(pc_W1(P[[cn$pair[1]]], P[[cn$pair[2]]], a))
    }, 0)
  }
  out <- t(vapply(seq_len(draws$S), one, numeric(length(contrasts))))
  if (length(contrasts) == 1) out <- t(out)
  colnames(out) <- names(contrasts); out
}

#' Unit-period W1 and M between the fitted composition and the unit's reference distribution (posterior draws)
pc_regionyear <- function(draws, obj) {
  if (is.null(draws$predictor)) stop("draws were taken with predictor = FALSE")
  Qm <- obj$q[obj$unit, , drop = FALSE]; W <- M <- matrix(NA_real_, obj$N, draws$S)
  for (s in seq_len(draws$S)) { P <- pc_softmax(matrix(draws$predictor[, s], obj$N, obj$K, byrow = TRUE))
    W[, s] <- pc_W1(P, Qm, obj$a); M[, s] <- pc_M(P, Qm, obj$a) }
  list(W1 = W, M = M)
}

#' Tolerance (near-null) summaries of posterior draws
#' For a signed contrast: P(|x| < delta); for a distance (distance = TRUE): P(x < delta). Also the median, the 95% equal-tailed
#' interval, P(x > 0) and the 97.5% posterior quantile of |x|. No test of an exact zero is made.
pc_tolerance <- function(x, delta = c(1, 2, 5), distance = FALSE) {
  x <- as.matrix(x); if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
  do.call(rbind, lapply(colnames(x), function(nm) { v <- x[, nm]
    out <- data.frame(quantity = nm, mean = mean(v), sd = sd(v), median = median(v), q025 = quantile(v, .025, names = FALSE),
                      q975 = quantile(v, .975, names = FALSE), p_gt0 = mean(v > 0), ub975_abs = quantile(abs(v), .975, names = FALSE))
    for (d in delta) out[[sprintf("p_within_%g", d)]] <- if (distance) mean(v < d) else mean(abs(v) < d)
    out }))
}
