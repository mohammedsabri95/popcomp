# Model set-up and fitting: multinomial counts over K ordered categories through the multinomial-Poisson transformation,
#   C_itk ~ Poisson(lambda_itk),
#   log lambda_itk = alpha_it + log q_ik + log n_it + sum_{j=1,2} B_j(a_k) {gamma_jt + x_it' beta_j + u_ij}
#                    + sum_{j>=3} B_j(a_k) z_ij + eps_itk,
# with a diffuse Gaussian intercept alpha_it per unit-period (an approximation to conditional inference),
# BYM2 (or iid) u_i1, iid u_i2, iid z_ij with one shared
# precision, and an iid cell-level effect eps_itk (logistic-normal overdispersion).

#' Orthonormal polynomial basis on the category values
pc_basis <- function(a, degree = 2) {
  K <- length(a); stopifnot(K >= degree + 2, !is.unsorted(a, strictly = TRUE))
  Bf <- unclass(poly(a, K - 1))
  B <- Bf[, seq_len(degree), drop = FALSE]; dimnames(B) <- NULL
  Bh <- Bf[, -seq_len(degree), drop = FALSE]; dimnames(Bh) <- NULL
  list(a = a, B = B, Bh = Bh)
}

#' Long-format data for the count model
pc_data <- function(Y, unit, time, q, X, a, basis = pc_basis(a), time_labels = NULL, offset_log_n = TRUE) {
  Y <- as.matrix(Y); X <- as.matrix(X); q <- as.matrix(q)
  if (!is.numeric(Y) || is.complex(Y) || any(!is.finite(Y)) ||
      any(Y < 0) || any(Y != floor(Y)))
    stop("Y must contain finite non-negative integer counts", call. = FALSE)
  if (!is.numeric(X) || is.complex(X) || any(!is.finite(X)))
    stop("X must contain finite numeric covariates", call. = FALSE)
  if (!is.numeric(q) || is.complex(q) || any(!is.finite(q)) || any(q < 0))
    stop("q must contain finite non-negative probabilities", call. = FALSE)
  if (!is.numeric(time) || is.complex(time) || any(!is.finite(time)) ||
      any(time < 1) || any(time != floor(time)))
    stop("time must contain positive integer indices", call. = FALSE)
  N <- nrow(Y); K <- ncol(Y); n <- rowSums(Y)
  stopifnot(length(unit) == N, length(time) == N, nrow(X) == N, ncol(q) == K, length(a) == K, !is.null(colnames(X)),
            all(abs(rowSums(q) - 1) < 1e-8), all(unit %in% seq_len(nrow(q))), ncol(basis$B) == 2)
  if (any(Y > 0 & q[unit, , drop = FALSE] == 0))
    stop("positive case count at a zero reference probability is incompatible with the model", call. = FALSE)
  if (any(n == 0)) stop("rows with n = 0 carry no information on composition; remove them before calling pc_data()")
  I <- nrow(q); TT <- max(time); if (is.null(time_labels)) time_labels <- seq_len(TT)
  covs <- colnames(X); B <- basis$B
  it <- rep(seq_len(N), each = K); k <- rep(seq_len(K), N); i <- unit[it]; t <- time[it]
  d <- data.frame(y = as.vector(t(Y)), it = it, off = log(q[cbind(i, k)]) + if (offset_log_n) log(n[it]) else 0,
                  r1 = i, r2 = i, w1 = B[k, 1], w2 = B[k, 2], rz = seq_along(it), cell = seq_along(it))
  for (j in 1:2) {
    for (tt in seq_len(TT)) d[[sprintf("B%d_y%s", j, time_labels[tt])]] <- B[k, j] * (t == tt)
    for (v in covs) d[[sprintf("B%d_%s", j, v)]] <- B[k, j] * X[it, v]
  }
  fixed <- c(outer(c(sprintf("y%s", time_labels), covs), 1:2, function(x, j) sprintf("B%d_%s", j, x)))
  structure(list(data = d, I = I, TT = TT, K = K, N = N, covs = covs, a = a, B = B, Bh = basis$Bh, q = q, unit = unit, time = time,
                 n = n, fixed = fixed, time_labels = time_labels), class = "pc_data")
}

.region_Z <- function(obj) {
  J <- ncol(obj$Bh); d <- obj$data; k <- (d$cell - 1L) %% obj$K + 1L
  Matrix::sparseMatrix(i = rep(seq_len(nrow(d)), J), j = (rep(d$r1, J) - 1L) * J + rep(seq_len(J), each = nrow(d)),
                       x = as.vector(obj$Bh[k, ]), dims = c(nrow(d), obj$I * J))
}

#' INLA formula of the count model
pc_formula <- function(obj, graph = NULL, spatial = c("bym2", "iid"), higher = TRUE, eps = TRUE,
                       prior_prec = c(1, 0.01), prior_phi = c(0.5, 2 / 3)) {
  spatial <- match.arg(spatial)
  if (spatial == "bym2" && is.null(graph)) stop("a BYM2 slope effect needs an adjacency matrix 'graph'")
  PC1 <- list(prec = list(prior = "pc.prec", param = prior_prec)); PCPHI <- list(prior = "pc", param = prior_phi)
  r1 <- if (spatial == "bym2") 'f(r1, w1, model = "bym2", graph = W, scale.model = TRUE, hyper = list(prec = PC1$prec, phi = PCPHI))'
        else 'f(r1, w1, model = "iid", hyper = PC1, values = seq_len(I))'
  terms <- c("-1", paste(obj$fixed, collapse = " + "), 'f(it, model = "iid", hyper = list(prec = list(initial = -10, fixed = TRUE)))', r1,
             'f(r2, w2, model = "iid", hyper = PC1, values = seq_len(I))', if (higher) 'f(rz, model = "z", Z = Z, hyper = PC1)',
             if (eps) 'f(cell, model = "iid", hyper = PC1)')
  f <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
  environment(f) <- list2env(list(W = graph, PC1 = PC1, PCPHI = PCPHI, I = obj$I, Z = if (higher) .region_Z(obj)), parent = globalenv())
  attr(f, "pc") <- list(higher = higher, eps = eps, spatial = spatial)
  f
}

#' Fit the count model with R-INLA
pc_fit <- function(obj, formula = pc_formula(obj, ...), ..., control.fixed = list(mean = 0, prec = 0.1),
                   control.compute = list(config = TRUE), num.threads = "1:1") {
  if (!requireNamespace("INLA", quietly = TRUE)) stop("pc_fit() needs the INLA package (https://www.r-inla.org)")
  d <- obj$data; d$.pc_offset <- d$off
  go <- function(ci) INLA::inla(formula, data = d, family = "poisson", offset = .pc_offset, control.fixed = control.fixed,
                                control.compute = control.compute, control.inla = ci, num.threads = num.threads)
  fit <- tryCatch(go(list()), error = function(e) NULL)
  if (is.null(fit) || anyNA(fit$summary.fixed$mean)) {             # one retry with empirical-Bayes integration
    message("pc_fit: default fit failed; retrying with int.strategy = 'eb'")
    fit <- go(list(int.strategy = "eb"))
  }
  attr(fit, "pc") <- attr(formula, "pc")
  fit
}
