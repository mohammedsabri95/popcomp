# Distances between distributions on ordered categories with representative values a_1 < ... < a_K.

#' Row-wise softmax (numerically stable)
pc_softmax <- function(eta) {
  eta <- as.matrix(eta)
  m <- eta[, 1]; for (k in seq_len(ncol(eta))[-1]) m <- pmax(m, eta[, k])
  e <- exp(eta - m); e / rowSums(e)
}

.as_rows <- function(Q, P) if (is.null(dim(Q))) matrix(Q, nrow(P), length(Q), byrow = TRUE) else Q

#' Signed mean displacement M(p, q) = sum_k a_k (p_k - q_k), row-wise
pc_M <- function(P, Q, a) {
  P <- as.matrix(P); Q <- .as_rows(Q, P)
  as.vector((P - Q) %*% a)
}

#' Wasserstein-1 distance W1(p, q) = sum_{k<K} |P_k - Q_k| (a_{k+1} - a_k), row-wise (P_k, Q_k cumulative)
pc_W1 <- function(P, Q, a) {
  P <- as.matrix(P)
  K <- ncol(P)
  if (!is.numeric(a) || is.complex(a) || length(a) != K || K < 1L ||
      any(!is.finite(a)) || any(diff(a) <= 0))
    stop("a must contain one finite, strictly increasing value per probability column", call. = FALSE)
  if (is.null(dim(Q)) && length(Q) != K)
    stop("Q must have one probability per support value", call. = FALSE)
  Q <- .as_rows(Q, P)
  if (!is.numeric(P) || is.complex(P) || !is.numeric(Q) || is.complex(Q) ||
      !identical(dim(P), dim(Q)) || any(!is.finite(P)) || any(!is.finite(Q)) ||
      any(P < 0) || any(Q < 0) || any(abs(rowSums(P) - 1) > 1e-8) ||
      any(abs(rowSums(Q) - 1) > 1e-8))
    stop("P and Q must be matching finite, non-negative probability rows summing to one", call. = FALSE)
  C <- ((P - Q) %*% upper.tri(diag(K), diag = TRUE))[, -K, drop = FALSE]
  as.vector(abs(C) %*% diff(a))
}
