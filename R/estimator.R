# ============================================================
# GMM Estimator for Endogenous Peer Effects
# ============================================================
# Two-step GMM estimator for the endogenous peer-effect parameter in a
# peer-group random-assignment (urn) design with optional covariates.
# Heteroskedasticity-robust variance estimator with finite-sample bias
# correction. Pure base R; no external dependencies.
#
# Public API: peer_cra(). All other functions are internal helpers.
#
# The estimator exploits the block-diagonal structure of the spatial
# operators (W, I - lambda W, and its inverse) to run in O(sum_r n_r)
# time per moment evaluation. No n x n matrices are formed.

# ------------------------------------------------------------
# Closed-form p I^* + q J^* operator
# ------------------------------------------------------------

#' Apply the group-diagonal operator diag_g{p_g I_g^* + q_g J_g^*} to a
#' vector (or column-wise to a matrix).
#'
#' For each group block of length m_g, split v_g into the group mean and
#' the demeaned residual; scale the mean by q_g and the demeaned residual
#' by p_g; recombine. O(length(v)) total.
#'
#' @noRd
apply_pq_group_diag <- function(v, group_sizes, p_by_group, q_by_group) {
  if (length(p_by_group) == 1L) p_by_group <- rep(p_by_group, length(group_sizes))
  if (length(q_by_group) == 1L) q_by_group <- rep(q_by_group, length(group_sizes))
  if (is.null(dim(v))) {
    out <- numeric(length(v))
    off <- 0L
    for (g in seq_along(group_sizes)) {
      m_g <- group_sizes[g]
      idx <- (off + 1L):(off + m_g)
      v_g <- v[idx]
      mean_g <- sum(v_g) / m_g
      out[idx] <- p_by_group[g] * (v_g - mean_g) + q_by_group[g] * mean_g
      off <- off + m_g
    }
    out
  } else {
    out <- matrix(0, nrow(v), ncol(v))
    off <- 0L
    for (g in seq_along(group_sizes)) {
      m_g <- group_sizes[g]
      idx <- (off + 1L):(off + m_g)
      v_g <- v[idx, , drop = FALSE]
      mean_g <- .colSums(v_g, m_g, ncol(v_g)) / m_g
      dem <- v_g - rep(mean_g, each = m_g)
      out[idx, ] <- p_by_group[g] * dem +
        rep(q_by_group[g] * mean_g, each = m_g)
      off <- off + m_g
    }
    out
  }
}

#' Apply I_r^* (subtract urn mean) to a vector or matrix of urn-length.
#' @noRd
apply_Istar_urn <- function(v) {
  if (is.null(dim(v))) {
    v - sum(v) / length(v)
  } else {
    nr <- nrow(v); nc <- ncol(v)
    v - rep(.colSums(v, nr, nc) / nr, each = nr)
  }
}

#' Apply W_r in closed form: per group, p = -1/(m_g - 1), q = 1.
#' @noRd
apply_W_urn <- function(v, group_sizes) {
  apply_pq_group_diag(v, group_sizes, -1 / (group_sizes - 1), 1)
}

#' Apply (I_r - lambda W_r)^{-1} in closed form: per group,
#' p = (m_g - 1) / (m_g - 1 + lambda), q = 1 / (1 - lambda).
#' @noRd
apply_IlambdaW_inv_urn <- function(lambda, v, group_sizes) {
  apply_pq_group_diag(v, group_sizes,
    (group_sizes - 1) / (group_sizes - 1 + lambda),
    rep(1 / (1 - lambda), length(group_sizes)))
}

# ------------------------------------------------------------
# Whole-sample group and urn operators
# ------------------------------------------------------------

#' Sum a vector or the columns of a matrix by a consecutive integer index.
#'
#' `build_urn_info()` constructs indices in first-occurrence order, so
#' `rowsum(reorder = FALSE)` returns rows in the same order as the associated
#' size vectors. Keeping this primitive in one place prevents accidental
#' reordering in the matrix-free formulas below.
#' @noRd
sum_by_index <- function(x, index) {
  out <- rowsum(x, group = index, reorder = FALSE)
  if (is.null(dim(x))) {
    as.numeric(out)
  } else {
    dimnames(out) <- list(NULL, colnames(x))
    out
  }
}

#' Retrieve the whole-sample index and size mappings from `urn_info`.
#' @noRd
sample_map <- function(urn_info) {
  map <- attr(urn_info, "sample_map", exact = TRUE)
  if (is.null(map))
    stop("Internal error: urn_info is missing its sample map.", call. = FALSE)
  map
}

#' Apply a group-block operator to the full estimation sample.
#'
#' Each block is `p_g I_g^* + q_g J_g^*`. Unlike
#' `apply_pq_group_diag()`, this version uses the global group index and
#' therefore avoids one R loop per urn and group.
#' @noRd
apply_group_operator_sample <- function(v, urn_info,
                                        p_by_group, q_by_group) {
  map <- sample_map(urn_info)
  n_groups <- length(map$group_sizes)
  if (length(p_by_group) == 1L) p_by_group <- rep(p_by_group, n_groups)
  if (length(q_by_group) == 1L) q_by_group <- rep(q_by_group, n_groups)

  group_mean <- sum_by_index(v, map$group_index) / map$group_sizes
  if (is.null(dim(v))) {
    mean_obs <- group_mean[map$group_index]
    p_by_group[map$group_index] * (v - mean_obs) +
      q_by_group[map$group_index] * mean_obs
  } else {
    mean_obs <- group_mean[map$group_index, , drop = FALSE]
    p_obs <- p_by_group[map$group_index]
    q_obs <- q_by_group[map$group_index]
    out <- p_obs * (v - mean_obs) + q_obs * mean_obs
    dimnames(out) <- dimnames(v)
    out
  }
}

#' Apply W to a vector or matrix over the full estimation sample.
#' @noRd
apply_W_sample <- function(v, urn_info) {
  map <- sample_map(urn_info)
  apply_group_operator_sample(
    v, urn_info,
    p_by_group = -1 / (map$group_sizes - 1),
    q_by_group = 1
  )
}

#' Apply I^* to a vector or matrix over the full estimation sample.
#' @noRd
apply_Istar_sample <- function(v, urn_info) {
  map <- sample_map(urn_info)
  urn_mean <- sum_by_index(v, map$urn_index) / map$urn_sizes
  if (is.null(dim(v))) {
    v - urn_mean[map$urn_index]
  } else {
    out <- v - urn_mean[map$urn_index, , drop = FALSE]
    dimnames(out) <- dimnames(v)
    out
  }
}

# ------------------------------------------------------------
# Numerical helpers
# ------------------------------------------------------------

#' Solve `M x = b` (or invert `M` if `b` is NULL) with a truncated-SVD
#' pseudo-inverse fallback when `M` is singular or near-singular.
#' @noRd
solve_robust <- function(M, b = NULL, tol = sqrt(.Machine$double.eps)) {
  tryCatch({
    if (is.null(b)) solve(M) else solve(M, b)
  }, error = function(e) {
    s <- svd(M)
    pos <- s$d > max(tol * s$d[1], 0)
    if (is.null(b)) {
      s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
    } else {
      UTb <- crossprod(s$u[, pos, drop = FALSE], b)
      out <- s$v[, pos, drop = FALSE] %*% (UTb / s$d[pos])
      if (is.null(dim(b))) as.numeric(out) else out
    }
  })
}

#' Repair a symmetric matrix to be positive definite by flooring only its
#' eigenvalues. This preserves the raw moment-variance matrix whenever it is
#' already well behaved, and changes only the problematic eigendirections.
#' @noRd
repair_pd_matrix <- function(B, eig_floor_rel = 1e-10) {
  if (nrow(B) == 0L || ncol(B) == 0L) {
    return(list(matrix = B, projected = FALSE, min_eig_raw = NA_real_,
                n_eig_repaired = 0L, eig_floor = NA_real_))
  }
  B_sym <- (B + t(B)) / 2
  ev <- eigen(B_sym, symmetric = TRUE)
  d <- ev$values
  scale <- max(mean(abs(diag(B_sym))), 1)
  eig_floor <- eig_floor_rel * scale
  min_eig_raw <- min(d)
  bad <- d < eig_floor
  if (any(bad)) {
    d[bad] <- eig_floor
    B_sym <- ev$vectors %*% (d * t(ev$vectors))
    B_sym <- (B_sym + t(B_sym)) / 2
  }
  list(matrix = B_sym, projected = any(bad), min_eig_raw = min_eig_raw,
       n_eig_repaired = sum(bad), eig_floor = eig_floor)
}

#' Use a candidate moment-variance matrix only when it is positive definite.
#'
#' The raw bias-corrected estimator can be singular in a finite sample even
#' though it converges to a positive-definite limit. On that event, use the
#' positive-definite first-step matrix instead of a pseudoinverse that would
#' assign zero weight to a singular moment direction.
#' @noRd
factor_pd_matrix <- function(V) {
  V_sym <- (V + t(V)) / 2
  R <- tryCatch(chol(V_sym), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  list(matrix = V_sym, inverse = chol2inv(R))
}

#' @noRd
select_pd_or_first <- function(V_candidate, V_first) {
  if (!identical(dim(V_candidate), dim(V_first)))
    stop("Candidate and first-step variance matrices must have equal dimensions.")

  selected <- factor_pd_matrix(V_candidate)
  used_first <- is.null(selected)
  if (used_first) selected <- factor_pd_matrix(V_first)
  if (is.null(selected))
    stop("The first-step moment-variance matrix is not positive definite.")

  selected$used_first <- used_first
  selected
}

#' Select the reported covariance under the four fallback branches.
#'
#' Point estimation always uses the stage-two weight supplied through `Xi2`.
#' A positive-definite final residual covariance uses the maintained Option B
#' calculation. If only the final residual covariance fails, `Xi2` is the
#' inverse of the positive-definite stage-one residual covariance and supplies
#' the inverse-bread repair. If both residual covariances fail, inference is
#' unavailable and no covariance is manufactured from the initial matrix.
#' @noRd
compute_final_covariance <- function(D_hat, V_hat_s1, V_hat_final, Xi2,
                                     stage2_weight_used_first, n,
                                     report_transform, coef_labels = NULL) {
  k_theta <- ncol(D_hat)
  if (k_theta < 1L || nrow(D_hat) != nrow(V_hat_final) ||
      !identical(dim(V_hat_final), dim(V_hat_s1)) ||
      !identical(dim(V_hat_final), dim(Xi2)) ||
      !identical(dim(report_transform), c(k_theta, k_theta))) {
    stop("Final-covariance inputs have incompatible dimensions.")
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0)
    stop("n must be one positive finite number.")
  if (!is.logical(stage2_weight_used_first) ||
      length(stage2_weight_used_first) != 1L ||
      is.na(stage2_weight_used_first)) {
    stop("stage2_weight_used_first must be TRUE or FALSE.")
  }

  final_factor <- factor_pd_matrix(V_hat_final)
  final_vcov_nonpd <- is.null(final_factor)

  if (!final_vcov_nonpd) {
    bread <- crossprod(D_hat, final_factor$inverse %*% D_hat)
    Sigma_internal <- solve_robust(bread) / n
    vcov_status <- if (stage2_weight_used_first) {
      "stage2_weight_fallback_final_option_b"
    } else {
      "final_option_b"
    }
    inference_available <- TRUE
  } else if (!stage2_weight_used_first) {
    if (is.null(factor_pd_matrix(V_hat_s1))) {
      stop("The final-only covariance branch requires a positive-definite stage-one residual covariance.")
    }
    bread <- crossprod(D_hat, Xi2 %*% D_hat)
    Sigma_internal <- solve_robust(bread) / n
    vcov_status <- "stage1_residual_final_fallback"
    inference_available <- TRUE
  } else {
    Sigma_internal <- matrix(NA_real_, k_theta, k_theta)
    vcov_status <- "unavailable_both_fallback"
    inference_available <- FALSE
  }

  if (inference_available)
    Sigma_internal <- (Sigma_internal + t(Sigma_internal)) / 2
  Sigma_reported <- report_transform %*% Sigma_internal %*%
    t(report_transform)
  if (inference_available)
    Sigma_reported <- (Sigma_reported + t(Sigma_reported)) / 2
  if (!is.null(coef_labels))
    dimnames(Sigma_reported) <- list(coef_labels, coef_labels)

  list(
    Sigma_internal = Sigma_internal,
    Sigma_reported = Sigma_reported,
    se = sqrt(diag(Sigma_reported)),
    final_vcov_nonpd = final_vcov_nonpd,
    inference_available = inference_available,
    vcov_status = vcov_status
  )
}

#' Stable Euclidean norm that avoids overflow and underflow.
#' @noRd
vector_l2_norm <- function(x) {
  scale <- max(abs(x))
  if (scale == 0 || !is.finite(scale)) return(scale)
  scale * sqrt(sum((x / scale)^2))
}

#' Normalize nonzero matrix columns to unit Euclidean norm.
#' @noRd
normalize_columns_for_rank <- function(M) {
  if (ncol(M) == 0L) return(M)
  out <- M
  norms <- vapply(seq_len(ncol(M)),
                  function(j) vector_l2_norm(M[, j]), numeric(1))
  nonzero <- is.finite(norms) & norms > 0
  if (any(nonzero)) {
    out[, nonzero] <- sweep(out[, nonzero, drop = FALSE], 2L,
                            norms[nonzero], "/")
  }
  out
}

#' Relative-tolerance numerical rank after column normalization.
#' @noRd
matrix_rank_relative <- function(M, tol = 1e-7) {
  if (min(dim(M)) == 0L) return(0L)
  d <- svd(normalize_columns_for_rank(M), nu = 0, nv = 0)$d
  if (length(d) == 0L || d[1] == 0) return(0L)
  sum(d > tol * max(dim(M)) * d[1])
}

#' Keep transformed columns that are not absorbed and add numerical rank.
#'
#' Absorption is measured relative to the corresponding level column, so
#' changing a covariate's units cannot alter the decision. Rank is evaluated
#' after normalizing every surviving transformed column to unit norm.
#' @noRd
select_rank_columns <- function(transformed, reference = transformed,
                                absorption_tol = 1e-7,
                                rank_tol = 1e-7) {
  if (!identical(dim(transformed), dim(reference)))
    stop("transformed and reference matrices must have equal dimensions.")
  if (ncol(transformed) == 0L) return(integer(0))

  kept <- integer(0)
  normalized_kept <- matrix(0, nrow(transformed), 0L)
  for (j in seq_len(ncol(transformed))) {
    candidate <- transformed[, j]
    candidate_norm <- vector_l2_norm(candidate)
    reference_norm <- vector_l2_norm(reference[, j])
    absorbed <- candidate_norm == 0 ||
      (reference_norm > 0 &&
       candidate_norm <= absorption_tol * reference_norm)
    if (!absorbed) {
      trial <- cbind(normalized_kept, candidate / candidate_norm)
      if (matrix_rank_relative(trial, tol = rank_tol) ==
          length(kept) + 1L) {
        kept <- c(kept, j)
        normalized_kept <- trial
      }
    }
  }
  kept
}

#' Efficient quadratic form `eps' A_r eps` in O(n_r), using the identity
#' eps' A_r eps = sum_g (S_g^2 - SS_g)/(m_g - 1) + SS/(n_r - 1),
#' where S_g and SS_g are the per-group sum and sum-of-squares, and SS is
#' the urn-level sum-of-squares. No matrix materialization.
#' @noRd
quad_A_r <- function(eps, group_sizes) {
  n_r <- length(eps)
  eps2 <- eps * eps
  SS <- sum(eps2)
  val <- SS / (n_r - 1)
  off <- 0L
  for (g in seq_along(group_sizes)) {
    m_g <- group_sizes[g]
    idx <- (off + 1L):(off + m_g)
    S_g <- sum(eps[idx])
    SS_g <- sum(eps2[idx])
    val <- val + (S_g * S_g - SS_g) / (m_g - 1)
    off <- off + m_g
  }
  val
}

# ------------------------------------------------------------
# Per-urn small dense matrices
# ------------------------------------------------------------

#' Build A_r^* = I_r^* A_r I_r^* in closed form:
#'   A_r^* = -diag_g{m_g / (m_g - 1) I_g^*} + (n_r / (n_r - 1)) I_r^*.
#' @noRd
make_A_star_r <- function(group_sizes) {
  n_r <- sum(group_sizes)
  I_nr <- diag(n_r)
  g_part <- apply_pq_group_diag(I_nr, group_sizes,
    -group_sizes / (group_sizes - 1), 0)
  u_part <- (n_r / (n_r - 1)) * apply_Istar_urn(I_nr)
  g_part + u_part
}

#' Paper's f-map: f(M_r) such that E(tilde(varsigma)' M_r tilde(varsigma)) =
#' varsigma' f(M_r) varsigma under Lemma rxr(ii). Used by make_A_dagger_r
#' below to fit (beta, gamma) in the three-term closed form for A_r^dagger.
#'
#' This is eq. (f(x)) of bias_correction1.tex:
#'   (n_r - 2)^2 f(M_r) = (a_r/2)(M_r + M_r')
#'                     - (4/(n_r-1)) (MJ - Diag(MJ) + JM - Diag(JM))
#'                     + (2 * 1'M1 / (n_r-1)^2) (J - I)
#' with a_r = 4 + (n_r - 2)^2. MJ has constant columns equal to the row
#' sums of M, so (MJ - Diag(MJ)) is "row sums broadcast off-diagonal".
#' Likewise (JM - Diag(JM)) is "column sums broadcast off-diagonal".
#' @noRd
paper_f <- function(M, n_r) {
  J_minus_I <- 1 - diag(n_r)
  Msym     <- (M + t(M)) / 2
  r        <- rowSums(M)
  c        <- colSums(M)
  off_term <- outer(r, rep(1, n_r)) + outer(rep(1, n_r), c)
  diag(off_term) <- 0

  ((4 + (n_r - 2)^2) * Msym
   + 2 * sum(M) / (n_r - 1)^2 * J_minus_I
   - 4 / (n_r - 1) * off_term) /
    (n_r - 2)^2
}

#' Build A_r^{dagger} for bias-corrected variance.
#'
#' Uses the proof's unified closed form
#'   A_r^dagger = alpha * M + beta * C + gamma * (J - I)
#' where M = A_r^* (hadamard) A_r^*, C = Ddot (J - I) + (J - I) Ddot,
#' Ddot = D_q - tr(D_q) I / (4(n_r - 1)), and
#' q_i = (n_r - m_g) / ((m_g - 1)(n_r - 1)). This covers the minimum
#' admissible urn size without a separate branch.
#' @noRd
make_A_dagger_r <- function(group_sizes) {
  n_r <- sum(group_sizes)
  if (n_r < 4) stop("n_r < 4 not supported (need G_r >= 2 and m_g >= 2).")

  A_star    <- make_A_star_r(group_sizes)
  M         <- A_star * A_star
  J_minus_I <- 1 - diag(n_r)
  q         <- rep((n_r - group_sizes) / ((group_sizes - 1) * (n_r - 1)),
                   times = group_sizes)
  d_dot     <- q - sum(q) / (4 * (n_r - 1))
  C         <- outer(d_dot, rep(1, n_r)) + outer(rep(1, n_r), d_dot)
  diag(C)   <- 0

  a <- (n_r - 2)^2 + 4
  b <- (n_r - 2)^2 * (n_r - 1) + 4
  tau <- sum(d_dot)

  alpha <- (n_r - 2)^2 / a
  beta  <- 4 * (n_r - 2)^2 / (a * b)
  gamma <- 16 * (n_r - 2) * tau / (a * b * n_r * (n_r - 3))

  alpha * M + beta * C + gamma * J_minus_I
}

# ------------------------------------------------------------
# Sample cleaning & urn info assembly
# ------------------------------------------------------------

#' Iteratively drop singleton groups (m_g = 1) and urns with G_r = 1.
#'
#' Dropping a group can leave an urn with fewer than 2 groups, which then
#' needs to be dropped too; dropping that urn is a no-op for other urns but
#' is required for correctness. Loop until the dataset is fixed.
#' @noRd
clean_sample <- function(data) {
  repeat {
    n_before <- nrow(data)
    group_keys <- paste(data$urn_id, data$group_id, sep = "__")
    group_sizes <- table(group_keys)
    keep_groups <- names(group_sizes)[group_sizes >= 2]
    data <- data[group_keys %in% keep_groups, , drop = FALSE]
    if (nrow(data) == 0) break
    urn_group_counts <- tapply(data$group_id, data$urn_id,
                               function(g) length(unique(g)))
    keep_urns <- names(urn_group_counts)[urn_group_counts >= 2]
    data <- data[as.character(data$urn_id) %in% keep_urns, , drop = FALSE]
    if (nrow(data) == n_before) break
  }
  data
}

#' Build one urn_info entry per urn. Assumes `data` is sorted by
#' (urn_id, group_id, i_id) and that rows within an urn are contiguous.
#' @noRd
build_urn_info <- function(data) {
  urn_ids <- unique(data$urn_id)
  group_ids <- unique(data$group_id)
  urn_index <- match(data$urn_id, urn_ids)
  group_index <- match(data$group_id, group_ids)
  urn_sizes <- tabulate(urn_index, nbins = length(urn_ids))
  group_sizes <- tabulate(group_index, nbins = length(group_ids))
  first_group_row <- match(seq_along(group_ids), group_index)
  group_urn_index <- urn_index[first_group_row]
  rows_by_urn <- split(seq_len(nrow(data)), urn_index)

  urn_info <- lapply(seq_along(urn_ids), function(k) {
    idx <- rows_by_urn[[k]]
    groups_k <- unique(group_index[idx])
    list(
      urn_id = urn_ids[k],
      row_idx = idx,
      group_sizes = as.integer(group_sizes[groups_k]),
      n_r = as.integer(urn_sizes[k])
    )
  })
  attr(urn_info, "sample_map") <- list(
    urn_index = urn_index,
    group_index = group_index,
    urn_sizes = as.integer(urn_sizes),
    group_sizes = as.integer(group_sizes),
    group_urn_index = group_urn_index,
    urn_size_by_obs = urn_sizes[urn_index],
    group_size_by_obs = group_sizes[group_index]
  )
  urn_info
}

# ------------------------------------------------------------
# Model-matrix construction
# ------------------------------------------------------------

#' Apply W column-wise to a data-length matrix, urn-by-urn.
#' @noRd
W_times_mat <- function(X, urn_info) {
  if (is.null(X) || ncol(X) == 0) return(X)
  apply_W_sample(X, urn_info)
}

#' Apply I^* column-wise to a data-length matrix, urn-by-urn.
#' @noRd
Istar_times_mat <- function(X, urn_info) {
  if (is.null(X) || ncol(X) == 0) return(X)
  apply_Istar_sample(X, urn_info)
}

#' Build Z, H candidate set, and H^* = I^* H. Z is the concatenation
#' of X_own with W applied to X_peer.
#'
#' The LD drop is performed on H^* rather than H because H^* is what enters
#' the moment-covariance matrix V_n. Candidate columns of H that remain
#' linearly independent in H but collapse after I^* demeaning (e.g.,
#' W X_peer vs W^2 X_peer when m_g is large) are dropped in H^*-space,
#' with the same column indices dropped in H for consistency.
#' @noRd
build_ZH <- function(X_own, X_peer, urn_info, rank_tol = 1e-7) {
  k_own  <- if (is.null(X_own))  0L else ncol(X_own)
  k_peer <- if (is.null(X_peer)) 0L else ncol(X_peer)
  n <- if (length(urn_info) > 0) sum(vapply(urn_info, function(u) u$n_r, integer(1))) else 0L

  WX_peer  <- if (k_peer > 0) W_times_mat(X_peer, urn_info) else NULL
  WX_own   <- if (k_own  > 0) W_times_mat(X_own,  urn_info) else NULL
  W2X_peer <- if (k_peer > 0) W_times_mat(WX_peer, urn_info) else NULL

  Z_parts <- list()
  if (k_own  > 0) Z_parts[[length(Z_parts) + 1]] <- X_own
  if (k_peer > 0) Z_parts[[length(Z_parts) + 1]] <- WX_peer
  Z <- if (length(Z_parts) == 0) matrix(0, n, 0) else do.call(cbind, Z_parts)

  H_parts <- list()
  if (k_own  > 0) H_parts[[length(H_parts) + 1]] <- X_own
  if (k_peer > 0) H_parts[[length(H_parts) + 1]] <- WX_peer
  if (k_own  > 0) H_parts[[length(H_parts) + 1]] <- WX_own
  if (k_peer > 0) H_parts[[length(H_parts) + 1]] <- W2X_peer
  H_full <- if (length(H_parts) == 0) matrix(0, n, 0) else do.call(cbind, H_parts)
  H_full_star <- Istar_times_mat(H_full, urn_info)

  kept <- select_rank_columns(H_full_star, H_full,
                              absorption_tol = rank_tol,
                              rank_tol = rank_tol)
  H     <- H_full     [, kept, drop = FALSE]
  Hstar <- H_full_star[, kept, drop = FALSE]

  list(Z = Z, H = H, Hstar = Hstar, k_own = k_own, k_peer = k_peer)
}

#' Check that Z has full column rank. H is not checked because
#' `build_ZH` already guarantees H^* = I^* H is full column rank, which
#' implies H itself is full column rank (H a = 0 => I^* H a = 0, so
#' H* rank-deficient follows from H rank-deficient; contrapositive).
#' @noRd
rank_check <- function(Z, rank_tol = 1e-7) {
  rank_Z <- matrix_rank_relative(Z, tol = rank_tol)
  if (ncol(Z) > 0 && rank_Z < ncol(Z))
    stop("Z is rank-deficient (ncol = ", ncol(Z),
         ", rank = ", rank_Z, ").")
  invisible(TRUE)
}

#' Drop own/peer covariates that are absorbed by urn demeaning or collinear.
#'
#' The effective beta system uses `Z* = I* [X_own, W X_peer]`. Keep columns in
#' the user's original order whenever each addition increases the rank of Z*.
#' @noRd
drop_collinear_covariates <- function(X_own, X_peer, urn_info,
                                      rank_tol = 1e-7) {
  zh <- build_ZH(X_own, X_peer, urn_info, rank_tol = rank_tol)
  Z <- zh$Z
  if (ncol(Z) == 0L) {
    return(list(X_own = X_own, X_peer = X_peer, zh = zh,
                dropped = character(0), n_dropped = 0L))
  }

  Zstar <- Istar_times_mat(Z, urn_info)
  kept <- select_rank_columns(Zstar, Z,
                              absorption_tol = rank_tol,
                              rank_tol = rank_tol)

  if (length(kept) == ncol(Z)) {
    return(list(X_own = X_own, X_peer = X_peer, zh = zh,
                dropped = character(0), n_dropped = 0L))
  }

  k_own_full <- if (is.null(X_own)) 0L else ncol(X_own)
  k_peer_full <- if (is.null(X_peer)) 0L else ncol(X_peer)
  own_names <- if (k_own_full > 0L) colnames(X_own) else character(0)
  peer_names <- if (k_peer_full > 0L) colnames(X_peer) else character(0)

  kept_own <- kept[kept <= k_own_full]
  kept_peer <- kept[kept > k_own_full] - k_own_full
  dropped_own <- setdiff(seq_len(k_own_full), kept_own)
  dropped_peer <- setdiff(seq_len(k_peer_full), kept_peer)

  X_own_eff <- if (length(kept_own) > 0L) {
    X_own[, kept_own, drop = FALSE]
  } else {
    NULL
  }
  X_peer_eff <- if (length(kept_peer) > 0L) {
    X_peer[, kept_peer, drop = FALSE]
  } else {
    NULL
  }

  dropped <- c(own_names[dropped_own],
               if (length(dropped_peer) > 0L)
                 paste0("ave_", peer_names[dropped_peer]) else character(0))
  zh_eff <- build_ZH(X_own_eff, X_peer_eff, urn_info,
                     rank_tol = rank_tol)
  list(X_own = X_own_eff, X_peer = X_peer_eff, zh = zh_eff,
       dropped = dropped, n_dropped = length(dropped))
}

#' Check beta identification after urn demeaning.
#'
#' `rank_check(Z)` catches collinearity in the level covariates. This catches
#' covariates that survive in levels but are absorbed by urn fixed effects, or
#' become collinear in the usable moments after applying I^*.
#' @noRd
identification_check <- function(Z, Zstar, Hstar, coef_labels,
                                 rank_tol = 1e-7) {
  k_beta <- ncol(Zstar)
  if (k_beta == 0L) return(invisible(TRUE))

  Hstar_scaled <- normalize_columns_for_rank(Hstar)
  Zstar_scaled <- normalize_columns_for_rank(Zstar)
  HIZ_scaled <- crossprod(Hstar_scaled, Zstar_scaled)
  rank_HIZ <- matrix_rank_relative(HIZ_scaled, tol = rank_tol)
  if (rank_HIZ < k_beta) {
    z_norm <- vapply(seq_len(ncol(Zstar)),
                     function(j) vector_l2_norm(Zstar[, j]), numeric(1))
    z_level_norm <- vapply(seq_len(ncol(Z)),
                           function(j) vector_l2_norm(Z[, j]), numeric(1))
    absorbed <- coef_labels[-1][
      z_norm == 0 |
        (z_level_norm > 0 & z_norm <= rank_tol * z_level_norm)
    ]
    detail <- if (length(absorbed) > 0L) {
      paste0(" Absorbed covariate(s): ", paste(absorbed, collapse = ", "), ".")
    } else {
      ""
    }
    stop("Model is not identified after urn demeaning: H' I* Z has rank ",
         rank_HIZ, " but ", k_beta, " beta coefficient(s) were requested.",
         detail,
         " Drop covariates that are constant within urns or collinear after ",
         "the within-urn transformation.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------
# Moment function, V_n, D_n (all urn-wise)
# ------------------------------------------------------------

#' epsilon^*(theta) from precomputed stacks: I*Y, I*WY, Z* = I*Z.
#'
#' By linearity of I^*: eps^*(theta) = I^*(I - lambda W) Y - I^* Z beta
#'                                   = I*Y - lambda * I*WY - (I*Z) beta.
#' Using the precomputed stacked vectors avoids rebuilding W*Y urn-by-urn
#' on every BFGS evaluation.
#' @noRd
epsilon_star <- function(theta, IstarY, IstarWY, Zstar) {
  eps <- IstarY - theta[1] * IstarWY
  beta <- theta[-1]
  if (length(beta) > 0) eps <- eps - as.numeric(Zstar %*% beta)
  eps
}

#' h_n(theta) — length (1 + ncol(H)). Uses the O(n_r) quadratic-form
#' shortcut over global group and urn sufficient statistics.
#' @noRd
h_n <- function(theta, IstarY, IstarWY, Zstar, H, urn_info) {
  eps <- epsilon_star(theta, IstarY, IstarWY, Zstar)
  n <- length(eps)
  q_moment <- quadratic_moment_sum(eps, urn_info)
  if (ncol(H) > 0) {
    c(q_moment, as.numeric(crossprod(H, eps))) / n
  } else {
    q_moment / n
  }
}

#' Sum eps_r' A_r eps_r over urns without materializing any A_r.
#' @noRd
quadratic_moment_sum <- function(eps, urn_info) {
  map <- sample_map(urn_info)
  eps_sq <- eps * eps
  group_sum <- sum_by_index(eps, map$group_index)
  group_sum_sq <- sum_by_index(eps_sq, map$group_index)
  urn_sum_sq <- sum_by_index(eps_sq, map$urn_index)
  sum((group_sum * group_sum - group_sum_sq) /
        (map$group_sizes - 1)) +
    sum(urn_sum_sq / (map$urn_sizes - 1))
}

#' Apply the block-diagonal A operator needed by the quadratic-moment
#' derivative. For urn-demeaned inputs this is also A^* eps.
#' @noRd
apply_A_sample <- function(eps, urn_info) {
  map <- sample_map(urn_info)
  apply_W_sample(eps, urn_info) + eps / (map$urn_size_by_obs - 1)
}

#' Return s_r' (A_r^* o A_r^*) s_r for every urn.
#'
#' A_r^* has zero diagonal, within-group off-diagonal value
#' `1 / (m_g - 1) - 1 / (n_r - 1)`, and between-group value
#' `-1 / (n_r - 1)`. Group and urn sums therefore evaluate the quadratic
#' forms in O(n), without storing dense per-urn matrices.
#' @noRd
quad_Astar_sq_by_urn <- function(s, urn_info) {
  map <- sample_map(urn_info)
  group_sum <- sum_by_index(s, map$group_index)
  group_sum_sq <- sum_by_index(s * s, map$group_index)
  urn_sum <- sum_by_index(s, map$urn_index)
  urn_group_sum_sq <- sum_by_index(
    group_sum * group_sum, map$group_urn_index
  )
  urn_sizes_by_group <- map$urn_sizes[map$group_urn_index]
  within_coef <- 1 / (map$group_sizes - 1) -
    1 / (urn_sizes_by_group - 1)
  within <- sum_by_index(
    within_coef * within_coef *
      (group_sum * group_sum - group_sum_sq),
    map$group_urn_index
  )
  within +
    (urn_sum * urn_sum - urn_group_sum_sq) /
      (map$urn_sizes - 1)^2
}

#' Sum s_r' A_r^dagger s_r over urns using the paper's closed form.
#' @noRd
quad_Adagger_sum <- function(s, urn_info) {
  map <- sample_map(urn_info)
  M_by_urn <- quad_Astar_sq_by_urn(s, urn_info)
  urn_sizes_by_group <- map$urn_sizes[map$group_urn_index]
  q_by_group <- (urn_sizes_by_group - map$group_sizes) /
    ((map$group_sizes - 1) * (urn_sizes_by_group - 1))
  q_sum_by_urn <- sum_by_index(
    q_by_group * map$group_sizes, map$group_urn_index
  )
  d_dot <- q_by_group[map$group_index] -
    q_sum_by_urn[map$urn_index] /
      (4 * (map$urn_size_by_obs - 1))

  urn_sum <- sum_by_index(s, map$urn_index)
  urn_sum_sq <- sum_by_index(s * s, map$urn_index)
  C_by_urn <- 2 * (
    sum_by_index(d_dot * s, map$urn_index) * urn_sum -
      sum_by_index(d_dot * s * s, map$urn_index)
  )
  J_by_urn <- urn_sum * urn_sum - urn_sum_sq

  a <- (map$urn_sizes - 2)^2 + 4
  b <- (map$urn_sizes - 2)^2 * (map$urn_sizes - 1) + 4
  tau <- sum_by_index(d_dot, map$urn_index)
  alpha <- (map$urn_sizes - 2)^2 / a
  beta <- 4 * (map$urn_sizes - 2)^2 / (a * b)
  gamma <- 16 * (map$urn_sizes - 2) * tau /
    (a * b * map$urn_sizes * (map$urn_sizes - 3))

  sum(alpha * M_by_urn + beta * C_by_urn + gamma * J_by_urn)
}

#' Construct exact objective and gradient callbacks for joint GMM.
#'
#' The stage weight Xi is fixed inside each callback pair. The moment
#' Jacobian has quadratic-moment row `-2 C' A eps / n` and linear-moment
#' block `-H' C / n`, where `C = [I* W Y, I* Z]`. The objective derivative
#' is `J' (Xi + Xi') h`, which remains exact if Xi differs from its transpose
#' by floating-point noise.
#' @noRd
make_gmm_objective <- function(Xi, IstarY, IstarWY, Zstar, H, urn_info) {
  n <- length(IstarY)
  Ctheta <- cbind(IstarWY, Zstar)
  linear_grad <- if (ncol(H) > 0L) {
    -crossprod(H, Ctheta) / n
  } else {
    matrix(0, nrow = 0L, ncol = ncol(Ctheta))
  }
  Xi_gradient <- Xi + t(Xi)

  last_theta <- NULL
  last_result <- NULL
  evaluate <- function(theta) {
    theta <- as.numeric(theta)
    if (!is.null(last_theta) && identical(theta, last_theta))
      return(last_result)

    eps <- epsilon_star(theta, IstarY, IstarWY, Zstar)
    q_moment <- quadratic_moment_sum(eps, urn_info)
    hv <- if (ncol(H) > 0L) {
      c(q_moment, as.numeric(crossprod(H, eps))) / n
    } else {
      q_moment / n
    }
    q_grad <- -2 * as.numeric(crossprod(
      Ctheta, apply_A_sample(eps, urn_info)
    )) / n
    J <- rbind(q_grad, linear_grad)
    dimnames(J) <- NULL
    value <- as.numeric(crossprod(hv, Xi %*% hv))
    gradient <- as.numeric(crossprod(J, Xi_gradient %*% hv))

    last_theta <<- theta
    last_result <<- list(
      value = value,
      gradient = gradient,
      moments = hv,
      jacobian = J
    )
    last_result
  }

  list(
    fn = function(theta) evaluate(theta)$value,
    gr = function(theta) evaluate(theta)$gradient,
    moments = function(theta) evaluate(theta)$moments,
    jacobian = function(theta) evaluate(theta)$jacobian
  )
}

#' Assemble the block-diagonal V matrix from its top-left scalar and
#' bottom-right H^{*'} Omega H^* block, then divide by n.
#' @noRd
assemble_V_block <- function(top_left, sigma_sq_vec, Hstar, n,
                             repair_lower = TRUE,
                             eig_floor_rel = 1e-10) {
  if (ncol(Hstar) == 0) {
    V <- matrix(top_left / n, 1, 1)
    attr(V, "lower_projected") <- FALSE
    attr(V, "lower_min_eig_raw") <- NA_real_
    attr(V, "lower_n_eig_repaired") <- 0L
    attr(V, "lower_eig_floor") <- NA_real_
    return(V)
  }
  bottom_right <- crossprod(Hstar * sigma_sq_vec, Hstar)
  lower_diag <- repair_pd_matrix(bottom_right, eig_floor_rel)
  if (isTRUE(repair_lower)) bottom_right <- lower_diag$matrix
  dim_V <- 1 + ncol(Hstar)
  V <- matrix(0, dim_V, dim_V)
  V[1, 1]   <- top_left
  V[-1, -1] <- bottom_right
  V <- V / n
  attr(V, "lower_projected") <- isTRUE(repair_lower) &&
    isTRUE(lower_diag$projected)
  attr(V, "lower_min_eig_raw") <- lower_diag$min_eig_raw / n
  attr(V, "lower_n_eig_repaired") <- if (isTRUE(repair_lower))
    lower_diag$n_eig_repaired else 0L
  attr(V, "lower_eig_floor") <- lower_diag$eig_floor / n
  V
}

#' V_n(Omega) for a diagonal Omega passed as a length-n vector of sigma^2_i.
#'   V_n(1,1) = (2/n) sum_r sigma_r' (A_r^* (.) A_r^*) sigma_r
#'   V_n(lower-right) = (1/n) Hstar' Omega Hstar
#' Uses tr(A^* Omega A^* Omega) = sum_{ij} (A^*_{ij})^2 sigma^2_i sigma^2_j on diagonal Omega.
#' @noRd
V_n_from_Omega <- function(sigma_sq_vec, urn_info, Hstar) {
  top_left <- sum(quad_Astar_sq_by_urn(sigma_sq_vec, urn_info))
  assemble_V_block(2 * top_left, sigma_sq_vec, Hstar, length(sigma_sq_vec))
}

#' D_n at theta, using closed-form (I - lambda W)^{-1} per urn.
#'   D_n(1,1)    = -(2/n) sum_r tr(A_r^* W_r (I-lambda W)^{-1} Omega_r)
#'   D_n(linear, lambda) = -(1/n) H' I^* B Z beta,  B = W(I-lambda W)^{-1}
#'   D_n(linear, beta)   = -(1/n) H' I^* Z  (precomputed outside; passed in as HIZ)
#' @noRd
D_n_at <- function(theta, sigma_sq_vec, urn_info, Z, H, HIZ) {
  n <- length(sigma_sq_vec)
  lambda <- theta[1]
  beta   <- theta[-1]

  map <- sample_map(urn_info)
  urn_sizes_by_group <- map$urn_sizes[map$group_urn_index]
  sigma_sum_by_group <- sum_by_index(sigma_sq_vec, map$group_index)
  diagonal_by_group <- (urn_sizes_by_group - map$group_sizes) /
    (map$group_sizes * (urn_sizes_by_group - 1)) *
    (1 / (1 - lambda) +
       1 / (map$group_sizes - 1 + lambda))
  top_left <- -2 * sum(diagonal_by_group * sigma_sum_by_group) / n

  if (ncol(H) == 0) return(matrix(top_left, 1, 1))

  # Bottom-left: H' I^* W(I-lambda W)^{-1} Z beta.
  bottom_left <- numeric(ncol(H))
  if (length(beta) > 0) {
    Zbeta <- as.numeric(Z %*% beta)
    B_times_Zbeta <- apply_group_operator_sample(
      Zbeta, urn_info,
      p_by_group = -1 / (map$group_sizes - 1 + lambda),
      q_by_group = 1 / (1 - lambda)
    )
    demeaned <- apply_Istar_sample(B_times_Zbeta, urn_info)
    bottom_left <- as.numeric(crossprod(H, demeaned))
  }
  bottom_left <- -bottom_left / n

  # Bottom-right is -(1/n) H' I^* Z = -HIZ / n.
  dim_r <- 1 + ncol(H)
  dim_c <- 1 + ncol(Z)
  D <- matrix(0, dim_r, dim_c)
  D[1, 1] <- top_left
  D[-1, 1] <- bottom_left
  if (ncol(Z) > 0) D[-1, -1] <- -HIZ / n
  D
}

#' Bias-corrected estimator of V_n using the per-urn kernel A_r^{dagger}.
#'
#' Implements the top-left (quadratic-moment) block of \eqn{\widehat V} from
#' the paper's convergence-of-Sigma theorem. A_r^{dagger} is designed so
#' that, under Assumption (innovations) -- independent errors with finite
#' 6+c_eps moments and the zero-diagonal property of A_r^* -- the quadratic
#' form `varsigma' A_r^{dagger} varsigma` is an unbiased estimator of
#' `tr(A_r^* Omega_r A_r^* Omega_r)` on diagonal Omega. No symmetry or
#' normality assumption on the error distribution is required; the
#' distribution-free property follows from diag(A_r^*) = 0 killing all
#' fourth-cumulant contributions. See thm:convergence_Sigma in the paper
#' for the derivation and make_A_dagger_r above for the case-split on n_r.
#'
#' Takes the raw per-urn varsigma vector for the top-left quadratic-moment
#' block and applies a single aggregate `max(0, .)` at the scalar level. The
#' true top-left `2 * sum_r tr(A_r^* Omega_r A_r^* Omega_r) >= 0` is positive
#' by construction, so the raw estimator is unbiased for a non-negative
#' quantity; the aggregate floor only activates on realizations where the sum
#' goes negative, and is asymptotically unbiased (unlike a per-entry floor,
#' which has O(1/n_r) asymptotic bias in the R->infty, bounded-n_r regime).
#' The bottom-right block uses raw sigma_sq_vec by default, then applies a
#' matrix-level eigenvalue repair only when `H*' diag(sigma_sq) H*` is not
#' positive definite.
#' @noRd
V_hat_bias_corrected <- function(urn_info, sigma_sq_vec, Hstar,
                                 repair_lower = TRUE,
                                 eig_floor_rel = 1e-10) {
  top_left_raw <- quad_Adagger_sum(sigma_sq_vec, urn_info)
  top_left <- max(0, top_left_raw)
  assemble_V_block(2 * top_left, sigma_sq_vec, Hstar, length(sigma_sq_vec),
                   repair_lower = repair_lower,
                   eig_floor_rel = eig_floor_rel)
}

#' Structural residuals, fitted values, and urn fixed effects at theta_hat.
#'
#' residual_i = Y_i - lambda_hat * (W Y)_i - (Z beta_hat)_i - alpha_hat_{r(i)},
#' where alpha_hat_r is the urn mean of (Y - lambda_hat W Y - Z beta_hat).
#' @noRd
compute_urn_fe <- function(Y, theta_hat, Z, urn_info) {
  n <- length(Y)
  beta_hat <- theta_hat[-1]
  map <- sample_map(urn_info)
  WY <- apply_W_sample(Y, urn_info)
  Zbeta <- if (length(beta_hat) > 0) as.numeric(Z %*% beta_hat) else rep(0, n)
  pre_alpha <- Y - theta_hat[1] * WY - Zbeta
  alpha_hat <- sum_by_index(pre_alpha, map$urn_index) / map$urn_sizes
  names(alpha_hat) <- vapply(urn_info, function(u) as.character(u$urn_id),
                             character(1))
  residuals_vec <- pre_alpha - alpha_hat[map$urn_index]
  list(alpha_hat = alpha_hat,
       residuals = residuals_vec,
       fitted    = Y - residuals_vec)
}

# ------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------

#' Estimation and Inference for Peer Effects under Conditional Random Assignment
#'
#' Fits the paper's two-step GMM estimator for endogenous peer effects. Each
#' row of `data` is an individual; `urn_var` identifies assignment urns and
#' `group_var` identifies peer groups nested within urns. Use `own_vars` for
#' individual covariates and `peer_vars` for covariates entering through peer
#' averages.
#'
#' @param data a data.frame containing the outcome, identifiers, and
#'   covariates. Rows index individuals.
#' @param y_var character; name of the outcome column in `data`.
#' @param urn_var character; name of the urn (outer cluster) column.
#' @param group_var character; name of the peer-group column (nested in urn).
#' @param own_vars character vector; covariate column names entering the
#'   structural equation as own effects (coefficients named as the column
#'   name itself). `NULL` or empty => no own-effect covariates.
#' @param peer_vars character vector; covariate column names entering as peer
#'   averages (W X). Coefficients are named `ave_<colname>`. `NULL` or empty
#'   => no peer-effect covariates.
#' @param lambda_bounds length-2 numeric, the search interval for lambda in
#'   the 1D optimizer. Default `c(-0.99, 0.99)`.
#' @param lower_eig_floor_rel non-negative scalar controlling the
#'   matrix-level eigenvalue floor for the lower-right block of \eqn{\hat V}.
#'   The floor is `lower_eig_floor_rel * max(mean(abs(diag(B))), 1)` for
#'   `B = H*' diag(sigma_sq) H*` before division by n. Default `1e-10`.
#' @param verbose if TRUE, print stage progress.
#' @param optim_control list passed as `control` to `optim(method = "L-BFGS-B")`
#'   in the joint refinement step. Default `list()` uses R's defaults. Tighten
#'   e.g. `list(factr = 1e1, pgtol = 1e-14)` for cross-validation against
#'   alternative implementations where bit-level agreement is required. The
#'   joint objective uses an exact analytic gradient, so finite-difference
#'   controls such as `ndeps` have no effect.
#' @param rank_tol positive relative tolerance used for absorbed-column and
#'   numerical-rank decisions after temporary column normalization. Columns
#'   are considered in the order supplied, so the first admissible column is
#'   retained when later columns are numerically dependent. Default `1e-7`.
#'
#' @return a `peer_cra` object. The main results are:
#'   - `theta_hat`: named numeric vector of point estimates (lambda, beta).
#'   - `Sigma_hat`: asymptotic variance (dim(theta) x dim(theta)); a named
#'     `NA_real_` matrix when residual-based covariance inference is unavailable.
#'   - `se`: standard errors (length dim(theta)); `NA_real_` when inference is
#'     unavailable.
#'   - `alpha_hat`: estimated urn fixed effects (length n_urns).
#'   - `residuals`: structural-equation residuals (length n).
#'   - `fitted`: fitted values \eqn{\hat Y_i = Y_i - \mathrm{residual}_i} (length n).
#'     Both vectors follow the original input order among observations used
#'     for estimation.
#'   - `sample_index`: integer input-row positions used for estimation, in the
#'     same order as `residuals` and `fitted`.
#'   - `estimation_sample`: logical vector with one entry per input row.
#'   - `n`, `n_urns`, `n_groups`: sample sizes.
#'
#' Standard R methods are provided: `print`, `summary`, `coef`, `vcov`,
#' `confint`, `residuals`, `fitted`, `nobs`.
#'
#' @section Coefficient interpretation:
#' `lambda` is the endogenous peer-effect coefficient. Coefficients named for
#' variables in `own_vars` are own effects. Coefficients prefixed by `ave_` are
#' effects of the corresponding peer averages from `peer_vars`.
#'
#' @section Advanced diagnostics:
#' The fitted object also retains the following diagnostics:
#'   - `n_missing_dropped`: observations dropped for missing estimator variables.
#'   - `n_structural_dropped`: observations removed because their peer group
#'     was a singleton or their urn retained fewer than two peer groups.
#'   - `n_collinear_dropped`: own/peer covariates dropped because they are
#'     absorbed by urn demeaning or collinear with earlier covariates.
#'   - `collinear_dropped`: names of the dropped own/peer covariates.
#'   - `rank_tol`: numerical-rank tolerance used for this fit.
#'   - `internal_covariate_scales`: named within-urn RMS scales used only in
#'     the internal calculation; reported coefficients and covariance are in
#'     the input variables' original units.
#'   - `converged`, `optimizer_status`, `convergence_code`: optimizer status.
#'   - `stage2_weight_used_first`: indicator that the initial matrix supplied
#'     the stage-two estimation weight.
#'   - `final_vcov_nonpd`: indicator that the final residual-based covariance
#'     was not positive definite. `final_vcov_used_first` is a deprecated alias
#'     with the same event value.
#'   - `inference_available`: whether covariance-based inference is available.
#'   - `vcov_status`: one of `final_option_b`,
#'     `stage2_weight_fallback_final_option_b`,
#'     `stage1_residual_final_fallback`, or
#'     `unavailable_both_fallback`.
#'   - final lower-right variance repair diagnostics:
#'     `vhat_lower_projected`, `vhat_lower_min_eig_raw`,
#'     `vhat_lower_n_eig_repaired`, `vhat_lower_eig_floor`.
#'   - stage-one lower-right variance repair diagnostics:
#'     `vhat_s1_lower_projected`, `vhat_s1_lower_min_eig_raw`,
#'     `vhat_s1_lower_n_eig_repaired`, `vhat_s1_lower_eig_floor`.
#'   - `stage1`, `stage2`: per-stage estimates, criteria, convergence codes.
#'
#' @section Asymptotic regime:
#' The estimator's consistency and asymptotic normality require
#'   (i) the average group size n/G to be uniformly bounded, AND
#'   (ii) the total number of groups G to go to infinity.
#' In finite samples, designs with small G or large n/G may give unreliable
#' estimates and standard errors. Warnings emitted by `peer_cra` fire only
#' on observable symptoms of a bad fit (lambda at the boundary of
#' `lambda_bounds`, optim non-convergence), not on design characteristics
#' alone.
#'
#' @references
#' Ying Zeng. "Estimation and Inference for Peer Effects under Conditional
#' Random Assignment." Manuscript.
#'
#' @examples
#' \donttest{
#' # Minimal example using a small synthetic dataset:
#' set.seed(1)
#' R <- 10; G_r <- 5; m_g <- 10; n <- R * G_r * m_g
#' df <- data.frame(
#'   Y      = rnorm(n),
#'   urn    = rep(seq_len(R), each = G_r * m_g),
#'   grp    = rep(rep(seq_len(G_r), each = m_g), times = R),
#'   x1     = rnorm(n),
#'   x2     = rnorm(n)
#' )
#' fit <- peer_cra(df, y_var = "Y", urn_var = "urn", group_var = "grp",
#'                own_vars = c("x1"), peer_vars = c("x2"))
#' print(fit)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#' confint(fit)
#' head(residuals(fit))
#' head(fitted(fit))
#' }
#' @export
peer_cra <- function(data, y_var, urn_var, group_var,
                    own_vars = NULL, peer_vars = NULL,
                    lambda_bounds = c(-0.99, 0.99),
                    lower_eig_floor_rel = 1e-10,
                    verbose = FALSE,
                    optim_control = list(),
                    rank_tol = 1e-7) {

  # --- Input validation ---
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  scalar_name <- function(x, arg) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
      stop("`", arg, "` must be one non-empty column name.")
  }
  scalar_name(y_var, "y_var")
  scalar_name(urn_var, "urn_var")
  scalar_name(group_var, "group_var")
  for (x in list(own_vars = own_vars, peer_vars = peer_vars)) {
    if (!is.null(x) && (!is.character(x) || anyNA(x) || any(!nzchar(x))))
      stop("`own_vars` and `peer_vars` must contain valid column names.")
  }
  for (v in c(y_var, urn_var, group_var))
    if (!v %in% names(data)) stop("Column '", v, "' not found in `data`.")
  if (length(own_vars)  > 0 && !all(own_vars  %in% names(data)))
    stop("own_vars contain names not in `data`: ",
         paste(setdiff(own_vars, names(data)), collapse = ", "))
  if (length(peer_vars) > 0 && !all(peer_vars %in% names(data)))
    stop("peer_vars contain names not in `data`: ",
         paste(setdiff(peer_vars, names(data)), collapse = ", "))
  if (!is.numeric(lambda_bounds) || length(lambda_bounds) != 2L ||
      any(!is.finite(lambda_bounds)) || lambda_bounds[1] >= lambda_bounds[2] ||
      lambda_bounds[1] <= -1 || lambda_bounds[2] >= 1)
    stop("lambda_bounds must be a length-2 numeric inside (-1, 1) with first < second.")
  if (!is.numeric(rank_tol) || length(rank_tol) != 1L ||
      !is.finite(rank_tol) || rank_tol <= 0 || rank_tol >= 1)
    stop("rank_tol must be a finite scalar strictly between 0 and 1.")
  if (!is.numeric(lower_eig_floor_rel) || length(lower_eig_floor_rel) != 1L ||
      !is.finite(lower_eig_floor_rel) || lower_eig_floor_rel < 0)
    stop("lower_eig_floor_rel must be a non-negative scalar.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.")
  if (!is.list(optim_control)) stop("optim_control must be a list.")
  if (length(optim_control) > 0L &&
      (is.null(names(optim_control)) || any(!nzchar(names(optim_control)))))
    stop("Every optim_control element must be named.")
  for (nm in intersect(names(optim_control), c("factr", "pgtol"))) {
    value <- optim_control[[nm]]
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 0)
      stop("optim_control$", nm, " must be a non-negative finite scalar.")
  }
  if ("maxit" %in% names(optim_control)) {
    value <- optim_control$maxit
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 1 || value != as.integer(value))
      stop("optim_control$maxit must be a positive integer.")
  }

  numeric_column <- function(v, role, allow_logical = FALSE) {
    x <- data[[v]]
    if (is.logical(x) && allow_logical) x <- as.numeric(x)
    if (!is.numeric(x))
      stop(role, " '", v, "' must be numeric; received ",
           paste(class(data[[v]]), collapse = "/"), ".")
    if (any(!is.finite(x) & !is.na(x)))
      stop(role, " '", v, "' contains non-finite values.")
    x
  }
  y_values <- numeric_column(y_var, "Outcome")
  extra_cols <- union(own_vars, peer_vars)
  extra_values <- setNames(vector("list", length(extra_cols)), extra_cols)
  for (v in extra_cols)
    extra_values[[v]] <- numeric_column(v, "Covariate", allow_logical = TRUE)

  # --- Assemble working frame + iterative sample cleaning ---
  # Store ids as character so `table()` doesn't generate spurious zero
  # counts for factor levels that were dropped during cleaning.
  work <- data.frame(
    .row_id  = seq_len(nrow(data)),
    urn_id   = as.character(data[[urn_var]]),
    group_id = paste(as.character(data[[urn_var]]),
                      as.character(data[[group_var]]), sep = "__"),
    Y        = y_values,
    stringsAsFactors = FALSE
  )
  for (v in extra_cols) work[[v]] <- extra_values[[v]]

  n_before_missing <- nrow(work)
  complete_rows <- stats::complete.cases(work)
  n_missing_dropped <- n_before_missing - sum(complete_rows)
  if (n_missing_dropped > 0L) {
    warning("peer_cra: dropped ", n_missing_dropped,
            " observation(s) with missing values in y_var, urn_var, ",
            "group_var, own_vars, or peer_vars.",
            call. = FALSE)
    work <- work[complete_rows, , drop = FALSE]
  }

  n_after_missing <- nrow(work)
  work <- clean_sample(work)
  n_structural_dropped <- n_after_missing - nrow(work)
  if (nrow(work) == 0) stop("No observations remain after cleaning.")
  work <- work[order(work$urn_id, work$group_id), , drop = FALSE]
  rownames(work) <- NULL

  Y      <- work$Y
  n      <- length(Y)
  X_own  <- if (length(own_vars)  > 0) as.matrix(work[, own_vars,  drop = FALSE]) else NULL
  X_peer <- if (length(peer_vars) > 0) as.matrix(work[, peer_vars, drop = FALSE]) else NULL
  if (!is.null(X_own))  colnames(X_own)  <- own_vars
  if (!is.null(X_peer)) colnames(X_peer) <- peer_vars

  urn_info <- build_urn_info(work)

  # --- Build Z, H, Z*, H* (with LD drop on H*) ---
  dropped_covariates <- character(0)
  drop_info <- drop_collinear_covariates(
    X_own, X_peer, urn_info, rank_tol = rank_tol
  )
  X_own <- drop_info$X_own
  X_peer <- drop_info$X_peer
  dropped_covariates <- drop_info$dropped
  if (drop_info$n_dropped > 0L) {
    warning("peer_cra: dropped ", drop_info$n_dropped,
            " collinear/absorbed covariate(s): ",
            paste(dropped_covariates, collapse = ", "), ".",
            call. = FALSE)
  }

  # Keep the retained design in standardized coordinates throughout every
  # numerical calculation. The scale of each structural covariate is its
  # within-urn RMS in Z* = I* [X_own, W X_peer]. Dividing the input column once
  # also scales every W-transform constructed from it, e.g. W(X/s) = WX/s.
  # The raw retained matrices remain available only for public-unit fitted
  # values, residuals, and coefficient/covariance back-transformation.
  zh_report <- drop_info$zh
  Z_report <- zh_report$Z
  Zstar_report <- Istar_times_mat(Z_report, urn_info)
  covariate_scales <- if (ncol(Zstar_report) > 0L) {
    vapply(seq_len(ncol(Zstar_report)), function(j) {
      vector_l2_norm(Zstar_report[, j]) / sqrt(n)
    }, numeric(1))
  } else {
    numeric(0)
  }
  if (any(!is.finite(covariate_scales) | covariate_scales <= 0))
    stop("Internal covariate scaling failed after rank selection.")

  k_own_report <- if (is.null(X_own)) 0L else ncol(X_own)
  k_peer_report <- if (is.null(X_peer)) 0L else ncol(X_peer)
  X_own_internal <- if (k_own_report > 0L) {
    sweep(X_own, 2L, covariate_scales[seq_len(k_own_report)], "/")
  } else {
    NULL
  }
  X_peer_internal <- if (k_peer_report > 0L) {
    peer_idx <- k_own_report + seq_len(k_peer_report)
    sweep(X_peer, 2L, covariate_scales[peer_idx], "/")
  } else {
    NULL
  }
  if (!is.null(X_own_internal)) colnames(X_own_internal) <- colnames(X_own)
  if (!is.null(X_peer_internal)) colnames(X_peer_internal) <- colnames(X_peer)

  zh <- build_ZH(X_own_internal, X_peer_internal, urn_info,
                 rank_tol = rank_tol)
  Z <- zh$Z; H <- zh$H; Hstar <- zh$Hstar
  rank_check(Z, rank_tol = rank_tol)
  Zstar <- Istar_times_mat(Z, urn_info)
  HIZ <- if (ncol(H) > 0 && ncol(Z) > 0) crossprod(H, Zstar) else matrix(0, ncol(H), ncol(Z))

  k_own  <- zh$k_own
  k_peer <- zh$k_peer
  k_beta <- ncol(Z)
  coef_labels <- c("lambda",
                   if (k_own  > 0) colnames(X_own)                else character(0),
                   if (k_peer > 0) paste0("ave_", colnames(X_peer)) else character(0))
  identification_check(Z, Zstar, Hstar, coef_labels,
                       rank_tol = rank_tol)

  names(covariate_scales) <- coef_labels[-1]
  report_transform <- diag(c(1, 1 / covariate_scales),
                           nrow = length(coef_labels))
  dimnames(report_transform) <- list(coef_labels, coef_labels)
  to_reported_theta <- function(theta_internal) {
    out <- as.numeric(report_transform %*% theta_internal)
    names(out) <- coef_labels
    out
  }

  # --- Precompute urn-stacked stars ---
  IstarY <- apply_Istar_sample(Y, urn_info)
  IstarWY <- apply_Istar_sample(apply_W_sample(Y, urn_info), urn_info)

  # --- Stage 1 weight: Xi1 = V_n^{-1} | Omega = I ---
  sigma_sq_init <- rep(1, n)
  V1 <- V_n_from_Omega(sigma_sq_init, urn_info, Hstar)
  stage1_weight <- select_pd_or_first(V1, V1)
  Xi1 <- stage1_weight$inverse

  # --- Precompute lambda-independent crossproducts for the concentrated
  # step. HZ, HIY, HIWY are all n-scale and constant across optimize
  # evaluations; within the concentrated objective, Hy(lambda) is assembled
  # as the linear combination HIY - lambda * HIWY.
  HZ_cache <- crossprod(Hstar, Zstar)                 # k_H x k_Z
  HIY      <- as.numeric(crossprod(Hstar, IstarY))    # k_H
  HIWY     <- as.numeric(crossprod(Hstar, IstarWY))   # k_H

  # --- Closure factory: given Xi22, factor (Z'H) Xi22 (H'Z) ONCE per weight
  # and return a lambda-only callback for the concentrated search.
  make_compute_beta <- function(Xi22) {
    if (k_beta == 0) return(function(lambda) numeric(0))
    ZHXi    <- crossprod(HZ_cache, Xi22)       # k_Z x k_H
    lhs_inv <- solve_robust(ZHXi %*% HZ_cache) # k_Z x k_Z, factored once
    function(lambda) {
      Hy <- HIY - lambda * HIWY                # k_H
      as.numeric(lhs_inv %*% (ZHXi %*% Hy))    # k_Z
    }
  }

  extract_Xi22 <- function(Xi) {
    if (ncol(H) == 0) matrix(0, 0, 0) else Xi[-1, -1, drop = FALSE]
  }

  # Hybrid optimizer: concentrated 1D `optimize` gets a warm start, then
  # joint L-BFGS-B over the full (lambda, beta) refines to the exact joint
  # minimizer of Q. The warm start makes BFGS converge in a handful of
  # iterations and avoids local-minimum / starting-value fragility.
  joint_refine <- function(theta_start, Xi) {
    k_b <- length(theta_start) - 1L
    callbacks <- make_gmm_objective(
      Xi, IstarY, IstarWY, Zstar, H, urn_info
    )
    q_start <- callbacks$fn(theta_start)

    # With no covariates, the concentrated 1D problem is the full joint
    # problem, so its minimizer can be returned directly. When covariates are
    # present, always run the joint refinement: an absolute criterion cutoff
    # is unit-dependent and does not establish exact identification or a
    # first-order condition.
    if (k_b == 0L) {
      return(list(par = theta_start, value = q_start, convergence = 0L,
                  status = "short_circuit_no_covariates"))
    }

    # Wrap optim in tryCatch so that hard errors ("L-BFGS-B needs finite
    # values of 'fn'", non-finite gradient at the warm start, etc.) surface
    # as status = "failed" and the warm start is returned, matching the
    # Mata side's _optimize-with-rc-recovery path. Without this, a single
    # non-finite objective evaluation aborts the entire fit with no e(b)
    # analogue on the R side. convergence = 99 is a sentinel outside optim's
    # documented codes (0-52, mostly); see ?optim.
    res <- tryCatch(
      optim(theta_start, callbacks$fn, gr = callbacks$gr,
            method = "L-BFGS-B",
            lower = c(lambda_bounds[1], rep(-Inf, k_b)),
            upper = c(lambda_bounds[2], rep( Inf, k_b)),
            control = optim_control),
      error = function(e) list(par = theta_start, value = q_start,
                               convergence = 99L,
                               message = conditionMessage(e))
    )
    res$status <- if (res$convergence == 0L) "converged" else "failed"
    res
  }

  # --- One stage of the two-step: concentrated 1D warm start, then joint
  # L-BFGS-B refinement to the exact joint arg-min for weight matrix Xi.
  run_stage <- function(Xi, label) {
    if (verbose) cat(label, "1D concentrated warm start, then joint refinement\n")
    compute_beta <- make_compute_beta(extract_Xi22(Xi))
    obj <- function(lambda) {
      beta <- compute_beta(lambda)
      hv <- h_n(c(lambda, beta), IstarY, IstarWY, Zstar, H, urn_info)
      as.numeric(crossprod(hv, Xi %*% hv))
    }
    conc <- optimize(obj, interval = lambda_bounds, tol = 1e-8)
    theta_conc <- c(conc$minimum, compute_beta(conc$minimum))
    joint_refine(theta_conc, Xi)
  }

  # --- Bias-corrected per-observation variance estimator, used by both the
  # stage-2 weight and the final sandwich. No per-entry floor is applied:
  # the theory targets the raw bias-corrected varsigma vector. Only the
  # top-left quadratic-moment sum in V_hat_bias_corrected is aggregate-floored
  # at zero because its population target is non-negative.
  compute_varsigma <- function(eps) {
    map <- sample_map(urn_info)
    eps_sq <- eps * eps
    urn_mean <- sum_by_index(eps_sq, map$urn_index) / map$urn_sizes
    mean_obs <- urn_mean[map$urn_index]
    p_obs <- map$urn_size_by_obs / (map$urn_size_by_obs - 2)
    q_obs <- map$urn_size_by_obs / (map$urn_size_by_obs - 1)
    p_obs * (eps_sq - mean_obs) + q_obs * mean_obs
  }

  # --- Stage 1 ---
  opt1   <- run_stage(Xi1, "Stage 1:")
  theta1_internal <- opt1$par

  # --- Stage-2 weight: Xi2 = V_hat(theta1)^{-1}, bias-corrected (two-step) ---
  vs1      <- compute_varsigma(epsilon_star(theta1_internal, IstarY,
                                            IstarWY, Zstar))
  V_hat_s1 <- V_hat_bias_corrected(urn_info, vs1, Hstar,
                                   repair_lower = TRUE,
                                   eig_floor_rel = lower_eig_floor_rel)
  stage2_weight <- select_pd_or_first(V_hat_s1, V1)
  Xi2 <- stage2_weight$inverse

  # --- Stage 2 ---
  opt2               <- run_stage(Xi2, "Stage 2:")
  theta_hat_internal <- opt2$par
  names(theta_hat_internal) <- coef_labels
  theta1 <- to_reported_theta(theta1_internal)
  theta_hat <- to_reported_theta(theta_hat_internal)

  stage1 <- list(par = theta1,    value = opt1$value, convergence = opt1$convergence,
                 status = opt1$status)
  stage2 <- list(par = theta_hat, value = opt2$value, convergence = opt2$convergence,
                 status = opt2$status)

  # --- Final variance: re-compute V_hat at theta_hat (Option B) ---
  # Matches R `gmm::gmm` (vcov = "MDS"/"iid") and Stata `gmm` conventions.
  # Stage-2 weight above stays at V_hat(theta1) per textbook two-step.
  vs_final <- compute_varsigma(epsilon_star(theta_hat_internal, IstarY,
                                            IstarWY, Zstar))

  V_hat <- V_hat_bias_corrected(urn_info, vs_final, Hstar,
                                repair_lower = TRUE,
                                eig_floor_rel = lower_eig_floor_rel)
  D_hat <- D_n_at(theta_hat_internal, vs_final, urn_info, Z, H, HIZ)

  covariance <- compute_final_covariance(
    D_hat = D_hat,
    V_hat_s1 = V_hat_s1,
    V_hat_final = V_hat,
    Xi2 = Xi2,
    stage2_weight_used_first = stage2_weight$used_first,
    n = n,
    report_transform = report_transform,
    coef_labels = coef_labels
  )
  Sigma_hat <- covariance$Sigma_reported
  se <- covariance$se

  if (!covariance$inference_available) {
    warning(
      "peer_cra: the stage-one and final residual-based moment-covariance matrices are not positive definite; point estimates are returned, but covariance inference is unavailable.",
      call. = FALSE
    )
  }

  # ---- Diagnostic warnings -------------------------------------------------
  # Fire only on observable symptoms of a problematic fit. Theoretical
  # concerns about the design (small G, large n/G) are documented in the
  # function doc but not warned about unless they manifest in a bad fit.

  lambda_hat <- theta_hat[1]
  bound_tol <- 1e-3
  if (abs(lambda_hat - lambda_bounds[1]) < bound_tol ||
      abs(lambda_hat - lambda_bounds[2]) < bound_tol) {
    warning("peer_cra: lambda_hat = ", formatC(lambda_hat, digits = 4),
            " is at the boundary of lambda_bounds = (",
            lambda_bounds[1], ", ", lambda_bounds[2], "). ",
            "This typically indicates weak identification; point estimate and ",
            "standard errors may be unreliable.",
            call. = FALSE)
  }

  # Warn only on stage-2 non-convergence: stage 2 is what determines the
  # final theta_hat, so its success is the operational "converged" state.
  # A stage-1 non-zero code with stage-2 success usually indicates an
  # intermediate line-search issue; stage 2 is the operational fit.
  if (stage2$convergence != 0) {
    warning("peer_cra: stage-2 optimizer did not converge (code ",
            stage2$convergence, "). Point estimates are returned but ",
            "may be unreliable. Stage-1 code was ", stage1$convergence,
            ". See ?optim for code meanings.",
            call. = FALSE)
  }

  vhat_lower_projected <- isTRUE(attr(V_hat, "lower_projected"))
  vhat_lower_min_eig_raw <- attr(V_hat, "lower_min_eig_raw")
  vhat_lower_n_eig_repaired <- attr(V_hat, "lower_n_eig_repaired")
  vhat_lower_eig_floor <- attr(V_hat, "lower_eig_floor")
  vhat_s1_lower_projected <- isTRUE(attr(V_hat_s1, "lower_projected"))
  vhat_s1_lower_min_eig_raw <- attr(V_hat_s1, "lower_min_eig_raw")
  vhat_s1_lower_n_eig_repaired <- attr(V_hat_s1, "lower_n_eig_repaired")
  vhat_s1_lower_eig_floor <- attr(V_hat_s1, "lower_eig_floor")

  urn_sizes    <- vapply(urn_info, function(u) u$n_r,                integer(1))
  urn_ngroups  <- vapply(urn_info, function(u) length(u$group_sizes), integer(1))
  total_groups <- sum(urn_ngroups)

  fe <- compute_urn_fe(Y, theta_hat, Z_report, urn_info)
  alpha_hat     <- fe$alpha_hat
  residuals_vec <- fe$residuals
  fitted_vec    <- fe$fitted
  output_order <- order(work$.row_id)
  sample_index <- work$.row_id[output_order]
  estimation_sample <- rep(FALSE, nrow(data))
  estimation_sample[sample_index] <- TRUE
  residuals_vec <- residuals_vec[output_order]
  fitted_vec <- fitted_vec[output_order]
  names(residuals_vec) <- rownames(data)[sample_index]
  names(fitted_vec) <- rownames(data)[sample_index]

  # optimizer_status is an enum that distinguishes the three paths
  # `joint_refine` can take on stage 2:
  #   - "converged": L-BFGS-B ran and reported convergence.
  #   - "failed":    L-BFGS-B ran but did not converge.
  #   - "short_circuit_no_covariates": k_b == 0, BFGS skipped by design.
  # `converged` below is the derived boolean kept for backward compat;
  # it is TRUE for "converged" and the structural no-covariate shortcut.
  out <- list(
    theta_hat = theta_hat,
    Sigma_hat = Sigma_hat,
    se = setNames(se, coef_labels),
    coef_labels = coef_labels,
    alpha_hat = alpha_hat,
    residuals = residuals_vec,
    fitted = fitted_vec,
    sample_index = sample_index,
    estimation_sample = estimation_sample,
    n = n,
    n_urns = length(urn_info),
    n_groups = total_groups,
    n_missing_dropped = n_missing_dropped,
    n_structural_dropped = n_structural_dropped,
    n_collinear_dropped = length(dropped_covariates),
    collinear_dropped = dropped_covariates,
    rank_tol = rank_tol,
    internal_covariate_scales = covariate_scales,
    urn_sizes = urn_sizes,
    urn_ngroups = urn_ngroups,
    converged = (stage2$convergence == 0),
    optimizer_status = stage2$status,
    convergence_code = c(stage1 = stage1$convergence,
                         stage2 = stage2$convergence),
    stage2_weight_used_first = stage2_weight$used_first,
    final_vcov_nonpd = covariance$final_vcov_nonpd,
    final_vcov_used_first = covariance$final_vcov_nonpd,
    inference_available = covariance$inference_available,
    vcov_status = covariance$vcov_status,
    vhat_lower_projected = vhat_lower_projected,
    vhat_lower_min_eig_raw = vhat_lower_min_eig_raw,
    vhat_lower_n_eig_repaired = vhat_lower_n_eig_repaired,
    vhat_lower_eig_floor = vhat_lower_eig_floor,
    vhat_s1_lower_projected = vhat_s1_lower_projected,
    vhat_s1_lower_min_eig_raw = vhat_s1_lower_min_eig_raw,
    vhat_s1_lower_n_eig_repaired = vhat_s1_lower_n_eig_repaired,
    vhat_s1_lower_eig_floor = vhat_s1_lower_eig_floor,
    stage1 = list(theta = theta1,
                  criterion = stage1$value,
                  convergence = stage1$convergence,
                  status = stage1$status),
    stage2 = list(theta = theta_hat,
                  criterion = stage2$value,
                  convergence = stage2$convergence,
                  status = stage2$status)
  )
  class(out) <- "peer_cra"
  out
}

# ------------------------------------------------------------
# S3 methods for the "peer_cra" class
# ------------------------------------------------------------

#' Print a fitted peer-GMM object.
#'
#' Prints sample sizes, actionable notices, and a coefficient table with point
#' estimates, standard errors, z-statistics, two-sided p-values, and 95% Wald
#' confidence intervals. Detailed optimizer and covariance diagnostics remain
#' available in the fitted object.
#'
#' @param x a `peer_cra` object returned by `peer_cra`.
#' @param digits number of significant digits for numeric columns (default 4).
#' @param ... ignored.
#' @export
print.peer_cra <- function(x, digits = 4, ...) {
  cat("PeerCRA\n")
  cat(sprintf("Observations: %d   Urns: %d   Peer groups: %d\n",
              x$n, x$n_urns, x$n_groups))
  if (!is.null(x$n_missing_dropped) && x$n_missing_dropped > 0L) {
    cat(sprintf("Data cleaning removed %d observation(s) with missing values.\n",
                x$n_missing_dropped))
  }
  if (!is.null(x$n_structural_dropped) && x$n_structural_dropped > 0L) {
    cat(sprintf(paste0("Structural cleaning removed %d observation(s) in ",
                       "singleton groups or urns with fewer than two groups.\n"),
                x$n_structural_dropped))
  }
  if (!is.null(x$n_collinear_dropped) && x$n_collinear_dropped > 0L) {
    dropped_names <- if (length(x$collinear_dropped) > 0L) {
      paste0(": ", paste(x$collinear_dropped, collapse = ", "))
    } else {
      ""
    }
    cat(sprintf("Dropped %d collinear or absorbed covariate(s)%s.\n",
                x$n_collinear_dropped, dropped_names))
  }
  if (!isTRUE(x$converged))
    cat("Warning: optimization did not converge; estimates may be unreliable.\n")
  if (isTRUE(x$stage2_weight_used_first))
    cat("Note: stage two retained the first-step weighting matrix.\n")
  if (identical(x$vcov_status, "stage1_residual_final_fallback"))
    cat("Note: inference uses the stage-one residual covariance.\n")
  if (isFALSE(x$inference_available))
    cat("Warning: inference is unavailable; point estimates are shown without standard errors.\n")
  if (isTRUE(x$vhat_lower_projected) || isTRUE(x$vhat_s1_lower_projected))
    cat("Note: numerical covariance repair was applied.\n")
  cat("\nCoefficients:\n")

  est <- x$theta_hat
  se  <- x$se
  zstat <- est / se
  pval  <- 2 * pnorm(-abs(zstat))
  q95   <- qnorm(0.975)
  ci_lo <- est - q95 * se
  ci_hi <- est + q95 * se

  tab <- data.frame(
    Estimate = signif(est,   digits),
    `Std. Err` = signif(se,  digits),
    `z value` = round(zstat, 2),
    `Pr(>|z|)` = format.pval(pval, digits = 3, eps = .Machine$double.eps),
    `95% CI low`  = signif(ci_lo, digits),
    `95% CI high` = signif(ci_hi, digits),
    check.names = FALSE, row.names = x$coef_labels
  )
  print(tab)
  invisible(x)
}

#' Extract coefficients from a fitted peer-GMM object.
#' @param object a `peer_cra` object.
#' @param ... ignored.
#' @export
coef.peer_cra <- function(object, ...) object$theta_hat

#' Extract the asymptotic variance-covariance matrix.
#' @param object a `peer_cra` object.
#' @param ... ignored.
#' @export
vcov.peer_cra <- function(object, ...) object$Sigma_hat

#' Wald confidence intervals for a fitted peer-GMM object.
#' @param object a `peer_cra` object.
#' @param parm character or integer selecting which parameters; default all.
#' @param level confidence level (default 0.95).
#' @param ... ignored.
#' @return a matrix with two columns (lower, upper) and rows for each parameter.
#' @export
confint.peer_cra <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be one finite number strictly between 0 and 1.")
  est <- object$theta_hat
  se  <- object$se
  labels <- object$coef_labels
  if (missing(parm)) parm <- seq_along(est)
  else if (is.character(parm)) {
    unknown <- setdiff(parm, labels)
    if (length(unknown) > 0L)
      stop("Unknown parameter(s): ", paste(unknown, collapse = ", "), ".")
    parm <- match(parm, labels)
  } else if (!is.numeric(parm) || anyNA(parm) ||
             any(parm != as.integer(parm)) ||
             any(parm < 1L | parm > length(est))) {
    stop("parm must contain valid parameter names or integer positions.")
  }
  q <- qnorm(1 - (1 - level) / 2)
  lo <- est[parm] - q * se[parm]
  hi <- est[parm] + q * se[parm]
  ci <- cbind(lo, hi)
  colnames(ci) <- sprintf("%.1f %%", c((1 - level) / 2, 1 - (1 - level) / 2) * 100)
  rownames(ci) <- labels[parm]
  ci
}

#' Residuals from a fitted peer-GMM object.
#'
#' Structural-equation residuals: \eqn{r_{i} = Y_i - \hat\lambda (W Y)_i -
#' Z_i' \hat\beta - \hat\alpha_{r(i)}}, where \eqn{\hat\alpha_r} is the urn
#' mean of \eqn{Y - \hat\lambda W Y - Z \hat\beta}.
#'
#' @param object a `peer_cra` object.
#' @param ... ignored.
#' @return numeric vector of length n.
#' @export
residuals.peer_cra <- function(object, ...) object$residuals

#' Fitted values from a fitted peer-GMM object.
#' @param object a `peer_cra` object.
#' @param ... ignored.
#' @return numeric vector \eqn{\hat Y_i = Y_i - r_i}.
#' @export
fitted.peer_cra <- function(object, ...) object$fitted

#' Number of observations in a fitted peer-GMM object.
#' @param object a `peer_cra` object.
#' @param ... ignored.
#' @export
nobs.peer_cra <- function(object, ...) object$n

#' Summary of a fitted peer-GMM object.
#'
#' Displays the ordinary coefficient table. Residuals and urn fixed effects
#' remain available through \code{residuals()} and the fitted object's
#' \code{alpha_hat} component.
#'
#' @param object a `peer_cra` object.
#' @param ... ignored.
#' @export
summary.peer_cra <- function(object, ...) {
  print(object)
  invisible(object)
}
