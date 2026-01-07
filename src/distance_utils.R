# =========================================================
# Distance Calculation and Matrix Utilities
# =========================================================
# This module contains functions for computing Mahalanobis distances
# and matrix manipulation utilities used in clustering.
# =========================================================

#' Matrix Trace
#'
#' Computes the trace (sum of diagonal elements) of a matrix.
#'
#' @param M Matrix
#' @return Numeric scalar, the trace of M
tr  <- function(M) sum(diag(M))


#' Symmetrize Matrix
#'
#' Symmetrizes a matrix by averaging with its transpose.
#'
#' @param M Matrix
#' @return Symmetric matrix, (M + t(M))/2
sym <- function(M) (M + t(M))/2


#' Safe Determinant-Powered Matrix Normalization
#'
#' Normalizes a matrix by scaling by the p-th root of its determinant
#' to ensure numerical stability.
#'
#' @param M Matrix to normalize
#' @param eps Numeric, small positive value for numerical stability (default: 1e-8)
#' @param min_eig Numeric, minimum eigenvalue threshold (default: 1e-6)
#' @return Normalized matrix with determinant-adjusted scaling
#'
#' @details
#' The matrix is symmetrized, then scaled by the p-th root of its determinant
#' to prevent numerical issues with very small or very large determinants.
safe_det_powm <- function(M, eps = 1e-8, min_eig = 1e-6) {
  if (!is.matrix(M)) {
    M <- as.matrix(M)
  }
  
  M <- sym(M)
  p <- nrow(M)
  det_val <- det(M + eps*diag(p))
  scale_fac <- abs(det_val)^(1/p)
  if (scale_fac < eps) scale_fac <- eps
  M / scale_fac
}


#' Strictly Elimination Matrix
#'
#' Constructs the elimination matrix Lbar for converting a matrix
#' to its strictly half-vectorized form.
#'
#' @param n Integer, dimension of the square matrix
#' @return Matrix of dimension (n*(n-1)/2) x n^2
#'
#' @details
#' This matrix is used in the conversion between orthogonal matrices
#' and their parameterized forms for optimization.
Lbar <- function(n){
  Lb      <- matrix(0,(n*(n-1)/2),n^2)
  indcol  <- sapply(0:(n-2), function(i) 
    sapply((i*n+(i+2)):(n+i*n), function(j) j )
  )
  
  c1      <- 1:(n*(n-1)/2)
  c2      <- unlist(indcol)
  ind     <- cbind(c1,c2)
  Lb[ind] <- 1
  
  return(Lb)
}


#' Strictly Half-Vectorization Operator
#'
#' Converts a square matrix to its strictly half-vectorized form
#' (excluding diagonal elements).
#'
#' @param A Square matrix
#' @return Numeric vector of length n*(n-1)/2 containing off-diagonal elements
#'
#' @details
#' This extracts only the strictly lower triangular elements of the matrix.
bvech <- function(A){   # A: square matrix
  n   <- ncol(A)
  res <- Lbar(n)%*%c(A)
  return(res)
}


#' Convert Orthogonal Matrix to Permutation-Lower Form
#'
#' Decomposes an orthogonal matrix Q into a permutation matrix P
#' and a lower triangular parameter vector l.
#'
#' @param Q Orthogonal matrix
#' @return List containing:
#'   \describe{
#'     \item{P}{Permutation matrix}
#'     \item{l}{Numeric vector of n(n-1)/2 elements representing lower triangular part}
#'   }
#'
#' @details
#' Uses LU decomposition to extract the permutation and lower triangular structure
#' for parameterization in optimization routines.
QPl <- function(Q){   # Orthogonal matrix
  luQ    <- Matrix::lu(Q)
  eluQ   <- Matrix::expand(luQ)
  l      <- bvech(as.matrix(eluQ$L))
  U      <- eluQ$U
  P      <- eluQ$P
  return(list(      # It returns:
    P=P,  # Permutation matrix
    l=l   # Vector of n(n-1)/2 elements
  ))
}


#' Convert Permutation-Lower Form to Orthogonal Matrix
#'
#' Reconstructs an orthogonal matrix Q from a permutation matrix P
#' and a lower triangular parameter vector l.
#'
#' @param P Permutation matrix
#' @param l Numeric vector of n(n-1)/2 elements representing lower triangular part
#' @return Orthogonal matrix Q
#'
#' @details
#' This is the inverse operation of QPl, reconstructing the orthogonal matrix
#' from its parameterized form.
PlQ <- function(P,  # Permutation matrix
                l   # Vector of n(n-1)/2 elements
){
  n  <- (1+sqrt(1+8*length(l)))/2
  L  <- diag(n) + matrix(t(Lbar(n))%*%l,n,n) 
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(as.matrix(Q))         # It returns an orthogonal matrix
}


#' Calculate Mahalanobis Distance
#'
#' Computes the Mahalanobis distance between a data point and cluster center
#' based on the specified model's covariance structure.
#'
#' @param model Character string, model identifier (F1-F14)
#' @param diff Numeric vector, difference between data point and center (x - mu)
#' @param A Matrix or NULL, shared A matrix (for models F3-F4, F7-F8, F11-F12)
#' @param A_k Matrix or NULL, cluster-specific A_k matrix (for models F5-F6, F9-F10, F13-F14)
#' @param D Matrix or NULL, shared D matrix (for models F7-F10)
#' @param D_k Matrix or NULL, cluster-specific D_k matrix (for models F11-F14)
#' @param eps Numeric, small positive value for numerical stability (default: 1e-8)
#' @param min_eig Numeric, minimum eigenvalue threshold (default: 1e-6)
#'
#' @return Numeric scalar, Mahalanobis distance
#'
#' @details
#' The distance is computed based on the model's covariance structure:
#' \itemize{
#'   \item F1-F2: Euclidean distance (spherical covariance)
#'   \item F3-F4: Uses shared A matrix
#'   \item F5-F6: Uses cluster-specific A_k matrix
#'   \item F7-F8: Uses D * A * D'
#'   \item F9-F10: Uses D * A_k * D'
#'   \item F11-F12: Uses D_k * A * D_k'
#'   \item F13-F14: Uses D_k * A_k * D_k'
#' }
calculate_distance <- function(model, diff, A, A_k, D, D_k, eps = 1e-8, min_eig = 1e-6) {
  p <- length(diff)
  if (model %in% c("F1","F2")) return(sqrt(sum(diff^2)) + eps)
  
  Sigma <-
    if (model %in% c("F3","F4")) A else
      if (model %in% c("F5","F6")) A_k else
        if (model %in% c("F7","F8")) D %*% A %*% t(D) else
          if (model %in% c("F9","F10")) D %*% A_k %*% t(D) else
            if (model %in% c("F11","F12")) D_k %*% A %*% t(D_k) else
              if (model %in% c("F13","F14")) D_k %*% A_k %*% t(D_k) else stop("Unknown model")
  
  Sigma <- sym(Sigma)
  eig <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(eig$values, min_eig)
  vecs <- eig$vectors
  Sigma_psd <- vecs %*% diag(vals) %*% t(vecs) + diag(eps, p)
  inv_sig <- pseudoinverse(Sigma_psd)
  sqrt(drop(t(diff) %*% inv_sig %*% diff)) + eps
}

