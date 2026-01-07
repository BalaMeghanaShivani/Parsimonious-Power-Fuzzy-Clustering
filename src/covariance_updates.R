# =========================================================
# Covariance Matrix Update Functions
# =========================================================
# This module contains functions for updating D and D_k matrices
# for models F7-F14 that involve orthogonal transformations.
# =========================================================
# Note: This module requires functions from distance_utils.R
# which should be sourced before this module.


#' Shared D Matrix Update (for models F7-F10)
#'
#' Updates the shared D matrix by optimizing over the space of orthogonal matrices.
#' Uses parameterization via permutation and lower triangular matrices.
#'
#' @param X Matrix, data matrix (n x p)
#' @param mu_list List of length G, each element is a p-vector (cluster centers)
#' @param p_m_matrix Matrix, membership probabilities raised to power m (n x G)
#' @param A_array Array, A matrices for each cluster (p x p x G)
#' @param Gamma0 Matrix, initial D matrix (orthogonal, p x p)
#' @param rho Numeric vector, rho values for each cluster (length G)
#' @param q Numeric, power parameter
#' @param method Character, optimization method (default: "BFGS")
#' @param iter.max Integer, maximum iterations for optimization (default: 500)
#' @param reltol Numeric, relative tolerance for convergence (default: 1e-15)
#' @param trace Integer, trace level for optimization (default: 0)
#'
#' @return Matrix, updated orthogonal D matrix (p x p)
#'
#' @details
#' The optimization minimizes the objective function over orthogonal matrices
#' by parameterizing them through their LU decomposition.
D.update.Red <- function(X, mu_list, p_m_matrix, A_array, Gamma0, rho, q, 
                             method = "BFGS", iter.max = 500, reltol = 1e-15, trace = 0) {
  
  tempfac <- QPl(Gamma0)
  P <- tempfac$P
  initial.values <- tempfac$l
  
  G <- length(mu_list)
  n <- nrow(X)
  p <- ncol(X)
  
  # Precompute A inverses
  A_inv_array <- array(0, dim = dim(A_array))
  constants <- rho^-q
  for (k in 1:G) {
    A_inv_array[,,k] <- corpcor::pseudoinverse(A_array[,,k] + diag(1e-8, p))
  }
  
  f <- function(par, P, X, mu_list, p_m_matrix, A_inv_array, constants, q) {
    Gamma <- PlQ(P = P, l = par)
    total_objective <- 0
    
    for (k in 1:G) {
      A_k_inv <- A_inv_array[,,k]
      mu_k <- mu_list[[k]]
      constant_k <- constants[k]
      
      for (i in 1:n) {
        if (p_m_matrix[i, k] > 1e-8) { 
          diff <- X[i, ] - mu_k
          quad_form <- t(diff) %*% Gamma %*% A_k_inv %*% t(Gamma) %*% diff
          term <- constant_k * p_m_matrix[i, k] * (pmax(drop(quad_form), 1e-8)^(q/2))
          total_objective <- total_objective + term
        }
      }
    }
    return(total_objective)
  }
  
  res <- optim(par = initial.values, fn = f,
               P = P, X = X, mu_list = mu_list, p_m_matrix = p_m_matrix,
               A_inv_array = A_inv_array, constants = constants, q = q,
               method = method,
               control = list(maxit = iter.max, reltol = reltol, trace = trace))
  
  Gamma_opt <- PlQ(P = P, l = res$par)
  return(Gamma_opt)
}


#' Cluster-Specific D_k Matrix Update (for models F11-F14)
#'
#' Updates a cluster-specific D_k matrix by optimizing over orthogonal matrices.
#'
#' @param X Matrix, data matrix (n x p)
#' @param mu_k Numeric vector, cluster center (length p)
#' @param p_ik_m Numeric vector, membership probabilities for cluster k raised to power m (length n)
#' @param A_k Matrix, A matrix for this cluster (p x p)
#' @param Dk0 Matrix, initial D_k matrix (orthogonal, p x p)
#' @param rho_k Numeric, rho value for this cluster
#' @param q Numeric, power parameter
#' @param method Character, optimization method (default: "BFGS")
#' @param iter.max Integer, maximum iterations for optimization (default: 500)
#' @param reltol Numeric, relative tolerance for convergence (default: 1e-15)
#' @param trace Integer, trace level for optimization (default: 0)
#'
#' @return Matrix, updated orthogonal D_k matrix (p x p)
Dk.update.Red <- function(X, mu_k, p_ik_m, A_k, Dk0, rho_k, q,
                          method = "BFGS", iter.max = 500, reltol = 1e-15, trace = 0) {
  
  tempfac <- QPl(Dk0)
  P <- tempfac$P
  initial.values <- tempfac$l
  
  p <- ncol(X)
  A_k_inv <- corpcor::pseudoinverse(A_k + diag(1e-8, p))
  constant <- (rho_k^-q)
  
  f <- function(par, P, X, mu_k, p_ik_m, A_k_inv, constant, q) {
    Dk <- PlQ(P = P, l = par)
    objective_sum <- 0
    
    for (i in 1:length(p_ik_m)) {
      if (p_ik_m[i] > 1e-8) {  # Skip if membership is negligible
        diff <- X[i, ] - mu_k
        quad_form <- t(diff) %*% Dk %*% A_k_inv %*% t(Dk) %*% diff
        term <- constant * p_ik_m[i] * (pmax(drop(quad_form), 1e-8)^(q/2))
        objective_sum <- objective_sum + term
      }
    }
    return(objective_sum)
  }
  
  res <- optim(par = initial.values, fn = f,
               P = P, X = X, mu_k = mu_k, p_ik_m = p_ik_m, 
               A_k_inv = A_k_inv, constant = constant, q = q,
               method = method,
               control = list(maxit = iter.max, reltol = reltol, trace = trace))
  
  Dk_opt <- PlQ(P = P, l = res$par)
  return(Dk_opt)
}


# =========================
# MODEL-SPECIFIC WRAPPERS
# =========================

#' D Update for Model F7
#'
#' Updates shared D matrix for model F7: rho * D * A * D'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A Matrix, shared A matrix (p x p)
#' @param D0 Matrix, initial D matrix (orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return Matrix, updated D matrix (p x p)
D.update.F7.Red <- function(X, C, P_m, A, D0, q, m = 2) {
  # F7: rho * D * A * D'
  G <- nrow(C)
  p <- ncol(X)
  
  mu_list <- lapply(1:G, function(k) C[k, ])
  A_array <- array(A, dim = c(p, p, G))
  rho <- rep(1, G)  # Shared rho
  
  D_shared <- D.update.Red(
    X = X,
    mu_list = mu_list,
    p_m_matrix = P_m,
    A_array = A_array,
    Gamma0 = D0,
    rho = rho,
    q = q,
    iter.max = 500
  )
  
  return(D_shared)
}


#' D Update for Model F8
#'
#' Updates shared D matrix for model F8: rho_k * D * A * D'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A Matrix, shared A matrix (p x p)
#' @param rho_k Numeric vector, rho values for each cluster (length G)
#' @param D0 Matrix, initial D matrix (orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return Matrix, updated D matrix (p x p)
D.update.F8.Red <- function(X, C, P_m, A, rho_k, D0, q, m = 2) {
  # F8: rho_k * D * A * D'
  G <- nrow(C)
  p <- ncol(X)
  
  mu_list <- lapply(1:G, function(k) C[k, ])
  A_array <- array(A, dim = c(p, p, G))
  
  D_shared <- D.update.Red(
    X = X,
    mu_list = mu_list,
    p_m_matrix = P_m,
    A_array = A_array,
    Gamma0 = D0,
    rho = rho_k,
    q = q,
    iter.max = 500
  )
  
  return(D_shared)
}


#' D Update for Model F9
#'
#' Updates shared D matrix for model F9: rho * D * A_k * D'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A_k_list List, cluster-specific A_k matrices (each p x p)
#' @param D0 Matrix, initial D matrix (orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return Matrix, updated D matrix (p x p)
D.update.F9.Red <- function(X, C, P_m, A_k_list, D0, q, m = 2) {
  # F9: rho * D * A_k * D'
  G <- nrow(C)
  p <- ncol(X)
  
  mu_list <- lapply(1:G, function(k) C[k, ])
  A_array <- array(0, dim = c(p, p, G))
  for (k in 1:G) {
    A_array[,,k] <- A_k_list[[k]]
  }
  rho <- rep(1, G)  # Shared rho
  
  D_shared <- D.update.Red(
    X = X,
    mu_list = mu_list,
    p_m_matrix = P_m,
    A_array = A_array,
    Gamma0 = D0,
    rho = rho,
    q = q,
    iter.max = 500
  )
  
  return(D_shared)
}


#' D Update for Model F10
#'
#' Updates shared D matrix for model F10: rho_k * D * A_k * D'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A_k_list List, cluster-specific A_k matrices (each p x p)
#' @param rho_k Numeric vector, rho values for each cluster (length G)
#' @param D0 Matrix, initial D matrix (orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return Matrix, updated D matrix (p x p)
D.update.F10.Red <- function(X, C, P_m, A_k_list, rho_k, D0, q, m = 2) {
  # F10: rho_k * D * A_k * D'
  G <- nrow(C)
  p <- ncol(X)
  
  mu_list <- lapply(1:G, function(k) C[k, ])
  A_array <- array(0, dim = c(p, p, G))
  for (k in 1:G) {
    A_array[,,k] <- A_k_list[[k]]
  }
  
  D_shared <- D.update.Red(
    X = X,
    mu_list = mu_list,
    p_m_matrix = P_m,
    A_array = A_array,
    Gamma0 = D0,
    rho = rho_k,
    q = q,
    iter.max = 500
  )
  
  return(D_shared)
}


#' D_k Update for Model F11
#'
#' Updates cluster-specific D_k matrices for model F11: rho * D_k * A * D_k'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A Matrix, shared A matrix (p x p)
#' @param Dk_list List, initial D_k matrices (each orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return List, updated D_k matrices (each p x p)
Dk.update.F11.Red <- function(X, C, P_m, A, Dk_list, q, m = 2) {
  # F11: rho * D_k * A * D_k'
  G <- nrow(C)
  
  Dk_updated <- list()
  for (k in 1:G) {
    Dk_updated[[k]] <- Dk.update.Red(
      X = X,
      mu_k = C[k, ],
      p_ik_m = P_m[, k],
      A_k = A,
      Dk0 = Dk_list[[k]],
      rho_k = 1,  # Shared rho
      q = q,
      iter.max = 500
    )
  }
  return(Dk_updated)
}


#' D_k Update for Model F12
#'
#' Updates cluster-specific D_k matrices for model F12: rho_k * D_k * A * D_k'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A Matrix, shared A matrix (p x p)
#' @param rho_k Numeric vector, rho values for each cluster (length G)
#' @param Dk_list List, initial D_k matrices (each orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return List, updated D_k matrices (each p x p)
Dk.update.F12.Red <- function(X, C, P_m, A, rho_k, Dk_list, q, m = 2) {
  # F12: rho_k * D_k * A * D_k'
  G <- nrow(C)
  
  Dk_updated <- list()
  for (k in 1:G) {
    Dk_updated[[k]] <- Dk.update.Red(
      X = X,
      mu_k = C[k, ],
      p_ik_m = P_m[, k],
      A_k = A,
      Dk0 = Dk_list[[k]],
      rho_k = rho_k[k],
      q = q,
      iter.max = 500
    )
  }
  return(Dk_updated)
}


#' D_k Update for Model F13
#'
#' Updates cluster-specific D_k matrices for model F13: rho * D_k * A_k * D_k'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A_k_list List, cluster-specific A_k matrices (each p x p)
#' @param Dk_list List, initial D_k matrices (each orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return List, updated D_k matrices (each p x p)
Dk.update.F13.Red <- function(X, C, P_m, A_k_list, Dk_list, q, m = 2) {
  # F13: rho * D_k * A_k * D_k'
  G <- nrow(C)
  
  Dk_updated <- list()
  for (k in 1:G) {
    Dk_updated[[k]] <- Dk.update.Red(
      X = X,
      mu_k = C[k, ],
      p_ik_m = P_m[, k],
      A_k = A_k_list[[k]],
      Dk0 = Dk_list[[k]],
      rho_k = 1,  # Shared rho
      q = q,
      iter.max = 500
    )
  }
  return(Dk_updated)
}


#' D_k Update for Model F14
#'
#' Updates cluster-specific D_k matrices for model F14: rho_k * D_k * A_k * D_k'
#'
#' @param X Matrix, data matrix (n x p)
#' @param C Matrix, cluster centers (G x p)
#' @param P_m Matrix, membership probabilities raised to power m (n x G)
#' @param A_k_list List, cluster-specific A_k matrices (each p x p)
#' @param rho_k Numeric vector, rho values for each cluster (length G)
#' @param Dk_list List, initial D_k matrices (each orthogonal, p x p)
#' @param q Numeric, power parameter
#' @param m Numeric, fuzziness parameter (default: 2, not used in this function)
#' @return List, updated D_k matrices (each p x p)
Dk.update.F14.Red <- function(X, C, P_m, A_k_list, rho_k, Dk_list, q, m = 2) {
  # F14: rho_k * D_k * A_k * D_k'
  G <- nrow(C)
  
  Dk_updated <- list()
  for (k in 1:G) {
    Dk_updated[[k]] <- Dk.update.Red(
      X = X,
      mu_k = C[k, ],
      p_ik_m = P_m[, k],
      A_k = A_k_list[[k]],
      Dk0 = Dk_list[[k]],
      rho_k = rho_k[k],
      q = q,
      iter.max = 500
    )
  }
  return(Dk_updated)
}

