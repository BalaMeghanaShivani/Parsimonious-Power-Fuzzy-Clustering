# =========================================================
# Data Generation for Simulations
# =========================================================
# This module contains functions for generating synthetic datasets
# with various covariance structures for simulation studies.
# =========================================================


#' Build Base Covariance Matrices
#'
#' Creates three base covariance matrices with different correlation structures
#' for use in data generation.
#'
#' @param p Integer, number of dimensions (default: 4)
#' @param scale Numeric, scaling factor for off-diagonal elements (default: 0.9)
#'
#' @return List containing three covariance matrices:
#'   \describe{
#'     \item{sig1}{Matrix with correlation 0.5}
#'     \item{sig2}{Matrix with correlation 0.9}
#'     \item{sig3}{Matrix with correlation -0.5}
#'   }
#'
#' @details
#' All matrices have unit diagonal elements and scaled off-diagonal elements.
build_base_sigmas <- function(p = 4, scale = 0.9) {
  sig1 <- matrix(0.5,  p, p) * scale
  sig2 <- matrix(0.9,  p, p) * scale
  sig3 <- matrix(-0.5, p, p) * scale
  diag(sig1) <- 1; diag(sig2) <- 1; diag(sig3) <- 1
  list(sig1 = sig1, sig2 = sig2, sig3 = sig3)
}


#' Generate Covariance Matrices for Data Generation
#'
#' Creates covariance matrices based on a generation code that specifies
#' the covariance structure type.
#'
#' @param gen_code Character string, generation code specifying structure:
#'   \itemize{
#'     \item "rhoI": Spherical (rho * I)
#'     \item "rhoA": Shared A matrix (rho * A)
#'     \item "rhoAk": Cluster-specific A_k matrices (rho * A_k)
#'     \item "rhoDADt": Shared D and A (rho * D * A * D')
#'     \item "rhoDAkDt": Shared D, cluster A_k (rho * D * A_k * D')
#'     \item "rhoDkADk": Cluster D_k, shared A (rho * D_k * A * D_k')
#'     \item "rhoDkAkDk": Fully general (rho * D_k * A_k * D_k')
#'   }
#' @param p Integer, number of dimensions (default: 4)
#' @param scale_A Numeric, scaling factor for A matrices (default: 0.9)
#'
#' @return List of K covariance matrices (K=3), each p x p
#'
#' @details
#' The function decomposes base covariance matrices into eigenstructures
#' and reconstructs them according to the specified generation code.
make_sigmas_for_gen <- function(gen_code, p = 4, scale_A = 0.9) {
  symM <- function(S) (S + t(S)) / 2
  bases <- build_base_sigmas(p = p, scale = 0.9)
  sig1 <- bases$sig1; sig2 <- bases$sig2; sig3 <- bases$sig3
  K <- 3

  e1 <- eigen(sig1, symmetric = TRUE); D1 <- e1$vectors; A1 <- diag(e1$values)
  e2 <- eigen(sig2, symmetric = TRUE); D2 <- e2$vectors; A2 <- diag(e2$values)
  e3 <- eigen(sig3, symmetric = TRUE); D3 <- e3$vectors; A3 <- diag(e3$values)

  I_p <- diag(p)

  if (gen_code == "rhoI") {
    S <- 0.9 * I_p
    return(rep(list(symM(S)), K))
  }
  if (gen_code == "rhoDkAkDk") {
    return(lapply(list(sig1, sig2, sig3), symM))
  }
  if (gen_code == "rhoDADt") {
    S <- 0.9 * (D2 %*% A2 %*% t(D2))
    return(rep(list(symM(S)), K))
  }
  if (gen_code == "rhoDkADk") {
    Ak <- A2
    S1 <- 0.9 * (D1 %*% Ak %*% t(D1))
    S2 <- 0.9 * (D2 %*% Ak %*% t(D2))
    S3 <- 0.9 * (D3 %*% Ak %*% t(D3))
    return(lapply(list(S1, S2, S3), symM))
  }
  if (gen_code == "rhoDAkDt") {
    S1 <- 0.9 * (D2 %*% A1 %*% t(D2))
    S2 <- 0.9 * (D2 %*% A2 %*% t(D2))
    S3 <- 0.9 * (D2 %*% A3 %*% t(D2))
    return(lapply(list(S1, S2, S3), symM))
  }
  if (gen_code == "rhoAk") {
    S1 <- 0.9 * diag(diag(A1))
    S2 <- 0.9 * diag(diag(A2))
    S3 <- 0.9 * diag(diag(A3))
    return(lapply(list(S1, S2, S3), symM))
  }
  if (gen_code == "rhoA") {
    S <- 0.9 * diag(diag(A2))
    return(rep(list(symM(S)), K))
  }
  stop("Unknown generator: ", gen_code)
}


#' Simulate Dataset with Specified Covariance Structure
#'
#' Generates a synthetic dataset with three clusters using multivariate t-distributions
#' with specified covariance structures.
#'
#' @param gen_code Character string, generation code (see \code{\link{make_sigmas_for_gen}})
#' @param n Integer, total number of observations
#' @param p Integer, number of dimensions (default: 4)
#' @param seed Integer or NULL, random seed for reproducibility
#'
#' @return List containing:
#'   \describe{
#'     \item{X}{Matrix, data matrix (n x p)}
#'     \item{y}{Integer vector, true cluster labels (length n)}
#'   }
#'
#' @details
#' The dataset consists of three clusters with sizes approximately 50%, 30%, and 20%
#' of the total sample size. Each cluster is generated from a multivariate t-distribution
#' with 3 degrees of freedom and cluster-specific mean and covariance.
simulate_dataset <- function(gen_code, n, p = 4, seed = NULL) {
  set.seed(seed)
  Sigmas <- make_sigmas_for_gen(gen_code, p = p, scale_A = 0.9)
  mu_list <- list(c(0,0,0,0), c(0,3,0,3), c(3,0,3,0))
  n1 <- as.integer(round(n * 0.5))
  n2 <- as.integer(round(n * 0.3))
  n3 <- n - n1 - n2
  X1 <- mvtnorm::rmvt(n1, sigma = Sigmas[[1]], df = 3, delta = mu_list[[1]])
  X2 <- mvtnorm::rmvt(n2, sigma = Sigmas[[2]], df = 3, delta = mu_list[[2]])
  X3 <- mvtnorm::rmvt(n3, sigma = Sigmas[[3]], df = 3, delta = mu_list[[3]])


  X <- rbind(X1, X2, X3)
  y <- c(rep(1, n1), rep(2, n2), rep(3, n3))
  list(X = X, y = y)
}


#' Xie-Beni Index
#'
#' Computes the Xie-Beni index for cluster validation.
#'
#' @param prob Matrix, membership probabilities (n x G)
#' @param dist Matrix, distances from data points to cluster centers (n x G)
#' @param centers Matrix, cluster centers (G x p)
#' @param m Numeric, fuzziness parameter (default: 2)
#' @param q Numeric, power parameter (default: 2)
#'
#' @return Numeric scalar, Xie-Beni index value
#'
#' @details
#' The Xie-Beni index measures cluster compactness and separation.
#' Lower values indicate better clustering.
GXB <- function(prob, dist, centers, m = 2, q = 2) {
  n <- nrow(prob)
  minD <- min(stats::dist(centers))
  NUM  <- sum((prob ^ m) * (dist ^ q))
  NUM / (n * (minD ^ q))
}

