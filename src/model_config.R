# =========================================================
# Model Configuration and Parameter Initialization
# =========================================================
# This module contains model configuration definitions and
# parameter initialization functions for Probabilistic Distance Clustering.
# =========================================================

#' Model Configuration for Probabilistic Distance Clustering
#'
#' Defines 14 model variants (F1-F14) with different covariance structures.
#' Each model specifies which parameters are shared or cluster-specific.
#'
#' @format A named list with 14 elements (F1-F14), each containing:
#' \describe{
#'   \item{description}{Character description of the model}
#'   \item{has_A}{Logical: whether model uses shared A matrix}
#'   \item{has_A_k}{Logical: whether model uses cluster-specific A_k matrices}
#'   \item{has_D}{Logical: whether model uses shared D matrix}
#'   \item{has_D_k}{Logical: whether model uses cluster-specific D_k matrices}
#'   \item{has_rho}{Logical: whether model uses rho parameter}
#'   \item{shared_rho}{Logical: whether rho is shared across clusters}
#' }
#'
#' @details
#' Models F1-F14 represent different covariance structures:
#' \itemize{
#'   \item F1: Spherical, shared rho (Sigma_k = rho * I)
#'   \item F2: Spherical, cluster-specific rho (Sigma_k = rho_k * I)
#'   \item F3: Shared A, shared rho (Sigma_k = rho * A)
#'   \item F4: Shared A, cluster rho (Sigma_k = rho_k * A)
#'   \item F5: Cluster A_k, shared rho (Sigma_k = rho * A_k)
#'   \item F6: Cluster A_k, cluster rho (Sigma_k = rho_k * A_k)
#'   \item F7: Shared D, A, rho (Sigma_k = rho * D * A * D')
#'   \item F8: Shared D, A; cluster rho (Sigma_k = rho_k * D * A * D')
#'   \item F9: Shared D; cluster A_k; shared rho (Sigma_k = rho * D * A_k * D')
#'   \item F10: Shared D; cluster A_k; cluster rho (Sigma_k = rho_k * D * A_k * D')
#'   \item F11: Cluster D_k; shared A, rho (Sigma_k = rho * D_k * A * D_k')
#'   \item F12: Cluster D_k; shared A; cluster rho (Sigma_k = rho_k * D_k * A * D_k')
#'   \item F13: Cluster D_k, A_k; shared rho (Sigma_k = rho * D_k * A_k * D_k')
#'   \item F14: Fully general (Sigma_k = rho_k * D_k * A_k * D_k')
#' }
ModelConfig <- list(
  F1  = list(description="Spherical, shared rho",
             has_A=FALSE, has_A_k=FALSE, has_D=FALSE, has_D_k=FALSE, has_rho=TRUE,  shared_rho=TRUE),
  F2  = list(description="Spherical, cluster-specific rho",
             has_A=FALSE, has_A_k=FALSE, has_D=FALSE, has_D_k=FALSE, has_rho=TRUE,  shared_rho=FALSE),
  F3  = list(description="Shared A, shared rho",
             has_A=TRUE,  has_A_k=FALSE, has_D=FALSE, has_D_k=FALSE, has_rho=TRUE,  shared_rho=TRUE),
  F4  = list(description="Shared A, cluster rho",
             has_A=TRUE,  has_A_k=FALSE, has_D=FALSE, has_D_k=FALSE, has_rho=TRUE,  shared_rho=FALSE),
  F5  = list(description="Cluster A_k, shared rho",
             has_A=FALSE, has_A_k=TRUE,  has_D=FALSE, has_D_k=FALSE, has_rho=TRUE,  shared_rho=TRUE),
  F6  = list(description="Cluster A_k, cluster rho",
             has_A=FALSE, has_A_k=TRUE,  has_D=FALSE, has_D_k=FALSE, has_rho=TRUE,  shared_rho=FALSE),
  F7  = list(description="Shared D, A, rho",
             has_A=TRUE,  has_A_k=FALSE, has_D=TRUE,  has_D_k=FALSE, has_rho=TRUE,  shared_rho=TRUE),
  F8  = list(description="Shared D, A; cluster rho",
             has_A=TRUE,  has_A_k=FALSE, has_D=TRUE,  has_D_k=FALSE, has_rho=TRUE,  shared_rho=FALSE),
  F9  = list(description="Shared D; cluster A_k; shared rho",
             has_A=FALSE, has_A_k=TRUE,  has_D=TRUE,  has_D_k=FALSE, has_rho=TRUE,  shared_rho=TRUE),
  F10 = list(description="Shared D; cluster A_k; cluster rho",
             has_A=FALSE, has_A_k=TRUE,  has_D=TRUE,  has_D_k=FALSE, has_rho=TRUE,  shared_rho=FALSE),
  F11 = list(description="Cluster D_k; shared A, rho",
             has_A=TRUE,  has_A_k=FALSE, has_D=FALSE, has_D_k=TRUE,  has_rho=TRUE,  shared_rho=TRUE),
  F12 = list(description="Cluster D_k; shared A; cluster rho",
             has_A=TRUE,  has_A_k=FALSE, has_D=FALSE, has_D_k=TRUE,  has_rho=TRUE,  shared_rho=FALSE),
  F13 = list(description="Cluster D_k, A_k; shared rho",
             has_A=FALSE, has_A_k=TRUE,  has_D=FALSE, has_D_k=TRUE,  has_rho=TRUE,  shared_rho=TRUE),
  F14 = list(description="Fully general (cluster D_k, A_k, rho)",
             has_A=FALSE, has_A_k=TRUE,  has_D=FALSE, has_D_k=TRUE,  has_rho=TRUE,  shared_rho=FALSE)
)


#' Initialize Parameters for Probabilistic Distance Clustering
#'
#' Initializes model parameters based on the specified model configuration.
#' Parameters are either initialized from provided values or estimated from data.
#'
#' @param model Character string specifying the model (F1-F14)
#' @param p Integer, number of dimensions
#' @param G Integer, number of clusters
#' @param c Numeric, scaling parameter for rho initialization
#' @param eps Numeric, small positive value for numerical stability
#' @param init_params Optional list with pre-initialized parameters:
#'   \itemize{
#'     \item A_shared: Shared A matrix (if applicable)
#'     \item A_k: List of cluster-specific A_k matrices (if applicable)
#'     \item D_shared: Shared D matrix (if applicable)
#'     \item D_k: List of cluster-specific D_k matrices (if applicable)
#'   }
#' @param X_data Matrix, data matrix used for initialization if init_params is NULL
#'
#' @return List containing initialized parameters:
#' \describe{
#'   \item{A}{Shared A matrix or NULL}
#'   \item{A_k}{List of cluster-specific A_k matrices or NULL}
#'   \item{D}{Shared D matrix or NULL}
#'   \item{D_k}{List of cluster-specific D_k matrices or NULL}
#'   \item{rho_k}{Numeric vector of rho values for each cluster}
#'   \item{Sigma}{List of covariance matrices for each cluster}
#' }
#'
#' @details
#' If init_params is NULL, parameters are initialized from the data using
#' shrinkage covariance estimation. Otherwise, provided parameters are used.
initialize_parameters <- function(model, p, G, c, eps, init_params = NULL, X_data) {
  config <- ModelConfig[[model]]
  params <- list()
  
  if (!is.null(init_params)) {
    A_shared <- init_params$A_shared
    A_k_in   <- init_params$A_k
    D_shared <- init_params$D_shared
    D_k_in   <- init_params$D_k
  } else {
    baseSig <- cov.shrink(X_data, verbose = FALSE)
    eig0    <- eigen(baseSig, symmetric = TRUE)
    evals   <- pmax(eig0$values, 1e-6)
    baseA   <- diag(evals / (abs(prod(evals))^(1/p)))
    A_shared <- baseA
    A_k_in   <- replicate(G, baseA, simplify = FALSE)
    D_shared <- eig0$vectors
    D_k_in   <- replicate(G, eig0$vectors, simplify = FALSE)
  }
  
  params$A     <- if (config$has_A)    A_shared else NULL
  params$A_k   <- if (config$has_A_k)  A_k_in   else NULL
  params$D     <- if (config$has_D)    D_shared else NULL
  params$D_k   <- if (config$has_D_k)  D_k_in   else NULL
  params$rho_k <- if (config$has_rho)  rep(c/G, G) else rep(1, G)
  
  params$Sigma <- lapply(1:G, function(i) {
    if (model %in% c("F1","F2")) return(diag(p))
    if (model %in% c("F3","F4"))  return(params$A)
    if (model %in% c("F5","F6"))  return(params$A_k[[i]])
    if (model %in% c("F7","F8"))  return(params$D %*% params$A %*% t(params$D))
    if (model %in% c("F9","F10")) return(params$D %*% params$A_k[[i]] %*% t(params$D))
    if (model %in% c("F11","F12")) return(params$D_k[[i]] %*% params$A %*% t(params$D_k[[i]]))
    if (model %in% c("F13","F14")) return(params$D_k[[i]] %*% params$A_k[[i]] %*% t(params$D_k[[i]]))
    diag(p)
  })
  
  params
}


#' Create PDClustResult Object
#'
#' Constructs a result object for Probabilistic Distance Clustering.
#'
#' @param model Character string, model identifier (F1-F14)
#' @param C Matrix, cluster centers (G x p)
#' @param P Matrix, membership probabilities (n x G)
#' @param rho_k Numeric vector, rho values for each cluster
#' @param Sigma List, covariance matrices for each cluster
#' @param jdf_values Numeric vector, objective function values across iterations
#' @param diagnostics List, diagnostic information about the fit
#' @param extras List, additional results (distances, parameter matrices, etc.)
#'
#' @return Object of class "PDClustResult" containing all clustering results
PDClustResult <- function(model, C, P, rho_k, Sigma, jdf_values, diagnostics, extras) {
  structure(
    list(model=model, C=C, P=P, rho_k=rho_k, Sigma=Sigma,
         jdf_values=jdf_values, diagnostics=diagnostics, extras=extras),
    class="PDClustResult"
  )
}

