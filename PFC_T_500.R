# =========================================================
# Probabilistic Distance Clustering 
# =========================================================
#   F1:  Sigma_k = rho * I                       
#   F2:  Sigma_k = rho_k * I                     
#   F3:  Sigma_k = rho * A                     
#   F4:  Sigma_k = rho_k * A                     
#   F5:  Sigma_k = rho * A_k                    
#   F6:  Sigma_k = rho_k * A_k                   
#   F7:  Sigma_k = rho * D * A * D'              
#   F8:  Sigma_k = rho_k * D * A * D'            
#   F9:  Sigma_k = rho * D * A_k * D'            
#   F10: Sigma_k = rho_k * D * A_k * D'          
#   F11: Sigma_k = rho * D_k * A * D_k'          
#   F12: Sigma_k = rho_k * D_k * A * D_k'        
#   F13: Sigma_k = rho * D_k * A_k * D_k'        
#   F14: Sigma_k = rho_k * D_k * A_k * D_k'      
# =========================================================



library(mvtnorm)
library(mclust)
library(Matrix)
library(cluster)
library(corpcor)
library(openxlsx)
library(future.apply)
library(ggplot2)
library(future)
library(dplyr)
library(tibble)


tr  <- function(M) sum(diag(M))
sym <- function(M) (M + t(M))/2

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

#-----------------------------#
# Strictly elimination matrix #
#-----------------------------#
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

#-------------------------------------#
# Stricly half-vectorization operator #
#-------------------------------------#
bvech <- function(A){   # A: square matrix
  n   <- ncol(A)
  res <- Lbar(n)%*%c(A)
  return(res)
}

#--------------#
# From Q to Pl #
#--------------#
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

#--------------#
# From Pl to Q #
#--------------#
PlQ <- function(P,  # Permutation matrix
                l   # Vector of n(n-1)/2 elements
){
  n  <- (1+sqrt(1+8*length(l)))/2
  L  <- diag(n) + matrix(t(Lbar(n))%*%l,n,n) 
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(as.matrix(Q))         # It returns an orthogonal matrix
}


# =========================
#  D / D_k UPDATE FUNCTIONS
# =========================

# Shared D update (for models F7-F10)
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

# Cluster-specific D_k update (for models F11-F14)
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
# MODEL-SPECIFIC UPDATES
# =========================

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

# =========================
# Model configuration 
# =========================
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


initialize_parameters <- function(model, p, G, c, eps, init_params = NULL, X_data ) {
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

PDClustResult <- function(model, C, P, rho_k, Sigma, jdf_values, diagnostics, extras) {
  structure(
    list(model=model, C=C, P=P, rho_k=rho_k, Sigma=Sigma,
         jdf_values=jdf_values, diagnostics=diagnostics, extras=extras),
    class="PDClustResult"
  )
}


pd_clust <- function(
    X, G,
    model = c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14"),
    iterations = 50, q = 2, m = 2, c = 1,
    eps = 1e-8, plot_jdf = FALSE, verbose = FALSE,
    skip_on_increase = TRUE, conv_tol = 1e-8,
    conv_patience = 3, centers_init = NULL, init_params = NULL
) {
  model <- match.arg(model)
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  
  params <- initialize_parameters(model, p, G, c, eps, init_params, X_data = X)
  A <- params$A; A_k <- params$A_k; D <- params$D; D_k <- params$D_k
  rho_k <- params$rho_k; Sigma <- params$Sigma
  
  C <- if (!is.null(centers_init)) centers_init else X[sample(n, G), , drop = FALSE]
  P <- matrix(1/G, n, G)
  jdf_values <- rep(NA_real_, iterations)
  d <- matrix(0, n, G); d_star <- matrix(0, n, G)
  stable_runs <- 0
  min_eig <- 1e-6
  
  for (iter in 1:iterations) {
    
    # Distances (Mahalanobis)
    for (i in 1:n) for (k in 1:G) {
      diff <- X[i,] - C[k,]
      d[i,k] <- calculate_distance(model, diff, A,
                                   if (!is.null(A_k)) A_k[[k]] else NULL,
                                   D,
                                   if (!is.null(D_k)) D_k[[k]] else NULL,
                                   eps, min_eig)
    }
    
    # Membership update
    d_star <- d / matrix(pmax(rho_k, eps), nrow = n, ncol = G, byrow = TRUE)
    recip_d <- 1 / pmax(d_star, eps)
    dp <- recip_d^(q / pmax(m - 1, eps))
    row_sums <- pmax(rowSums(dp), eps)
    P <- dp / row_sums
    P[!is.finite(P)] <- 1/G; P <- P / rowSums(P)
    P_m <- P^m
    
    # Centers
    for (k in 1:G) {
      dc <- (pmax(d[,k], eps))^(q - 2); dc[!is.finite(dc)] <- 1
      rho_q <- (pmax(rho_k[k], eps))^q 
      w_ik <- (P_m[,k] * dc) / rho_q  
      w_ik[!is.finite(w_ik)] <- 0
      num <- colSums(X * w_ik); den <- sum(w_ik) + eps
      C[k,] <- num / den
    }
    
   
    SigmaInv_list <- lapply(Sigma, function(S) corpcor::pseudoinverse(sym(S) + diag(eps, p)))
    Wk <- vector("list", G)
    for (k in 1:G) {
      W <- matrix(0, p, p)
      SigmaInv_k <- SigmaInv_list[[k]]
      for (i in 1:n) {
        diff  <- X[i,] - C[k,]
        quad  <- drop(t(diff) %*% SigmaInv_k %*% diff)
        weight <- (P_m[i,k] / (pmax(rho_k[k], eps)^q)) * (pmax(quad, eps)^(q/2 - 1))
        W <- W + weight * tcrossprod(diff)
      }
      Wk[[k]] <- sym(W)
    }
    
    if (model == "F1") {
      # F1: rho * I 
      rho_k <- rep(c/G, G)
      
    } else if (model == "F2") {
      # F2: rho_k * I 
      Tvec <- sapply(Wk, function(W) sum(diag(W)))
      rho_unscaled <- Tvec^(q/(2*(q+1))); rho_k <- (rho_unscaled / sum(rho_unscaled)) * c
      
    } else if (model == "F3") {
      # F3: rho * A 
      numA <- matrix(0, p, p)
      for (k in 1:G) {
        Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
        Tk <- sum(diag(Ainv %*% Wk[[k]]))
        numA <- numA + (Tk^(q/2 - 1)) * Wk[[k]]
      }
      A <- safe_det_powm(sym(numA), eps, min_eig)
      for (k in 1:G) Sigma[[k]] <- A
      rho_k <- rep(c/G, G)
      
    } else if (model == "F4") {
      # F4: rho_k * A 
      numA <- matrix(0, p, p)
      for (k in 1:G) {
        Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
        Tk <- sum(diag(Ainv %*% Wk[[k]]))
        numA <- numA + (rho_k[k]^(-q)) * (Tk^(q/2 - 1)) * Wk[[k]]
      }
      A <- safe_det_powm(sym(numA), eps, min_eig)
      Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
      Tvec <- sapply(Wk, function(W) sum(diag(Ainv %*% W)))
      rho_unscaled <- Tvec^(q/(2*(q+1))); rho_k <- (rho_unscaled / sum(rho_unscaled)) * c
      for (k in 1:G) Sigma[[k]] <- A
      
    } else if (model == "F5") {
      # F5: rho * A_k 
      for (k in 1:G) {
        Ak_new <- safe_det_powm(Wk[[k]], eps, min_eig)
        A_k[[k]] <- Ak_new
        Sigma[[k]] <- Ak_new + diag(eps, p)
      }
      rho_k <- rep(c/G, G)
      
    } else if (model == "F6") {
      # F6: rho_k * A_k 
      for (k in 1:G) {
        Ak_new <- safe_det_powm(Wk[[k]], eps, min_eig)
        A_k[[k]] <- Ak_new
        Sigma[[k]] <- Ak_new + diag(eps, p)
      }
      Tvec <- sapply(1:G, function(k) {
        Ainv_k <- corpcor::pseudoinverse(A_k[[k]] + diag(eps, p))
        sum(diag(Ainv_k %*% Wk[[k]]))
      })
      rho_unscaled <- pmax(Tvec, eps)^(q/(2*(q+1)))
      rho_k <- (rho_unscaled / sum(rho_unscaled)) * c
      
    } else if (model == "F7") {
      # F7: rho * D * A * D' 
      numA <- matrix(0, p, p)
      for (k in 1:G) {
        Mk <- t(D) %*% Wk[[k]] %*% D
        Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
        Tk <- sum(diag(Ainv %*% Mk))
        numA <- numA + (Tk^(q/2 - 1)) * Mk
      }
      A <- safe_det_powm(sym(numA), eps, min_eig)
      D <- D.update.F7.Red(X = X, C = C, P_m = P_m, A = A, D0 = D, q = q, m = m)
      for (k in 1:G) Sigma[[k]] <- D %*% A %*% t(D)
      rho_k <- rep(c/G, G)
      
    } else if (model == "F8") {
      # F8: rho_k * D * A * D' 
      numA <- matrix(0, p, p)
      for (k in 1:G) {
        Mk <- t(D) %*% Wk[[k]] %*% D
        Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
        Tk <- sum(diag(Ainv %*% Mk))
        numA <- numA + (rho_k[k]^(-q)) * (Tk^(q/2 - 1)) * Mk
      }
      A <- safe_det_powm(sym(numA), eps, min_eig)
      D <- D.update.F8.Red(X = X, C = C, P_m = P_m, A = A, rho_k = rho_k, D0 = D, q = q, m = m)
      Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
      Tvec <- sapply(Wk, function(W) sum(diag(Ainv %*% (t(D) %*% W %*% D))))
      rho_unscaled <- Tvec^(q/(2*(q+1))); rho_k <- (rho_unscaled / sum(rho_unscaled)) * c
      for (k in 1:G) Sigma[[k]] <- D %*% A %*% t(D)
      
    } else if (model == "F9") {
      # Update A_k for each cluster
      for (k in 1:G) {
        A_k[[k]] <- safe_det_powm(Wk[[k]], eps, min_eig)
      }
      
      # Then update shared D using new A_k
      D <- D.update.F9.Red(X = X, C = C, P_m = P_m,
                           A_k_list = A_k, D0 = D, q = q, m = m)
      
      for (k in 1:G) Sigma[[k]] <- D %*% A_k[[k]] %*% t(D)
      rho_k <- rep(c/G, G)
    }
    else if (model == "F10") {
      # Update A_k for each cluster
      for (k in 1:G) {
        A_k[[k]] <- safe_det_powm(Wk[[k]], eps, min_eig)
      }
      
      # Update shared D
      D <- D.update.F10.Red(X = X, C = C, P_m = P_m,
                            A_k_list = A_k, rho_k = rho_k, D0 = D, q = q, m = m)
      
      # Update rho_k 
      Tvec <- sapply(1:G, function(k) {
        Ainv <- corpcor::pseudoinverse(A_k[[k]] + diag(eps, p))
        sum(diag(Ainv %*% (t(D) %*% Wk[[k]] %*% D)))
      })
      rho_unscaled <- Tvec^(q/(2*(q+1)))
      rho_k <- (rho_unscaled / sum(rho_unscaled)) * c
      
      for (k in 1:G) Sigma[[k]] <- D %*% A_k[[k]] %*% t(D)
      
    } else if (model == "F11") {
      # F11: rho * D_k * A * D_k' 
      
      numA <- matrix(0, p, p)
      for (k in 1:G) {
        Mk <- t(D_k[[k]]) %*% Wk[[k]] %*% D_k[[k]]  
        Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
        Tk <- sum(diag(Ainv %*% Mk))
        numA <- numA + (Tk^(q/2 - 1)) * Mk
      }
      A <- safe_det_powm(sym(numA), eps, min_eig)  
      
    
      D_k_new <- Dk.update.F11.Red(X = X, C = C, P_m = P_m, A = A, Dk_list = D_k, q = q, m = m)
      for (k in 1:G) {
        D_k[[k]] <- D_k_new[[k]]
        Sigma[[k]] <- D_k[[k]] %*% A %*% t(D_k[[k]])
      }
      rho_k <- rep(c/G, G)
      
    } else if (model == "F12") {
      # F12: rho_k * D_k * A * D_k' 
      
      numA <- matrix(0, p, p)
      for (k in 1:G) {
        Mk <- t(D_k[[k]]) %*% Wk[[k]] %*% D_k[[k]]  
        Ainv <- corpcor::pseudoinverse(A + diag(eps, p))
        Tk <- sum(diag(Ainv %*% Mk))
        numA <- numA + (rho_k[k]^(-q)) * (Tk^(q/2 - 1)) * Mk
      }
      A <- safe_det_powm(sym(numA), eps, min_eig)  
      
    
      D_k_new <- Dk.update.F12.Red(X = X, C = C, P_m = P_m, A = A, rho_k = rho_k, Dk_list = D_k, q = q, m = m)
      for (k in 1:G) {
        D_k[[k]] <- D_k_new[[k]]
        Sigma[[k]] <- D_k[[k]] %*% A %*% t(D_k[[k]])
      }
      
     
      Tvec <- sapply(1:G, function(k) {
        Sigma_k <- D_k[[k]] %*% A %*% t(D_k[[k]])
        sum(diag(corpcor::pseudoinverse(Sigma_k + diag(eps, p)) %*% Wk[[k]]))
      })
      rho_unscaled <- Tvec^(q/(2*(q+1))); rho_k <- (rho_unscaled / sum(rho_unscaled)) * c
      
    } else if (model == "F13") {
     
      for (k in 1:G) {
        Ak_new <- safe_det_powm(Wk[[k]], eps, min_eig)
        A_k[[k]] <- Ak_new
      }
      
     
      D_k_new <- Dk.update.F13.Red(X = X, C = C, P_m = P_m, 
                                   A_k_list = A_k, Dk_list = D_k, q = q, m = m)
      for (k in 1:G) {
        D_k[[k]] <- D_k_new[[k]]
        Sigma[[k]] <- D_k[[k]] %*% A_k[[k]] %*% t(D_k[[k]])
      }
      rho_k <- rep(c/G, G)
    } else if (model == "F14") {
      # Update A_k for each cluster
      for (k in 1:G) {
        Ak_new <- safe_det_powm(Wk[[k]], eps, min_eig)
        A_k[[k]] <- Ak_new
      }
      
      # Update D_k using the new A_k
      D_k_new <- Dk.update.F14.Red(X = X, C = C, P_m = P_m, 
                                   A_k_list = A_k, rho_k = rho_k, Dk_list = D_k, q = q, m = m)
      for (k in 1:G) {
        D_k[[k]] <- D_k_new[[k]]
        Sigma[[k]] <- D_k[[k]] %*% A_k[[k]] %*% t(D_k[[k]])
      }
      
      # Update rho_k (same as before)
      Tvec <- sapply(1:G, function(k) {
        Ainv <- corpcor::pseudoinverse(A_k[[k]] + diag(eps, p))
        sum(diag(Ainv %*% (t(D_k[[k]]) %*% Wk[[k]] %*% D_k[[k]])))
      })
      rho_unscaled <- Tvec^(q/(2*(q+1)))
      rho_k <- (rho_unscaled / sum(rho_unscaled)) * c
    }
    
    
    # JDF calculation 
    jdf <- 0
    for (i in 1:n) for (k in 1:G) {
      diff <- X[i,] - C[k,]
      d_val <- calculate_distance(model, diff, A,
                                  if (!is.null(A_k)) A_k[[k]] else NULL,
                                  D,
                                  if (!is.null(D_k)) D_k[[k]] else NULL,
                                  eps, min_eig)
      term <- (P_m[i,k] / (pmax(rho_k[k], eps)^q)) * (pmax(d_val, eps)^q)
      if (is.finite(term)) jdf <- jdf + term
    }
    jdf_values[iter] <- jdf
    
    # Convergence check (
    if (verbose && iter %% 10 == 0) cat(sprintf("Iteration %d: JDF = %.6f\n", iter, jdf))
    if (iter > 1 && is.finite(jdf_values[iter-1]) && is.finite(jdf) && iter >= 5) {
      rel_drop <- abs(jdf - jdf_values[iter-1]) / max(abs(jdf_values[iter-1]), eps)
      if (rel_drop < 1e-6) {
        jdf_values <- jdf_values[1:iter]
        if (verbose) message("Converged after ", iter, " iterations")
        break
      }
    }
  }
  
  # Return results 
  jdf_values <- jdf_values[is.finite(jdf_values)]
  diagnostics <- list(
    final_jdf = if (length(jdf_values)) tail(jdf_values, 1) else NA_real_,
    iterations_run = length(jdf_values),
    converged = (length(jdf_values) < iterations)
  )
  extras <- list(d = d, d_star = d_star)
  if (!is.null(A)) extras$A <- A
  if (!is.null(A_k)) extras$A_k_list <- A_k
  if (!is.null(D)) extras$D <- D
  if (!is.null(D_k)) extras$D_k_list <- D_k
  
  result <- PDClustResult(model, C, P, rho_k, Sigma, jdf_values, diagnostics, extras)
  
  if (plot_jdf && length(jdf_values) >= 2) {
    plot(jdf_values, type="l", lwd=2, main=paste0("JDF path (", model, ", m=", m, ", q=", q,")"),
         xlab="Iteration", ylab="JDF")
    lines(cummin(jdf_values), col="red", lty=2)
    legend("topright", legend=c("JDF","Cumulative min"), col=c("black","red"), lty=c(1,2), bty="n")
  }
  if (verbose) print(result)
  result
}
# Example usage 
# X <- as.matrix(iris[,1:4]); G <- 3
# 
# fit <- pd_clust(
#   X, G,
#   model = "F12", iterations = 50, m = 3, q = 1,
#   centers_init = pam(X, k = G)$medoids,
#   verbose = TRUE, plot_jdf = TRUE
# )
# 
# cat("Final JDF:", tail(fit$jdf_values, 1), "\n")
# print(fit$C)
# print(head(fit$P))

# =========================
# SIMULATION FUNCTIONS 
# =========================

build_base_sigmas <- function(p = 4, scale = 0.9) {
  sig1 <- matrix(0.5,  p, p) * scale
  sig2 <- matrix(0.9,  p, p) * scale
  sig3 <- matrix(-0.5, p, p) * scale
  diag(sig1) <- 1; diag(sig2) <- 1; diag(sig3) <- 1
  list(sig1 = sig1, sig2 = sig2, sig3 = sig3)
}

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

GXB <- function(prob, dist, centers, m = 2, q = 2) {
  n <- nrow(prob)
  minD <- min(stats::dist(centers))
  NUM  <- sum((prob ^ m) * (dist ^ q))
  NUM / (n * (minD ^ q))
}

fit_once <- function(
    X, y, model_id, m, q, iterations,
    centers_init, G_sim,
    max_restarts = 5,
    debug = FALSE
) {

  best_result <- NULL
  best_jdf <- Inf
  n <- nrow(X)

  iter_use <- if (q >= 3) iterations * 2 else iterations

  for (restart in 1:max_restarts) {

    if (debug)
      cat(sprintf(" Restart %d/%d...\n", restart, max_restarts))

    if (restart == 1 && !is.null(centers_init)) {

      current_centers <- centers_init

    }  else {

      current_centers <- X[sample(n, G_sim), ]
    }

   
    result <- tryCatch({

      res <- pd_clust(
        X, G = G_sim,
        model = model_id,
        iterations = iter_use,
        m = m, q = q,
        verbose = FALSE,
        plot_jdf = FALSE,
        skip_on_increase = TRUE,
        conv_tol = 1e-8, conv_patience = 3,
        centers_init = current_centers,
        init_params = NULL
      )

      jdf_final <- if (length(res$jdf_values) > 0)
        tail(res$jdf_values, 1)
      else Inf

      U <- res$P
      C <- res$C
      d <- if (!is.null(res$extras) && !is.null(res$extras$d)) res$extras$d else NULL

      hard <- tryCatch(max.col(U, ties.method = "first"),
                       error = function(e) rep(1, length(y)))

      ARI_val <- tryCatch(
        mclust::adjustedRandIndex(hard, y),
        error = function(e) NA_real_
      )

      XB_val <- tryCatch(
        if (!is.null(d)) GXB(U, d, C, m = m, q = q) else NA_real_,
        error = function(e) NA_real_
      )

      list(
        ARI = ARI_val,
        XB = XB_val,
        error = NULL,
        jdf = jdf_final
      )

    }, error = function(e) {

      if (debug)
        cat("        ERROR:", e$message, "\n")

      list(ARI = NA, XB = NA, error = e$message, jdf = Inf)
    })

    if (is.finite(result$jdf) && result$jdf < best_jdf) {
      best_jdf <- result$jdf
      best_result <- result
    }

    if (!is.na(result$ARI)) {
      if (debug)
        cat(" Early stopping — high ARI\n")
      break
    }

  }

  if (!is.null(best_result)) {
    return(list(
      ARI = best_result$ARI,
      XB = best_result$XB,
      error = NULL
    ))
  } else {
    return(list(ARI = NA, XB = NA, error = "All restarts failed"))
  }
}



run_simulation <- function(
    run_tag = "tdist_n500",
    n_values = c(500),
    m_values = c(2,3),
    q_values = c(1,2,3),
    datasets_per_n = 50,
    iterations = 200
) {
  
  gens <- c("rhoI","rhoA","rhoAk","rhoDADt","rhoDAkDt","rhoDkADk","rhoDkAkDk")
  models <- paste0("F", 1:14)
  
  cat("RUNNING FULL SIMULATION for MODELS:", paste(models, collapse=", "), "\n")
  

  cat("[INFO]", format(Sys.time()), "- Running sequentially \n")
  
  # Backup folder
  if (!dir.exists("backup_results")) dir.create("backup_results")
  
  all_detailed_results <- data.frame()
  
  # LOOP over sample sizes
  for (n_val in n_values) {
    cat("", format(Sys.time()), "- Starting simulations for n =", n_val, "\n")
    
    all_combinations <- expand.grid(
      Model = models,
      m = m_values,
      q = q_values,
      stringsAsFactors = FALSE
    )
    
    #all_combinations = data.frame with 84 rows × 3 columns
    
    
    for (gen_code in gens) {
      cat("\n", format(Sys.time()), "- Starting Generator:", gen_code, "\n")
      
      start_time_gen <- Sys.time()
      seeds <- 100 + seq_len(datasets_per_n) 
      
  
      dataset_results_list <- lapply(seeds, function(seed) {
        
        
        dataset <- simulate_dataset(gen_code, n = n_val, p = 4, seed = seed)
        X <- dataset$X      
        y_true <- dataset$y   
        centers_init <- pam(X, k = 3)$medoids  
        
        dataset_df <- data.frame(stringsAsFactors = FALSE)  
        
        # Loop through all model combinations
        for (combo_idx in 1:nrow(all_combinations)) {
          model <- all_combinations$Model[combo_idx]  
          m_val <- all_combinations$m[combo_idx]      
          q_val <- all_combinations$q[combo_idx]     
          
         
          result <- fit_once(
            X, y_true,
            model_id = model,
            m = m_val, q = q_val,
            iterations = iterations,
            centers_init = centers_init,
            G_sim = 3
          )
          
          # Append result to dataset_df (84*7)
          dataset_df <- rbind(dataset_df, data.frame(
            Model = model,
            m = m_val,
            q = q_val,
            Dataset = paste0("Seed_", seed),
            Generator = gen_code,
            n = n_val,
            ARI = result$ARI,  
            XB = result$XB,    
            stringsAsFactors = FALSE
          ))
        }
        
        
        dataset_df
      })
      
      
      gen_results <- do.call(rbind, dataset_results_list)
      all_detailed_results <- rbind(all_detailed_results, gen_results)
      
      #  gen_results = data.frame with 4,200 rows × 7 columns
      # (50 datasets × 84 combinations = 4,200 rows)

      
      # Save backup after each generator
      backup_file <- paste0("backup_results/Backup_Results_", run_tag, "_", gen_code, ".xlsx")
      openxlsx::write.xlsx(all_detailed_results, backup_file, overwrite = TRUE)
      
      end_time_gen <- Sys.time()
      cat("Completed Generator:", gen_code, "in", 
          round(difftime(end_time_gen, start_time_gen, units = "mins"), 2), "minutes\n")
    }
  }
  
  cat("\n", format(Sys.time()), "- Simulation completed for", run_tag, "\n")
  cat("Backup files created in 'backup_results/' folder\n")
  
  invisible(list(
    detailed_results = all_detailed_results
  ))
}


final_partial <- run_simulation(
  run_tag = "tdist_n500",
  n_values = c(500),
  m_values = c(2,3),
  q_values = c(1,2,3),
  datasets_per_n = 50,
  iterations = 200
)
