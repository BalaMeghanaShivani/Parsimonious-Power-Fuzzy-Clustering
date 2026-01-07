# =========================================================
# Probabilistic Distance Clustering - Main Function
# =========================================================
# This module contains the main pd_clust function that performs
# fuzzy clustering with various covariance structures.
# =========================================================

# Source required modules
source("src/model_config.R")
source("src/distance_utils.R")
source("src/covariance_updates.R")

# Note: Ensure required packages are loaded (see src/setup.R)


#' Probabilistic Distance Clustering
#'
#' Performs fuzzy clustering using probabilistic distance-based objective function
#' with various covariance structures (14 model variants F1-F14).
#'
#' @param X Matrix or data frame, data matrix (n x p)
#' @param G Integer, number of clusters
#' @param model Character string, model identifier. One of:
#'   "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", "F13", "F14"
#' @param iterations Integer, maximum number of iterations (default: 50)
#' @param q Numeric, power parameter for distance (default: 2)
#' @param m Numeric, fuzziness parameter (default: 2)
#' @param c Numeric, scaling parameter for rho initialization (default: 1)
#' @param eps Numeric, small positive value for numerical stability (default: 1e-8)
#' @param plot_jdf Logical, whether to plot objective function values (default: FALSE)
#' @param verbose Logical, whether to print progress information (default: FALSE)
#' @param skip_on_increase Logical, whether to stop if objective increases (default: TRUE)
#' @param conv_tol Numeric, convergence tolerance (default: 1e-8)
#' @param conv_patience Integer, patience for convergence (default: 3)
#' @param centers_init Matrix or NULL, initial cluster centers (G x p). If NULL, randomly sampled.
#' @param init_params List or NULL, initial parameter values. See \code{\link{initialize_parameters}}
#'
#' @return Object of class "PDClustResult" containing:
#' \describe{
#'   \item{model}{Character, model identifier}
#'   \item{C}{Matrix, cluster centers (G x p)}
#'   \item{P}{Matrix, membership probabilities (n x G)}
#'   \item{rho_k}{Numeric vector, rho values for each cluster}
#'   \item{Sigma}{List, covariance matrices for each cluster}
#'   \item{jdf_values}{Numeric vector, objective function values across iterations}
#'   \item{diagnostics}{List, diagnostic information}
#'   \item{extras}{List, additional results (distances, parameter matrices)}
#' }
#'
#' @details
#' The algorithm alternates between updating:
#' \itemize{
#'   \item Membership probabilities (P)
#'   \item Cluster centers (C)
#'   \item Covariance parameters (A, A_k, D, D_k, rho_k)
#' }
#'
#' The objective function (JDF) is minimized iteratively until convergence
#' or maximum iterations are reached.
#'
#' @examples
#' \dontrun{
#' # Example with iris data
#' X <- as.matrix(iris[,1:4])
#' G <- 3
#' 
#' fit <- pd_clust(
#'   X, G,
#'   model = "F12", iterations = 50, m = 3, q = 1,
#'   centers_init = pam(X, k = G)$medoids,
#'   verbose = TRUE, plot_jdf = TRUE
#' )
#' 
#' cat("Final JDF:", tail(fit$jdf_values, 1), "\n")
#' print(fit$C)
#' print(head(fit$P))
#' }
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

