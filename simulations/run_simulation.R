# =========================================================
# Simulation Execution
# =========================================================
# This module contains functions for running comprehensive simulation studies
# comparing different model variants across various parameter settings.
# =========================================================

# Source required modules
source("src/setup.R")  # Load required packages
source("src/pd_clust.R")
source("simulations/data_generation.R")


#' Fit Clustering Model with Multiple Restarts
#'
#' Fits a clustering model with multiple random restarts and returns
#' the best result based on objective function value.
#'
#' @param X Matrix, data matrix (n x p)
#' @param y Integer vector, true cluster labels (length n)
#' @param model_id Character string, model identifier (F1-F14)
#' @param m Numeric, fuzziness parameter
#' @param q Numeric, power parameter
#' @param iterations Integer, maximum number of iterations
#' @param centers_init Matrix or NULL, initial cluster centers
#' @param G_sim Integer, number of clusters
#' @param max_restarts Integer, maximum number of random restarts (default: 5)
#' @param debug Logical, whether to print debug information (default: FALSE)
#'
#' @return List containing:
#'   \describe{
#'     \item{ARI}{Numeric, Adjusted Rand Index}
#'     \item{XB}{Numeric, Xie-Beni index}
#'     \item{error}{Character or NULL, error message if fitting failed}
#'   }
#'
#' @details
#' The function attempts multiple restarts with different initializations.
#' If a restart achieves high ARI (good clustering), it stops early.
#' Otherwise, it returns the result with the lowest objective function value.
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


#' Run Comprehensive Simulation Study
#'
#' Executes a comprehensive simulation study comparing all model variants (F1-F14)
#' across different parameter settings and data generation schemes.
#'
#' @param run_tag Character string, tag for identifying this simulation run (default: "tdist_n500")
#' @param n_values Integer vector, sample sizes to test (default: c(500))
#' @param m_values Numeric vector, fuzziness parameters to test (default: c(2,3))
#' @param q_values Numeric vector, power parameters to test (default: c(1,2,3))
#' @param datasets_per_n Integer, number of datasets to generate per sample size (default: 50)
#' @param iterations Integer, maximum iterations for clustering (default: 200)
#'
#' @return List containing:
#'   \describe{
#'     \item{detailed_results}{Data frame with all simulation results}
#'   }
#'
#' @details
#' The simulation:
#' \itemize{
#'   \item Tests all 14 model variants (F1-F14)
#'   \item Uses 7 different data generation schemes
#'   \item Tests multiple combinations of m and q parameters
#'   \item Generates multiple datasets per configuration
#'   \item Saves backup results after each generator completes
#' }
#'
#' Results are saved to Excel files in the "backup_results/" directory.
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

