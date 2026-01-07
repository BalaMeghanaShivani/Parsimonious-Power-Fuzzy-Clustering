# =========================================================
# Package Setup and Loading
# =========================================================
# This file loads all required R packages for the clustering framework.
# Source this file before using any clustering functions.
# =========================================================

# Load required libraries
library(mvtnorm)      # Multivariate t-distribution for data generation
library(mclust)       # Adjusted Rand Index for evaluation
library(Matrix)       # Sparse matrix operations
library(cluster)      # PAM clustering for initialization
library(corpcor)      # Shrinkage covariance estimation and pseudoinverse
library(openxlsx)     # Excel file writing for simulation results
library(future.apply) # Parallel processing (optional)
library(ggplot2)      # Plotting (optional)
library(future)       # Parallel processing backend (optional)
library(dplyr)        # Data manipulation (optional)
library(tibble)       # Data frames (optional)

