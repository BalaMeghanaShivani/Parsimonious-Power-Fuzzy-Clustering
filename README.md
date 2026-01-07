# Parsimonious-Power-Fuzzy-Clustering


This repository contains an implementation of **Parsimonious Power Fuzzy Clustering** with a
systematic exploration of **parsimonious covariance structures** and **fuzziness–distance tradeoffs**.

The goal of this project is to study how different covariance parameterizations and distance scaling
affect clustering stability and accuracy under noisy, heavy-tailed data.

---

## 🔍 Problem

Most clustering algorithms assume:
- fixed distance metrics,
- limited covariance flexibility,
- or clean, well-behaved data.

In practice, real datasets often exhibit:
- overlapping clusters,
- anisotropic covariance structures,
- heavy tails and outliers.

This project builds a **configurable clustering framework** that:
- supports **14 covariance models** (from spherical to fully general),
- decouples **membership fuzziness (m)** from **distance scaling (q)**,
- and evaluates behavior under controlled simulation settings.

---

##  What This System Does

- Implements **14 parsimonious covariance models (F1–F14)** inspired by eigen-decomposition
- Supports **shared vs cluster-specific**:
  - scale (ρ),
  - shape (A, Aₖ),
  - orientation (D, Dₖ)
- Uses **probabilistic distance weighting** instead of likelihood-based EM
- Allows independent control of:
  - fuzziness parameter `m`
  - distance exponent `q`
- Provides full **simulation + evaluation pipeline**

---

##  Core Components

### 1. Clustering Engine
- `pd_clust()` — main iterative clustering routine
- Modular update steps for:
  - cluster centers
  - memberships
  - covariance components
- Objective tracked via **Joint Distance Function (JDF)**

### 2. Covariance Model System
Each model specifies a different structure for Σₖ, including:
- spherical,
- diagonal,
- shared orientation,
- cluster-specific orientation,
- fully general covariance.

Model selection is handled via a unified configuration interface.

### 3. Numerical Stability
To ensure robustness:
- eigenvalue regularization
- determinant normalization
- pseudoinverse-based Mahalanobis distances
- convergence guards and restart logic

---

##  Simulation Framework

The simulation pipeline:
- generates multivariate **Student-t datasets**
- varies:
  - sample size,
  - fuzzifier `m`,
  - distance exponent `q`,
  - covariance generator
- evaluates performance using:
  - Adjusted Rand Index (ARI)
  - Generalized Xie–Beni index (GXB)

This allows controlled comparison across models and parameter regimes.

---

## 📊Key Takeaways (High Level)

- Increasing covariance flexibility improves clustering only when paired with appropriate `m` and `q`
- Fully general models are **not always superior** under heavy noise
- Distance scaling (`q`) has a stronger impact than fuzziness (`m`) in several regimes
- Parsimony often yields more stable solutions than maximum flexibility

---

##  Tech Stack

- **Language:** R
- **Libraries:** mvtnorm, mclust, Matrix, corpcor, cluster, ggplot2
- **Methods:** numerical optimization (BFGS), eigen-decomposition, probabilistic distance metrics

---

## ▶️ Example Usage

```r
X <- as.matrix(iris[,1:4])

fit <- pd_clust(
  X,
  G = 3,
  model = "F12",
  m = 3,
  q = 1,
  iterations = 50,
  centers_init = pam(X, k = 3)$medoids
)

fit$C      # cluster centers
fit$P      # membership matrix
