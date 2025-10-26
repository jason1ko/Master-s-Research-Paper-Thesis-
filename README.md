# Asymptotic Error and Regularization Path in High-Dimensional LDA: A Random Matrix Approach

This repository accompanies my master's research paper, which studies **regularized discriminant analysis (rLDA)** in the **high-dimensional, low-sample-size (HDLSS)** regime using tools from **random matrix theory (RMT)**.  
It contains selected **R codes** that reproduce the theoretical derivations, numerical experiments, and figures presented in the paper.

---

## 📘 Paper Overview

**Title:** *Regularized Discriminant Analysis in the HDLSS Regime: A Random Matrix Perspective*  
**Author:** Han Jun Ko  
**Date:** January 2025  
**Institution:** Yonsei University, Department of Statistics  

### Abstract (Summary)

This work analyzes the asymptotic classification error of ridge-regularized LDA in the HDLSS regime.  
We establish explicit limiting forms for the misclassification risk:

$$
\mathrm{Err}(\hat{w}_\lambda) \to \Phi(-\Theta(\lambda))
$$

and provide theoretical and numerical characterization of the optimal regularization parameter λ as a function of the aspect ratio  
γ = p/n and correlation strength ρ in Toeplitz covariance structures (AR(1), Matérn).

---

## 🧮 Key Contributions

1. **Closed-form limit under Σ = I**  
   - Derived a piecewise closed-form expression for the MDP limit Θ(0) across regimes γ < 1, γ = 1, γ > 1.  
   - Theoretical error curve Φ(−Θ(0)) peaks at γ = 1, reconciling classical HDLSS observations.

2. **Toeplitz Covariance Validation**  
   - Verified via **Szegő’s theorem** that AR(1) and Matérn covariances satisfy RMT spectral assumptions.

3. **Stable Approximation Pipeline**  
   - Implemented a numerically stable computation of Θ(λ) using empirical spectral sums.  
   - Mapped **λ_opt(γ, ρ)** and analyzed the transition between MDP (λ=0) and MD (λ=∞).

---

## 🧠 Theoretical Framework

Regularized LDA forms a continuum between the **Maximal Data Piling (MDP)** and **Mean Difference (MD)** rules:

$$
v(\lambda) = (CC^\top + \lambda I)^{-1} w, \quad \lambda > 0.
$$

As λ varies:
- λ = 0 → Maximal Data Piling (MDP)
- λ = ∞ → Mean Difference (MD)
- Intermediate λ → Ridge-regularized direction connecting both ends

The limiting misclassification risk depends on the **Stieltjes transform** of the sample covariance spectral distribution (Marchenko–Pastur law).

---

## 📊 Repository Contents

| File | Description |
|------|--------------|
| `ar1_example.R` | plot of $\Theta(\lambda)$, from Dobriban, E. & Wager, S. (2018) |
| `rda_functions.R` | calculate Stieltjes functions, from Dobriban, E. & Wager, S. (2018) |
| `rda_rcpp2.cpp` | generate multivariate Normal using Rcpp |
| `results.R` | includes all codes and functions needed for the paper |

---

## 🔍 Numerical Patterns Observed

- $\lambda_{opt}$ **increases** with γ and **decreases** with ρ.  
- $\Theta(\lambda_{opt})$ is close to $\Theta(0)$ at small γ and approaches $\Theta(\infty)$ as γ grows.  
- Under Matérn covariance, $\Phi(−\Theta(\lambda_{opt}))$ and $\Phi(−\Theta(\infty))$ nearly coincide, indicating smoother effective regularization.

---

## 🔗 Key References

- Dobriban, E. & Wager, S. (2018). High-dimensional asymptotics of prediction: Ridge regression and classification. Annals of Statistics.
- Ahn, J. & Marron, J. S. (2010). The maximal data piling direction for discrimination. Biometrika.
- Lee, M. H., Ahn, J., & Jeon, Y. (2013). HDLSS discrimination with adaptive data piling. JCGS.
- Bai, Z. & Silverstein, J. (2010). Spectral Analysis of Large Dimensional Random Matrices. Springer.

---
