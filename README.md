# Regularized LDA in the HDLSS Regime: Random-Matrix Limits & Toeplitz Covariances

> R code accompanying my master’s research paper on regularized discriminant analysis (rLDA) in high-dimension, low-sample-size (HDLSS) settings. The project connects the MDP–MD continuum to asymptotic classification risk using random matrix theory, and provides a stable pipeline to compute the margin $\Theta(\lambda)$ and the optimal ridge $\lambda_{\text{opt}}(\gamma,\rho)$ for Toeplitz covariances (AR(1), Matérn).

**Paper**: Ko Han Jun (Jan 2025), *Regularized Discriminant Analysis in HDLSS via Random Matrix Limits* (working title). :contentReference[oaicite:0]{index=0}

---

## 🚀 What’s in this repo?

- **Reusable R functions** to compute Stieltjes-based quantities and the margin
  $\Theta(\lambda)$, then trace the error $\Phi(-\Theta(\lambda))$ along the regularization path.
- **Experiments for AR(1) and Matérn** Toeplitz covariances: maps of $\lambda_{\text{opt}}(\gamma,\rho)$
  and comparisons with the MDP ($\lambda=0$) and MD ($\lambda=\infty$) ends.
- **Figure scripts** to reproduce the main plots (risk curves, optimal-$\lambda$ trends).

> Key takeaways from the paper:
> - Under $\Sigma=I$, we derive a **piecewise closed form** for the MDP limit $\Theta(0)$ across $\gamma<1$, $\gamma=1$, $\gamma>1$, producing an error curve peaking at $\gamma=1$.
> - **Toeplitz families (AR(1), Matérn)** satisfy the spectral assumptions for the limiting risk (via Szegő arguments), justifying their use.
> - We provide a **numerically stable pipeline** to evaluate $\Theta(\lambda)$ and **map $\lambda_{\text{opt}}(\gamma,\rho)$**:
>   - $\lambda_{\text{opt}}$ ↑ as aspect ratio $\gamma$ ↑; $\lambda_{\text{opt}}$ ↓ as correlation strength $\rho$ ↑.
>   - $\Theta(\lambda_{\text{opt}})$ is close to $\Theta(0)$ for small $\gamma$ and close to $\Theta(\infty)$ for large $\gamma$.
>   - Under Matérn, the gap to MD is uniformly smaller than under AR(1).
