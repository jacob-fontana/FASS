# fass: Fast Approximate Stochastic Search

`fass` implements the methods of *Massive Samples, Misspecified Models: Robust and Tractable Spike and Slab Variable Selection* (Fontana and Sansó). It performs Bayesian variable selection in linear regression when the sample size is very large and none of the candidate models is exactly correct (the M-open setting).

Two problems appear in that regime, and the package addresses both:

- **Model superinduction.** When the model is misspecified, the posterior odds between nested models favour the larger model at a rate that grows exponentially in *n*. Searches drift toward models with hundreds of tiny, near-spurious effects. This happens with point-mass spikes and with Gaussian or Laplace continuous spikes. **Student-t spikes and slabs are robust to it**: the penalty for an extra parameter becomes polynomial instead of exponential.
- **Computational cost.** Gibbs samplers cost O(np + p³) per iteration, and skinny Gibbs costs O(np + |γ|³). Both slow down further once superinduction inflates the model size. **FASS** scores each candidate move with a first-order surrogate for the posterior odds. After a one-time pre-computation of (XᵀX)⁻¹, β̂ and the residual variance, every search iteration costs **O(p), whatever the value of *n***.

## Installation

The package contains C++ code. Windows users need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to install it.

```r
# from the folder that contains fass/
install.packages("fass", repos = NULL, type = "source")
# or
devtools::install("fass")
```

## Quick start

```r
library(fass)

set.seed(1)
n <- 5000; p <- 50
X <- matrix(rnorm(n * p), n, p)
y <- drop(2 + X[, 1:5] %*% c(3, -2, 1.5, -1, 0.5) + rnorm(n))

# Student-t spike and slab (nu = 5), paper defaults tau0 = n, tau1 = 1
fit <- FASS.sampler(y, X, niter = 2000,
                    family = "Student", nu = 5, tau1 = 1, tau0 = n,
                    full_posterior = "Block-VB")

# Highest-odds model visited (odds[1] is the null model; odds[i + 1] <-> gamma[i, ])
best <- which.max(fit$odds)
selected <- if (best == 1) integer(0) else which(fit$gamma[best - 1, ] == 1)
selected

# Post-hoc posterior for the selected coefficients
fit$posterior$mu
fit$posterior$alpha     # intercept
fit$posterior$sigma2    # residual variance

# Out-of-sample prediction and scoring
X_new <- matrix(rnorm(200 * p), 200, p)
y_new <- drop(2 + X_new[, 1:5] %*% c(3, -2, 1.5, -1, 0.5) + rnorm(200))
pred <- predict_approx(fit, X_new, y_new)
pred$RMSE; mean(pred$CRPS); mean(pred$coverage)
```

## Prior parameterisation

All families are continuous spike-and-slab mixtures, π(βⱼ | γⱼ) = γⱼ π₁(βⱼ) + (1 − γⱼ) π₀(βⱼ), with γⱼ ~ Bernoulli(θ) and θ ~ Beta(`a_pi`, `b_pi`). The coefficient priors are not scaled by σ². In `FASS.sampler` and `fit_approx_posterior`, `tau0` (spike) and `tau1` (slab) are inverse-scale parameters:

| `family` | Prior on βⱼ | `tau` is | Paper setting (τ₀ = √n, τ₁ = 1) |
|---|---|---|---|
| `"Gaussian"` | N(0, 1/`tau`) | precision (τ² in the paper) | `tau0 = n`, `tau1 = 1` |
| `"Student"` | t_ν with scale 1/√`tau` | squared inverse scale (τ²) | `tau0 = n`, `tau1 = 1`, `nu` ≥ 5 |
| `"Laplace"` | (`tau`/2) exp(−`tau`·\|βⱼ\|) | rate (τ) | `tau0 = sqrt(n)`, `tau1 = 1` |

The benchmark Gibbs samplers (`Gibbs`, `Gibbs_PC`, `skinnyGibbs`) take squared **scales** instead: `tau0_2 = 1/n, tau1_2 = 1` gives the same priors as the paper settings above.

Small degrees of freedom (ν = 2) flatten the difference between signal and noise and give poor recovery; ν ≥ 5 is stable in our experiments.

## Main functions

| Function | Purpose |
|---|---|
| `FASS.sampler()` | FASS model search with Gaussian, Laplace or Student-t spike-and-slab surrogates, with optional post-hoc inference on the best model |
| `fit_approx_posterior()` | Posterior approximation for a given model: `"Block-VB"` (recommended), `"MFVB-CAVI"`, `"MMFVB"`, `"LA"` (Laplace approximation), `"EP"` (expectation propagation). σ² is estimated by EM |
| `VB()`, `VB_MF()`, `VB_Direct()`, `Laplace_Approx()`, `EP()` | The individual approximation methods (they expect centred `X` when `intercept = TRUE`; `fit_approx_posterior()` handles the centring) |
| `Gibbs()`, `Gibbs_PC()`, `skinnyGibbs()` | Benchmark samplers: standard (uncollapsed), partially collapsed, and skinny Gibbs |
| `predict_approx()` / `predict_FASS()` | Gaussian posterior predictive, CRPS, RMSE and interval coverage for FASS or approximate-posterior fits |
| `predict_gibbs()` | Posterior predictive from Gibbs draws using the median probability model |
| `compute_gaussian_crps()`, `compute_empirical_crps()`, `compute_gibbs_metrics_cpp()` | CRPS utilities |

### How selection works

`FASS.sampler()` begins at the null model. At each iteration it evaluates the approximate posterior odds of all *p* single-variable add/remove moves in O(p), then samples a move (or stays put) in proportion to those odds. `fit$odds` records the approximate log posterior odds of each visited model against the null model, and the model with the highest odds is selected. The Gibbs samplers select the median probability model, `colMeans(fit$outgamma) > 0.5`.

When the first-order surrogate is non-positive (possible for the Gaussian and Student-t families, though rare for large *n*), that term is replaced by its zeroth-order version, which is the prior evaluated at β̂. The paper shows that strong consistency still holds with this fallback.

### Requirements on the data

- `n > p`, and `X` must have full column rank. FASS inverts the full Gram matrix once.
- Do not include a column of ones in `X`. Use `intercept = TRUE`, which is the default. The intercept has a flat prior, and the returned intercepts refer to the original, uncentred covariates.
- Inputs are validated up front, with clear errors for missing or non-finite values, dimension mismatches, collinear designs and invalid hyperparameters.

## Reproducing the paper's simulations

The simulation scripts are in the sibling `simulations/` folder. See `simulations/README.md`.

## Citation

Fontana, J. and Sansó, B. *Massive Samples, Misspecified Models: Robust and Tractable Spike and Slab Variable Selection.*
