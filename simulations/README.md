# Simulations

These scripts reproduce the synthetic experiments in Section 5 of *Massive Samples, Misspecified Models: Robust and Tractable Spike and Slab Variable Selection* (Fontana and Sansó). They use the `fass` package in the sibling folder.

## Layout

```
simulations/
├── R/
│   ├── dgp_generators.R     # Data-generating processes (M-closed, superinduction, post-hoc designs)
│   └── sim_helpers.R        # Metrics: TPR/TNR and model size; RMSE, coverage and sigma^2 error
├── 01_superinduction.R      # Section 5.1: misspecified model, FASS vs Gibbs / PC-Gibbs / skinny Gibbs
├── 02_m_closed.R            # Section 5: M-closed recovery rates of FASS
├── 03_posthoc_inference.R   # Post-hoc inference: EP, Block-VB, MFVB-CAVI, MMFVB, LA
├── 04_plot_results.R        # Figures from the saved results
├── results/                 # .rds data frames written by 01-03
├── figures/                 # .png figures written by 04
└── legacy/                  # Original one-off scripts, kept for reference (not maintained)
```

## Running

1. Install the package. From the folder that contains both `fass/` and `simulations/`:
   ```r
   install.packages("fass", repos = NULL, type = "source")
   ```
2. **Run every script from that same parent folder**, because the scripts use paths such as `simulations/R/...`. Each script accepts `test` (the default, a few minutes in total) or `full` (the paper grid):
   ```bash
   Rscript simulations/01_superinduction.R          # quick check
   Rscript simulations/01_superinduction.R full     # paper grid
   Rscript simulations/02_m_closed.R full
   Rscript simulations/03_posthoc_inference.R full
   Rscript simulations/04_plot_results.R
   ```
   From an interactive R session, use `setwd()` to move to the parent folder, then `source("simulations/01_superinduction.R")`. This runs in test mode.

Results are saved to `simulations/results/` and figures to `simulations/figures/`. Each run overwrites the results of the previous run.

## Experiments

All prior settings follow the paper: a Beta(1, 1) prior on the inclusion probability, spike inverse scale τ₀ = √n and slab τ₁ = 1. In package arguments that is `tau0 = n, tau1 = 1` for the Gaussian and Student-t families and `tau0 = sqrt(n), tau1 = 1` for the Laplace family in `FASS.sampler()`, and `tau0_2 = 1/n, tau1_2 = 1` in the Gibbs samplers. Student-t priors use ν ∈ {2, 5, 10}. FASS selects the highest-odds model it visits, and the Gibbs samplers select the median probability model.

### 1. M-closed recovery (`02_m_closed.R`, Figure 1)

The data come from y = α + Kβ + ε, where K holds the true covariates. Pure-noise covariates J ~ N(0, 10) are added, so the true model is among the candidates. The full grid uses (p, p\*) ∈ {(200, 100), (500, 10)} and n from 5,000 to 10⁶, with 10,000 FASS iterations. The script records TPR, TNR, model size and runtime for FASS under each prior.

### 2. Superinduction (`01_superinduction.R`, Figures 2–3)

The same design gets an extra latent effect μ\* ~ N(0, 10) in the mean that no candidate covariate captures. Each nuisance covariate is Jₖ = (1 − η)Z + ημ\* with η = 0.01, so it is very weakly correlated with the missing effect. The full grid uses p = 500, p\* = 10 and n from 5,000 to 10⁶. FASS runs with every prior and is compared against the standard Gibbs, partially collapsed Gibbs and skinny Gibbs samplers under the same priors, each with 2,000 burn-in and 8,000 retained iterations.

The Gibbs baselines dominate the runtime of the full grid and take days at n = 10⁶. Set `include_mcmc_baselines <- FALSE` near the top of the script to run only the FASS part.

### 3. Post-hoc inference (`03_posthoc_inference.R`)

For a fixed model, the script compares EP, Block-VB, MFVB-CAVI, MMFVB and the Laplace approximation under Laplace and Student-t (ν ∈ {2, 5, 10}) slab priors. It uses uncorrelated designs and AR(1) designs with ρ = 0.8, p ∈ {10, 100} and n up to 10⁵. The metrics are coefficient RMSE and MAE, 95% interval coverage and length, and the absolute error of σ².

### Figures (`04_plot_results.R`)

| File | Content |
|---|---|
| `m_closed_accuracy.png` | TPR/TNR of FASS by prior in the M-closed setting (Figure 1) |
| `m_closed_runtime.png` | FASS runtime in the M-closed setting |
| `superinduction_accuracy.png` | TPR/TNR by algorithm and prior under superinduction (Figures 2–3) |
| `superinduction_runtime.png` | Runtime of FASS and the Gibbs samplers (Figure 2) |
| `posthoc_rmse.png` | Coefficient RMSE of the post-hoc approximations |

## Notes

- Test mode uses small n and p and short chains. Its numbers only check that the pipeline runs end to end and do not reproduce the paper's figures.
- A fixed seed (4747) is used for data generation, so each script's datasets are reproducible.
- `legacy/` contains the original scripts behind the paper's first draft. They source files that are no longer in this repository and are kept only as a record.
