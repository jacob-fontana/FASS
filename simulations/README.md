# Simulations

These scripts reproduce the synthetic experiments in Section 5 of *Massive Samples, Misspecified Models: Robust and Tractable Spike and Slab Variable Selection* (Fontana and Sansó). The scripts require use of the `fass` package.

## Layout

```
simulations/
├── R/
│   ├── dgp_generators.R     # Data-generating processes (M-closed, superinduction, post-hoc designs)
│   └── sim_helpers.R        # Metrics: TPR/TNR and model size; RMSE, coverage and sigma^2 error
├── 01_superinduction.R      # Section 5.1: Experiment with M-open misspecified model --- compares FASS vs Gibbs, partially collapsed Gibbs and skinny Gibbs
├── 02_m_closed.R            # Section 5: Experiment with M-closed DGP to benchmark FASS performance against ground truth (recovery rates, runtimes, etc.)
├── 03_posthoc_inference.R   # Experiment to compare several post-hoc inference methods with FASS, including: Expectation propagation, Block-VB, Mean-field VB-CAVI, MMFVB, and Laplace approximations
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
2. **Run every script from that same parent folder**, because the scripts use paths such as `simulations/R/...`. Each script accepts `test` (the default, small sample runs that take a few minutes in total) or `full` (the full sample sequence runs done in the paper, can take several hours):
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

All prior settings follow the paper: a Beta(1, 1) prior on the inclusion probability, spike inverse scale τ₀ = √n and slab τ₁ = 1. In package arguments that is `tau0 = n, tau1 = 1` for the Gaussian and Student-t families and `tau0 = sqrt(n), tau1 = 1` for the Laplace family in `FASS.sampler()`, and `tau0_2 = 1/n, tau1_2 = 1` in the Gibbs samplers. For Student-t priors, the degrees of freedom considered are ν ∈ {2, 5, 10}. FASS selects the maximum a-posteriori model that the search has visited. The Gibbs samplers select the median probability model.

### 1. M-closed recovery (`02_m_closed.R`, Figure 1)

The data are generated from y = α + Kβ + ε, where K holds the p\* true covariates. Additionally, p nuissance covariates are added and are generated as J ~ N(0, 10). The experiment considered two settings: (p, p\*) ∈ {(200, 100), (500, 10)} and a sequence of n values ranging from 5,000 to 10⁶. For each setting and sample size, FASS is iterated 10000 times. The script records the TPR, TNR, selected model size and runtime for FASS under each prior (Gaussian, laplace, Student).

### 2. Superinduction (`01_superinduction.R`, Figures 2–3)

The data are now generated from y = α + Kβ + μ\* + ε, where the missing effect μ\* is not captured by any covariate (μ\* ~ N(0, 10)). Each nuisance covariate is now generated as Jₖ = (1 − η)Z + ημ\* with η = 0.01, so each is very weakly correlated with the missing effect, but is mostly Gaussian noise. The experiment considers the setting p = 500, p\* = 10 and a sequence of n from 5,000 to 10⁶. FASS is again run under each prior (Gaussian, Laplace, Student) and its performance compared against the standard Gibbs sampler, a partially collapsed Gibbs sampler, and a skinny Gibbs sampler (all under the same priors). FASS is again run for 10000 iterations, and each Gibbs sampler is run for 10000 iterations, with the first 2,000 discarded for burn-in.

The Gibbs samplers are very slow at large sample sizes, and can take days to run to the full sequence depending on the machine. Set `include_mcmc_baselines <- FALSE` near the top of the script to run only the FASS part.

### 3. Post-hoc inference (`03_posthoc_inference.R`)

This experiment assesses the post-hoc inferences methods available for FASS. The script compares EP, Block-VB, MFVB-CAVI, MMFVB and the Laplace approximation under Laplace and Student-t (ν ∈ {2, 5, 10}) slab priors (Gaussian priors don't require approximate inference on the posterior of β). Two experiments are performed, considering both an uncorrelated design and an AR(1) design with ρ = 0.8. Each are assessed under settings: p ∈ {10, 100} and n sequences ranging from 10^3 to 10⁵. The metrics are coefficient-level RMSE and MAE, 95% interval coverage and length, and the absolute error of σ².

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
