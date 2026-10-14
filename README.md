# FASS

Code for *Massive Samples, Misspecified Models: Robust and Tractable Spike and Slab Variable Selection* (Fontana and Sansó).

- [`fass/`](fass): the R package, which implements Fast Approximate Stochastic Search (FASS), post-hoc posterior approximations and the benchmark Gibbs samplers. See [`fass/README.md`](fass/README.md).
- [`simulations/`](simulations): scripts that reproduce the simulation studies in the paper. See [`simulations/README.md`](simulations/README.md).

## Installation

```r
# install.packages("remotes")
remotes::install_github("jacob-fontana/FASS", subdir = "fass")
```

The package contains C++ code (Rcpp). Windows users need [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
