# Nonlinear Correlation

This subproject accompanies the article ["An Undeservedly Forgotten Correlation Coefficient"](https://medium.com/data-science/an-undeservedly-forgotten-correlation-coefficient-86245ccb774c), also mirrored on [Towards Data Science](https://towardsdatascience.com/an-undeservedly-forgotten-correlation-coefficient-86245ccb774c/).

It contains a compact R script that compares three association measures:

- Pearson correlation `rho`, via `cor()`
- a mutual-information-based coefficient `R`, implemented here as `sqrt(1 - exp(-2 * max(0, mutinfo(x, y, k))))`
- Chatterjee's rank correlation `xi`, via `xicor()`

## What is in this folder

- `mic.R`: the full experiment script, including synthetic data generation, coefficient computation, benchmarking, plotting, and PNG export.

## Main idea

The article argues that a useful nonlinear correlation measure should satisfy several practical properties:

- nonlinearity: detect dependence beyond linear or monotonic structure
- symmetry: give the same value for `A(x, y)` and `A(y, x)`
- consistency: agree with Pearson correlation under a bivariate normal distribution
- scalability: remain computable on moderately large samples
- precision: have low estimator variance

The script evaluates those ideas empirically by comparing `rho`, `R`, and `xi` on several synthetic datasets.

## What the script implements

`mic.R` runs four experiment blocks.

1. Nonlinearity:
   generates donut-shaped data with varying thickness and shows that `R` captures dependence even when Pearson correlation is near zero.
2. Symmetry:
   generates noisy parabolic data and compares `xi(x, y)` with `xi(y, x)` to highlight that `xi` is directional, while `rho` and `R` are symmetric.
3. Consistency:
   samples bivariate normal data with target correlations `0.4`, `0.7`, and `1.0`, then checks whether each coefficient tracks the corresponding linear correlation.
4. Scalability:
   benchmarks runtime on independent uniform data for sample sizes `100`, `1000`, `10000`, and `50000`.

Each block saves a figure:

- `nonlinearity_plot.png`
- `symmetry_plot.png`
- `consistency_plot.png`
- `scalability_plot.png`

## How to run it

From this directory:

```r
install.packages(c("XICOR", "FNN", "MASS", "ggplot2", "gridExtra", "microbenchmark"))
source("mic.R")
```

Or from the shell:

```bash
Rscript mic.R
```

The script writes the four PNG files listed above into the current working directory.

## Notes

- Despite the helper function name `mic()`, this code is implementing the mutual-information-based coefficient `R` described in the article, not the maximal information coefficient usually abbreviated as MIC.
- Exact coefficient values and benchmark timings vary slightly with platform, package versions, and random sampling, though the script fixes seeds for the main experiment blocks.
