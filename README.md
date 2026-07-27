# BayesMAP

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue.svg)]()
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)]()

## Overview

**BayesMAP** is an R package for **Bayesian Multi-layer Association Mapping**, designed to jointly identify trait-associated **SNPs, genes, and cell/pathway annotations** within a unified hierarchical Bayesian framework.

Unlike conventional GWAS methods that analyse variants independently, BayesMAP integrates multiple layers of biological information to improve the identification of causal variants and biologically relevant functional units.

The package implements an efficient **Markov chain Monte Carlo (MCMC)** algorithm, with computationally intensive updates accelerated using **RcppArmadillo**, providing substantial improvements in computational efficiency while maintaining the flexibility and usability of R.

---

## Methodological Framework

BayesMAP models three hierarchical layers simultaneously:

```
SNP  →  Gene  →  Cell / Pathway
```

The model estimates

- SNP effects (β)
- SNP inclusion indicators (δ)
- Gene inclusion indicators (Δ)
- Cell/pathway inclusion indicators (γ)

while borrowing information across biological layers through hierarchical priors.

Posterior inference is performed using Gibbs sampling.

---

## Features

- Bayesian multi-layer association mapping
- Joint modelling of SNPs, genes and cell/pathway annotations
- Hierarchical Bayesian variable selection
- Posterior inclusion probability (PIP) estimation
- Bayesian estimation of SNP effects
- Efficient Gibbs sampler
- Computationally intensive updates implemented in C++ using Rcpp
- Example datasets included
- Fully documented R package

---

## Installation

Install the latest development version from GitHub

```r
install.packages("remotes")

remotes::install_github("mommy003/BayesMAP")

library(BayesMAP)
```

---

## Quick Start
### Individual-level data
```r
library(BayesMAP)

data(X)
data(y)
data(L)
data(A)
data(B)

fit <- bayesmap(
    X = X,
    y = y,
    L = L,
    A = A,
    B = B,
    niter = 1000,
    verbose = TRUE
)

summary(fit)
```

### Summary statistics with block-wise LD matrices

```r
library(BayesMAP)

data(X)
data(y)
data(L)
data(A)
data(B)

N <- nrow(X)

## Standardise genotypes using denominator N
X_std <- scale(
  X,
  center = TRUE,
  scale = FALSE
)

X_std <- sweep(
  X_std,
  MARGIN = 2,
  STATS = sqrt(colMeans(X_std^2)),
  FUN = "/"
)

## Standardise phenotype using denominator N
y_std <- y - mean(y)
y_std <- y_std / sqrt(mean(y_std^2))

## Construct marginal SNP effects and dense LD matrix
bhat <- as.numeric(
  crossprod(X_std, y_std) / N
)

LD <- crossprod(X_std) / N

## Divide SNPs into four blocks
block_indices <- split(
  seq_len(ncol(LD)),
  ceiling(seq_len(ncol(LD)) / 50)
)

block_indices <- unname(block_indices)

## Extract block-specific LD matrices
LD_blocks <- lapply(
  block_indices,
  function(index) {
    LD[index, index, drop = FALSE]
  }
)

set.seed(123)

fit_block <- bayesmap_summary_blocks(
  bhat = bhat,
  LD_blocks = LD_blocks,
  block_indices = block_indices,
  N = N,
  L = L,
  A = A,
  B = B,
  phenotype_variance = 1,
  niter = 5000,
  burnin = 1000,
  thin = 2,
  verbose = TRUE
)

summary(fit_block)
```


```markdown
### Summary-statistics analyses

Dense and block-wise summary-statistics implementations are currently under
development and are not yet available through the exported package interface.
They will be added in a future release.
```


## Inputs

BayesMAP requires

| Object | Description |
|---------|-------------|
| `X` | Genotype matrix (individuals × SNPs) |
| `y` | Phenotype vector |
| `L` | SNP-to-gene annotation matrix |
| `A` | SNP-to-cell/pathway annotation matrix |
| `B` | Gene-to-cell/pathway annotation matrix |

Optional inputs include

- baseline SNP priors
- pathway priors
- starting values
- MCMC settings

---

## Outputs

The fitted model returns

- Posterior SNP effects
- SNP posterior inclusion probabilities (PIP)
- Gene posterior inclusion probabilities
- Cell/pathway posterior inclusion probabilities
- Estimated SNP heritability
- Genetic and residual variance estimates
- Posterior samples of model parameters

---

## Current Development Status

The current version includes

- Hierarchical Bayesian model
- Gibbs sampler
- Rcpp implementation of computational bottlenecks
- Unit testing
- Simulated example dataset

The following developments are currently underway

- Full C++ implementation
- Parallel MCMC
- Multi-threaded computation (OpenMP)
- Prediction of independent datasets
- Advanced visualisation functions
- Vignettes
- CRAN submission

---

## Software Architecture

```
R
│
├── bayesmap()
├── MCMC Driver
├── Parameter Updates
│
└── Rcpp
      ├── updateBetaDelta.cpp
      ├── updateDelta.cpp
      └── updateGamma.cpp
```

---

## Citation

If you use BayesMAP in your research, please cite

> ####
> BayesMAP: #######.
> 

---

## License

GPL-3

---

## Contact

**Md. Moksedul Momin**
Email:
m.momin@uq.edu.au

**Jian Zeng**
Email:
j.zeng@imb.uq.edu.au

Institute for Molecular Bioscience  
The University of Queensland

GitHub:
https://github.com/mommy003

---

## Acknowledgements

BayesMAP has been developed at the Institute for Molecular Bioscience, The University of Queensland, Australia.
