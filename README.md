# BayesMAP

**BayesMAP** is an R package for Bayesian multi-layer association mapping that integrates SNPs, genes, and cell/pathway annotations within a hierarchical Bayesian framework.

The package implements an efficient Markov chain Monte Carlo (MCMC) algorithm with computationally intensive components accelerated using **Rcpp/C++**, providing improved performance while maintaining the flexibility of R.

## Features

* Bayesian variable selection for SNPs, genes, and cell/pathway annotations
* Hierarchical probabilistic modeling
* Posterior inclusion probability (PIP) estimation
* C++-accelerated MCMC updates using Rcpp
* Fully documented R package with automated testing
* Example datasets included

## Installation

### Install the development version

```r
# install.packages("remotes")
remotes::install_github("mommy003/BayesMAP")
```

### Load the package

```r
library(BayesMAP)
```

## Example

```r
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

## Output

The fitted object contains posterior samples and summary statistics, including:

* Posterior samples of model parameters
* SNP posterior inclusion probabilities
* Gene posterior inclusion probabilities
* Cell/pathway posterior inclusion probabilities
* Estimated heritability

## Development Status

BayesMAP is under active development. Planned improvements include:

* Additional C++ optimization
* Parallel MCMC implementation
* Enhanced visualization functions
* Prediction utilities
* Comprehensive vignettes

## Citation

If you use BayesMAP in your research, please cite the associated methodological publication (to be added).

## License

GPL-3 (or the license you selected for the package).
