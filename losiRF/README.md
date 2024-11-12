
# losiRF

<!-- badges: start -->
<!-- badges: end -->

The losiRF R package provides tools to compute the local stability importances and p-values for features and interactions from a fitted random forest, as described in "Epistasis regulates genetic control of cardiac hypertrophy" by [Wang et al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2).

## Installation

You can install the development version of losiRF from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Yu-Group/epistasis-cardiac-hypertrophy", subdir = "losiRF")
```

## Example

This is a basic example which shows you how to run the lo-siRF permutation test on simulated data:

``` r
library(losiRF)

# simulate data
n <- 1500
n_train <- 1000
p <- 10
X <- as.data.frame(matrix(rnorm(n * p), nrow = n, ncol = p))
y_prob <- 0.8 * ((X[, 1] > 0) & (X[, 2] > 0))
y <- rbinom(n, size = 1, prob = y_prob)

train_idx <- sample(1:n, n_train, replace = FALSE)
X_train <- X[train_idx, , drop = FALSE]
y_train <- y[train_idx]
X_test <- X[-train_idx, , drop = FALSE]
y_test <- y[-train_idx]

# # fit (iterative) random forest
# rf_fit <- iRF::iRF(
#   x = X_train,
#   y = y_train,
#   n.iter = 3,
#   int.return = 3,
#   type = "ranger",
#   n.bootstrap = 50,
#   ncores = 1
# )
rf_fit <- ranger::ranger(
  formula = y ~ .,
  data = data.frame(y = as.factor(y_train), X_train)
)

results <- local_rf_stability_importance(
  rf_fit = rf_fit,
  X = X_test,
  y = y_test,
  features = colnames(X_test),
  ints = c("X1_X2", "X3_X4")
)

results$feature_stability_pvals
results$int_stability_pvals
```

Note that in [Wang et al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2), the lo-siRF importances and permutation tests involved additional steps and preprocessing that are not included in this basic example. These additional steps include data preprocessing and dimension reduction prior to fitting the iterative random forest, a thorough prediction check to ensure that the fitted forest is indeed an appropriate fit for the given data, the grouping of SNV *features* into gene *groups* (which can be specified using the `feature_groups` argument in `local_rf_stability_importance()`), and rigorous screening/filtering of genes and interactions prior to conducting the permutation test. We refer interested readers to the original paper ([Wang et al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2)) for more details and the [01_lo-siRF/scripts/](./01_lo-siRF/scripts) folder for scripts to reproduce our original lo-siRF analysis.

