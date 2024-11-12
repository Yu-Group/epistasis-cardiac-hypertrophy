
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

## Examples

This is a basic example which shows you how to run the lo-siRF permutation test on simulated data:

``` r
library(losiRF)
set.seed(331)

# simulate data
n <- 1000
p <- 100
X <- as.data.frame(
  matrix(sample(0:2, n * p, replace = TRUE), nrow = n, ncol = p)
)
y <- rbinom(n, size = 1, prob = 0.5)

# split data
n_train <- 750
train_idx <- sample(1:n, n_train, replace = FALSE)
X_train <- X[train_idx, , drop = FALSE]
y_train <- y[train_idx]
X_test <- X[-train_idx, , drop = FALSE]
y_test <- y[-train_idx]

# fit signed iterative random forest (siRF)
irf_fit <- iRF::iRF(
  x = X_train,
  y = as.factor(y_train),
  n.iter = 3,
  int.return = 3,
  type = "ranger",
  n.bootstrap = 50,
  n.core = 1
)

# interaction candidates from siRF output
irf_interactions <- tibble::as_tibble(irf_fit$interaction)
irf_interactions

# compute local stability importances and lo-siRF p-values
results <- local_rf_stability_importance(
  rf_fit = irf_fit$rf.list,
  X = X_test,
  y = y_test,
  ints = irf_interactions$int,
)

# lo-siRF feature importance results
results$feature_stability_pvals

# lo-siRF interaction importance results
results$int_stability_pvals
```

Note that in [Wang et al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2), the lo-siRF importances and permutation tests involved additional steps and preprocessing that are not included in this basic example. These additional steps include data preprocessing and dimension reduction prior to fitting the iterative random forest, a thorough prediction check to ensure that the fitted forest is indeed an appropriate fit for the given data, the grouping of SNV *features* into gene *groups* (which can be specified using the `feature_groups` argument in `local_rf_stability_importance()`), and careful screening/filtering of the features (i.e., genes) and interactions prior to conducting the permutation test. We refer interested readers to the original paper ([Wang et al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2)) for more details and the [01_lo-siRF/scripts/](./01_lo-siRF/scripts) folder for scripts to reproduce our original lo-siRF analysis.

We also provide an example which includes screening/filtering of features and interactions prior to conducting the permutation test below:

``` r
library(losiRF)
set.seed(331)

# simulate data
n <- 1000
p <- 100
X <- as.data.frame(
  matrix(sample(0:2, n * p, replace = TRUE), nrow = n, ncol = p)
)
y_prob <- 0.4 * ((X[, 1] > 0) & (X[, 2] > 0)) + 0.4 * ((X[, 3] > 0) & (X[, 4] > 0))
y <- rbinom(n, size = 1, prob = y_prob)

# split data
n_train1 <- 600
n_train2 <- 150
train_idx1 <- sample(1:n, n_train1, replace = FALSE)
train_idx2 <- sample(setdiff(1:n, train_idx1), n_train2, replace = FALSE)
X_train1 <- X[train_idx1, , drop = FALSE]
y_train1 <- y[train_idx1]
X_train2 <- X[train_idx2, , drop = FALSE]
y_train2 <- y[train_idx2]
X_test <- X[-c(train_idx1, train_idx2), , drop = FALSE]
y_test <- y[-c(train_idx1, train_idx2)]

# fit signed iterative random forest (siRF)
irf_fit <- iRF::iRF(
  x = X_train1,
  y = as.factor(y_train1),
  n.iter = 3,
  int.return = 3,
  type = "ranger",
  n.bootstrap = 50,
  n.core = 1
)

# interaction candidates from siRF output
irf_interactions <- tibble::as_tibble(irf_fit$interaction)
irf_interactions

# get features to test (can modify the choice of k depending on the data problem)
top_features <- get_top_stable_features(
  rf_fit = irf_fit$rf.list, X = X_train2, y = y_train2, min_stability = 0.1
)

# get interactions to test (can modify these thresholds, depending on data problem)
keep_ints <- irf_interactions |>
  dplyr::filter(
    cpe > 0, fsd > 0, mip > 0,
    stability > 0.5, sta.cpe > 0.5, sta.fsd > 0.5, sta.mip > 0.5
  )

# compute local stability importances and lo-siRF p-values
results <- local_rf_stability_importance(
  rf_fit = irf_fit$rf.list,
  X = X_test,
  y = y_test,
  features = top_features,
  ints = keep_ints$int,
)

# lo-siRF feature importance results
results$feature_stability_pvals

# lo-siRF interaction importance results
results$int_stability_pvals
```



