
<!-- README.md is generated from README.Rmd. Please edit that file -->

# losiRF

<!-- badges: start -->
<!-- badges: end -->

The losiRF R package provides tools to compute the local stability
importances and p-values for features and interactions from a fitted
random forest, as described in “Epistasis regulates genetic control of
cardiac hypertrophy” by [Wang et
al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2).

## Installation

You can install the development version of losiRF from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Yu-Group/epistasis-cardiac-hypertrophy", subdir = "losiRF")
```

## Examples

This is a basic example which shows you how to run the lo-siRF
permutation test on simulated data:

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
#> Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display
#> Warning: 'rgl.init' failed, running with 'rgl.useNULL = TRUE'.
#> [1] "iteration =  1"
#> [1] "iteration =  2"
#> [1] "iteration =  3"
#> finding interactions...
#> evaluating interactions...

# interaction candidates from siRF output
irf_interactions <- tibble::as_tibble(irf_fit$interaction)
irf_interactions
#> # A tibble: 17 × 10
#>    int   prevalence precision      cpe sta.cpe      fsd sta.fsd      mip sta.mip
#>    <chr>      <dbl>     <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>
#>  1 V86+…   0.00245      0.540  0.264      0.96  0.190      0.82  0.0567     0.76
#>  2 V9+_…   0.00316      0.525  0.202      0.98  0.0932     0.74  0.0285     0.74
#>  3 V28+…   0.00185      0.506  0.128      0.86 -0.210      0.32  0.00774    0.54
#>  4 V16-…   0.00144      0.501  0.106      0.74  0.148      0.64  0.0168     0.5 
#>  5 V100…   0.00132      0.501  0.105      0.74 -0.169      0.42  0.00942    0.6 
#>  6 V81+…   0.00187      0.499  0.0999     0.7  -0.163      0.3  -0.0224     0.46
#>  7 V17+…   0.00326      0.491  0.0675     0.76 -0.107      0.42 -0.0361     0.18
#>  8 V37-…   0.00168      0.489  0.0604     0.68 -0.00179    0.5  -0.0130     0.52
#>  9 V21+…   0.000593     0.489  0.0590     0.6  -0.0742     0.5   0.00745    0.54
#> 10 V14+…   0.000654     0.478  0.0143     0.5  -0.215      0.32 -0.0164     0.46
#> 11 V50-…   0.00103      0.476  0.00334    0.5  -0.0606     0.42 -0.00990    0.5 
#> 12 V54-…   0.000965     0.468 -0.0264     0.48 -0.0530     0.5  -0.0222     0.5 
#> 13 V1-_…   0.00127      0.461 -0.0555     0.42 -0.202      0.3  -0.0501     0.34
#> 14 V28-…   0.00181      0.458 -0.0666     0.32 -0.0820     0.38 -0.0601     0.24
#> 15 V12-…   0.00126      0.457 -0.0722     0.4  -0.107      0.42 -0.0378     0.34
#> 16 V39-…   0.00356      0.387 -0.356      0     0.206      0.86 -0.164      0   
#> 17 V30-…   0.00415      0.368 -0.437      0    -0.186      0.18 -0.163      0   
#> # ℹ 1 more variable: stability <dbl>

# compute local stability importances and lo-siRF p-values
results <- local_rf_stability_importance(
  rf_fit = irf_fit$rf.list,
  X = X_test,
  y = y_test,
  ints = irf_interactions$int,
)

# lo-siRF feature importance results
results$feature_stability_pvals
#> # A tibble: 200 × 4
#>    feature   pval     T_obs perm_dist     
#>    <chr>    <dbl>     <dbl> <list>        
#>  1 V1+     0.553  -0.00289  <dbl [10,000]>
#>  2 V2+     0.644  -0.00286  <dbl [10,000]>
#>  3 V3+     0.0917 -0.00615  <dbl [10,000]>
#>  4 V4+     0.817   0.000786 <dbl [10,000]>
#>  5 V5+     0.800   0.00125  <dbl [10,000]>
#>  6 V6+     0.0479 -0.00577  <dbl [10,000]>
#>  7 V7+     0.0464 -0.0100   <dbl [10,000]>
#>  8 V8+     0.126  -0.00383  <dbl [10,000]>
#>  9 V9+     0.745   0.00316  <dbl [10,000]>
#> 10 V10+    0.350   0.00332  <dbl [10,000]>
#> # ℹ 190 more rows

# lo-siRF interaction importance results
results$int_stability_pvals
#> # A tibble: 17 × 4
#>    int         pval     T_obs perm_dist     
#>    <chr>      <dbl>     <dbl> <list>        
#>  1 V86+_V91+ 0.702   0.000155 <dbl [10,000]>
#>  2 V9+_V92-  0.588   0.000258 <dbl [10,000]>
#>  3 V28+_V53+ 0.865  -0.000216 <dbl [10,000]>
#>  4 V16-_V2-  0.165   0.000486 <dbl [10,000]>
#>  5 V100+_V5- 0.442  -0.000288 <dbl [10,000]>
#>  6 V81+_V9+  0.0366  0.000715 <dbl [10,000]>
#>  7 V17+_V28- 0.552   0.000305 <dbl [10,000]>
#>  8 V37-_V70+ 0.171  -0.00101  <dbl [10,000]>
#>  9 V21+_V78+ 0.330   0.000144 <dbl [10,000]>
#> 10 V14+_V16+ 0.292  -0.000275 <dbl [10,000]>
#> 11 V50-_V65- 0.520  -0.000188 <dbl [10,000]>
#> 12 V54-_V80- 0.339  -0.000307 <dbl [10,000]>
#> 13 V1-_V74-  0.397   0.000398 <dbl [10,000]>
#> 14 V28-_V92- 0.661  -0.000200 <dbl [10,000]>
#> 15 V12-_V91- 0.518  -0.000116 <dbl [10,000]>
#> 16 V39-_V47- 0.284   0.000875 <dbl [10,000]>
#> 17 V30-_V47- 0.273  -0.00116  <dbl [10,000]>
```

Note that in [Wang et
al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2),
the lo-siRF importances and permutation tests involved additional steps
and preprocessing that are not included in this basic example. These
additional steps include data preprocessing and dimension reduction
prior to fitting the iterative random forest, a thorough prediction
check to ensure that the fitted forest is indeed an appropriate fit for
the given data, the grouping of SNV *features* into gene *groups* (which
can be specified using the `feature_groups` argument in
`local_rf_stability_importance()`), and careful screening/filtering of
the features (i.e., genes) and interactions prior to conducting the
permutation test. We refer interested readers to the original paper
([Wang et
al. (2024)](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v2))
for more details and the [01_lo-siRF/scripts/](./01_lo-siRF/scripts)
folder for scripts to reproduce our original lo-siRF analysis.

We also provide an example which includes screening/filtering of
features and interactions prior to conducting the permutation test
below:

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
#> [1] "iteration =  1"
#> [1] "iteration =  2"
#> [1] "iteration =  3"
#> finding interactions...
#> evaluating interactions...

# interaction candidates from siRF output
irf_interactions <- tibble::as_tibble(irf_fit$interaction)
irf_interactions
#> # A tibble: 27 × 10
#>    int         prevalence precision   cpe sta.cpe     fsd sta.fsd    mip sta.mip
#>    <chr>            <dbl>     <dbl> <dbl>   <dbl>   <dbl>   <dbl>  <dbl>   <dbl>
#>  1 V1+_V2+_V3…     0.0788     0.779 1.98        1  0.429     1    0.0568       1
#>  2 V1+_V2+_V4+     0.0806     0.736 1.75        1  0.270     1    0.157        1
#>  3 V1+_V3+_V4+     0.174      0.677 1.46        1  0.173     1    0.0729       1
#>  4 V1+_V2+_V3+     0.134      0.630 1.26        1  0.339     1    0.125        1
#>  5 V1+_V4+         0.177      0.629 1.25        1  0.0118    0.64 0.302        1
#>  6 V2+_V3+_V4+     0.152      0.626 1.24        1 -0.0374    0.12 0.116        1
#>  7 V2+_V46-        0.0165     0.568 0.999       1  0.307     1    0.226        1
#>  8 V3+_V4+         0.515      0.557 0.951       1  0.138     1    0.218        1
#>  9 V1+_V3+         0.281      0.556 0.947       1  0.0325    1    0.178        1
#> 10 V3+_V46-        0.0316     0.552 0.932       1  0.0909    1    0.187        1
#> # ℹ 17 more rows
#> # ℹ 1 more variable: stability <dbl>

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
#> # A tibble: 8 × 4
#>   feature   pval   T_obs perm_dist     
#>   <chr>    <dbl>   <dbl> <list>        
#> 1 V3+     0       0.261  <dbl [10,000]>
#> 2 V2+     0       0.226  <dbl [10,000]>
#> 3 V4+     0       0.267  <dbl [10,000]>
#> 4 V3-     0      -0.195  <dbl [10,000]>
#> 5 V1+     0       0.171  <dbl [10,000]>
#> 6 V1-     0.180  -0.0439 <dbl [10,000]>
#> 7 V2-     0.0001 -0.169  <dbl [10,000]>
#> 8 V4-     0.0003 -0.137  <dbl [10,000]>

# lo-siRF interaction importance results
results$int_stability_pvals
#> # A tibble: 6 × 4
#>   int          pval T_obs perm_dist     
#>   <chr>       <dbl> <dbl> <list>        
#> 1 V1+_V3+_V4+     0 0.159 <dbl [10,000]>
#> 2 V1+_V2+_V3+     0 0.167 <dbl [10,000]>
#> 3 V1+_V4+         0 0.163 <dbl [10,000]>
#> 4 V3+_V4+         0 0.298 <dbl [10,000]>
#> 5 V1+_V3+         0 0.186 <dbl [10,000]>
#> 6 V1+_V2+         0 0.197 <dbl [10,000]>
```
