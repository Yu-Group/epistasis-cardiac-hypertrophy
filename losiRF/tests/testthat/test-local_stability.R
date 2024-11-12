test_that("local_rf_stability_importance works", {

  # simulate data
  n <- 500
  p <- 10
  X <- as.data.frame(
    matrix(sample(0:2, n * p, replace = TRUE), nrow = n, ncol = p)
  )
  y <- rbinom(n, size = 1, prob = 0.5)
  data <- cbind(y, X)

  # fit random forest
  rf_fit <- ranger::ranger(
    formula = y ~ .,
    data = data,
    num.trees = 10,
    probability = TRUE
  )

  expect_no_error(
    local_rf_stability_importance(
      rf_fit = rf_fit,
      X = X,
      y = y,
      ints = c("V1+_V2+", "V3+_V4-")
    )
  )
})
