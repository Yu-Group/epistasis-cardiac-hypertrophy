rm(list = ls())
library(magrittr)
library(simChef)
library(future)

RESULTS_DIR <- here::here(
  file.path("01_lo-siRF", "results", "simulations_efficiency")
)
SAVE <- TRUE
USE_CACHED <- TRUE

options(simChef.plot_theme = "vthemes")

set.seed(331)
# plan(multisession, workers = 5)

#### Create experiment parts ####

efficiency_dgp_fun <- function(n, beta0, beta1, beta2, beta12,
                               eff1 = 1, eff2 = 1, int_eff = 1, 
                               err_fun = rnorm, ...) {
  eff1_vec <- rbinom(n, 1, eff1)
  eff2_vec <- rbinom(n, 1, eff2)
  int_eff_vec <- rbinom(n, 1, int_eff)
  
  y0 <- beta0 + err_fun(n, ...)
  y1 <- beta0 + beta1 * eff1_vec + err_fun(n, ...)
  y2 <- beta0 + beta2 * eff2_vec + err_fun(n, ...)
  y12 <- beta0 + (beta1 + beta2 + beta12) * int_eff_vec + err_fun(n, ...)
  
  x <- matrix(0, nrow = 4 * n, ncol = 2)
  x[(n + 1):(2 * n), 1] <- 1
  x[(2 * n + 1):(3 * n), 2] <- 1
  x[(3 * n + 1):(4 * n), 1:2] <- 1
  y <- c(y0, y1, y2, y12)
  
  out <- list(
    x = x,
    y = y,
    eff1_vec = eff1_vec,
    eff2_vec = eff2_vec,
    int_eff_vec = int_eff_vec
  )
}

quant_reg_fun <- function(x, y, tau = 0.5, min_thr = 1e-16, ...) {
  fit_df <- data.frame(.y = y, x)
  fit <- quantreg::rq(formula = .y ~ X1 * X2, tau = tau, data = fit_df)
  fit_summary <- summary(fit, se = "boot")$coefficients %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::rename(Coefficient = Value) %>%
    dplyr::filter(Term == "X1:X2") %>%
    dplyr::mutate(`Pr(>|t|)` = max(min_thr, `Pr(>|t|)`))
  return(fit_summary)
  # return(broom::tidy(fit_summary))
}

quant_reg <- create_method(
  .method_fun = quant_reg_fun, .name = "Quantile Regression (tau = 0.5)"
)

nested_feature_cols <- NULL
feature_col <- "Term"
pval_col <- "Pr(>|t|)"
coef_col <- "Coefficient"

fi_pval <- create_evaluator(
  .eval_fun = summarize_feature_importance,
  .name = 'P-values Summary',
  eval_id = 'pval',
  nested_cols = nested_feature_cols,
  feature_col = feature_col,
  imp_col = pval_col
)

plot_results <- function(fit_results = NULL, eval_results, vary_params = NULL,
                         eval_name, eval_id, add_ggplot_layers = NULL) {
  vary_params <- unique(vary_params)
  plt_df <- eval_results[[eval_name]] %>%
    dplyr::filter(!is.na(!!rlang::sym(vary_params)))
  plt <- plt_df %>%
    tidyr::unnest(cols = tidyselect::all_of(sprintf("raw_%s", eval_id))) %>%
    vdocs::plot_boxplot(
      x_str = vary_params, y_str = sprintf("raw_%s", eval_id)
    ) +
    ggplot2::labs(y = "P-value")
  if (vary_params == "sd") {
    plt <- plt + ggplot2::labs(x = bquote("SD of Noise" ~ (sigma)))
  } else if (vary_params == "beta12") {
    plt <- plt + ggplot2::labs(x = bquote("Interaction Effect" ~ (beta[12])))
  }
  n_dgps <- length(unique(plt_df$.dgp_name))
  if (n_dgps > 1) {
    plt <- plt + ggplot2::facet_wrap(~ .dgp_name)
  }
  if (!is.null(add_ggplot_layers)) {
    for (ggplot_layer in add_ggplot_layers) {
      plt <- plt + ggplot_layer
    }
  }
  return(plt)
}

fi_pval_plot <- create_visualizer(
  .viz_fun = plot_results,
  .name = 'P-values Plot',
  eval_name = 'P-values Summary',
  eval_id = 'pval',
  add_ggplot_layers = list(
    ggplot2::scale_y_continuous(trans = "log10"),
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "red"),
    ggplot2::coord_cartesian(ylim = c(NA, 1))
  )
)


#### Experiment configurations ####
n <- 500
beta0 <- 0
beta1 <- -1
beta2 <- -1
beta12 <- 0
err_fun <- rnorm
noise_sd <- 1

experiment_config_ls <- list(
  `Efficiency = 1` = list(
    dgp_args = list(
      n = n, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta12 = beta12,
      eff1 = 1, eff2 = 1, int_eff = 1,
      err_fun = err_fun, sd = noise_sd
    )
  ),
  `CCDC141-IGF1R, Healthy` = list(
    dgp_args = list(
      n = n, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta12 = beta12,
      eff1 = 0.66, eff2 = 0.9, int_eff = 0.725,
      err_fun = err_fun, sd = noise_sd
    )
  ),
  `CCDC141-TTN, Healthy` = list(
    dgp_args = list(
      n = n, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta12 = beta12,
      eff1 = 0.66, eff2 = 0.8, int_eff = 0.93,
      err_fun = err_fun, sd = noise_sd
    )
  ),
  `CCDC141-IGF1R, Diseased` = list(
    dgp_args = list(
      n = n, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta12 = beta12,
      eff1 = 0.92, eff2 = 0.83, int_eff = 0.695,
      err_fun = err_fun, sd = noise_sd
    )
  ),
  `CCDC141-TTN, Diseased` = list(
    dgp_args = list(
      n = n, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta12 = beta12,
      eff1 = 0.92, eff2 = 0.99, int_eff = 0.92,
      err_fun = err_fun, sd = noise_sd
    )
  )
)

#### Run experiments ####

n_reps <- 100
beta12s <- c(0, -0.1, -0.2, -0.5, -1)
noise_sds <- c(0.1, 0.2, 0.5, 1, 1.5, 2)

for (exp_name in names(experiment_config_ls)) {
  exp_args <- experiment_config_ls[[exp_name]]
  
  dgp <- do.call(
    create_dgp,
    args = c(
      list(
        .dgp_fun = efficiency_dgp_fun,
        .name = exp_name
      ),
      exp_args$dgp_args
    )
  )
  
  experiment_name <- sprintf("%s Simulation", exp_name)
  experiment <- create_experiment(
    name = experiment_name,
    save_dir = file.path(RESULTS_DIR, experiment_name)
  ) %>%
    add_dgp(dgp) %>%
    add_method(quant_reg) %>%
    add_evaluator(fi_pval) %>%
    add_visualizer(fi_pval_plot)
  
  # vary across beta12 signal
  results <- experiment %>%
    add_vary_across(.dgp = dgp$name, beta12 = beta12s) %>%
    run_experiment(
      n_reps = n_reps, save = SAVE, use_cached = USE_CACHED
    )
  
  # vary across noise levels
  results <- experiment %>%
    remove_vary_across() %>%
    add_vary_across(.dgp = dgp$name, sd = noise_sds) %>%
    run_experiment(
      n_reps = n_reps, save = SAVE, use_cached = USE_CACHED
    )
}

render_docs(
  save_dir = RESULTS_DIR, show_eval = FALSE,
  title = "Examining the effect of varying gene-silencing efficiencies on epistasis testing: A simulation study"
)
