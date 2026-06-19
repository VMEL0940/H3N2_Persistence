# ==============================================================================
# Script Name: 11_ML_DL_RollingWidth10_Comparison_WithDNN_and_Figures.R
# Purpose:
#   Compare ML/DL models against the final mixed-effects logistic regression model
#   under the same rolling-width-10 temporal validation splits.
#
# Final regression model to compare against:
#   dom_binary ~ HA * PB1 + PB2 + NA_seg + (1 | vaccine_code)
#
# Design:
#   - Same dataset and period order as regression scripts.
#   - Same rolling width = 10, last CV test period = TH22.
#   - CR23 is excluded from CV and kept for future deployment.
#   - Scaling is fit on each outer training split only, then applied to test split.
#   - Inner temporal CV is used for hyperparameter selection within each training split.
#   - Feature schemes are aligned with regression analysis:
#       1) Selected_dN_final: HA, PB1, PB2, NA_seg, HA_x_PB1
#       2) HA_dN_only
#       3) All_dN
#       4) All_dN_All_dS
#   - Outputs include train/test AUROC, deltas versus selected regression model,
#     publication-ready summary tables, and manuscript-style A/B/C figures including DNN.
# ===============================================================================

library(reticulate)
reticulate::py_require("tensorflow")

library(tensorflow)
library(keras)

tensorflow::tf_config()

# ---- 0A. TensorFlow/Keras backend initialization for DNN ----
# This block must run before any keras/tensorflow model call.
# New reticulate/keras backends may otherwise repeatedly print:
#   Hint: To use tensorflow with py_require(), call py_require("tensorflow")
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "2")

initialize_dnn_backend <- function() {
  if (!requireNamespace("keras", quietly = TRUE) || !requireNamespace("tensorflow", quietly = TRUE)) {
    message("DNN disabled: R packages 'keras' and/or 'tensorflow' are not installed.")
    return(FALSE)
  }

  # For recent reticulate-managed Python environments.
  if (requireNamespace("reticulate", quietly = TRUE)) {
    if ("py_require" %in% getNamespaceExports("reticulate")) {
      try(reticulate::py_require("tensorflow"), silent = TRUE)
    }
  }

  # For legacy-keras installations. Safe to skip if the function is unavailable.
  if ("py_require_legacy_keras" %in% getNamespaceExports("keras")) {
    try(keras::py_require_legacy_keras(), silent = TRUE)
  }

  ok <- tryCatch({
    invisible(tensorflow::tf$constant(1L))
    TRUE
  }, error = function(e) {
    message("DNN disabled: TensorFlow backend could not be initialized: ", e$message)
    FALSE
  })

  ok
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(pROC)
  library(glmnet)
  library(e1071)
  library(randomForest)
})

# ---- 0. User settings ----
base_path <- "~/"
file_name <- "1_Data/0_Metadata/H3N2_777strains_Metadata.csv"
data_path <- file.path(base_path, file_name)

final_window_width <- 10
last_test_period <- "TH22"

# Regression output from Script 10. Used only for side-by-side comparison.
# The first path is the original directory. The second path is accepted if you renamed
# the final regression output directory to explicitly include RollingWidth10.
regression_result_dir_candidates <- c(
  file.path(base_path, "3_Result/8_Final_Rolling_TrainTest_AUROC"),
  file.path(base_path, "3_Result/8_Final_RollingWidth10_TrainTest_AUROC")
)
regression_table_f2_candidates <- file.path(
  regression_result_dir_candidates,
  "Table_F2_Train_Test_AUROC_By_Rolling_Split.csv"
)
regression_table_f2 <- regression_table_f2_candidates[file.exists(regression_table_f2_candidates)][1]
if (is.na(regression_table_f2)) {
  regression_table_f2 <- regression_table_f2_candidates[1]
}
regression_selected_model_name <- "1_Final_Rolling_Best_dN_HA_x_PB1_plus_PB2_NA"

# ML output directory
result_dir <- file.path(base_path, "3_Result/9_ML_DL_RollingWidth10_Comparison")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# Computational settings
set.seed(777)
run_xgboost <- requireNamespace("xgboost", quietly = TRUE)
run_dnn <- initialize_dnn_backend()

# Toggle manually if needed
# run_xgboost <- FALSE
# run_dnn <- FALSE

# DNN settings. Kept intentionally small because outer rolling splits and inner temporal CV are small.
dnn_units_grid <- c(4, 8)
dnn_dropout_grid <- c(0.0, 0.2)
dnn_lr_grid <- c(0.001)
dnn_epochs <- 120
dnn_batch_size <- 32
dnn_verbose <- 0

# ---- 1. Load data ----
if (!file.exists(data_path)) stop(paste("Target dataset metadata not found at:", data_path))
full_df <- read.csv(data_path, stringsAsFactors = FALSE)

base_levels <- c(
  "Mos99", "Fuj02", "Cal04", "Wis05", "Bris07", "Prth09", "Vic11",
  "Swtz13", "HK15", "Kan17", "HK19", "Cam20", "Dar21", "TH22", "CR23"
)
unique_codes <- unique(full_df$vaccine_code)
all_periods <- unique(c(base_levels, unique_codes))

if (!(last_test_period %in% all_periods)) stop(paste("last_test_period not found:", last_test_period))

raw_df <- full_df %>%
  mutate(
    date = as_date(date),
    decimal_year = decimal_date(date),
    Period = factor(vaccine_code, levels = all_periods),
    dom_binary = as.integer(dom_binary)
  )

# ---- 2. Feature preprocessing ----
dn_raw <- c("PB2_Nonsyn", "PB1_Nonsyn", "PA_Nonsyn", "HA_Nonsyn", "NP_Nonsyn", "NA_Nonsyn", "M_Nonsyn", "NS_Nonsyn")
ds_raw <- c("PB2_Syn",    "PB1_Syn",    "PA_Syn",    "HA_Syn",    "NP_Syn",    "NA_Syn",    "M_Syn",    "NS_Syn")
all_raw_features <- c(dn_raw, ds_raw)

fit_scaler <- function(df, cols) {
  stats <- lapply(cols, function(col) {
    x <- df[[col]]
    mu <- mean(x, na.rm = TRUE)
    sig <- sd(x, na.rm = TRUE)
    if (is.na(mu)) mu <- 0
    if (is.na(sig) || sig == 0) sig <- 1
    list(mean = mu, sd = sig)
  })
  names(stats) <- cols
  stats
}

apply_scaler <- function(df, scaler) {
  for (col in names(scaler)) {
    out_col <- paste0(col, "_std")
    df[[out_col]] <- (df[[col]] - scaler[[col]]$mean) / scaler[[col]]$sd
    df[[out_col]][is.na(df[[out_col]])] <- 0
  }

  out <- df %>%
    transmute(
      dom_binary, dom_prob, dom_category, StrainName, vaccine_code, decimal_year, Period,
      PB2 = PB2_Nonsyn_std,
      PB1 = PB1_Nonsyn_std,
      PA = PA_Nonsyn_std,
      HA = HA_Nonsyn_std,
      NP = NP_Nonsyn_std,
      NA_seg = NA_Nonsyn_std,
      M = M_Nonsyn_std,
      NS = NS_Nonsyn_std,
      PB2_syn = PB2_Syn_std,
      PB1_syn = PB1_Syn_std,
      PA_syn = PA_Syn_std,
      HA_syn = HA_Syn_std,
      NP_syn = NP_Syn_std,
      NA_syn = NA_Syn_std,
      M_syn = M_Syn_std,
      NS_syn = NS_Syn_std
    ) %>%
    mutate(
      HA_x_PB1 = HA * PB1
    )

  out
}

recalc_weights <- function(df) {
  zeros <- sum(df$dom_binary == 0, na.rm = TRUE)
  ones <- sum(df$dom_binary == 1, na.rm = TRUE)
  if (ones <= 0 || zeros <= 0) return(mutate(df, norm_wts = 1))
  df %>%
    mutate(
      norm_wts = ifelse(dom_binary == 1, zeros / ones, 1),
      norm_wts = norm_wts / mean(norm_wts, na.rm = TRUE)
    )
}

safe_auc <- function(y, p) {
  keep <- !is.na(y) & !is.na(p)
  y <- y[keep]
  p <- p[keep]
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
}

clip_prob <- function(p) {
  p <- as.numeric(p)
  p[p < 1e-6] <- 1e-6
  p[p > 1 - 1e-6] <- 1 - 1e-6
  p
}

# ---- 3. Feature schemes ----
feature_schemes <- list(
  "1_Selected_dN_Final_HA_x_PB1_plus_PB2_NA" = c("HA", "PB1", "PB2", "NA_seg", "HA_x_PB1"),
  "2_HA_dN_Only" = c("HA"),
  "3_All_dN_PB2_to_NS" = c("PB2", "PB1", "PA", "HA", "NP", "NA_seg", "M", "NS"),
  "4_All_dN_All_dS_PB2_to_NS" = c(
    "PB2", "PB1", "PA", "HA", "NP", "NA_seg", "M", "NS",
    "PB2_syn", "PB1_syn", "PA_syn", "HA_syn", "NP_syn", "NA_syn", "M_syn", "NS_syn"
  )
)

feature_scheme_defs <- tibble(
  Feature_Scheme = names(feature_schemes),
  Predictors = sapply(feature_schemes, paste, collapse = ","),
  N_Predictors = sapply(feature_schemes, length)
)

# ---- 4. Rolling window grid ----
make_rolling_windows <- function(all_periods, last_test_period, rolling_width) {
  max_idx <- which(all_periods == last_test_period)
  if (rolling_width >= max_idx) stop("rolling_width must be smaller than last_test_period index.")

  windows <- list()
  for (start_i in 1:(max_idx - rolling_width)) {
    end_i <- start_i + rolling_width - 1
    windows[[length(windows) + 1]] <- tibble(
      Window_Scheme = paste0("rolling_width_", rolling_width),
      Window_Type = "rolling",
      Window_Size = rolling_width,
      Train_Start_Index = start_i,
      Train_End_Index = end_i,
      Test_Index = end_i + 1,
      Train_Periods = paste(all_periods[start_i:end_i], collapse = ","),
      Test_Period = all_periods[end_i + 1],
      Train_Range_Label = paste0(all_periods[start_i], "-", all_periods[end_i], " -> ", all_periods[end_i + 1])
    )
  }

  bind_rows(windows) %>%
    mutate(Test_Period = as.character(Test_Period)) %>%
    filter(Test_Period != "CR23")
}

window_grid <- make_rolling_windows(all_periods, last_test_period, final_window_width)

# Inner temporal folds for hyperparameter tuning within each outer train split.
# Uses expanding inner training to avoid future-to-past leakage in validation.
make_inner_folds <- function(train_periods, min_train_periods = 6) {
  n <- length(train_periods)
  if (n <= min_train_periods) return(list())
  folds <- list()
  for (val_i in (min_train_periods + 1):n) {
    folds[[length(folds) + 1]] <- list(
      train_periods = train_periods[1:(val_i - 1)],
      val_period = train_periods[val_i]
    )
  }
  folds
}

get_scaled_xy_from_periods <- function(train_periods, eval_periods, predictors) {
  raw_train <- raw_df %>% filter(as.character(Period) %in% train_periods)
  raw_eval <- raw_df %>% filter(as.character(Period) %in% eval_periods)
  scaler <- fit_scaler(raw_train, all_raw_features)
  train_df <- apply_scaler(raw_train, scaler) %>% recalc_weights()
  eval_df <- apply_scaler(raw_eval, scaler)

  list(
    train_df = train_df,
    eval_df = eval_df,
    X_train = as.matrix(train_df %>% select(all_of(predictors))),
    y_train = train_df$dom_binary,
    w_train = train_df$norm_wts,
    X_eval = as.matrix(eval_df %>% select(all_of(predictors))),
    y_eval = eval_df$dom_binary
  )
}

# ---- 5. Model training / prediction functions ----
fit_logistic <- function(X, y, w) {
  dat <- as.data.frame(X)
  dat$y <- y
  # Use the vector w directly as case weights. Do not store .weights in dat,
  # because some glm calls can accidentally look for .weights in newdata/formula env.
  suppressWarnings(glm(y ~ ., data = dat, family = binomial(), weights = w))
}

pred_logistic <- function(fit, X) {
  dat <- as.data.frame(X)
  clip_prob(predict(fit, newdata = dat, type = "response"))
}

fit_glmnet_model <- function(X, y, w, alpha, lambda) {
  glmnet::glmnet(X, y, family = "binomial", alpha = alpha, lambda = lambda, weights = w, standardize = FALSE)
}

pred_glmnet_model <- function(fit, X, lambda) {
  clip_prob(as.numeric(predict(fit, newx = X, s = lambda, type = "response")))
}

fit_svm_model <- function(X, y, C, gamma) {
  dat <- as.data.frame(X)
  dat$y <- factor(y, levels = c(0, 1))
  zeros <- sum(y == 0, na.rm = TRUE)
  ones <- sum(y == 1, na.rm = TRUE)
  cw <- if (ones > 0 && zeros > 0) c("0" = 1, "1" = zeros / ones) else NULL
  e1071::svm(
    y ~ ., data = dat,
    kernel = "radial",
    cost = C,
    gamma = gamma,
    probability = TRUE,
    class.weights = cw
  )
}

pred_svm_model <- function(fit, X) {
  dat <- as.data.frame(X)
  pp <- predict(fit, newdata = dat, probability = TRUE)
  prob <- attr(pp, "probabilities")
  if ("1" %in% colnames(prob)) return(clip_prob(prob[, "1"]))
  clip_prob(prob[, ncol(prob)])
}

fit_rf_model <- function(X, y, mtry, ntree) {
  dat <- as.data.frame(X)
  yfac <- factor(y, levels = c(0, 1))
  zeros <- sum(y == 0, na.rm = TRUE)
  ones <- sum(y == 1, na.rm = TRUE)
  cw <- if (ones > 0 && zeros > 0) c("0" = 1, "1" = zeros / ones) else NULL
  randomForest::randomForest(
    x = dat,
    y = yfac,
    mtry = mtry,
    ntree = ntree,
    classwt = cw
  )
}

pred_rf_model <- function(fit, X) {
  dat <- as.data.frame(X)
  prob <- predict(fit, newdata = dat, type = "prob")
  if ("1" %in% colnames(prob)) return(clip_prob(prob[, "1"]))
  clip_prob(prob[, ncol(prob)])
}

fit_xgb_model <- function(X, y, w, params, nrounds) {
  dtrain <- xgboost::xgb.DMatrix(data = X, label = y, weight = w)
  xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )
}

pred_xgb_model <- function(fit, X) {
  dtest <- xgboost::xgb.DMatrix(data = X)
  clip_prob(predict(fit, dtest))
}

# DNN/Keras functions
# Notes:
# - Uses the same training-split-only scaling as all other ML models.
# - sample_weight uses the same class balancing weights as the classical ML models.
# - All calls are wrapped by tryCatch later, so a local TensorFlow/Keras issue will fail only DNN rows.
# ============================================================
# Keras 3 / TensorFlow 2.16 compatible DNN model
# Convert R matrix -> TensorFlow tensor explicitly
# ============================================================

fit_dnn_model <- function(X, y, w,
                          units = 8,
                          dropout = 0.1,
                          lr = 0.001,
                          epochs = 120,
                          batch_size = 32) {
  if (!isTRUE(run_dnn)) {
    stop("DNN disabled or keras/tensorflow not available")
  }
  
  keras::k_clear_session()
  
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  y <- as.numeric(y)
  w <- as.numeric(w)
  
  input_dim <- ncol(X)
  
  # Explicit TensorFlow tensors
  X_tf <- tensorflow::tf$convert_to_tensor(X, dtype = tensorflow::tf$float32)
  y_tf <- tensorflow::tf$convert_to_tensor(matrix(y, ncol = 1), dtype = tensorflow::tf$float32)
  w_tf <- tensorflow::tf$convert_to_tensor(w, dtype = tensorflow::tf$float32)
  
  input <- keras::layer_input(
    shape = c(input_dim),
    dtype = "float32",
    name = "input"
  )
  
  output <- input %>%
    keras::layer_dense(
      units = units,
      activation = "relu",
      name = "hidden_dense_1"
    ) %>%
    keras::layer_dropout(
      rate = dropout,
      name = "dropout_1"
    ) %>%
    keras::layer_dense(
      units = 1,
      activation = "sigmoid",
      name = "output"
    )
  
  model <- keras::keras_model(
    inputs = input,
    outputs = output
  )
  
  do.call(
    model$compile,
    list(
      optimizer = keras::optimizer_adam(learning_rate = lr),
      loss = "binary_crossentropy"
    )
  )
  
  do.call(
    model$fit,
    list(
      x = X_tf,
      y = y_tf,
      sample_weight = w_tf,
      epochs = as.integer(epochs),
      batch_size = as.integer(min(batch_size, nrow(X))),
      verbose = as.integer(dnn_verbose)
    )
  )
  
  model
}


pred_dnn_model <- function(fit, X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  X_tf <- tensorflow::tf$convert_to_tensor(X, dtype = tensorflow::tf$float32)
  
  # Avoid keras predict() wrapper; call model directly
  p_tf <- fit(X_tf, training = FALSE)
  
  p <- as.numeric(as.array(p_tf))
  
  clip_prob(p)
}

# ---- 6. Hyperparameter tuning ----
mean_auc_cv <- function(vals) {
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

tune_glmnet <- function(train_periods, predictors, alpha) {
  folds <- make_inner_folds(train_periods)
  lambda_grid <- 10 ^ seq(-4, 1, length.out = 25)
  if (length(folds) == 0) return(list(lambda = 0.01, CV_AUC = NA_real_))

  res <- map_dfr(lambda_grid, function(lambda) {
    aucs <- map_dbl(folds, function(fold) {
      dat <- get_scaled_xy_from_periods(fold$train_periods, fold$val_period, predictors)
      if (length(unique(dat$y_train)) < 2 || length(unique(dat$y_eval)) < 2) return(NA_real_)
      fit <- tryCatch(fit_glmnet_model(dat$X_train, dat$y_train, dat$w_train, alpha, lambda), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      p <- pred_glmnet_model(fit, dat$X_eval, lambda)
      safe_auc(dat$y_eval, p)
    })
    tibble(lambda = lambda, CV_AUC = mean_auc_cv(aucs))
  })

  best <- res %>% arrange(desc(CV_AUC), lambda) %>% slice(1)
  list(lambda = best$lambda, CV_AUC = best$CV_AUC)
}

tune_svm <- function(train_periods, predictors) {
  folds <- make_inner_folds(train_periods)
  p <- length(predictors)
  gamma_default <- 1 / max(p, 1)
  grid <- expand_grid(
    C = c(0.25, 1, 4, 16),
    gamma = gamma_default * c(0.25, 1, 4)
  )
  if (length(folds) == 0) return(list(C = 1, gamma = gamma_default, CV_AUC = NA_real_))

  res <- pmap_dfr(grid, function(C, gamma) {
    aucs <- map_dbl(folds, function(fold) {
      dat <- get_scaled_xy_from_periods(fold$train_periods, fold$val_period, predictors)
      if (length(unique(dat$y_train)) < 2 || length(unique(dat$y_eval)) < 2) return(NA_real_)
      fit <- tryCatch(fit_svm_model(dat$X_train, dat$y_train, C, gamma), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      p <- pred_svm_model(fit, dat$X_eval)
      safe_auc(dat$y_eval, p)
    })
    tibble(C = C, gamma = gamma, CV_AUC = mean_auc_cv(aucs))
  })

  best <- res %>% arrange(desc(CV_AUC), C, gamma) %>% slice(1)
  list(C = best$C, gamma = best$gamma, CV_AUC = best$CV_AUC)
}

tune_rf <- function(train_periods, predictors) {
  folds <- make_inner_folds(train_periods)
  p <- length(predictors)
  grid <- expand_grid(
    mtry = unique(pmax(1, pmin(p, c(1, floor(sqrt(p)), floor(p / 2), p)))),
    ntree = c(500, 1000)
  )
  if (length(folds) == 0) return(list(mtry = max(1, floor(sqrt(p))), ntree = 500, CV_AUC = NA_real_))

  res <- pmap_dfr(grid, function(mtry, ntree) {
    aucs <- map_dbl(folds, function(fold) {
      dat <- get_scaled_xy_from_periods(fold$train_periods, fold$val_period, predictors)
      if (length(unique(dat$y_train)) < 2 || length(unique(dat$y_eval)) < 2) return(NA_real_)
      fit <- tryCatch(fit_rf_model(dat$X_train, dat$y_train, mtry, ntree), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      p <- pred_rf_model(fit, dat$X_eval)
      safe_auc(dat$y_eval, p)
    })
    tibble(mtry = mtry, ntree = ntree, CV_AUC = mean_auc_cv(aucs))
  })

  best <- res %>% arrange(desc(CV_AUC), mtry, ntree) %>% slice(1)
  list(mtry = best$mtry, ntree = best$ntree, CV_AUC = best$CV_AUC)
}

tune_xgb <- function(train_periods, predictors) {
  folds <- make_inner_folds(train_periods)
  grid <- expand_grid(
    max_depth = c(1, 2, 3),
    eta = c(0.03, 0.1),
    nrounds = c(50, 100, 200),
    subsample = c(0.8),
    colsample_bytree = c(0.8)
  )
  if (length(folds) == 0) {
    return(list(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 2, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8), nrounds = 100, CV_AUC = NA_real_))
  }

  res <- pmap_dfr(grid, function(max_depth, eta, nrounds, subsample, colsample_bytree) {
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = max_depth,
      eta = eta,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      min_child_weight = 1
    )
    aucs <- map_dbl(folds, function(fold) {
      dat <- get_scaled_xy_from_periods(fold$train_periods, fold$val_period, predictors)
      if (length(unique(dat$y_train)) < 2 || length(unique(dat$y_eval)) < 2) return(NA_real_)
      fit <- tryCatch(fit_xgb_model(dat$X_train, dat$y_train, dat$w_train, params, nrounds), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      p <- pred_xgb_model(fit, dat$X_eval)
      safe_auc(dat$y_eval, p)
    })
    tibble(max_depth = max_depth, eta = eta, nrounds = nrounds, subsample = subsample, colsample_bytree = colsample_bytree, CV_AUC = mean_auc_cv(aucs))
  })

  best <- res %>% arrange(desc(CV_AUC), max_depth, eta, nrounds) %>% slice(1)
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = best$max_depth,
    eta = best$eta,
    subsample = best$subsample,
    colsample_bytree = best$colsample_bytree,
    min_child_weight = 1
  )
  list(params = params, nrounds = best$nrounds, CV_AUC = best$CV_AUC)
}

tune_dnn <- function(train_periods, predictors) {
  if (!isTRUE(run_dnn)) {
    return(list(units = NA_integer_, dropout = NA_real_, lr = NA_real_, epochs = NA_integer_, CV_AUC = NA_real_))
  }

  folds <- make_inner_folds(train_periods)
  grid <- expand_grid(
    units = dnn_units_grid,
    dropout = dnn_dropout_grid,
    lr = dnn_lr_grid
  )

  if (length(folds) == 0) {
    return(list(units = dnn_units_grid[1], dropout = dnn_dropout_grid[1], lr = dnn_lr_grid[1], epochs = dnn_epochs, CV_AUC = NA_real_))
  }

  res <- pmap_dfr(grid, function(units, dropout, lr) {
    aucs <- map_dbl(folds, function(fold) {
      dat <- get_scaled_xy_from_periods(fold$train_periods, fold$val_period, predictors)
      if (length(unique(dat$y_train)) < 2 || length(unique(dat$y_eval)) < 2) return(NA_real_)
      fit <- tryCatch(
        fit_dnn_model(dat$X_train, dat$y_train, dat$w_train, units = units, dropout = dropout, lr = lr, epochs = dnn_epochs, batch_size = dnn_batch_size),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      p <- pred_dnn_model(fit, dat$X_eval)
      safe_auc(dat$y_eval, p)
    })
    tibble(units = units, dropout = dropout, lr = lr, epochs = dnn_epochs, CV_AUC = mean_auc_cv(aucs))
  })

  if (nrow(res) == 0 || all(is.na(res$CV_AUC))) {
    return(list(units = dnn_units_grid[1], dropout = dnn_dropout_grid[1], lr = dnn_lr_grid[1], epochs = dnn_epochs, CV_AUC = NA_real_))
  }

  best <- res %>% arrange(desc(CV_AUC), units, dropout, lr) %>% slice(1)
  list(units = best$units, dropout = best$dropout, lr = best$lr, epochs = best$epochs, CV_AUC = best$CV_AUC)
}

# ---- 7. Fit/evaluate one outer split and one feature scheme ----
fit_ml_models_one_split_scheme <- function(win_row, feature_scheme, predictors) {
  train_periods <- strsplit(win_row$Train_Periods, ",")[[1]]
  test_period <- win_row$Test_Period

  raw_train <- raw_df %>% filter(as.character(Period) %in% train_periods)
  raw_test <- raw_df %>% filter(as.character(Period) == test_period)

  if (nrow(raw_train) < 10 || nrow(raw_test) < 5) return(NULL)
  if (length(unique(raw_train$dom_binary)) < 2 || length(unique(raw_test$dom_binary)) < 2) return(NULL)

  scaler <- fit_scaler(raw_train, all_raw_features)
  train_df <- apply_scaler(raw_train, scaler) %>% recalc_weights()
  test_df <- apply_scaler(raw_test, scaler)

  X_train <- as.matrix(train_df %>% select(all_of(predictors)))
  y_train <- train_df$dom_binary
  w_train <- train_df$norm_wts
  X_test <- as.matrix(test_df %>% select(all_of(predictors)))
  y_test <- test_df$dom_binary

  base_row <- tibble(
    Window_Scheme = win_row$Window_Scheme,
    Window_Type = win_row$Window_Type,
    Window_Size = win_row$Window_Size,
    Train_Start_Index = win_row$Train_Start_Index,
    Train_End_Index = win_row$Train_End_Index,
    Train_Periods = win_row$Train_Periods,
    Train_Range_Label = win_row$Train_Range_Label,
    Test_Index = win_row$Test_Index,
    Test_Period = test_period,
    Feature_Scheme = feature_scheme,
    Predictors = paste(predictors, collapse = ","),
    N_Predictors = length(predictors),
    Train_N = nrow(train_df),
    Train_Pos = sum(y_train == 1),
    Train_Neg = sum(y_train == 0),
    Test_N = nrow(test_df),
    Test_Pos = sum(y_test == 1),
    Test_Neg = sum(y_test == 0)
  )

  results <- list()
  selected_params <- list()

  # Logistic regression
  results[["Logistic"]] <- tryCatch({
    fit <- fit_logistic(X_train, y_train, w_train)
    p_train <- pred_logistic(fit, X_train)
    p_test <- pred_logistic(fit, X_test)
    bind_cols(base_row, tibble(Model = "Logistic", Tune_CV_AUC = NA_real_, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
  }, error = function(e) {
    bind_cols(base_row, tibble(Model = "Logistic", Tune_CV_AUC = NA_real_, Train_AUC = NA_real_, Test_AUC = NA_real_, Error = e$message))
  })

  # Ridge
  results[["Ridge"]] <- tryCatch({
    if (ncol(X_train) < 2) {
      # glmnet requires at least 2 columns. For HA-only, Ridge is effectively
      # redundant with weighted logistic regression, so use logistic fallback.
      fit <- fit_logistic(X_train, y_train, w_train)
      p_train <- pred_logistic(fit, X_train)
      p_test <- pred_logistic(fit, X_test)
      selected_params[["Ridge_fallback"]] <- bind_cols(base_row, tibble(Model = "Ridge", Param = "fallback", Value = "logistic_for_single_predictor", Tune_CV_AUC = NA_real_))
      bind_cols(base_row, tibble(Model = "Ridge", Tune_CV_AUC = NA_real_, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
    } else {
      tune <- tune_glmnet(train_periods, predictors, alpha = 0)
      fit <- fit_glmnet_model(X_train, y_train, w_train, alpha = 0, lambda = tune$lambda)
      p_train <- pred_glmnet_model(fit, X_train, tune$lambda)
      p_test <- pred_glmnet_model(fit, X_test, tune$lambda)
      selected_params[["Ridge"]] <- bind_cols(base_row, tibble(Model = "Ridge", Param = "lambda", Value = as.character(tune$lambda), Tune_CV_AUC = tune$CV_AUC))
      bind_cols(base_row, tibble(Model = "Ridge", Tune_CV_AUC = tune$CV_AUC, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
    }
  }, error = function(e) {
    bind_cols(base_row, tibble(Model = "Ridge", Tune_CV_AUC = NA_real_, Train_AUC = NA_real_, Test_AUC = NA_real_, Error = e$message))
  })

  # Lasso
  results[["Lasso"]] <- tryCatch({
    if (ncol(X_train) < 2) {
      # glmnet requires at least 2 columns. For HA-only, Lasso is effectively
      # redundant with weighted logistic regression, so use logistic fallback.
      fit <- fit_logistic(X_train, y_train, w_train)
      p_train <- pred_logistic(fit, X_train)
      p_test <- pred_logistic(fit, X_test)
      selected_params[["Lasso_fallback"]] <- bind_cols(base_row, tibble(Model = "Lasso", Param = "fallback", Value = "logistic_for_single_predictor", Tune_CV_AUC = NA_real_))
      bind_cols(base_row, tibble(Model = "Lasso", Tune_CV_AUC = NA_real_, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
    } else {
      tune <- tune_glmnet(train_periods, predictors, alpha = 1)
      fit <- fit_glmnet_model(X_train, y_train, w_train, alpha = 1, lambda = tune$lambda)
      p_train <- pred_glmnet_model(fit, X_train, tune$lambda)
      p_test <- pred_glmnet_model(fit, X_test, tune$lambda)
      selected_params[["Lasso"]] <- bind_cols(base_row, tibble(Model = "Lasso", Param = "lambda", Value = as.character(tune$lambda), Tune_CV_AUC = tune$CV_AUC))
      bind_cols(base_row, tibble(Model = "Lasso", Tune_CV_AUC = tune$CV_AUC, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
    }
  }, error = function(e) {
    bind_cols(base_row, tibble(Model = "Lasso", Tune_CV_AUC = NA_real_, Train_AUC = NA_real_, Test_AUC = NA_real_, Error = e$message))
  })

  # SVM radial
  results[["SVM_Radial"]] <- tryCatch({
    tune <- tune_svm(train_periods, predictors)
    fit <- fit_svm_model(X_train, y_train, tune$C, tune$gamma)
    p_train <- pred_svm_model(fit, X_train)
    p_test <- pred_svm_model(fit, X_test)
    selected_params[["SVM_Radial_C"]] <- bind_cols(base_row, tibble(Model = "SVM_Radial", Param = "C", Value = as.character(tune$C), Tune_CV_AUC = tune$CV_AUC))
    selected_params[["SVM_Radial_gamma"]] <- bind_cols(base_row, tibble(Model = "SVM_Radial", Param = "gamma", Value = as.character(tune$gamma), Tune_CV_AUC = tune$CV_AUC))
    bind_cols(base_row, tibble(Model = "SVM_Radial", Tune_CV_AUC = tune$CV_AUC, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
  }, error = function(e) {
    bind_cols(base_row, tibble(Model = "SVM_Radial", Tune_CV_AUC = NA_real_, Train_AUC = NA_real_, Test_AUC = NA_real_, Error = e$message))
  })

  # Random forest
  results[["RandomForest"]] <- tryCatch({
    tune <- tune_rf(train_periods, predictors)
    fit <- fit_rf_model(X_train, y_train, tune$mtry, tune$ntree)
    p_train <- pred_rf_model(fit, X_train)
    p_test <- pred_rf_model(fit, X_test)
    selected_params[["RandomForest_mtry"]] <- bind_cols(base_row, tibble(Model = "RandomForest", Param = "mtry", Value = as.character(tune$mtry), Tune_CV_AUC = tune$CV_AUC))
    selected_params[["RandomForest_ntree"]] <- bind_cols(base_row, tibble(Model = "RandomForest", Param = "ntree", Value = as.character(tune$ntree), Tune_CV_AUC = tune$CV_AUC))
    bind_cols(base_row, tibble(Model = "RandomForest", Tune_CV_AUC = tune$CV_AUC, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
  }, error = function(e) {
    bind_cols(base_row, tibble(Model = "RandomForest", Tune_CV_AUC = NA_real_, Train_AUC = NA_real_, Test_AUC = NA_real_, Error = e$message))
  })

  # XGBoost, if available
  if (isTRUE(run_xgboost)) {
    results[["XGBoost"]] <- tryCatch({
      tune <- tune_xgb(train_periods, predictors)
      fit <- fit_xgb_model(X_train, y_train, w_train, tune$params, tune$nrounds)
      p_train <- pred_xgb_model(fit, X_train)
      p_test <- pred_xgb_model(fit, X_test)
      selected_params[["XGBoost_nrounds"]] <- bind_cols(base_row, tibble(Model = "XGBoost", Param = "nrounds", Value = as.character(tune$nrounds), Tune_CV_AUC = tune$CV_AUC))
      selected_params[["XGBoost_max_depth"]] <- bind_cols(base_row, tibble(Model = "XGBoost", Param = "max_depth", Value = as.character(tune$params$max_depth), Tune_CV_AUC = tune$CV_AUC))
      selected_params[["XGBoost_eta"]] <- bind_cols(base_row, tibble(Model = "XGBoost", Param = "eta", Value = as.character(tune$params$eta), Tune_CV_AUC = tune$CV_AUC))
      bind_cols(base_row, tibble(Model = "XGBoost", Tune_CV_AUC = tune$CV_AUC, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
    }, error = function(e) {
      bind_cols(base_row, tibble(Model = "XGBoost", Tune_CV_AUC = NA_real_, Train_AUC = NA_real_, Test_AUC = NA_real_, Error = e$message))
    })
  }

  # DNN/Keras, if available
  if (isTRUE(run_dnn)) {
    results[["DNN"]] <- tryCatch({
      tune <- tune_dnn(train_periods, predictors)
      fit <- fit_dnn_model(
        X_train, y_train, w_train,
        units = tune$units,
        dropout = tune$dropout,
        lr = tune$lr,
        epochs = tune$epochs,
        batch_size = dnn_batch_size
      )
      p_train <- pred_dnn_model(fit, X_train)
      p_test <- pred_dnn_model(fit, X_test)
      selected_params[["DNN_units"]] <- bind_cols(base_row, tibble(Model = "DNN", Param = "units", Value = as.character(tune$units), Tune_CV_AUC = tune$CV_AUC))
      selected_params[["DNN_dropout"]] <- bind_cols(base_row, tibble(Model = "DNN", Param = "dropout", Value = as.character(tune$dropout), Tune_CV_AUC = tune$CV_AUC))
      selected_params[["DNN_lr"]] <- bind_cols(base_row, tibble(Model = "DNN", Param = "learning_rate", Value = as.character(tune$lr), Tune_CV_AUC = tune$CV_AUC))
      selected_params[["DNN_epochs"]] <- bind_cols(base_row, tibble(Model = "DNN", Param = "epochs", Value = as.character(tune$epochs), Tune_CV_AUC = tune$CV_AUC))
      bind_cols(base_row, tibble(Model = "DNN", Tune_CV_AUC = tune$CV_AUC, Train_AUC = safe_auc(y_train, p_train), Test_AUC = safe_auc(y_test, p_test), Error = NA_character_))
    }, error = function(e) {
      bind_cols(base_row, tibble(Model = "DNN", Tune_CV_AUC = NA_real_, Train_AUC = NA_real_, Test_AUC = NA_real_, Error = e$message))
    })
  }

  list(
    results = bind_rows(results),
    selected_params = bind_rows(selected_params)
  )
}

# ---- 8. Run all ML comparisons ----
message("Outer rolling splits: ", nrow(window_grid))
message("Feature schemes: ", length(feature_schemes))
message("XGBoost enabled: ", run_xgboost)
message("DNN enabled: ", run_dnn)

run_list <- list()
param_list <- list()
run_id <- 1

for (i in seq_len(nrow(window_grid))) {
  win_row <- window_grid[i, ]
  message("\nOuter split ", i, "/", nrow(window_grid), ": ", win_row$Train_Range_Label)

  for (scheme_name in names(feature_schemes)) {
    message("  Feature scheme: ", scheme_name)
    predictors <- feature_schemes[[scheme_name]]
    tmp <- fit_ml_models_one_split_scheme(win_row, scheme_name, predictors)
    run_list[[run_id]] <- tmp$results
    param_list[[run_id]] <- tmp$selected_params
    run_id <- run_id + 1
  }
}

ml_auc_by_split <- bind_rows(run_list) %>%
  mutate(
    Error = as.character(Error),
    Has_Test_AUC = !is.na(Test_AUC) & !is.infinite(Test_AUC),
    Run_Status = if_else(Has_Test_AUC, "OK", "FAILED")
  )
ml_selected_params <- bind_rows(param_list)

ml_auc_summary <- ml_auc_by_split %>%
  filter(Has_Test_AUC) %>%
  group_by(Feature_Scheme, Model, Predictors, N_Predictors) %>%
  summarise(
    N_Splits = n(),
    Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    Mean_Train_AUC = mean(Train_AUC, na.rm = TRUE),
    Median_Train_AUC = median(Train_AUC, na.rm = TRUE),
    Min_Train_AUC = ifelse(all(is.na(Train_AUC)), NA_real_, min(Train_AUC, na.rm = TRUE)),
    SD_Train_AUC = sd(Train_AUC, na.rm = TRUE),
    Mean_Test_AUC = mean(Test_AUC, na.rm = TRUE),
    Median_Test_AUC = median(Test_AUC, na.rm = TRUE),
    Min_Test_AUC = ifelse(all(is.na(Test_AUC)), NA_real_, min(Test_AUC, na.rm = TRUE)),
    SD_Test_AUC = sd(Test_AUC, na.rm = TRUE),
    Mean_Train_minus_Test_AUC = mean(Train_AUC - Test_AUC, na.rm = TRUE),
    Median_Train_minus_Test_AUC = median(Train_AUC - Test_AUC, na.rm = TRUE),
    Mean_Tune_CV_AUC = mean(Tune_CV_AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Test_AUC), desc(Median_Test_AUC))

# ---- 9. Compare against final selected regression model ----
regression_selected_by_split <- NULL

if (file.exists(regression_table_f2)) {
  regression_selected_by_split <- read_csv(regression_table_f2, show_col_types = FALSE) %>%
    filter(Model_Name == regression_selected_model_name) %>%
    select(
      Test_Period,
      Test_Index,
      Regression_Selected_Train_AUC_FixedOnly = Train_AUC_FixedOnly,
      Regression_Selected_Test_AUC_FixedOnly = Test_AUC_FixedOnly
    )

  ml_delta_vs_regression_by_split <- ml_auc_by_split %>%
    left_join(regression_selected_by_split, by = c("Test_Period", "Test_Index")) %>%
    mutate(
      Delta_Test_AUC_vs_Selected_Regression = Test_AUC - Regression_Selected_Test_AUC_FixedOnly,
      Delta_Train_AUC_vs_Selected_Regression = Train_AUC - Regression_Selected_Train_AUC_FixedOnly
    )

  ml_delta_vs_regression_summary <- ml_delta_vs_regression_by_split %>%
    filter(Has_Test_AUC, !is.na(Regression_Selected_Test_AUC_FixedOnly)) %>%
    group_by(Feature_Scheme, Model, Predictors, N_Predictors) %>%
    summarise(
      N_Splits = n(),
      Mean_Delta_Test_AUC_vs_Selected_Regression = mean(Delta_Test_AUC_vs_Selected_Regression, na.rm = TRUE),
      Median_Delta_Test_AUC_vs_Selected_Regression = median(Delta_Test_AUC_vs_Selected_Regression, na.rm = TRUE),
      N_Test_AUC_Wins_vs_Selected_Regression = sum(Delta_Test_AUC_vs_Selected_Regression > 0, na.rm = TRUE),
      N_Test_AUC_Losses_vs_Selected_Regression = sum(Delta_Test_AUC_vs_Selected_Regression < 0, na.rm = TRUE),
      Wilcoxon_p_Test_AUC_vs_Selected_Regression = tryCatch(
        wilcox.test(Test_AUC, Regression_Selected_Test_AUC_FixedOnly, paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Delta_Test_AUC_vs_Selected_Regression))
} else {
  warning("Regression Table_F2 not found. Delta-vs-regression tables will be skipped: ", regression_table_f2)
  ml_delta_vs_regression_by_split <- tibble()
  ml_delta_vs_regression_summary <- tibble()
}

# ---- 10. Export ----
write.csv(window_grid, file.path(result_dir, "Table_ML0_RollingWidth10_Window_Grid.csv"), row.names = FALSE)
write.csv(feature_scheme_defs, file.path(result_dir, "Table_ML1_Feature_Scheme_Definitions.csv"), row.names = FALSE)
write.csv(ml_auc_by_split, file.path(result_dir, "Table_ML2_ML_AUROC_By_Rolling_Split.csv"), row.names = FALSE)
write.csv(ml_auc_by_split %>% count(Model, Error, Run_Status, sort = TRUE), file.path(result_dir, "Table_ML2B_Model_Error_Diagnostics.csv"), row.names = FALSE)
write.csv(ml_auc_summary, file.path(result_dir, "Table_ML3_ML_AUROC_Summary.csv"), row.names = FALSE)
write.csv(ml_selected_params, file.path(result_dir, "Table_ML4_Selected_Hyperparameters.csv"), row.names = FALSE)
write.csv(ml_delta_vs_regression_by_split, file.path(result_dir, "Table_ML5_Delta_vs_Selected_Regression_By_Split.csv"), row.names = FALSE)
write.csv(ml_delta_vs_regression_summary, file.path(result_dir, "Table_ML6_Delta_vs_Selected_Regression_Summary.csv"), row.names = FALSE)

# ---- 11. Basic figures ----
fig_dir <- file.path(result_dir, "Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Summary barplot: mean test AUROC by feature scheme and model
plot_summary <- ml_auc_summary %>%
  mutate(
    Model = factor(Model, levels = c("Logistic", "Ridge", "Lasso", "SVM_Radial", "RandomForest", "XGBoost", "DNN")),
    Feature_Scheme_Label = case_when(
      Feature_Scheme == "1_Selected_dN_Final_HA_x_PB1_plus_PB2_NA" ~ "Selected dN\nHA×PB1+PB2+NA",
      Feature_Scheme == "2_HA_dN_Only" ~ "HA dN\nonly",
      Feature_Scheme == "3_All_dN_PB2_to_NS" ~ "All dN",
      Feature_Scheme == "4_All_dN_All_dS_PB2_to_NS" ~ "All dN\n+ all dS",
      TRUE ~ Feature_Scheme
    ),
    Feature_Scheme_Label = factor(
      Feature_Scheme_Label,
      levels = c("Selected dN\nHA×PB1+PB2+NA", "HA dN\nonly", "All dN", "All dN\n+ all dS")
    )
  )

p1 <- ggplot(plot_summary, aes(x = Model, y = Mean_Test_AUC, fill = Feature_Scheme_Label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4, color = "gray55") +
  geom_col(position = position_dodge(width = 0.78), width = 0.68, color = "black", linewidth = 0.2, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = Mean_Test_AUC - SD_Test_AUC, ymax = Mean_Test_AUC + SD_Test_AUC),
    position = position_dodge(width = 0.78),
    width = 0.18,
    linewidth = 0.4
  ) +
  coord_cartesian(ylim = c(0.35, 1.0)) +
  labs(
    title = "ML model performance under rolling-width-10 validation",
    subtitle = "Bars show mean test AUROC; error bars show ±1 SD across rolling splits",
    x = NULL,
    y = "Mean test AUROC",
    fill = "Feature scheme"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "Figure_ML1_Mean_Test_AUROC_By_Model_FeatureScheme.svg"), p1, width = 11, height = 6.5, dpi = 400)
ggsave(file.path(fig_dir, "Figure_ML1_Mean_Test_AUROC_By_Model_FeatureScheme.pdf"), p1, width = 11, height = 6.5)

# Delta vs selected regression
if (nrow(ml_delta_vs_regression_summary) > 0) {
  delta_plot <- ml_delta_vs_regression_summary %>%
    mutate(
      Label = paste(Model, Feature_Scheme, sep = " | ") %>% stringr::str_replace_all("_", " ") %>% stringr::str_wrap(width = 45),
      Label = fct_reorder(Label, Mean_Delta_Test_AUC_vs_Selected_Regression)
    ) %>%
    arrange(desc(abs(Mean_Delta_Test_AUC_vs_Selected_Regression))) %>%
    slice_head(n = 30)

  p2 <- ggplot(delta_plot, aes(x = Label, y = Mean_Delta_Test_AUC_vs_Selected_Regression)) +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "gray45") +
    geom_col(width = 0.72, color = "black", linewidth = 0.2, fill = "#D55E00", alpha = 0.9) +
    coord_flip() +
    labs(
      title = "ML performance relative to the selected regression model",
      subtitle = "Positive values indicate higher mean test AUROC than the selected mixed-effects regression model",
      x = NULL,
      y = expression(Delta~"mean test AUROC vs selected regression")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title.x = element_text(face = "bold"),
      axis.text.y = element_text(size = 8.5)
    )

  ggsave(file.path(fig_dir, "Figure_ML2_Delta_vs_Selected_Regression.svg"), p2, width = 10, height = 8, dpi = 400)
  ggsave(file.path(fig_dir, "Figure_ML2_Delta_vs_Selected_Regression.pdf"), p2, width = 10, height = 8)
}


# ---- 12. Manuscript-style ML/DL + regression figures: A/B/C ----
# A: train vs test AUROC divergence using the final selected feature scheme.
# B: mean predictive performance across feature schemes.
# C: model x feature-scheme heatmap.

message("\nBuilding manuscript-style ML/DL + regression figures...")

fig_dir2 <- file.path(result_dir, "Figures_Manuscript_With_Regression")
dir.create(fig_dir2, showWarnings = FALSE, recursive = TRUE)

feature_label_map <- tibble::tribble(
  ~Feature_Scheme, ~Feature_Display,
  "1_Selected_dN_Final_HA_x_PB1_plus_PB2_NA", "Selected dN\nHA × PB1 + PB2 + NA",
  "2_HA_dN_Only", "HA dN\nonly",
  "3_All_dN_PB2_to_NS", "All dN",
  "4_All_dN_All_dS_PB2_to_NS", "All dN\n+ all dS"
)

model_display_fun <- function(x) {
  dplyr::case_when(
    x == "Regression_GLMM" ~ "GLMM Baseline",
    x == "Logistic" ~ "Logistic Classifier",
    x == "Ridge" ~ "Ridge Classifier",
    x == "Lasso" ~ "Lasso Classifier",
    x == "SVM_Radial" ~ "Support Vector Machine (SVC)",
    x == "RandomForest" ~ "Random Forest (RFC)",
    x == "XGBoost" ~ "XGBoost Classifier",
    x == "DNN" ~ "Deep Neural Network (DNN)",
    TRUE ~ x
  )
}

reg_model_map <- tibble::tribble(
  ~Model_Name, ~Feature_Scheme, ~Predictors, ~N_Predictors,
  "1_Final_Rolling_Best_dN_HA_x_PB1_plus_PB2_NA", "1_Selected_dN_Final_HA_x_PB1_plus_PB2_NA", "HA,PB1,PB2,NA_seg,HA_x_PB1", 5L,
  "2_HA_dN_Only", "2_HA_dN_Only", "HA", 1L,
  "3_Control1_All_dN_PB2_to_NS", "3_All_dN_PB2_to_NS", "PB2,PB1,PA,HA,NP,NA_seg,M,NS", 8L,
  "4_Control2_All_dN_All_dS_PB2_to_NS", "4_All_dN_All_dS_PB2_to_NS", "PB2,PB1,PA,HA,NP,NA_seg,M,NS,PB2_syn,PB1_syn,PA_syn,HA_syn,NP_syn,NA_syn,M_syn,NS_syn", 16L
)

# ML/DL split table
ml_split_for_fig <- ml_auc_by_split %>%
  filter(Has_Test_AUC) %>%
  left_join(feature_label_map, by = "Feature_Scheme") %>%
  mutate(
    Model_Group = "ML/DL",
    Model_Display = model_display_fun(Model)
  )

# Regression split table, if available
if (file.exists(regression_table_f2)) {
  reg_split_for_fig <- readr::read_csv(regression_table_f2, show_col_types = FALSE) %>%
    inner_join(reg_model_map, by = "Model_Name") %>%
    left_join(feature_label_map, by = "Feature_Scheme") %>%
    transmute(
      Window_Scheme,
      Window_Type,
      Window_Size,
      Train_Start_Index,
      Train_End_Index,
      Train_Periods,
      Train_Range_Label,
      Test_Index,
      Test_Period,
      Feature_Scheme,
      Predictors,
      N_Predictors,
      Train_N,
      Train_Pos,
      Train_Neg,
      Test_N,
      Test_Pos,
      Test_Neg,
      Model = "Regression_GLMM",
      Tune_CV_AUC = NA_real_,
      Train_AUC = Train_AUC_FixedOnly,
      Test_AUC = Test_AUC_FixedOnly,
      Error = NA_character_,
      Has_Test_AUC = !is.na(Test_AUC) & !is.infinite(Test_AUC),
      Run_Status = if_else(Has_Test_AUC, "OK", "FAILED"),
      Feature_Display,
      Model_Group = "Regression",
      Model_Display = "GLMM Baseline"
    ) %>%
    filter(Has_Test_AUC)
} else {
  reg_split_for_fig <- tibble()
  warning("Regression Table_F2 not found. A/B/C manuscript figures will use ML/DL models only.")
}

combined_split <- bind_rows(ml_split_for_fig, reg_split_for_fig) %>%
  mutate(
    Feature_Display = factor(
      Feature_Display,
      levels = c("Selected dN\nHA × PB1 + PB2 + NA", "HA dN\nonly", "All dN", "All dN\n+ all dS")
    ),
    Model_Group = factor(Model_Group, levels = c("Regression", "ML/DL"))
  )

combined_summary <- combined_split %>%
  filter(!is.na(Test_AUC), !is.infinite(Test_AUC)) %>%
  group_by(Feature_Scheme, Feature_Display, Model, Model_Display, Model_Group, Predictors, N_Predictors) %>%
  summarise(
    N_Splits = n(),
    Mean_Train_AUC = mean(Train_AUC, na.rm = TRUE),
    Median_Train_AUC = median(Train_AUC, na.rm = TRUE),
    SD_Train_AUC = sd(Train_AUC, na.rm = TRUE),
    Mean_Test_AUC = mean(Test_AUC, na.rm = TRUE),
    Median_Test_AUC = median(Test_AUC, na.rm = TRUE),
    SD_Test_AUC = sd(Test_AUC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(combined_split, file.path(result_dir, "Table_ML7_ML_DL_plus_Regression_By_Rolling_Split.csv"))
readr::write_csv(combined_summary, file.path(result_dir, "Table_ML8_ML_DL_plus_Regression_Summary.csv"))

# -----------------------------
# Panel A: Train vs test AUROC
# -----------------------------
selected_feature_scheme <- "1_Selected_dN_Final_HA_x_PB1_plus_PB2_NA"

panel_a_df <- combined_split %>%
  filter(Feature_Scheme == selected_feature_scheme) %>%
  select(Test_Period, Test_Index, Model, Model_Display, Model_Group, Train_AUC, Test_AUC) %>%
  pivot_longer(cols = c(Train_AUC, Test_AUC), names_to = "Phase", values_to = "AUROC") %>%
  mutate(
    Phase = recode(Phase, "Train_AUC" = "Train (In-Period)", "Test_AUC" = "Test (Out-of-Period)"),
    Phase = factor(Phase, levels = c("Train (In-Period)", "Test (Out-of-Period)"))
  ) %>%
  filter(!is.na(AUROC), !is.infinite(AUROC))

model_order_a <- panel_a_df %>%
  filter(Phase == "Test (Out-of-Period)") %>%
  group_by(Model_Display) %>%
  summarise(Mean_Test_AUC = mean(AUROC, na.rm = TRUE), .groups = "drop") %>%
  arrange(Mean_Test_AUC) %>%
  pull(Model_Display)

panel_a_df <- panel_a_df %>%
  mutate(Model_Display = factor(Model_Display, levels = model_order_a))

phase_colors <- c("Train (In-Period)" = "#BDBDBD", "Test (Out-of-Period)" = "#67A9CF")

pA <- ggplot(panel_a_df, aes(x = AUROC, y = Model_Display, fill = Phase)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.4, color = "gray50") +
  geom_boxplot(
    position = position_dodge2(width = 0.72, preserve = "single"),
    width = 0.58,
    alpha = 0.95,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.35
  ) +
  geom_point(
    aes(color = Phase),
    position = position_jitterdodge(jitter.width = 0.015, dodge.width = 0.72),
    size = 1.7,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = phase_colors) +
  scale_color_manual(values = phase_colors) +
  coord_cartesian(xlim = c(0.4, 1.0)) +
  labs(
    title = "Generalization capacity: train vs test AUROC divergence",
    subtitle = "Final selected dN feature scheme: HA × PB1 + PB2 + NA",
    x = "Predictive AUROC",
    y = NULL,
    fill = "Phase"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 10.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir2, "Figure_MLDL_A_Train_vs_Test_AUROC_SelectedFeatureScheme.svg"), pA, width = 8.4, height = 6.8, dpi = 400)
ggsave(file.path(fig_dir2, "Figure_MLDL_A_Train_vs_Test_AUROC_SelectedFeatureScheme.pdf"), pA, width = 8.4, height = 6.8)

# -----------------------------
# Panel B: performance scaling across feature schemes
# Mean across model architectures based on each model's mean test AUROC.
# -----------------------------
panel_b_df <- combined_summary %>%
  filter(!is.na(Mean_Test_AUC), !is.infinite(Mean_Test_AUC)) %>%
  group_by(Feature_Display) %>%
  summarise(
    N_Models = n(),
    Mean_AUC = mean(Mean_Test_AUC, na.rm = TRUE),
    SD_AUC = sd(Mean_Test_AUC, na.rm = TRUE),
    SE_AUC = SD_AUC / sqrt(N_Models),
    CI95_Low = Mean_AUC - 1.96 * SE_AUC,
    CI95_High = Mean_AUC + 1.96 * SE_AUC,
    .groups = "drop"
  ) %>%
  mutate(
    Feature_Display = factor(Feature_Display, levels = levels(combined_split$Feature_Display)),
    CI95_Low = pmax(CI95_Low, 0),
    CI95_High = pmin(CI95_High, 1)
  )

pB <- ggplot(panel_b_df, aes(x = Feature_Display, y = Mean_AUC, group = 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4, color = "gray55") +
  geom_errorbar(aes(ymin = CI95_Low, ymax = CI95_High), width = 0.16, linewidth = 0.7, color = "#2C7FB8") +
  geom_line(linewidth = 0.8, color = "#2C7FB8") +
  geom_point(size = 2.8, color = "#2C7FB8") +
  geom_text(aes(label = sprintf("%.3f", Mean_AUC)), vjust = -1.1, size = 3.2, fontface = "bold") +
  coord_cartesian(ylim = c(0.45, 0.85)) +
  labs(
    title = "Performance scaling across genomic feature schemes",
    subtitle = "Mean test AUROC averaged across regression, ML, and DNN architectures; error bars show 95% CI across architectures",
    x = "Target predictor feature scheme",
    y = "Mean test AUROC"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", size = 10)
  )

ggsave(file.path(fig_dir2, "Figure_MLDL_B_Performance_Scaling_FeatureSchemes.svg"), pB, width = 8.2, height = 5.3, dpi = 400)
ggsave(file.path(fig_dir2, "Figure_MLDL_B_Performance_Scaling_FeatureSchemes.pdf"), pB, width = 8.2, height = 5.3)

# -----------------------------
# Panel C: architecture x feature scheme heatmap
# -----------------------------
model_levels_c <- c(
  "Ridge Classifier",
  "Lasso Classifier",
  "Logistic Classifier",
  "Support Vector Machine (SVC)",
  "Random Forest (RFC)",
  "XGBoost Classifier",
  "Deep Neural Network (DNN)",
  "GLMM Baseline"
)

panel_c_df <- combined_summary %>%
  mutate(
    Model_Display = factor(Model_Display, levels = rev(model_levels_c)),
    Feature_Display = factor(Feature_Display, levels = levels(combined_split$Feature_Display))
  ) %>%
  filter(!is.na(Model_Display), !is.na(Feature_Display))

pC <- ggplot(panel_c_df, aes(x = Feature_Display, y = Model_Display, fill = Mean_Test_AUC)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", Mean_Test_AUC)), size = 3.2, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#3B1F6D", "#2C7FB8", "#41AB5D", "#FDE725"),
    limits = c(0.45, 0.90),
    oob = scales::squish,
    name = "Mean test\nAUROC"
  ) +
  labs(
    title = "Architecture-feature interaction matrix",
    subtitle = "Cells show mean out-of-period test AUROC across rolling-width-10 splits",
    x = "Target predictor feature scheme",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 10.5),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir2, "Figure_MLDL_C_Architecture_FeatureScheme_Heatmap.svg"), pC, width = 8.7, height = 6.5, dpi = 400)
ggsave(file.path(fig_dir2, "Figure_MLDL_C_Architecture_FeatureScheme_Heatmap.pdf"), pC, width = 8.7, height = 6.5)

# Combined A/B/C panel if patchwork is available
if (requireNamespace("patchwork", quietly = TRUE)) {
  pABC <- (pA | (pB / pC)) + patchwork::plot_annotation(tag_levels = "A")
  ggsave(file.path(fig_dir2, "Figure_MLDL_ABC_Combined.svg"), pABC, width = 14.5, height = 10.5, dpi = 400)
  ggsave(file.path(fig_dir2, "Figure_MLDL_ABC_Combined.pdf"), pABC, width = 14.5, height = 10.5)
}
