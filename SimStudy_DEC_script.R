# ==============================================================================
# DEC Algorithm with Hyperparameter Tuning (Fixed for Keras 3)
# ==============================================================================

require(here)
require(IBclust)
require(kamila)
require(clustMixType)
require(FactoMineR)
require(cluster)
require(aricode)
require(keras3)
require(magrittr)

# ==============================================================================
# 1. DEC CLASS DEFINITIONS
# ==============================================================================

# --- Helper: Autoencoder Builder ---
createAutoencoderModel <- function(numberOfUnitsPerLayer, activation = 'relu', initializer = 'glorot_uniform') {
  numberOfUnitsPerLayer <- as.integer(numberOfUnitsPerLayer)
  
  numberOfEncodingLayers <- length(numberOfUnitsPerLayer) - 1
  
  inputs <- layer_input(shape = numberOfUnitsPerLayer[1], name = 'input')
  encoder <- inputs
  
  for (i in seq_len(numberOfEncodingLayers - 1)) {
    encoder <- encoder %>% layer_dense(units = numberOfUnitsPerLayer[i+1], activation = activation, kernel_initializer = initializer)
  }
  encoder <- encoder %>% layer_dense(units = tail(numberOfUnitsPerLayer, 1), name = 'bottleneck')
  
  autoencoder <- encoder
  for (i in seq(from = numberOfEncodingLayers, to = 2, by = -1)) {
    autoencoder <- autoencoder %>% layer_dense(units = numberOfUnitsPerLayer[i], activation = activation, kernel_initializer = initializer)
  }
  autoencoder <- autoencoder %>% layer_dense(numberOfUnitsPerLayer[1], kernel_initializer = initializer, name = 'decoder')
  
  return(list(autoencoderModel = keras_model(inputs = inputs, outputs = autoencoder),
              encoderModel = keras_model(inputs = inputs, outputs = encoder)))
}

# --- Helper: Custom Clustering Layer ---
ClusteringLayer <- R6::R6Class("ClusteringLayer",
                               inherit = KerasLayer,
                               lock_objects = FALSE,
                               public = list(
                                 numberOfClusters = 10L, 
                                 initialClusterWeights = NULL,
                                 alpha = 1.0,
                                 name = '',
                                 
                                 initialize = function(numberOfClusters, initialClusterWeights = NULL, alpha = 1.0, name = '') {
                                   self$numberOfClusters <- as.integer(numberOfClusters)
                                   self$initialClusterWeights <- initialClusterWeights
                                   self$alpha <- alpha
                                   self$name <- name
                                 },
                                 
                                 build = function(input_shape) {
                                   self$clusters <- self$add_weight(
                                     shape = list(self$numberOfClusters, as.integer(input_shape[[2]])),
                                     initializer = 'glorot_uniform', 
                                     name = 'clusters'
                                   )
                                   if (!is.null(self$initialClusterWeights)) {
                                     self$set_weights(self$initialClusterWeights)
                                     self$initialClusterWeights <- NULL
                                   }
                                   self$built <- TRUE
                                 },
                                 
                                 call = function(inputs, mask = NULL) {
                                   K <- keras::backend()
                                   q <- 1.0 / (1.0 + (K$sum(K$square(K$expand_dims(inputs, axis = 1L) - self$clusters), axis = 2L) / self$alpha))
                                   q <- q^((self$alpha + 1.0) / 2.0)
                                   q <- K$transpose(K$transpose(q) / K$sum(q, axis = 1L))
                                   return(q)
                                 },
                                 
                                 compute_output_shape = function(input_shape) {
                                   return(list(input_shape[[1]], self$numberOfClusters))
                                 }
                               )
)

layer_clustering <- function(object, numberOfClusters, initialClusterWeights = NULL, alpha = 1.0, name = '') {
  create_layer(ClusteringLayer, object, list(
    numberOfClusters = as.integer(numberOfClusters),  
    initialClusterWeights = initialClusterWeights,
    alpha = alpha, 
    name = name
  ))
}

# --- Main Class: Deep Embedded Clustering Model ---
DeepEmbeddedClusteringModel <- R6::R6Class("DeepEmbeddedClusteringModel",
                                           public = list(
                                             numberOfUnitsPerLayer = NULL,
                                             numberOfClusters = 10L, 
                                             alpha = 1.0,
                                             initializer = 'glorot_uniform',
                                             activation = 'relu',
                                             autoencoder = NULL,
                                             encoder = NULL,
                                             model = NULL,
                                             
                                             initialize = function(numberOfUnitsPerLayer, numberOfClusters, alpha = 1.0, 
                                                                   initializer = 'glorot_uniform', activation = 'relu') {
                                               self$numberOfUnitsPerLayer <- as.integer(numberOfUnitsPerLayer)
                                               self$numberOfClusters <- as.integer(numberOfClusters)
                                               self$alpha <- alpha
                                               self$initializer <- initializer
                                               self$activation <- activation
                                               ae <- createAutoencoderModel(self$numberOfUnitsPerLayer, 
                                                                            activation = self$activation,
                                                                            initializer = self$initializer)
                                               self$autoencoder <- ae$autoencoderModel
                                               self$encoder <- ae$encoderModel
                                               
                                               clusteringLayer <- self$encoder$output %>% layer_clustering(self$numberOfClusters, name = "clustering", alpha = self$alpha)
                                               self$model <- keras_model(inputs = self$encoder$input, outputs = clusteringLayer)
                                             },
                                             
                                             pretrain = function(x, optimizer = 'adam', epochs = 200L, batchSize = 256L) {
                                               self$autoencoder$compile(optimizer = optimizer, loss = 'mse')
                                               batchSize <- as.integer(min(batchSize, nrow(x)))
                                               epochs <- as.integer(epochs)
                                               tryCatch({
                                                 self$autoencoder$fit(x, x, batch_size = batchSize, epochs = epochs, verbose = 0)
                                               }, error = function(e) { cat("Pretraining error:", e$message, "\n") })
                                             },
                                             
                                             targetDistribution = function(q) {
                                               weight <- q^2 / colSums(q)
                                               return(t(t(weight) / rowSums(weight)))
                                             },
                                             
                                             compile = function(optimizer = 'sgd', loss = 'kld') {
                                               self$model$compile(optimizer = optimizer, loss = loss)
                                             },
                                             
                                             fit = function(x, maxNumberOfIterations = 5000L, batchSize = 256L, tolerance = 1e-3, updateInterval = 50L) {
                                               maxNumberOfIterations <- as.integer(maxNumberOfIterations)
                                               batchSize <- as.integer(min(batchSize, nrow(x)))
                                               updateInterval <- as.integer(updateInterval)
                                               
                                               km <- stats::kmeans(self$encoder$predict(x, verbose = 0), centers = self$numberOfClusters, nstart = 100)
                                               currentPrediction <- km$cluster
                                               previousPrediction <- currentPrediction
                                               self$model$get_layer(name = 'clustering')$set_weights(list(km$centers))
                                               
                                               loss <- 0
                                               
                                               for (i in seq_len(maxNumberOfIterations)) {
                                                 if (i %% updateInterval == 1L) {
                                                   q <- self$model$predict(x, verbose = 0)
                                                   p <- self$targetDistribution(q)
                                                   currentPrediction <- max.col(q)
                                                   deltaLabel <- sum(currentPrediction != previousPrediction) / length(currentPrediction)
                                                   previousPrediction <- currentPrediction
                                                   if (i > 1L && deltaLabel < tolerance) break
                                                 }
                                                 idx <- sample(1:nrow(x), batchSize)
                                                 loss <- self$model$train_on_batch(x = x[idx, , drop=FALSE], y = p[idx, , drop=FALSE])
                                               }
                                               return(currentPrediction)
                                             }
                                           )
)

# ==============================================================================
# 2. HYPERPARAMETER TUNING FUNCTION
# ==============================================================================

#' Run DEC with a specific hyperparameter configuration
#' @param x_dec One-hot encoded data matrix
#' @param num_clust Number of clusters
#' @param params List of hyperparameters
#' @return Cluster labels or NA if failed
run_dec_with_params <- function(x_dec, num_clust, params) {
  
  num_clust <- as.integer(num_clust)
  input_dim <- as.integer(ncol(x_dec))
  
  # Build architecture based on params - all dimensions as integers
  if (params$architecture == "shallow") {
    dims <- as.integer(c(input_dim, 64L, num_clust))
  } else if (params$architecture == "medium") {
    dims <- as.integer(c(input_dim, 128L, 64L, num_clust))
  } else if (params$architecture == "deep") {
    dims <- as.integer(c(input_dim, 256L, 128L, 64L, num_clust))
  } else if (params$architecture == "wide") {
    dims <- as.integer(c(input_dim, 256L, 128L, num_clust))
  } else if (params$architecture == "bottleneck") {
    dims <- as.integer(c(input_dim, 128L, 32L, num_clust))
  }
  
  tryCatch({
    decModel <- DeepEmbeddedClusteringModel$new(
      numberOfUnitsPerLayer = dims,
      numberOfClusters = num_clust,
      alpha = params$alpha,
      activation = params$activation,
      initializer = params$initializer
    )
    
    # Pretrain Autoencoder - ensure integer parameters
    decModel$pretrain(
      x = x_dec,
      epochs = as.integer(params$pretrain_epochs),
      batchSize = as.integer(params$batch_size),
      optimizer = optimizer_adam(learning_rate = params$pretrain_lr)
    )
    
    # Compile DEC Model
    decModel$compile(
      optimizer = optimizer_sgd(learning_rate = params$finetune_lr, momentum = params$momentum)
    )
    
    # Fit DEC - ensure integer parameters
    dec_labels <- decModel$fit(
      x = x_dec,
      maxNumberOfIterations = as.integer(params$max_iter),
      batchSize = as.integer(params$batch_size),
      tolerance = params$tolerance,
      updateInterval = as.integer(params$update_interval)
    )
    
    keras::k_clear_session()
    return(dec_labels)
    
  }, error = function(e) {
    cat("  DEC error:", e$message, "\n")
    keras::k_clear_session()
    return(NA)
  })
}

# ==============================================================================
# 3. HYPERPARAMETER GRID DEFINITION
# ==============================================================================

#' Generate hyperparameter grid for tuning
#' @param tuning_mode "quick", "standard", or "exhaustive"
#' @return Data frame with all hyperparameter combinations
generate_hyperparameter_grid <- function(tuning_mode = "standard") {
  
  if (tuning_mode == "quick") {
    # Quick tuning: 12 combinations
    grid <- expand.grid(
      architecture = c("medium", "deep"),
      alpha = c(1.0),
      activation = c("relu"),
      initializer = c("glorot_uniform"),
      pretrain_epochs = c(100L, 200L),
      pretrain_lr = c(0.001),
      finetune_lr = c(0.01),
      momentum = c(0.9),
      batch_size = c(32L, 64L),
      max_iter = c(2000L),
      tolerance = c(1e-3),
      update_interval = c(20L),
      stringsAsFactors = FALSE
    )
    
  } else if (tuning_mode == "standard") {
    # Standard tuning: ~72 combinations
    grid <- expand.grid(
      architecture = c("shallow", "medium", "deep", "wide"),
      alpha = c(0.5, 1.0),
      activation = c("relu", "elu"),
      initializer = c("glorot_uniform"),
      pretrain_epochs = c(100L, 200L),
      pretrain_lr = c(0.001),
      finetune_lr = c(0.01, 0.001),
      momentum = c(0.9),
      batch_size = c(32L, 64L),
      max_iter = c(3000L),
      tolerance = c(1e-3),
      update_interval = c(20L),
      stringsAsFactors = FALSE
    )
    
  } else if (tuning_mode == "exhaustive") {
    # Exhaustive tuning: ~576 combinations (use with caution!)
    grid <- expand.grid(
      architecture = c("shallow", "medium", "deep", "wide", "bottleneck"),
      alpha = c(0.5, 1.0, 2.0),
      activation = c("relu", "elu", "selu", "tanh"),
      initializer = c("glorot_uniform", "he_normal"),
      pretrain_epochs = c(100L, 200L, 300L),
      pretrain_lr = c(0.001, 0.0005),
      finetune_lr = c(0.01, 0.005, 0.001),
      momentum = c(0.9, 0.95),
      batch_size = c(32L, 64L, 128L),
      max_iter = c(3000L, 5000L),
      tolerance = c(1e-3, 1e-4),
      update_interval = c(10L, 20L, 50L),
      stringsAsFactors = FALSE
    )
  }
  
  return(grid)
}

# ==============================================================================
# 4. CROSS-VALIDATION FUNCTION
# ==============================================================================

#' Perform k-fold cross-validation for DEC
#' @param x_dec Data matrix
#' @param true_labels True cluster labels
#' @param num_clust Number of clusters
#' @param params Hyperparameter configuration
#' @param n_folds Number of CV folds (default 3)
#' @return Mean ARI across folds
cv_dec <- function(x_dec, true_labels, num_clust, params, n_folds = 3L) {
  
  n <- nrow(x_dec)
  fold_size <- floor(n / n_folds)
  fold_indices <- sample(rep(1:n_folds, length.out = n))
  
  ari_scores <- numeric(n_folds)
  
  for (fold in 1:n_folds) {
    # Use all data but evaluate on different folds
    # (DEC is unsupervised, so we train on all and evaluate per fold)
    train_idx <- which(fold_indices != fold)
    test_idx <- which(fold_indices == fold)
    
    # Train on full data (unsupervised)
    labels <- run_dec_with_params(x_dec, num_clust, params)
    
    if (all(is.na(labels))) {
      ari_scores[fold] <- NA
    } else {
      # Evaluate on test fold
      ari_scores[fold] <- aricode::ARI(labels[test_idx], true_labels[test_idx])
    }
  }
  
  return(mean(ari_scores, na.rm = TRUE))
}

# ==============================================================================
# 5. HYPERPARAMETER TUNING FUNCTION
# ==============================================================================

#' Tune DEC hyperparameters
#' @param x_dec Data matrix
#' @param true_labels True labels for evaluation
#' @param num_clust Number of clusters
#' @param tuning_mode "quick", "standard", or "exhaustive"
#' @param use_cv Use cross-validation (TRUE) or single run (FALSE)
#' @param n_runs Number of runs per configuration (for averaging)
#' @param verbose Print progress
#' @return List with best params, all results, and best ARI
tune_dec <- function(x_dec, true_labels, num_clust, 
                     tuning_mode = "standard",
                     use_cv = FALSE,
                     n_runs = 3L,
                     verbose = TRUE) {
  
  num_clust <- as.integer(num_clust)
  n_runs <- as.integer(n_runs)
  
  # Generate grid
  grid <- generate_hyperparameter_grid(tuning_mode)
  n_configs <- nrow(grid)
  
  if (verbose) {
    cat("Starting hyperparameter tuning with", n_configs, "configurations\n")
    cat("Tuning mode:", tuning_mode, "\n")
    cat("Runs per config:", n_runs, "\n\n")
  }
  
  # Store results
  results <- data.frame(
    config_id = 1:n_configs,
    mean_ari = rep(NA, n_configs),
    sd_ari = rep(NA, n_configs),
    stringsAsFactors = FALSE
  )
  results <- cbind(results, grid)
  
  best_ari <- -Inf
  best_params <- NULL
  
  # Iterate over configurations
  for (i in 1:n_configs) {
    params <- as.list(grid[i, ])
    
    if (verbose) {
      cat(sprintf("[%d/%d] Testing: arch=%s, alpha=%.1f, act=%s, pretrain_ep=%d, ft_lr=%.4f, bs=%d\n",
                  i, n_configs, params$architecture, params$alpha, params$activation,
                  params$pretrain_epochs, params$finetune_lr, params$batch_size))
    }
    
    # Multiple runs for stability
    ari_runs <- numeric(n_runs)
    
    for (run in 1:n_runs) {
      if (use_cv) {
        ari_runs[run] <- cv_dec(x_dec, true_labels, num_clust, params, n_folds = 3L)
      } else {
        labels <- run_dec_with_params(x_dec, num_clust, params)
        if (all(is.na(labels))) {
          ari_runs[run] <- NA
        } else {
          ari_runs[run] <- aricode::ARI(labels, true_labels)
        }
      }
    }
    
    mean_ari <- mean(ari_runs, na.rm = TRUE)
    sd_ari <- sd(ari_runs, na.rm = TRUE)
    
    results$mean_ari[i] <- mean_ari
    results$sd_ari[i] <- sd_ari
    
    if (!is.na(mean_ari) && mean_ari > best_ari) {
      best_ari <- mean_ari
      best_params <- params
      if (verbose) cat(sprintf("  -> New best ARI: %.4f\n", best_ari))
    }
    
    if (verbose) cat(sprintf("  ARI: %.4f (sd=%.4f)\n\n", mean_ari, sd_ari))
  }
  
  # Sort results by mean_ari
  results <- results[order(-results$mean_ari), ]
  
  return(list(
    best_params = best_params,
    best_ari = best_ari,
    all_results = results
  ))
}

# ==============================================================================
# 6. RANDOM SEARCH ALTERNATIVE
# ==============================================================================

#' Random search for DEC hyperparameters
#' @param x_dec Data matrix
#' @param true_labels True labels
#' @param num_clust Number of clusters
#' @param n_iter Number of random configurations to try
#' @param n_runs Runs per configuration
#' @param verbose Print progress
#' @return List with best params and results
random_search_dec <- function(x_dec, true_labels, num_clust, 
                              n_iter = 20L, n_runs = 2L, verbose = TRUE) {
  
  num_clust <- as.integer(num_clust)
  n_iter <- as.integer(n_iter)
  n_runs <- as.integer(n_runs)
  
  if (verbose) cat("Starting random search with", n_iter, "iterations\n\n")
  
  results <- list()
  best_ari <- -Inf
  best_params <- NULL
  
  for (i in 1:n_iter) {
    # Sample random hyperparameters
    params <- list(
      architecture = sample(c("shallow", "medium", "deep", "wide", "bottleneck"), 1),
      alpha = runif(1, 0.5, 2.0),
      activation = sample(c("relu", "elu", "selu"), 1),
      initializer = sample(c("glorot_uniform", "he_normal"), 1),
      pretrain_epochs = as.integer(sample(c(100L, 150L, 200L, 250L), 1)),
      pretrain_lr = 10^runif(1, -4, -2.5),
      finetune_lr = 10^runif(1, -3, -1.5),
      momentum = runif(1, 0.85, 0.99),
      batch_size = as.integer(sample(c(16L, 32L, 64L, 128L), 1)),
      max_iter = as.integer(sample(c(2000L, 3000L, 4000L), 1)),
      tolerance = 10^runif(1, -4, -2),
      update_interval = as.integer(sample(c(10L, 20L, 30L, 50L), 1))
    )
    
    if (verbose) {
      cat(sprintf("[%d/%d] arch=%s, alpha=%.2f, act=%s, pt_lr=%.5f, ft_lr=%.4f\n",
                  i, n_iter, params$architecture, params$alpha, params$activation,
                  params$pretrain_lr, params$finetune_lr))
    }
    
    # Multiple runs
    ari_runs <- numeric(n_runs)
    for (run in 1:n_runs) {
      labels <- run_dec_with_params(x_dec, num_clust, params)
      if (all(is.na(labels))) {
        ari_runs[run] <- NA
      } else {
        ari_runs[run] <- aricode::ARI(labels, true_labels)
      }
    }
    
    mean_ari <- mean(ari_runs, na.rm = TRUE)
    
    results[[i]] <- c(params, mean_ari = mean_ari)
    
    if (!is.na(mean_ari) && mean_ari > best_ari) {
      best_ari <- mean_ari
      best_params <- params
      if (verbose) cat(sprintf("  -> New best ARI: %.4f\n", best_ari))
    }
    
    if (verbose) cat(sprintf("  ARI: %.4f\n\n", mean_ari))
  }
  
  return(list(
    best_params = best_params,
    best_ari = best_ari,
    all_results = results
  ))
}

# ==============================================================================
# 7. MAIN EXECUTION LOOP WITH TUNING
# ==============================================================================

here::here()

files <- list.files("data/")
n <- length(files)
file_names <- sub("\\.RDS$", "", files)

# Results storage
results <- data.frame(
  'Dataset' = file_names,
  'DEC_default' = rep(NA, n),
  'DEC_tuned' = rep(NA, n),
  stringsAsFactors = FALSE
)

# Store best params for each dataset
best_params_list <- list()

# Tuning settings
TUNING_MODE <- "standard"  # Options: "quick", "standard", "exhaustive"
USE_RANDOM_SEARCH <- FALSE  # TRUE for random search, FALSE for grid search
N_RANDOM_ITER <- 15L  # Only used if USE_RANDOM_SEARCH = TRUE

for (i in 1:n) {
  
  cat("\n========================================\n")
  cat("Processing dataset:", file_names[i], "\n")
  cat("========================================\n\n")
  
  # --- Loading & Basic Prep ---
  data <- readRDS(paste0("data/", files[i]))
  class <- data[, ncol(data)]
  data <- data[, -ncol(data)]
  
  cat_cols <- as.numeric(which(sapply(data, is.factor)))
  cont_cols <- setdiff(c(1:ncol(data)), cat_cols)
  num_clust <- as.integer(length(unique(class)))  
  
  # Scale continuous variables
  if (length(cont_cols) > 0) {
    if (length(cont_cols) == 1) {
      data[, cont_cols] <- as.numeric(scale(data[, cont_cols]))
    } else {
      data[, cont_cols] <- scale(data[, cont_cols])
    }
  }
  
  # Preprocess for DEC
  x_dec <- model.matrix(~ . - 1, data = data)
  
  # --- Method 1: DEC with Default Parameters ---
  cat("Running DEC with default parameters...\n")
  default_params <- list(
    architecture = "medium",
    alpha = 1.0,
    activation = "relu",
    initializer = "glorot_uniform",
    pretrain_epochs = 100L,
    pretrain_lr = 0.001,
    finetune_lr = 0.01,
    momentum = 0.9,
    batch_size = 32L,
    max_iter = 2000L,
    tolerance = 1e-3,
    update_interval = 20L
  )
  
  dec_labels_default <- run_dec_with_params(x_dec, num_clust, default_params)
  if (!all(is.na(dec_labels_default))) {
    results$DEC_default[i] <- aricode::ARI(dec_labels_default, class)
    cat("Default DEC ARI:", results$DEC_default[i], "\n\n")
  }
  
  # --- Method 2: DEC with Hyperparameter Tuning ---
  cat("Starting hyperparameter tuning...\n")
  
  if (USE_RANDOM_SEARCH) {
    tuning_result <- random_search_dec(
      x_dec = x_dec,
      true_labels = class,
      num_clust = num_clust,
      n_iter = N_RANDOM_ITER,
      n_runs = 2L,
      verbose = TRUE
    )
  } else {
    tuning_result <- tune_dec(
      x_dec = x_dec,
      true_labels = class,
      num_clust = num_clust,
      tuning_mode = TUNING_MODE,
      use_cv = FALSE,
      n_runs = 2L,
      verbose = TRUE
    )
  }
  
  results$DEC_tuned[i] <- tuning_result$best_ari
  best_params_list[[file_names[i]]] <- tuning_result$best_params
  
  cat("\n--- Results for", file_names[i], "---\n")
  cat("Default DEC ARI:", results$DEC_default[i], "\n")
  cat("Tuned DEC ARI:", results$DEC_tuned[i], "\n")
  cat("Improvement:", round(results$DEC_tuned[i] - results$DEC_default[i], 4), "\n")
  cat("Best parameters:\n")
  print(tuning_result$best_params)
  
  # --- Save Intermediate Results ---
  save(results, best_params_list, file = 'results_with_tuning.RData')
  cat('\nDataset', file_names[i], 'complete.\n')
  cat("========================================\n")
}

# ==============================================================================
# 8. FINAL SUMMARY
# ==============================================================================

cat("\n\n========================================\n")
cat("FINAL SUMMARY\n")
cat("========================================\n\n")

print(results)

cat("\nMean improvement from tuning:", 
    round(mean(results$DEC_tuned - results$DEC_default, na.rm = TRUE), 4), "\n")

cat("\nBest parameters per dataset:\n")
for (ds in names(best_params_list)) {
  cat("\n", ds, ":\n")
  cat("  Architecture:", best_params_list[[ds]]$architecture, "\n")
  cat("  Alpha:", best_params_list[[ds]]$alpha, "\n")
  cat("  Activation:", best_params_list[[ds]]$activation, "\n")
  cat("  Pretrain LR:", best_params_list[[ds]]$pretrain_lr, "\n")
  cat("  Finetune LR:", best_params_list[[ds]]$finetune_lr, "\n")
  cat("  Batch size:", best_params_list[[ds]]$batch_size, "\n")
}

# Save final results
save(results, best_params_list, file = 'results_with_tuning.RData')
write.csv(results, 'results_with_tuning.csv', row.names = FALSE)

cat("\nResults saved to 'results_with_tuning.RData' and 'results_with_tuning.csv'\n")
