library(IBclust)
library(aricode)
library(parallel)

# Load data
X <- readRDS('data/Adult.RDS')
class <- X[, ncol(X)]
X <- X[, -ncol(X)]
load('data/prep_list.RData')

# Previously computed bandwidths for Adult data set
s <- rep(5.9, 6)
lambda <- c(0.8076923,
            0.9765021,
            0.8076923,
            0.9009901,
            0.7777778,
            0.7368421,
            0.4117647,
            0.9655172)

# m values to test
m_values <- c(50, 100, NA, 250, 500, 1000, 2500)

# Number of replications and cores
n_reps <- 50
n_cores <- 7

# Function to run single replication
run_single_rep <- function(rep_id, X, class, s, lambda, m) {
  start_time <- Sys.time()
  
  if (is.na(m)) {
    DIB_test <- DIBmix(X, ncl = 2, s = s, lambda = lambda, nystrom = TRUE, n_landmarks = NULL)
  } else {
    DIB_test <- DIBmix(X, ncl = 2, s = s, lambda = lambda, nystrom = TRUE, n_landmarks = m)
  }
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  ari <- ARI(DIB_test$Cluster, class)
  
  return(list(rep = rep_id, runtime = runtime, ARI = ari))
}

# Run benchmarks for each k
for (m in m_values) {
  m_label <- ifelse(is.na(m), "NULL", as.character(m))
  cat("Running m =", m_label, "\n")
  
  # Run in parallel using mclapply (forking, no cluster setup needed)
  results <- mclapply(1:n_reps, function(i) {
    run_single_rep(i, X, class, s, lambda, m)
  }, mc.cores = n_cores)
  
  # Convert to data frame
  results_df <- data.frame(
    rep = sapply(results, `[[`, "rep"),
    runtime = sapply(results, `[[`, "runtime"),
    ARI = sapply(results, `[[`, "ARI")
  )
  
  # Save results
  filename <- paste0("res/runtimes_ari_", m_label, ".RDS")
  saveRDS(results_df, filename)
}

cat("All benchmarks complete.\n")