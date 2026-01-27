require(here)
require(IBclust)
require(kamila)
require(clustMixType)
require(FactoMineR)
require(cluster)
require(parallel)

# Set working directory to path of script
base_dir <- here::here()
setwd(base_dir)

# Go into data directory and read files (only .RDS files, not directories)
all_files <- list.files("data/", full.names = FALSE)
files <- all_files[grepl("\\.RDS$", all_files, ignore.case = TRUE)]
n <- length(files)
file_names <- sub("\\.RDS$", "", files, ignore.case = TRUE)

# Function to process a single dataset
process_dataset <- function(i, files, file_names, base_dir) {
  
  # Set working directory on worker
  setwd(base_dir)
  
  data <- readRDS(paste0("data/", files[i]))
  # Extract class variable (last variable)
  class <- data[, ncol(data)]
  data <- data[, -ncol(data)]
  # Find categorical & continuous variables
  cat_cols <- as.numeric(which(sapply(data, is.factor)))
  cont_cols <- setdiff(c(1:ncol(data)), cat_cols)
  # Number of clusters
  num_clust <- length(unique(class))
  # Scale continuous variables
  if (length(cont_cols) == 1) {
    data[, cont_cols] <- as.numeric(scale(data[, cont_cols]))
  } else {
    data[, cont_cols] <- scale(data[, cont_cols])
  }
  
  # Initialize results for this dataset
  dataset_runtimes <- data.frame(
    "Dataset" = file_names[i],
    "DIBmix_search" = NA,
    "DIBmix_fixed" = NA,
    "KAMILA" = NA,
    "K-Prototypes" = NA,
    "FAMD/K-Means" = NA,
    "Gower/PAM" = NA,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  nystrom <- nrow(data) > 1000
  # DIBmix with hyperparameter search
  start_time <- Sys.time()
  DIBmix_res <- IBclust::DIBmix(X = data,
                                ncl = num_clust,
                                randinit = NULL,
                                s = -1,
                                lambda = -1,
                                scale = FALSE,
                                maxiter = 100,
                                nstart = 100,
                                nystrom = nystrom)
  end_time <- Sys.time()
  dataset_runtimes[1, 2] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # DIBmix with fixed hyperparameters from search
  start_time <- Sys.time()
  DIBmix_fixed <- IBclust::DIBmix(X = data,
                                  ncl = num_clust,
                                  randinit = NULL,
                                  s = DIBmix_res$s,
                                  lambda = DIBmix_res$lambda,
                                  scale = FALSE,
                                  maxiter = 100,
                                  nstart = 100,
                                  nystrom = nystrom)
  end_time <- Sys.time()
  dataset_runtimes[1, 3] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # KAMILA
  start_time <- Sys.time()
  kamila_res <- kamila::kamila(conVar = as.data.frame(data[, cont_cols]), 
                               catFactor = as.data.frame(data[, cat_cols]),
                               numClust = num_clust,
                               numInit = 100)
  end_time <- Sys.time()
  dataset_runtimes[1, 4] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # K-Prototypes
  start_time <- Sys.time()
  kproto_res <- clustMixType::kproto(data, num_clust, nstart = 100, verbose = FALSE)
  end_time <- Sys.time()
  dataset_runtimes[1, 5] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # FAMD & K-Means
  start_time <- Sys.time()
  outpcamix <- FactoMineR::FAMD(data, ncp = num_clust-1, graph = FALSE)
  famd_res <- kmeans(outpcamix$ind$coord, nstart=100, num_clust, iter.max=1000, algorithm = "MacQueen")
  end_time <- Sys.time()
  dataset_runtimes[1, 6] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # PAM with Gower's dissimilarity
  start_time <- Sys.time()
  gower_dist <- cluster::daisy(data, metric = "gower")
  pam_res <- cluster::pam(gower_dist, diss = TRUE, k = num_clust, do.swap = FALSE, cluster.only = TRUE, nstart=100)
  end_time <- Sys.time()
  dataset_runtimes[1, 7] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Save individual result to disk
  saveRDS(dataset_runtimes, file = paste0("runtime_temp/result_", i, ".RDS"))
  
  return(i)  # Just return the index, not the full data
}

# Run in parallel using mclapply
cat("Starting parallel processing on", n, "datasets...\n")
completed <- mclapply(1:n, 
                      process_dataset, 
                      files = files, 
                      file_names = file_names,
                      base_dir = base_dir,
                      mc.cores = 10)

# Read all individual results from disk
cat("Collecting results...\n")
results_list <- lapply(1:n, function(i) {
  readRDS(paste0("runtime_temp/result_", i, ".RDS"))
})

# Combine results into data frame
runtimes <- do.call(rbind, results_list)

# Clean up temporary files
unlink("runtime_temp", recursive = TRUE)

saveRDS(runtimes, file = 'runtimes.RDS')
cat("Done! Results saved to runtimes.RDS\n")
print(runtimes)
