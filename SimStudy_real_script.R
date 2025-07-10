require(here)
require(IBclust)
require(kamila)
require(clustMixType)
require(FactoMineR)
require(cluster)
require(aricode)

# Set working directory to path of script
here::here()

# Go into data directory and read files, perform clustering and save
# Store results in data frame
files <- list.files("data/")
n <- length(files)
file_names <- sub("\\.RDS$", "", files)
results <- data.frame('Dataset' = file_names,
                      'DIBmix' = rep(NA, n),
                      'KAMILA' = rep(NA, n),
                      'K-Prototypes' = rep(NA, n),
                      'FAMD/K-Means' = rep(NA, n),
                      'Gower/PAM' = rep(NA, n),
                      stringsAsFactors=FALSE)

for (i in 1:n){
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
  # Run methods and store ARI values
  # DIBmix
  DIBmix_res <- IBclust::DIBmix(X = data,
                                ncl = num_clust,
                                catcols = cat_cols,
                                contcols = cont_cols,
                                randinit = NULL,
                                s = -1,
                                lambda = -1,
                                scale = FALSE,
                                maxiter = 100,
                                nstart = 100)
  results[i, 2] <- aricode::ARI(DIBmix_res$Cluster, class)
  cat('DIBmix done on', file_names[i], 'dataset.\n')
  save(results, file = 'results.RData')
  # KAMILA
  kamila_res <- kamila::kamila(conVar = as.data.frame(data[, cont_cols]), 
                               catFactor = as.data.frame(data[, cat_cols]),
                               numClust = num_clust,
                               numInit = 100)
  results[i, 3] <- aricode::ARI(kamila_res$finalMemb, class)
  cat('KAMILA done on', file_names[i], 'dataset.\n')
  save(results, file = 'results.RData')
  # K-Prototypes
  kproto_res <- clustMixType::kproto(data, num_clust, nstart = 100, verbose = FALSE)
  results[i, 4] <- aricode::ARI(kproto_res$cluster, class)
  cat('K-Prototypes done on', file_names[i], 'dataset.\n')
  save(results, file = 'results.RData')
  # FAMD & K-Means
  outpcamix <- FactoMineR::FAMD(data, ncp = num_clust-1, graph = FALSE)
  famd_res <- kmeans(outpcamix$ind$coord, nstart=100, num_clust, iter.max=1000, algorithm = "MacQueen")
  results[i, 5] <- aricode::ARI(famd_res$cluster, class)
  cat('FAMD & K-Means done on', file_names[i], 'dataset.\n')
  save(results, file = 'results.RData')
  # PAM with Gower's dissimilarity
  gower_dist <- cluster::daisy(data, metric = "gower")
  pam_res <- cluster::pam(gower_dist, diss = TRUE, k = num_clust, do.swap = FALSE, cluster.only = TRUE, nstart=100)
  results[i, 6] <- aricode::ARI(pam_res, class)
  cat('PAM/Gower done on', file_names[i], 'dataset.\n')
  save(results, file = 'results.RData')
  cat('Dataset', file_names[i], 'complete.\n')
}
