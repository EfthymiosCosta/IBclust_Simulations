library(MixSim)
library(IBclust)
library(parallel)

# Set up parameters
clusters <- 4
overlap <- 0.01
rows <- 500
columns <- 8
pi_val <- 1.0
cat_ratio <- 0.5
nreps <- 100
B <- 100 
C_max <- 8
C_min <- 2  

# Detect number of cores
n_cores <- 100  # One core per reference dataset

# Initialize results dataframe
res_df_gap_sep <- data.frame(
  rep = integer(),
  C = integer(),
  MI_observed = numeric(),
  MI_ref_mean = numeric(),
  MI_ref_sd = numeric(),
  Gap = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nreps) {
  set.seed(i)
  cat(sprintf("Processing replicate %d/%d\n", i, nreps))
  
  # Generate artificial dataset
  mixsimaux <- MixSim(BarOmega = overlap, PiLow = pi_val,
                      K = clusters, p = columns, resN = 1000000, sph = TRUE)
  mixdtaux <- simdataset(n = rows, Pi = mixsimaux$Pi, 
                         Mu = mixsimaux$Mu, S = mixsimaux$S)
  
  # Discretize first half of attributes using 4 levels
  for (k in 1:(round(columns * cat_ratio))) {
    mixdtaux$X[, k] <- 
      as.factor(cut(mixdtaux$X[, k], intv(mixdtaux$X[, k], 4), labels = (1:4)))
  }
  
  mixdt1df <- as.data.frame(mixdtaux$X)
  for (k in 1:(round(columns * cat_ratio))) {
    mixdt1df[, k] <- as.factor(mixdt1df[, k])
  }
  colnames(mixdt1df) <- sprintf("a%d", 1:columns)
  
  # Get bandwidth parameters from C=2 run
  outdibmix2 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 2,
                                randinit = NULL,
                                lambda = -1,
                                s = -1,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  
  lambda_fixed <- outdibmix2$lambda
  s_fixed <- outdibmix2$s
  
  # Loop over different values of C
  for (C_val in C_min:C_max) {
    cat(sprintf("  Processing C = %d\n", C_val))
    
    # Compute MI on observed data
    out_observed <- IBclust::DIBmix(X = mixdt1df,
                                    ncl = C_val,
                                    randinit = NULL,
                                    lambda = lambda_fixed,
                                    s = s_fixed,
                                    maxiter = 100,
                                    nstart = 100,
                                    verbose = FALSE)
    MI_observed <- out_observed$MutualInfo
    
    # Function to compute MI on one reference dataset
    compute_ref_MI <- function(b) {
      # Independently permute each column
      ref_data <- mixdt1df
      for (j in 1:ncol(mixdt1df)) {
        ref_data[, j] <- sample(mixdt1df[, j], size = nrow(mixdt1df), replace = FALSE)
      }
      
      # Run DIBmix on reference dataset
      out_ref <- IBclust::DIBmix(X = ref_data,
                                 ncl = C_val,
                                 randinit = NULL,
                                 lambda = lambda_fixed,
                                 s = s_fixed,
                                 maxiter = 100,
                                 nstart = 100,
                                 verbose = FALSE)
      return(out_ref$MutualInfo)
    }
    
    # Parallelize over B reference datasets
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, library(IBclust))
    clusterExport(cl, varlist = c("mixdt1df", "C_val", "lambda_fixed", "s_fixed"),
                  envir = environment())
    
    MI_ref_values <- parSapply(cl, 1:B, compute_ref_MI)
    stopCluster(cl)
    
    # Compute Gap statistic and related quantities
    MI_ref_mean <- mean(MI_ref_values)
    MI_ref_sd <- sd(MI_ref_values)
    Gap <- MI_observed - MI_ref_mean
    
    # Store results
    res_df_gap_sep <- rbind(res_df_gap_sep, data.frame(
      rep = i,
      C = C_val,
      MI_observed = MI_observed,
      MI_ref_mean = MI_ref_mean,
      MI_ref_sd = MI_ref_sd,
      Gap = Gap,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save intermediate results after each replicate
  saveRDS(res_df_gap_sep, file = 'res_df_gap_sep.RDS')
}
