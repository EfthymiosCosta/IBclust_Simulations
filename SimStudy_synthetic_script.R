require(MixSim)
require(cluster)
require(clustrd)
require(clustMixType)
require(mclust)
require(kmed)
require(FactoMineR)
require(fpc)
require(kamila)
require(rdist)
require(aricode)
require(DIBclust)

# categorized numerical variable function
intv <- function(vec, class) {
  nbase <- (1:(class-1))/class
  nq <- numeric(length(nbase))
  for (i in 1:length(nq)) {
    nq[i] <- quantile(vec, nbase[i])
  }
  res <- c(min(vec), nq, max(vec)) 
  res[1] <- res[1]-1
  for (i in 2:length(res)){
    if (res[i-1]==res[i]){
      res[i] <- res[i]+2e-15
    }
  }
  return(res)
}

# Full factorial simulation study - NON-SPHERICAL CLUSTERS

# Empty data frame to store simulation results
fullfactorial_res <- data.frame(seed=numeric(),
                                nClust=numeric(),
                                overlap=numeric(),
                                nrows=numeric(),
                                ncols=numeric(),
                                catratio=numeric(),
                                pi=numeric(),
                                method=character(),
                                ARI=numeric(),
                                AMI=numeric(),
                                stringsAsFactors=FALSE)

cluster_sizes <- c(3, 5, 8)
overlap_levels <- c(0.01, seq(0.05, 0.15, by=0.05))
num_rows <- c(100, 600)
num_cols <- c(8, 16)
pi_values <- c(0.01, 1.0)
cat_ratios <-  c(0.2, 0.5, 0.8)

nreps <- 1
tot <- nreps*length(cluster_sizes)*length(overlap_levels)*length(num_rows)*length(num_cols)*length(pi_values)*length(cat_ratios)
run <- 1

for (l in 0:(nreps-1)){
  # Random seed
  set.seed(1234+l)
  # Number of clusters
  for (clusters in cluster_sizes){
    # Overlap level
    for (overlap in overlap_levels){
      # Number of rows
      for (rows in num_rows){
        # Number of variables/columns
        for (columns in num_cols){
          # Balanced/unbalanced design
          for (pi_val in pi_values){
            # Ratio of categorical variables
            for (cat_ratio in cat_ratios){
              # Construct artificial data set - sphericity set to FALSE by default
              mixsimaux <- MixSim(BarOmega = overlap, PiLow=pi_val,
                                  K = clusters, p = columns, resN = 1000000)
              mixdtaux <- simdataset(n = rows, Pi = mixsimaux$Pi, Mu = mixsimaux$Mu, S = mixsimaux$S)
              # Discretise first half attributes using 4 levels
              for (k in 1:(round(columns*cat_ratio))){
                mixdtaux$X[,k] <- 
                  as.factor(cut(mixdtaux$X[,k], intv(mixdtaux$X[,k], 4), labels = (1:4)))
              }
              mixdt1df <- as.data.frame(mixdtaux$X)
              for (k in 1:(round(columns*cat_ratio))){
                mixdt1df[,k] <- as.factor(mixdt1df[,k])
              }
              # this temporarily fixes a nasty bug of cluspcamix
              colnames(mixdt1df) <- sprintf("a%d", 1:columns)
              # Calculate Gower's dissimilarities and Ahmad distances
              gower_dist <- daisy(mixdt1df, metric = "gower")
              # Specify continuous & categorical attributes
              conDf <- data.frame(scale(mixdt1df[,((round(columns*cat_ratio)+1):columns)]))
              catDf <- dummyCodeFactorDf(data.frame(mixdt1df[,1:(round(columns*cat_ratio))]))
              # Create dummy-coded data frame
              mixdt1df_dummy <- as.data.frame(cat2bin(mixdt1df, categorical=1:(round(columns*cat_ratio)))$data)
              # above command makes all variables characters - convert to numeric
              for (i in 1:ncol(mixdt1df_dummy)){
                mixdt1df_dummy[,i] <- as.numeric(mixdt1df_dummy[,i])
              } 
              # Calculate distance factor for categorical vars standardisation
              for (i in 1:(round(columns*cat_ratio))){
                mixdt1df_dummy[,c((4*(i-1)+1):(4*i))] <- mixdt1df_dummy[,c((4*(i-1)+1):(4*i))]*distancefactor(cat=4, n=rows, type="categorical")
              }
              
              # Continuous columns variances
              columnvars <- lapply(mixdt1df_dummy[c((4*round(columns*cat_ratio)+1):ncol(mixdt1df_dummy))], FUN=var)
              inx <- 1
              # Standardise continuous vars by standard deviation
              for (i in c((4*(round(columns*cat_ratio))+1):ncol(mixdt1df_dummy))){
                mixdt1df_dummy[,i] <- mixdt1df_dummy[,i]/sqrt(columnvars[[inx]])
                inx <- inx + 1
              }
              
              # Create dissimilarity/distance matrix
              dissmat <- pdist(mixdt1df_dummy, metric="euclidean")
              # Extract principal components
              outpcamix <- FAMD(mixdt1df, ncp = clusters-1, graph=FALSE)
              # PAM with Gower's
              pam_fit <- pam(gower_dist, diss = TRUE, k = clusters, do.swap = FALSE, cluster.only = TRUE,
                             nstart=100)
              fullfactorial_res <- rbind(fullfactorial_res, data.frame(seed=l+1, nClust=clusters, overlap=overlap,
                                                               nrows=rows, ncols=columns, catratio=cat_ratio, pi=pi_val,
                                                               method='PAM/Gower', ARI=ARI(pam_fit, mixdtaux$id),
                                                               AMI=AMI(pam_fit, mixdtaux$id)))
              save(fullfactorial_res, file='fullfactorial_res.RData')
              cat('PAM/Gower done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                  'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
              # K-prototypes
              outk <- kproto(mixdt1df, clusters, nstart = 100,verbose = FALSE)
              fullfactorial_res <- rbind(fullfactorial_res, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                               ncols=columns, catratio=cat_ratio, pi=pi_val, method='K-Prot',
                                                               ARI=ARI(outk$cluster, mixdtaux$id),
                                                               AMI=AMI(outk$cluster, mixdtaux$id)))
              save(fullfactorial_res, file='fullfactorial_res.RData')
              cat('K-prototypes done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                  'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
              # KAMILA
              kamilafit <- kamila_fit <- kamila(conVar = conDf,
                                                catFactor = data.frame(mixdt1df[,1:(round(columns*cat_ratio))]),
                                                numClust = clusters,
                                                numInit = 100,
                                                maxIter = 10000)
              fullfactorial_res <- rbind(fullfactorial_res, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                               ncols=columns, catratio=cat_ratio, pi=pi_val, method='KAMILA',
                                                               ARI=ARI(kamilafit$finalMemb, mixdtaux$id),
                                                               AMI=AMI(kamilafit$finalMemb, mixdtaux$id)))
              save(fullfactorial_res, file='fullfactorial_res.RData')
              cat('KAMILA done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                  'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
              # FAMD + K-means
              outkm <- kmeans(outpcamix$ind$coord, nstart=100, clusters, iter.max=1000, algorithm = "MacQueen")
              fullfactorial_res <- rbind(fullfactorial_res, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                               ncols=columns, catratio=cat_ratio, pi=pi_val, method='FAMD',
                                                               ARI=ARI(outkm$cluster, mixdtaux$id),
                                                               AMI=AMI(outkm$cluster, mixdtaux$id)))
              save(fullfactorial_res, file='fullfactorial_res.RData')
              cat('FAMD done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                  'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
              # DIBmix
              outdibmix <- DIBclust::DIBmix(X = mixdt1df,
                                            ncl = clusters,
                                            catcols = c(1:(round(columns*cat_ratio))),
                                            contcols = c((round(columns*cat_ratio)+1):columns),
                                            randinit = NULL,
                                            lambda = -1,
                                            s = -1,
                                            maxiter = 100,
                                            nstart = 100,
                                            select_features = FALSE)
              fullfactorial_res <- rbind(fullfactorial_res, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                               ncols=columns, catratio=cat_ratio, pi=pi_val, method='DIBmix',
                                                               ARI=ARI(outdibmix[[1]], mixdtaux$id),
                                                               AMI=AMI(outdibmix[[1]], mixdtaux$id)))
              save(fullfactorial_res, file='fullfactorial_res.RData')
              cat('Simulations for data set',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                  'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
              cat('*****  Run',run,'/',tot,' *****','\n')
              run <- run + 1
            }
          }
        }
      }
    } 
  }
}
