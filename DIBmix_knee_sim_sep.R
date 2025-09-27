library(MixSim)
library(IBclust)
library(kneedle)

# Well-separated spherical clusters
clusters <- 4
overlap <- 0.01
rows <- 500
columns <- 8
pi_val <- 1.0
cat_ratio <-  0.5
nreps <- 100
res_df_knee_sep <- data.frame(mi2 = numeric(),
                              mi3 = numeric(),
                              mi4 = numeric(),
                              mi5 = numeric(),
                              mi6 = numeric(),
                              mi7 = numeric(),
                              mi8 = numeric(),
                              knee = numeric(),
                              stringsAsFactors = FALSE)

for (i in 1:nreps){
  set.seed(i)
  # Construct artificial data set - sphericity set to FALSE by default
  mixsimaux <- MixSim(BarOmega = overlap, PiLow=pi_val,
                      K = clusters, p = columns, resN = 1000000, sph = TRUE)
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
  outdibmix2 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 2,
                                randinit = NULL,
                                lambda = -1,
                                s = -1,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  res_df_knee_sep[i, 1] <- outdibmix2$MutualInfo
  outdibmix3 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 3,
                                randinit = NULL,
                                lambda = outdibmix2$lambda,
                                s = outdibmix2$s,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  res_df_knee_sep[i, 2] <- outdibmix3$MutualInfo
  outdibmix4 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 4,
                                randinit = NULL,
                                lambda = outdibmix2$lambda,
                                s = outdibmix2$s,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  res_df_knee_sep[i, 3] <- outdibmix4$MutualInfo
  outdibmix5 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 5,
                                randinit = NULL,
                                lambda = outdibmix2$lambda,
                                s = outdibmix2$s,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  res_df_knee_sep[i, 4] <- outdibmix5$MutualInfo
  outdibmix6 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 6,
                                randinit = NULL,
                                lambda = outdibmix2$lambda,
                                s = outdibmix2$s,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  res_df_knee_sep[i, 5] <- outdibmix6$MutualInfo
  outdibmix7 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 7,
                                randinit = NULL,
                                lambda = outdibmix2$lambda,
                                s = outdibmix2$s,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  res_df_knee_sep[i, 6] <- outdibmix7$MutualInfo
  outdibmix8 <- IBclust::DIBmix(X = mixdt1df,
                                ncl = 8,
                                randinit = NULL,
                                lambda = outdibmix2$lambda,
                                s = outdibmix2$s,
                                maxiter = 100,
                                nstart = 100,
                                verbose = FALSE)
  res_df_knee_sep[i, 7] <- outdibmix8$MutualInfo
  infos <- c(outdibmix2$MutualInfo, outdibmix3$MutualInfo,
             outdibmix4$MutualInfo, outdibmix5$MutualInfo,
             outdibmix6$MutualInfo, outdibmix7$MutualInfo,
             outdibmix8$MutualInfo)
  knee <- kneedle::kneedle(x = c(2:8), y = infos,
                           decreasing = FALSE)
  res_df_knee_sep[i, 8] <- knee[1]
  saveRDS(res_df_knee_sep, file = 'res_df_knee_sep.RDS')
}