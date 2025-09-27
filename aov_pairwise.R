library(afex)
afex_options(es_aov="pes")

load('res/fullfactorial_res.RData')
load('res/fullfactorial_res_sph.RData')

fullfactorial_full <- rbind(fullfactorial_res, fullfactorial_res_sph)
fullfactorial_full$sph <- c(rep(FALSE, nrow(fullfactorial_full)/2),
                            rep(TRUE, nrow(fullfactorial_full)/2))
fullfactorial_full$sph <- as.factor(fullfactorial_full$sph)

#### DIBmix VS. Kamila
fullfactorial_full <- fullfactorial_full[which(fullfactorial_full$method %in% c("DIBmix", "KAMILA")),]

idvec <- c()
for (i in 1:(nrow(fullfactorial_full)/2)){
  idvec <- c(idvec, rep(i, 2))
}

fullfactorial_full$ID <- idvec

# Repeated Measures Anova for ARI
aov_ari_dibmix_kamila <- aov_ez(id = 'ID', dv = 'ARI', fun_aggregate = mean, fullfactorial_full,
                                within ='method',
                                between = c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ari_dibmix_kamila, file = 'anova_res/aov_ari_dibmix_kamila.RDS')

# Repeated Measures Anova for AMI
aov_ami_dibmix_kamila <- aov_ez(id = 'ID', dv = 'AMI', fun_aggregate = mean, fullfactorial_full,
                                within = 'method',
                                between= c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ami_dibmix_kamila, file = 'anova_res/aov_ami_dibmix_kamila.RDS')

#### DIBmix VS. K-Proto
fullfactorial_full <- readRDS('fullfactorial_full.RDS')

fullfactorial_full <- fullfactorial_full[which(fullfactorial_full$method %in% c("DIBmix", "K-Prototypes")),]
fullfactorial_full$ID <- idvec

# Repeated Measures Anova for ARI
aov_ari_dibmix_kproto <- aov_ez(id = 'ID', dv = 'ARI', fun_aggregate = mean, fullfactorial_full,
                                within ='method',
                                between = c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ari_dibmix_kproto, file = 'anova_res/aov_ari_dibmix_kproto.RDS')

# Repeated Measures Anova for AMI
aov_ami_dibmix_kproto <- aov_ez(id = 'ID', dv = 'AMI', fun_aggregate = mean, fullfactorial_full,
                                within = 'method',
                                between= c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ami_dibmix_kproto, file = 'anova_res/aov_ami_dibmix_kproto.RDS')

#### DIBmix VS. FAMD
fullfactorial_full <- readRDS('fullfactorial_full.RDS')

fullfactorial_full <- fullfactorial_full[which(fullfactorial_full$method %in% c("DIBmix", "FAMD")),]
fullfactorial_full$ID <- idvec

# Repeated Measures Anova for ARI
aov_ari_dibmix_famd <- aov_ez(id = 'ID', dv = 'ARI', fun_aggregate = mean, fullfactorial_full,
                              within ='method',
                              between = c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ari_dibmix_famd, file = 'anova_res/aov_ari_dibmix_famd.RDS')

# Repeated Measures Anova for AMI
aov_ami_dibmix_famd <- aov_ez(id = 'ID', dv = 'AMI', fun_aggregate = mean, fullfactorial_full,
                              within = 'method',
                              between= c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ami_dibmix_famd, file = 'anova_res/aov_ami_dibmix_famd.RDS')

  #### DIBmix VS. Gower/PAM
fullfactorial_full <- readRDS('fullfactorial_full.RDS')

fullfactorial_full <- fullfactorial_full[which(fullfactorial_full$method %in% c("DIBmix", "PAM/Gower")),]
fullfactorial_full$ID <- idvec

# Repeated Measures Anova for ARI
aov_ari_dibmix_gower <- aov_ez(id = 'ID', dv = 'ARI', fun_aggregate = mean, fullfactorial_full,
                               within ='method',
                               between = c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ari_dibmix_gower, file = 'anova_res/aov_ari_dibmix_gower.RDS')

# Repeated Measures Anova for AMI
aov_ami_dibmix_gower <- aov_ez(id = 'ID', dv = 'AMI', fun_aggregate = mean, fullfactorial_full,
                               within = 'method',
                               between= c('nClust','overlap','nrows','ncols','catratio','pi','sph'))

saveRDS(aov_ami_dibmix_gower, file = 'anova_res/aov_ami_dibmix_gower.RDS')