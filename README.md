# Simulations using the IBclust package

This repository includes `R` scripts that can be used to reproduce the results of the simulation studies comparing the DIBmix algorithm to KAMILA, K-Prototypes, FAMD with K-Means and PAM with Gower's dissimilarity. Preliminary experiments on the use of the knee heuristic and the Gap statistic for estimating the number of clusters, a runtime analysis, and a script of repeated-measures analysis of variance (ANOVA) comparing DIBmix to each of the four competing methods are also included. The scripts can be run to reproduce the following:

- `DIBmix_GAP_sim_nonsep.R`: Random generation of 100 synthetic data sets including 4 non-spherical (elliptical) clusters with moderate overlap. The features for each data set are randomly independently permuted to produce 100 reference data sets, on which the DIBmix algorithm is run. The DIBmix algorithm is run with 2 up to 8 clusters for each data set and its reference data sets and the Gap statistic is computed.
- `DIBmix_GAP_sim_sep.R`: Random generation of 100 synthetic data sets including 4 spherical clusters with low overlap. The features for each data set are randomly independently permuted to produce 100 reference data sets, on which the DIBmix algorithm is run. The DIBmix algorithm is run with 2 up to 8 clusters for each data set and its reference data sets and the Gap statistic is computed.
- `DIBmix_knee_sim_nonsep.R`: Random generation of 100 synthetic data sets including 4 non-spherical (elliptical) clusters with moderate overlap. The DIBmix algorithm is run with 2 up to 8 clusters for each data set and the knee of the mutual information curve is detected.
- `DIBmix_knee_sim_sep.R`: Random generation of 100 synthetic data sets including 4 spherical clusters with low overlap. The DIBmix algorithm is run with 2 up to 8 clusters for each data set and the knee of the mutual information curve is detected.
- `DIBmix_runtimes.R`: Runtime analysis of DIBmix, KAMILA, K-Prototypes, FAMD/K-Means & Gower/PAM on publicly available data sets of mixed-type.
- `SimStudy_real_script.R`: Applications of DIBmix, KAMILA, K-Prototypes, FAMD/K-Means & Gower/PAM on publicy available data sets of mixed-type (included in the `data` directory).
- `SimStudy_DEC_script.R`: Applications of Deep Embedded Clustering (DEC) on publicly available data sets of mixed-type.
- `SimStudy_synthetic_script.R`: Benchmarking study on artificial data with non-spherical (elliptical) clusters.
- `SimStudy_synthetic_sph_script.R`: Benchmarking study on artificial data with spherical clusters.
- `aov_pairwise.R`: Script for performing a repeated-measures analysis of variance (ANOVA) comparing average ARI and AMI scores between DIBmix and each of the 4 competing methods. The ARI and AMI values are obtained from the simulation results on the synthetic data.
- `nystrom_benchmark.R`: Script for performing a benchmarking study on the effect of the number of landmark points used in Nyström approximation on the runtime and cluster recovery performance of the Adult/Census Income data set.

The `res` directory includes the results files for the 3 simulation studies (synthetic non-spherical, synthetic spherical, real-world data), the simulation results for the knee heuristic and the Gap statistic, the runtime analysis results and the results of the  benchmarking study on the number of landmark points used in Nyström approximation. The results for the pairwise repeated-measures ANOVAs are in the `anova_res` directory.

`DIBmix` is a function available in the `IBclust` package which you can install either from CRAN or directly from GitHub using `devtools`:

```r
install.packages("devtools")  # Install devtools if not already installed
devtools::install_github("amarkos/IBclust")  # Install IBclust from GitHub
install.packages("IBclust") # Install IBclust from CRAN
```
