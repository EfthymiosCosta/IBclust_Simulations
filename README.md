# Simulations using the IBclust package

This repository includes `R` scripts that can be used to reproduce the results of the simulation studies comparing the DIBmix algorithm to KAMILA, K-Prototypes, FAMD with K-Means and PAM with Gower's dissimilarity. Preliminary experiments on the use of the knee heuristic for estimating the number of clusters and a script of repeated-measures analysis of variance (ANOVA) comparing DIBmix to each of the four competing methods are also included. The six scripts can be run to reproduce the following:

- `DIBmix_knee_sim_nonsep.R`: Random generation of 100 synthetic data sets including 4 non-spherical (elliptical) clusters with moderate overlap. The DIBmix algorithm is run with 2 up to 8 clusters for each data set and the knee of the mutual information curve is detected.
- `DIBmix_knee_sim_sep.R`: Random generation of 100 synthetic data sets including 4 spherical clusters with low overlap. The DIBmix algorithm is run with 2 up to 8 clusters for each data set and the knee of the mutual information curve is detected.
- `SimStudy_real_script.R`: Simulations on ten real-world datasets of mixed-type (included in the `data` directory).
- `SimStudy_synthetic_script.R`: Benchmarking study on artificial data with non-spherical (elliptical) clusters.
- `SimStudy_synthetic_sph_script.R`: Benchmarking study on artificial data with spherical clusters.
- `aov_pairwise.R`: Script for performing a repeated-measures analysis of variance (ANOVA) comparing average ARI and AMI scores between DIBmix and each of the 4 competing methods. The ARI and AMI values are obtained from the simulation results on the synthetic data.

The `res` directory includes the results files for the 3 simulation studies (synthetic non-spherical, synthetic spherical, real-world data), and the simulation results for the knee heuristic. The results for the pairwise repeated-measures ANOVAs are in the `anova_res` director.

`DIBmix` is a function available in the `IBclust` package which you can install either from CRAN or directly from GitHub using `devtools`:

```r
install.packages("devtools")  # Install devtools if not already installed
devtools::install_github("amarkos/IBclust")  # Install IBclust from GitHub
install.packages("IBclust") # Install IBclust from CRAN
```
