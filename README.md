# Simulations using the IBclust package

This repository includes `R` scripts that can be used to reproduce the results for simulation studies comparing the DIBmix algorithm to KAMILA, K-Prototypes, FAMD with K-Means and PAM with Gower's dissimilarity. The three scripts can be run to reproduce the following:

- `SimStudy_real_script.R`: Simulations on ten real-world datasets of mixed-type (included in the `data` directory).
- `SimStudy_synthetic_script.R`: Benchmarking study on artificial data with non-spherical (elliptical) clusters.
- `SimStudy_synthetic_sph_script.R`: Benchmarking study on artificial data with spherical clusters.

The three scripts compare the DIBmix, KAMILA, K-Prototypes, FAMD with K-Means and PAM with Gower's dissimilarity clustering methods for mixed-type data. `DIBmix` is a function available in the `IBclust` package which you can install directly from GitHub using `devtools`:

```r
install.packages("devtools")  # Install devtools if not already installed
devtools::install_github("amarkos/IBclust")  # Install IBclust from GitHub
install.packages("IBclust") # Install IBclust from CRAN
```
