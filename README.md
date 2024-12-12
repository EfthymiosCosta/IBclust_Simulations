# Simulations using the DIBclust package

This repository includes `R` scripts that can be used to reproduce the results in (Costa, Papatsouma & Markos (2024))[https://arxiv.org/abs/2407.03389]. The three scripts can be run to reproduce the following:

- `SimStudy_real_script.R`: Simulations on ten real-world datasets of mixed-type (included in the `data` directory).
- `SimStudy_synthetic_script.R`: Benchmarking study on artificial data with non-spherical (elliptical) clusters.
- `SimStudy_synthetic_sph_script.R`: Benchmarking study on artificial data with spherical clusters.

The three scripts compare the DIBmix, KAMILA, K-Prototypes, FAMD with K-Means and PAM with Gower's dissimilarity clustering methods for mixed-type data. `DIBmix` is a function available in the DIBclust package which you can install directly from GitHub using `devtools`:

