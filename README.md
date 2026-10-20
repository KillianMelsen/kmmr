# kmmr: Hybrid GPU/CPU implementation of AI-REML for Multivariate Linear Mixed Models

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![arXiv](https://img.shields.io/badge/arXiv-2610.16671-b31b1b.svg)](https://arxiv.org/abs/2610.16671)
[![DOI](https://zenodo.org/badge/1178027814.svg)](https://doi.org/10.5281/zenodo.23398636)

`kmmr` is an experimental R package for variance component estimation in multivariate linear mixed models (LMMs).
It features an NVIDIA CUDA based hybrid GPU/CPU implementation of the Average Information Restricted Maximum Likelihood (AI-REML) algorithm.

---

## Installation
`kmmr` requires Linux or Windows Subsystem for Linux (WSL2) with the NVIDIA CUDA Toolkit installed.

You can install `kmmr` directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("KillianMelsen/kmmr")
```
## Citation
If you use (parts of) this software in your project, please cite the following pre-print:

Melsen, K.A.C., van Eeuwijk, F.A. & Peeters, C.F.W. (2026).
A hybrid GPU/CPU implementation of the average information restricted maximum likelihood algorithm for fitting multivariate linear mixed models.
*arXiv:2610.16671*
https://doi.org/10.48550/arXiv.2610.16671

## References
### Methodology
* Gilmour, A.R. (2019). Average information residual maximum likelihood in practice. Journal of Animal Breeding and Genetics 136, 262–272. https://doi.org/10.1111/jbg.12398
* Jensen, J., Mäntysaari, E.A., Madsen, P. & Thompson, R. (1997). Residual maximum likelihood estimation of (co)variance components in multivariate mixed linear models using average information. Journal of the Indian Society of Agricultural Statistics 49, 215–236.
* Gilmour, A.R., Thompson, R. & Cullis, B. (1995). Average information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51, 1440–1450. https://doi.org/10.2307/2533274

### Software
* Covarrubias-Pazaran, G. (2016). Genome-assisted prediction of quantitative traits using the R package sommer. PLoS One 11, e0156744. https://doi.org/10.1371/journal.pone.0156744
* Luttinen, J. (2012). CHOLMOD Extra Module. GitHub Repository. https://github.com/jluttine/cholmod-extra
* Chen, Y., Davis, T.A., Hager, W.W. & Rajamanickam, S. (2008). Algorithm 887: CHOLMOD, supernodal sparse Cholesky factorization and update/downdate. ACM Transactions on Mathematical Software 35, 22:1–22:14. https://doi.org/10.1145/1391989.1391995
* Luebke, D. (2008). CUDA: Scalable parallel programming for high-performance scientific computing. In 2008 5th IEEE International Symposium on Biomedical Imaging: From Nano to Macro. 461–468. https://doi.org/10.1109/ISBI.2008.4541126
