# kmmr: Hybrid GPU/CPU implementation of AI-REML for Multivariate Linear Mixed Models

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![arXiv](https://img.shields.io/badge/arXiv-2610.16671-b31b1b.svg)](https://arxiv.org/abs/2610.16671)

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

---
If you use (parts of) this software in your project, please cite the following pre-print:

Melsen, K.A.C., van Eeuwijk, F.A. & Peeters, C.F.W. (2026).
A hybrid GPU/CPU implementation of the average information restricted maximum likelihood algorithm for fitting multivariate linear mixed models.
*arXiv:2610.16671*
https://doi.org/10.48550/arXiv.2610.16671
