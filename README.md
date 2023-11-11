# GEMS for latent variable selection in MIRT models

## Introduction
The generalized expectation model selection (GEMS) algorithm is proposed for dealing with the model selection problem in presence of missing data.
For the latent variable selection in multidimension two-parameter logistic model (M2PLM),
we present an efficient implementation of GEMS to find
the optimal model (i.e., the structure of item-trait relationships) and the parameter estimates
(including the item discrimination and difficulty parameters) under optimal model
which results in the smallest BIC value.
The GEMS for M2PLM is more computationally efficient than the EMS proposed by Xu et al. (2022).

The M2plmGems_codes directory contains the following 3 files: 
- M2plmGems_algorithm.cpp includes the c++ implementations of both GEMS and StepGEMS algorithm for M2PLM.
- M2plmGems_fcns.R        comprises R functions for GEMS and StepGEMS algorithms that users can directly invoke.
- examples.R              provides illustrative examples for the aforementioned algorithms.

The GEMS implementation rely on Gauss-Hermite quadrature.
Please first install the R package 'mvQuad' to ensure that the codes run successfully.
The c++ implementation is based on R-packages 'Rcpp', 'RcppArmadillo' and 'RcppClock'.
To run the examples.R, the R-packages 'magrittr' and 'MASS' are required.

All codes are written in R 4.2.3 with the corresponding Rtools version 4.2.

## Citation
To cite these codes in publications, please use the following reference:\
Shang, L., Zheng, Q. Z., Xu, P. F., Shan, N., & Tang, M. L. (2023).
A generalized expectation model selection algorithm for latent
variable selection in multidimensional item response theory
models. Statistics and Computing. Submitted for publication. 

https://github.com/Laixu3107/Gems_for_MIRT.
