# Robust hierarchical modeling of counts using the rescaled beta distribution

This repository provides R code implementing robust hierarchical modeling of counts using the rescaled beta distribution, as proposed by the following paper.

Hamura, Y., Irie, K. and Sugasawa, S. (2020). Robust hierarchical modeling of counts under zero-inflation and outliers. [arXiv:2106.10503](https://arxiv.org/abs/2106.10503)

The repository includes the following 7 files.

- RSB-reg-function.R : Script implementing the proposed method for count regression
- RSB-TF-function.R : Script implementing the proposed method for locally adaptive smoothing (trend filtering)
- RSB-spatial-function.R : Script implementing the proposed method for spatial count regression 
- CrimeDataset.RData: Dataset used in ``Analysis-spatial.R``
- Analysis-spatial.R : Script for analyzing the crime data via robust spatial count regression 
- Sim-reg.R: Script for simulation study with count regression 
- Sim-TF.R: Script for simulation study with trend filtering 
