# Riemannian approach to CoDa regression

This repository contains the data and R scripts required to reproduce the results presented in the paper:
**"Riemannian approach to CoDa regression"** by N. Trendafilov, M. Gallo, V. Simonacci, and V. Todorov (2026).

## Overview
Standard compositional data (CoDa) regression problems often rely on log-ratio transformations (extrinsic approach) or standard Aitchison geometry. However, when closed-form solutions are unavailable, standard Euclidean optimization may fail to respect the simplex constraints. 

This repository provides an algorithmic framework that treats the CoDa sample space natively as a Riemannian manifold. We construct Interior Point Flows (IPFs) to strictly enforce the constant-sum and positivity constraints natively on the manifold, providing a geometrically coherent alternative to standard log-ratio approaches.

## Repository Structure

* `data/`: Contains the benchmark compositional datasets used in the paper.
    * `arctic_lake.csv`: Sediment samples (Aitchison, 1986).
    * `gemas_soils.csv`: Agricultural and grazing land soils (Reimann et al., 2014).
    * `coxite.csv`: Mineral compositions (Aitchison, 1986).
    * `probiotics.csv`: Probiotics Intervention dataset (Lahti et al., 2013).
* `R/`: Contains the core R functions.
    * `ipf_coda_reg.R`: The main function implementing the Riemannian IPF regression for both scalar and compositional predictors.

## Prerequisites

The code is written in R. The only external dependency required to run the differential equations (IPFs) is the `deSolve` package.

```R
install.packages("deSolve")