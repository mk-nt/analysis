# Mk-nt analysis

This repository contains scripts necessary to replicate the analyses of 

Smith & Yang (2026), "Directional Biases in Morphological Evolution:
Implications for Phylogenetic Models".

Usage details are provided in the package vignettes.


## Framework

Our analyses were conducted on the Hamilton HPC service at Durham University.

Each analysis (dataset × model) is conducted using an MPI installation of
RevBayes. 

Morpohological datasets are referred to using their MorphoBank project
identifiers.  Models are referred to using the name of their RevBayes
implementation file, found in `inst/rbScripts`.

Once available, results are committed to a GitHub repository named by 
concatenating the dataset and model identifiers.

This R package then fetches the results via git, and caches a sub-sample
of results locally.
Analyses performed by the package are also cached.


## Setting up

Get started by checking out this repository from GitHub. 
Instructions are in the [installation vignette](https://mk-nt.github.io/analysis/articles/install.html)

Details of how to configure matrices for analysis are given in the 
[mattrix processing vignette](https://mk-nt.github.io/analysis/articles/matrix-processing.html)

A workflow to begin analysis can then be followed using [`?EnqueueMC`](https://mk-nt.github.io/analysis/reference/EnqueueMC)

Completed analyses can be retrieved from the server using 
[`?Collect`](https://mk-nt.github.io/analysis/reference/Collect) or
[`?UpdateRecords`](https://mk-nt.github.io/analysis/reference/UpdateRecords).

Further analyses are then detailed in the individual [vignettes](https://mk-nt.github.io/analysis/articles/).


## Cached data

To reproduce our results from scratch, you will need to clear cached data.
This includes:

- `inst/matrices`: Matrices obtained from MorphoBank, and manual annotations.
  This directory is a git submodule linked to
  [neo-trans/matrices](https://github.com/neo-trans/matrices)
  
- `inst/data`: Cached tree similarity measures (obtained using
  `TreeSimilarities()`) 

- `inst/projects`: Nexus files corresponding to analysis partitions,
  produced using `PrepareMatrix()`

- `inst/results`: Results of MCMCMC analysis; convergence diagnostics
  (`-conv.txt`); and marginal likelihood estimates (`-stone.txt`).


Analyses were queued using `EnqueueML()` and `EnqueueMC()`, with slurm
configuration files in `inst/slurm`.  Metadata from completed slurm jobs is
recorded in `inst/sacct.csv`.
