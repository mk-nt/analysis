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




The analyses make use of several services in order to accomplish scale.

We assume that you have `git` and `gh` installed locally and in the system path.

To reproduce the analyses:

1. Fork this GitHub repository.

2. Open `R/_setup.R` and adjust configuration variables.
   You will need to replace the value of `githubAccount` with a (probably
   organizational) account.
   You will need to configure a personal access token (PAT) with permissions to
   create and write to repositories in that GitHub account.
   A PAT with write access should be stored on your cluster as the 
   environment variable `MKNT_WRITE`.
   A PAT with read access should be stored on your computer as the
   environment variable `MKNT_READ`.
   
3. Execute `R/1_PrepareMatrix.R`.
   This script uses the Excel spreadsheets in `matrices/*.xlsx` to identify each
   character in the corresponding nexus file `matrices/*.nex` as neomorphic
   or transformational.
   It creates two new nexus files, 
   `projects/project<projectID>.trans.nex` containing the transformational
   characters, `*.neo.nex`, the neomorphic characters.
   Parsimony-uninformative characters are removed.
   
   These matrices, along with RevBayes scripts for phylogenetic inference 
   (`mcmcmc.Rev`) and marginal likelihood calculation (`marginal.Rev`),
   created from templates in `rbScripts`,
   are commited to a new GitHub repository at
   `githubAccount/<projectID>_<modelID>`.
   
   A slurm script for each analysis is created in the `slurm` directory, using
   the templates `slurm/mcmcmc.sh` and `slurm/marginal.sh`.
   References to `mk-nt` in these templates should be updated to refer to
   your `githubAccount`.
   
4. Log into a computing cluster and configure RevBayes.
   We recommend building the development version, which fixes some bugs in
   v1.2.4 that caused some analyses to fail.
   
5. Check out your fork of the repository, navigate to its root directory,
   and execute `bash qsub.sh`.  This will submit the slurm jobs defined by
   `RevBayes()` to the slurm queue.  Once complete, each job will commit its
   results to its associated GitHub repository.
   
6. Back on your own computer, execute `R/2_CheckCompletion.R`.
   `UpdateRecords()` fetches RevBayes results from GitHub and performs initial
   analyses, such as computing convergence diagnostics and estimated sample
   size.

   Some datasets may require extended runtimes or additional memory to complete.
   You may need to modify `R/MakeSlurm.R` to reflect the available resources on
   your own compute cluster.
   
   Results once available are stored in `results/project<projectID>`:
   - Each `<projectID>_<modelID>.trees` file contains a sample of 256 topologies
     from the posterior.
   - `*.rds` stores the distances between each pair of trees.
   - `*-conv.txt` records convergence diagnostics and the estimated sample size
     after discarding a fraction of the results (given by `burnin`) as burnin.
     Note that RevBayes discards the 3000 generations in which move proposals
     were being optimised; in some cases these generations suffice as a burnin
     phase.
   - `*-stone.txt` records:
    - An estimate of the marginal likelihood obtained by stepping stone sampling
    - An estimate of the marginal likelihood obtained by path sampling
    - The number of generations discarded as burnin when deriving the power
      posterior
    - The number of generations used in the power posterior run
    - The number of steps used in the power posterior run

7. Once analyses are complete, execute `R/3_Analysis.R`.
   This script is admittedly inelegant.
   It contains code used to explore the properties of the results, as well as
   to create the manuscript figures and supplementary table.
   
   Exploring changes in tree topology between models requires distance
   computations to be cached. 
   Where a cache is unavailable, `.TreeStats()` will print a list of calls to
   `AnalysisPDF()` that should be conducted after running
   `source("R/AnalysisPDF.R")`.
   This will cache tree-to-tree distances in
   `rbPDF/project<projectID>-<modelID1>-<modelID2>.dist.rds`, 
   and visualize distances in a summary PDF `rbPDF/*.pdf`.


## Archiving analytical results on GitHub

`for repo in $(gh repo list $USERNAME --limit 300 --json name --jq '.[] | select(.name | endswith("_ki")) | .name'); do gh release create "2025-10" --repo "mk-nt/$repo" --title "2025-10" --notes "Archive $repo for Zenodo deposit."; done`