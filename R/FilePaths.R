#' List all projects
#'
#' Return identifiers for all available projects whose matrices have been
#' prepared for analysis.
#'
#' @param indices Optional integer vector specifying which projects to
#'   return (after sorting).
#' @param sort Logical; whether to sort project IDs numerically.
#' @return Character vector of project identifiers.
#' @export
AllProjects <- function(indices = NULL, sort = TRUE) {
  x <- gsub("^\\D*(\\d+)\\D*?$", "\\1", perl = TRUE,
            list.files(dirname(MatrixFile("*")),
                       gsub("project", "(project|syab)",
                            basename(MatrixFile(".*", "neo.nex$")))))
  x <- setdiff(x, "4173") # duplicate of 4285
  x <- setdiff(x, "4284") # Cannot get ML even for by_ki
  x <- x[order(as.numeric(x))]
  if (is.null(indices)) x else x[indices]
}

#' List projects suitable for _ki analyses
#'
#' Return identifiers for projects retained for downstream `_ki`
#' analyses, excluding those known to fail or duplicate others.
#'
#' @return Character vector of project identifiers.
#' @export
KiProjects <- function() {
  setdiff(AllProjects(), c(.config$syab[-3], 684, 691, 3345, 3670, 3707,
                           3708, 3763, 4103, 4173, 4284, 4285, 4817))
}

#' Path to directories in which files are stored locally.
#' 
#' Set the paths in which files should be stored using
#' `options("ntOutDir" = "path/to/analytical/outputs")`,
#' `options("ntSlurmDir" = "path/to/slurm/files")`
#' `options("ntRBScriptDir" = "path/to/RevBayes/scripts")`
#' and
#' `options("ntRepoDir" = "path/where/repos/will/be/checked/out/to")`
#' 
#' Set `options("ntRemoteDir" = "path/on/remote/server")` to the directory
#' on the remote server to which the `neotrans` repository is checked out.
#' 
#' These will default respectively to:
#' - `OutputDir()`: the package's installation location
#' - `SlurmDir()`: the folder "slurm" in the package's installation location.
#' - `RepoDir()`: a directory `revbayes-repos` placed within the directory
#'   into which the package is installed.
#' - `RemoteDir()`: `/nobackup/$USER`
#' 
#' @export
OutputDir <- function() {
  getOption("ntOutDir") %||% system.file(package = "neotrans")
}

#' @rdname OutputDir
#' @export
RepoDir <- function() {
  getOption("ntRepoDir") %||% {
    file.path(dirname(gsub("/inst$", "", perl = TRUE,
                           system.file(package = "neotrans"))), "revbayes-repos")
  }
}

#' @rdname OutputDir
#' @export
SlurmDir <- function() {
  getOption("ntSlurmDir") %||% system.file("slurm", package = "neotrans")
}

#' @rdname OutputDir
#' @export
RemoteDir <- function() {
  getOption("ntRemoteDir") %||% "/nobackup/$USER"
}

#' @rdname OutputDir
#' @export
RBScriptDir <- function() {
  getOption("ntRBScriptDir") %||% system.file("rbScripts", package = "neotrans")
}

#' Directory containing analysis output
#'
#' Returns the path to the directory containing RevBayes log and tree files
#' for a specific project and script.
#' @return Character string giving directory path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
AnalysisDir <- function(pID, scriptID = NULL) {
  dirname(ScriptFile(pID, scriptID))
}

#' Directory containing stored results
#'
#' Returns the path to where processed results (e.g. `.rds` files, summaries)
#' are written for a given project and script.  
#' Creates the directory if it does not yet exist.
#'
#' @return Character string giving directory path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
ResultsDir <- function(pID, scriptID = NULL) {
  d <- sprintf(file.path(OutputDir(), "results", "project%s"), pID)
  if (!dir.exists(d)) {
    dir.create(d)
  }
  d
}

#' Path to parsimony-evaluation results
#'
#' Constructs the path to the `.rda` file storing results of parsimony and
#' Bayesian evaluation for a given project, script, and parameter combination.
#'
#' @param nPars,nBayes Integer counts of parsimony and Bayesian replicates.
#' @return Character string giving `.rda` file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
ParsEvalFile <- function(pID, scriptID = NULL, nPars, nBayes) {
  d <- file.path(OutputDir(), "data", "pars-eval")
  if (!dir.exists(d)) {
    dir.create(d)
  }
  file.path(d, sprintf("%s-%dprs-%dbys.rda", ScriptBase(pID, scriptID),
                       nPars, nBayes))
}

#' Path to convergence-diagnostic file
#'
#' @return Character string giving `.txt` path summarizing convergence checks.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
ConvergenceFile <- function(pID, scriptID) {
  sprintf("%s/%s-conv.txt",
          ResultsDir(pID, scriptID), ScriptBase(pID, scriptID))
}

#' Path to tree-distance results
#'
#' @param model1,model2 Names of two models to compare.
#' @return Character string giving `.rds` path containing distance metrics.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
DistanceFile <- function(pID, model1, model2) {
  models <- sort(c(model1, model2))
  sprintf("%s/dist-%s-%s%s.rds", ResultsDir(pID), models[[1]], models[[2]],
          if (HasConverged(pID, model1) && HasConverged(pID, model2))
            "" else ".tmp")
}

#' Path to marginal-likelihood output
#' @return Character string giving file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
MarginalFile <- function(pID, scriptID) {
  ScriptFile(pID, scriptID, "marginal")
}

#' Path to posterior-predictive sample file
#' @return Character string giving file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
PPSamplerFile <- function(pID, scriptID) {
  ScriptFile(pID, scriptID, "ppsample")
}

#' Path to matrix file
#'
#' Returns the location of a project's morphological matrix file.
#'
#' @param ext File extension (default `"nex"`).
#' @return Character string giving path to matrix file.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
MatrixFile <- function(pID = "", ext = "nex") {
  sprintf(file.path(OutputDir(), "projects", "%s%s.%s"),
          if (startsWith(as.character(pID), "072")) "syab" else "project",
          pID, ext)
}

#' Path to stored parameter summary
#'
#' Returns path to the `.rds` file containing summarized parameter estimates.
#' @return Character string giving file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
ParameterFile <- function(pID, scriptID) {
  sprintf("%s/%s.rds", ResultsDir(pID, scriptID), ScriptBase(pID, scriptID))
}

#' List RevBayes log-file paths
#'
#' @return Character vector of full paths to `.log` files from MCMC runs.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
PFiles <- function(pID, scriptID) {
  list.files(
    path = AnalysisDir(pID, scriptID),
    pattern = sprintf("%s\\.p_run_.*\\.log$", scriptID),
    full.names = TRUE, ignore.case = TRUE
  )
}

#' Path to time-usage log
#' @return Character string giving path to timing log.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
TimeFile <- function(pID, scriptID) {
  file.path(AnalysisDir(pID, scriptID),
            paste0(ScriptBase(pID, scriptID), ".time.log"))
}

#' Path to memory-usage log
#' @return Character string giving path to memory-usage log.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
MemFile <- function(pID, scriptID) {
  file.path(AnalysisDir(pID, scriptID),
            paste0(ScriptBase(pID, scriptID), ".mem.log"))
}

#' Path to temporary disk-usage log
#' @return Character string giving path to disk-usage log.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
DiskFile <- function(pID, scriptID) {
  file.path(AnalysisDir(pID, scriptID),
            paste0(ScriptBase(pID, scriptID), ".tmp_usage.log"))
}

#' Path to PDF comparison file
#' @export
PDFFile <- function(pID, model1, model2) {
  file.path(OutputDir(), "pdf", sprintf("%s_%s-%s.pdf", pID, model1, model2))
}

#' Construct base name for a script
#'
#' @return Character string of the form `"pID_scriptID"`.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
ScriptBase <- function(pID, scriptID) {
  sprintf("%s_%s", pID, scriptID)
}

#' Path to dummy input matrix for testing
#'
#' Returns the path used when a placeholder matrix is needed to run RevBayes.
#' @return Character string giving path to dummy `.nex` file.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
DummyFile <- function(pID, scriptID) {
  sprintf(file.path(RepoDir(), "%s", "project%s.01.nex"),
          ScriptBase(pID, scriptID), pID)
}

#' Path to RevBayes script
#'
#' Constructs the path to a `.Rev` analysis script within the corresponding
#' RevBayes repository.
#'
#' @param analysis One of `"mcmcmc"`, `"marginal"`, or `"ppsample"`.
#' @return Character string giving file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
ScriptFile <- function(pID, scriptID, analysis = "mcmcmc") {
  sprintf(file.path(RepoDir(), "%s/%s.Rev"),
          ScriptBase(pID, scriptID), analysis)
}

#' Path to SLURM submission template
#'
#' Returns the template script appropriate for a given analysis mode.
#'
#' @param ml Logical; whether this is a marginal-likelihood run.
#' @param mpi Logical; whether the analysis uses MPI parallelization.
#' @return Character string giving path to SLURM template.
#' @family file-path helpers
#' @export
SlurmTemplate <- function(ml, mpi = TRUE) {
  file.path(
    SlurmDir(),
    if (ml) {
      "marginal.sh"
    } else if (mpi) {
      "mcmcmc.sh"
    } else {
      "mc3serial.sh"
    }
  )
}

#' Path to generated SLURM script
#'
#' Constructs the path to the SLURM job script for a given project and script.
#' @param ml Logical; whether this corresponds to a marginal-likelihood run.
#' @return Character string giving file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
SlurmFile <- function(pID, scriptID, ml = FALSE) {
  sprintf("%s%s_%s.sh", if (ml) "ml-" else "", pID, scriptID)
}

#' Path to stepping-stone results
#'
#' Returns the path to the file containing stepping-stone marginal-likelihood
#' estimates for a given project and script.
#'
#' @return Character string giving `.txt` file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
StoneFile <- function(pID, scriptID) {
  sprintf(file.path(ResultsDir(pID, scriptID), "%s-stone.txt"),
          ScriptBase(pID, scriptID))
}

#' Path to stepping-stone origin file
#'
#' Returns the path to the `.out.pp` file written by RevBayes during
#' stepping-stone sampling.
#' @return Character string giving path to `.out.pp` file.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
StoneOrigin <- function(pID, scriptID) {
  sprintf(file.path(AnalysisDir(pID, scriptID), "%s.out.pp"),
          ScriptBase(pID, scriptID))
}

#' List tree-file paths
#'
#' Returns a vector of paths to `.trees` (or compressed `.tar.gz`) files
#' generated by an analysis.
#'
#' @param compressed Logical; whether to list compressed files.
#' @return Character vector of tree-file paths.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
TreeFiles <- function(pID, scriptID, compressed = FALSE) {
  list.files(
    path = AnalysisDir(pID, scriptID),
    pattern = sprintf("%s\\_run_.*\\.%s$", scriptID,
                      if (compressed) "tar.gz" else "trees"),
    full.names = TRUE, ignore.case = TRUE
  )
}

#' Path to combined tree sample
#'
#' Returns the path to the combined posterior tree sample used in downstream
#' analyses.
#' @return Character string giving `.trees` file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
TreeSampleFile <- function(pID, scriptID) {
  sprintf(file.path(ResultsDir(pID, scriptID), "%s.trees"),
          ScriptBase(pID, scriptID))
}

#' Path to parsimony-score results
#'
#' Returns the path to the `.rds` file containing parsimony scores of sampled
#' trees for a given project and script.
#' @return Character string giving `.rds` file path.
#' @inheritParams MakeSlurm
#' @family file-path helpers
#' @export
ParsimonyFile <- function(pID, scriptID) {
  sprintf(file.path(ResultsDir(pID, scriptID), "%s.trees.pscore.rds"),
          ScriptBase(pID, scriptID))
}
