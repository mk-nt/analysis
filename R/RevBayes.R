#' Create a new RevBayes analysis repository
#'
#' `MakeRepo()` establishes a new private GitHub repository for a given
#' project–script combination.  
#' The repository will host all RevBayes scripts and associated matrix files
#' required for downstream analysis.
#'
#' @details
#' Before using this function, ensure that:
#' \enumerate{
#'   \item The [GitHub CLI](https://cli.github.com/) (`gh`) is installed and on
#'   your system `PATH`.
#'   \item You are authenticated with write access to the account specified in
#'   `.Renviron` as:
#'   \preformatted{
#'   ntGithubAccount=neo-trans
#'   }
#'   \item You have configured your `.Renviron` file with SSH details for your
#'   analysis cluster (see below).
#' }
#'
#' Typical `.Renviron` entries (Windows paths shown) are:
#' \preformatted{
#' sshLogin=username@remote-host
#' sshKey=C:/Users/you/.ssh/id_ed25519
#' sshPass=key-passphrase
#' rb.exe=C:/Program Files/revbayes-vX.Y.Z/bin/rb.exe
#' }
#'
#' You can edit this file by running `usethis::edit_r_environ()` and must
#' restart R for changes to take effect.
#'
#' On success, a new repository named
#' `<ntGithubAccount>/<project>_<script>.git` is created on GitHub and cloned
#' locally to the appropriate subdirectory.
#'
#' @inheritParams MakeSlurm
#'
#' @return
#' `MakeRepo()` invisibly returns `TRUE` when the repository has been created
#' and initialized successfully.  
#' The function is called primarily for its side-effect of establishing a
#' corresponding GitHub repository and local working copy.
#'
#' @seealso
#' [RevBayes()], [EnqueueMC()], [EnqueueML()],
#' and the `zzz.R` startup configuration notes for authentication setup.
#'
#' @export
MakeRepo <- function(pID, scriptID) {
  if (system2("gh", "--version",
              stdout = FALSE) == 127) {
    stop("gh not found; install it from https://cli.github.com/ and add to PATH")
  }
  suppressWarnings(
    system2("gh", sprintf("repo create %s/%s.git --private",
                          Sys.getenv("ntGithubAccount"),
                          ScriptBase(pID, scriptID)),
            stdout = NULL)
  )
  scriptDir <- dirname(ScriptFile(pID, scriptID))
  cloned <- .GitClone(pID, scriptID)
  
  wd <- setwd(scriptDir)
  on.exit(setwd(wd))
  init <- system2("git", "init", stdout = TRUE)
  branch <- system2("git", "branch -M main", stdout = TRUE)
}

#' Prepare and submit a new RevBayes analysis job
#'
#' `RevBayes()` constructs all necessary RevBayes input files for a given
#' project and model specification, synchronizes them with the appropriate
#' GitHub repository, and (optionally) prepares SLURM submission scripts.
#'
#' @details
#' The typical workflow is:
#' \enumerate{
#'   \item Use [PrepareMatrix()] to generate standardized `.nex` matrix files.
#'   \item Run [MakeRepo()] to create a GitHub repository for this
#'   project–model combination.
#'   \item Call `RevBayes(pID, scriptID)` to populate the repository with
#'   RevBayes model scripts.
#'   \item Use [MakeSlurm()] or `makeSlurm = TRUE` to create SLURM job files.
#' }
#'
#' The function draws its templates from the `rbScripts/` directory, performs
#' placeholder substitution (`%PID%`, `%SCRIPT%`, `%MATRIXBASE%`), and copies
#' the appropriate matrix subsets (`neo.nex`, `trans.nex`, or their randomised
#' variants) into the repository.
#'
#' When `makeSlurm = TRUE`, both MCMC and marginal-likelihood SLURM job files
#' are generated automatically.
#'
#' @param pID Character. Project identifier (e.g., `"1210"`).
#' @param scriptID Character. Identifier of the model script to be used.
#'   Prefix `"rm_"` denotes random partitioning.
#' @param makeSlurm Logical. Whether to generate SLURM job scripts immediately.
#'   Defaults to `FALSE`.
#'
#' @return
#' `RevBayes()` invisibly returns `TRUE` when script generation and Git
#' synchronization succeed.  
#' The function is called primarily for its side-effect of populating
#' and committing RevBayes model repositories.
#'
#' @seealso
#' [PrepareMatrix()], [MakeRepo()], [MakeSlurm()],
#' [EnqueueMC()], [EnqueueML()]
#'
#' @export
RevBayes <- function(pID, scriptID, makeSlurm = FALSE) {
  
  newRepo <- FALSE
  .FromTemplate <- function(templateFile, scriptFile) {
    template <- file.path(RBScriptDir(), templateFile)
    if (!dir.exists(dirname(scriptFile))) {
      clone <- .GitClone(pID, scriptID)
      if (!clone) {
        MakeRepo(pID, scriptID)
      }
    }
    
    bias <- switch(sub("^.*_k(.)\\.Rev$", "\\1", template, perl = TRUE),
                   "i" = "informative",
                   "v" = "variable",
                   "a" = "all", 
                   NA)
    if (!file.exists(template)) {
      template <- gsub("_k.\\.Rev$", "_ki.Rev", template)
    }
    if (!file.exists(template)) {
      stop(template, " does not exist")
    }
    
    scriptLines <- gsub(
      fixed = TRUE,
      "coding = \"informative\"", paste0("coding = \"", bias, "\""),
      gsub("%PID%", pID, fixed = TRUE,
           gsub("%SCRIPT%", scriptID, fixed = TRUE,
                gsub("%MATRIXBASE.%", basename(MatrixFile(pID, "")),
                     readLines(template), fixed = TRUE
                     )
                )
           )
    )
    
    write(scriptLines, scriptFile)
  }
  
  scriptDir <- dirname(ScriptFile(pID, scriptID))
  
  # Copy model definition file to github subproject
  .FromTemplate(sprintf("%s.Rev", gsub("^rm_", "", scriptID)),
                sprintf("%s/%s.Rev", scriptDir, scriptID))
  
  if (!dir.exists(scriptDir)) { # Created by .FromTemplate()
    stop("Could not create directory ", scriptDir)
  }
  .FromTemplate("mcmcmc.Rev", ScriptFile(pID, scriptID))
  .FromTemplate("marginal.Rev", MarginalFile(pID, scriptID))
  .FromTemplate("ppsample.Rev", PPSamplerFile(pID, scriptID))
  
  randomizePartitions <- startsWith(scriptID, "rm_")
  neo.nex <- if (randomizePartitions) "neo-rand.nex" else "neo.nex"
  trans.nex <- if (randomizePartitions) "trans-rand.nex" else "trans.nex"
  if (!file.exists(MatrixFile(pID, neo.nex))) {
    stop(MatrixFile(pID, neo.nex), " does not exist")
  }
  file.copy(MatrixFile(pID, neo.nex), overwrite = TRUE,
            sprintf("%s/%s", scriptDir, basename(MatrixFile(pID, "neo.nex"))))
  file.copy(MatrixFile(pID, trans.nex), overwrite = TRUE, 
            sprintf("%s/%s", scriptDir, basename(MatrixFile(pID, "trans.nex"))))
  
  # Commit changes
  wd <- setwd(scriptDir)
  on.exit(setwd(wd))
  
  # Ensure we are on main branch
  branch <- system2("git", "branch", stdout = TRUE)
  if (!length(branch) || branch[[1]] != "* main") {
    std <- suppressWarnings(
      system2("git", "checkout main", stdout = TRUE, stderr = TRUE)
    )
    if (length(std)) {
      if (any(grepl("branch is behind", std))) {
        stash <- system2("git", "stash", stdout = TRUE)
        rebase <- system2("git", "rebase origin/main",
                          stdout = TRUE, stderr = TRUE)
        if (!grepl("Success", rebase)) {
          warning(rebase)
        }
        if (stash[[1]] != "No local changes to save") {
          pop <- system2("git", "stash pop", stdout = TRUE)
        }
        
      } else if (any(grepl("pathspec 'main' did not match", std))) {
        newRepo <- TRUE
      } else if (!any(grepl("Already on 'main'", std))) {
        warning(std)
      }
    }
  }
  
  if (!newRepo) {
    # Pull if necessary
    stash <- system2("git", "stash", stdout = TRUE)
    pull <- system2("git", "pull --rebase", stdout = TRUE)
    if (pull[[1]] != "Already up to date.") {
      message(paste(pull, collapse = "\r\n"))
      system2("git", "status")
    }
    
    # Now we've stashed and pulled,
    # Create run files for this project and script
    gitWd <- setwd(wd)
    .FromTemplate("mcmcmc.Rev", ScriptFile(pID, scriptID))
    .FromTemplate("marginal.Rev", MarginalFile(pID, scriptID))
    setwd(gitWd)
    add <- system2("git", "add *.Rev", stdout = NULL, stderr = NULL)
    if (add != 0) {
      warning(add)
    }
    
    if (stash[[1]] != "No local changes to save" &&
        substr(stash[[1]], 1, 20) != "Successfully rebased") {
      pop <- suppressWarnings(system2("git", "stash pop", stdout = TRUE))
      popStatus <- attr(pop, "status") %||% 0
      if (popStatus != 0) {
        stop(paste0(
          c("Could not restore stash; check for a conflict.\r\n", pop),
          "\r\n  "))
      }
    }
  }
  
  add <- system2("git", "add *.Rev *.nex", stdout = NULL, stderr = NULL)
  if (add != 0) {
    warning(add)
  }
  
  # If anything changed, commit
  status <- system2("git", "status", stdout = TRUE)
  if (length(status) < 4 || !grepl("nothing to commit", status[[4]])) {
    commit <- suppressWarnings(system2("git", stdout = NULL, stderr = TRUE,
                   sprintf("commit -m \"Generate %s script files\"", scriptID)))
    if (length(commit)) {
      warning(commit, .immediate = TRUE)
    }
    .GitPush("origin HEAD:main")
  }
  
  # Return to main repo and stage change to submodule
  setwd(wd)
  
  if (makeSlurm) {
    MakeSlurm(pID, scriptID, ml = FALSE)
    MakeSlurm(pID, scriptID, ml = TRUE)
  }
  
  .GitPush()
}

#' Queue projects for MCMC analysis
#'
#' `EnqueueMC()` automates the complete setup of MCMC analyses across multiple
#' projects and model definitions.  
#' It constructs matrices, prepares RevBayes repositories and scripts, and
#' generates SLURM job files where needed.
#'
#' @details
#' For each project in `projects`, `EnqueueMC()`:
#' \enumerate{
#'   \item Calls [PrepareMatrix()] to ensure up-to-date NEXUS matrices.
#'   \item Uses [RevBayes()] to copy or generate RevBayes scripts.
#'   \item Checks for an existing convergence file via [`ConvergenceFile()`].
#'   \item If missing or incomplete, creates MCMC SLURM job files with
#'   `MakeSlurm(ml = FALSE)`.
#' }
#'
#' The function uses [cli::cli_progress_message()] to display progress and
#' finishes with [cli::cli_progress_done()].
#'
#' @param projects Character vector of project identifiers.
#' @param models Character vector of RevBayes model identifiers.
#' @param overwrite Logical. Whether to overwrite existing processed matrices.
#'   Defaults to `TRUE`.
#'
#' @return
#' `EnqueueMC()` invisibly returns `NULL`.  
#' The function is called for its side-effect of generating and enqueuing
#' RevBayes MCMC analyses on the configured cluster.
#'
#' @seealso
#' [EnqueueML()] for marginal likelihood estimation,
#' [RevBayes()], [MakeSlurm()], [PrepareMatrix()]
#'
#' @export
EnqueueMC <- function(projects, models, overwrite = TRUE) {
  for (pID in projects) {
    cli_progress_message(paste0("Project ", pID))
    if (PrepareMatrix(pID, overwrite = overwrite)) {
      for (scriptID in models) {
        RevBayes(pID, scriptID, makeSlurm = FALSE)
        if (!file.exists(ConvergenceFile(pID, scriptID))
            || is.na(read.table(ConvergenceFile(pID, scriptID))[["burnin"]])) {
          MakeSlurm(pID, scriptID, ml = FALSE)
        }
      }
    }
  }
  cli_progress_done()
}

#' Queue projects for marginal likelihood estimation
#'
#' `EnqueueML()` mirrors the workflow of [EnqueueMC()] but instead enqueues
#' jobs for marginal likelihood estimation using the stepping-stone or
#' path-sampling methods in RevBayes.
#'
#' @details
#' Each project’s matrix is prepared with [PrepareMatrix()] and paired with
#' model scripts generated via [RevBayes()].  
#' If the marginal likelihood output ([`StoneFile()`]) is absent, SLURM job files
#' are created using `MakeSlurm(ml = TRUE)`.
#'
#' This function requires that SSH credentials and GitHub authentication have
#' been configured in `.Renviron`, as described in [MakeRepo()].
#'
#' @param projects Character vector of project identifiers.
#' @param models Character vector of RevBayes model identifiers.
#' @param overwrite Logical. Whether to overwrite existing processed matrices.
#'   Defaults to `TRUE`.
#'
#' @return
#' `EnqueueML()` invisibly returns `NULL`.  
#' The function is called for its side-effect of preparing and queuing
#' marginal-likelihood analyses on the configured remote cluster.
#'
#' @seealso
#' [EnqueueMC()], [RevBayes()], [MakeSlurm()], [PrepareMatrix()]
#'
#' @importFrom cli cli_progress_done cli_progress_message
#' @export
EnqueueML <- function(projects, models, overwrite = TRUE) {
  for (pID in projects) {
    cli_progress_message(paste0("Project ", pID))
    if (PrepareMatrix(pID, overwrite = overwrite)) {
      for (scriptID in models) {
        RevBayes(pID, scriptID, makeSlurm = FALSE)
        if (!file.exists(StoneFile(pID, scriptID))) {
          MakeSlurm(pID, scriptID, ml = TRUE)
        }
      }
    }
  }
  cli_progress_done()
}
