source("R/GitFileExists.R")
source("R/Helpers.R")
source("R/MakeSlurm.R")

# Create a new repository on GitHub
# Uses GitHub command line (gh) - https://github.com/cli/cli/releases
MakeRepo <- function(pID, scriptID) {
  suppressWarnings(
    system2("gh", sprintf("repo create %s/%s.git --private",
                          githubAccount,
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

# Create nexus and .Rev files for RevBayes analysis, and push to GitHub
# if `makeSlurm`, also creates a slurm file in this repository;
# pull this to your cluster and execute qsub.sh to enqueue the RevBayes run.
RevBayes <- function(pID, scriptID, makeSlurm = TRUE) {
  
  newRepo <- FALSE
  .FromTemplate <- function(template, scriptFile) {
    if (!dir.exists(dirname(scriptFile))) {
      clone <- .GitClone(pID, scriptID)
      if (!clone) {
        MakeRepo(pID, scriptID)
      }
    }
    
    if (!file.exists(template)) {
      stop(template, " does not exist")
    }
    
    scriptLines <- gsub("%PID%", pID,
                        gsub("%SCRIPT%", scriptID,
                             readLines(template)
                        )
    )
    write(scriptLines, scriptFile)
  }
  
  # Create run files for this project and script
  .FromTemplate("rbScripts/mcmcmc.Rev", ScriptFile(pID, scriptID))
  .FromTemplate("rbScripts/marginal.Rev", MarginalFile(pID, scriptID))
  
  scriptDir <- dirname(ScriptFile(pID, scriptID))
  
  # Copy model definition file to github subproject
  file.copy(sprintf("rbScripts/%s.Rev", scriptID),
            scriptDir,
            overwrite = TRUE)
  
  file.copy(MatrixFile(pID, "neo.nex"),
            scriptDir,
            overwrite = TRUE)
  
  file.copy(MatrixFile(pID, "trans.nex"),
            scriptDir,
            overwrite = TRUE)
  
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
    if (stash[[1]] != "No local changes to save" &&
        substr(stash[[1]], 1, 20) != "Successfully rebased") {
      pop <- system2("git", "stash pop", stdout = TRUE)
    }
  }
  
  add <- system2("git", "add *.Rev *.nex", stdout = NULL, stderr = NULL)
  if (add != 0) {
    warning(add)
  }
  
  # If anything changed, commit
  status <- system2("git", "status", stdout = TRUE)
  if (length(status) < 4 || !grepl("nothing to commit", status[[4]])) {
    commit <- system2("git", stdout = NULL, stderr = TRUE,
                   sprintf("commit -m \"Generate %s script files\"", scriptID))
    if (length(commit)) {
      warning(commit)
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
