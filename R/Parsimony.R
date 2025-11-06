#' Parsimony lenths of trees in posterior distribution
#' @importFrom ape read.tree
#' @importFrom TreeTools MatrixToPhyDat ReadCharacters RootTree
#' @importFrom TreeSearch TreeLength
#' @export
PosteriorTreeSteps <- function(projects, models) {
  do.call(
    rbind,
    lapply(projects, function(pID) do.call(
      rbind,
      lapply(models, function(scriptID) {
        pFile <- ParsimonyFile(pID, scriptID)
        if (file.exists(pFile)) {
          return(readRDS(pFile))
        }
        if (!file.exists(TreeSampleFile(pID, scriptID))) {
          UpdateRecords(pID, scriptID)
          if (!file.exists(TreeSampleFile(pID, scriptID))) {
            MakeSlurm(pID, scriptID, ml = FALSE)
            return(NULL)
          }
        }
        tryCatch({
          nexusFiles <- list.files(AnalysisDir(pID, scriptID),
                                   pattern = "\\.(neo|trans)\\.nex$",
                                   full.names = TRUE)
          chars <- do.call(cbind, lapply(nexusFiles, ReadCharacters)) |> 
            MatrixToPhyDat()
          trees <- read.tree(TreeSampleFile(pID, scriptID)) |>
            RootTree(1)
          
          parsTrees <- list.files(AnalysisDir(pID, "by_ki"), full.names = TRUE,
                                  pattern = "(project|syab).*k(1|10|Inf)\\.trees$")
          if (length(parsTrees) < 3) {
            if (scriptID == models[[1]]) {
              # If not for 1, it won't be for 2, ...
              FetchResults(pID, "by_ki")
            }
            parsTrees <- list.files(AnalysisDir(pID, "by_ki"), full.names = TRUE,
                                    pattern = "(project|syab).*k(1|10|Inf)\\.trees$")
            
            if (length(parsTrees) < 3) {
              if (scriptID == models[[1]]) {
                TreeSearch(pID)
              }
              
              stop("Parsimony trees not yet available")
            }
          }
          if (length(parsTrees) > 3) {
            # This function is not resilient to unexpected input
            stop("Unexpected parsimony results found.")
          }
          mpts <- lapply(parsTrees, function(x)
            ape::read.tree(x, keep.multi = TRUE)[[1]])
          
          ret <- data.frame(
            ew = TreeLength(trees, chars, concavity = Inf) /
              TreeLength(mpts[[3]], chars, concavity = Inf),
            k10 = TreeLength(trees, chars, concavity = 10) /
              TreeLength(mpts[[2]], chars, concavity = 10),
            k1 = TreeLength(trees, chars, concavity = 1) /
              TreeLength(mpts[[1]], chars, concavity = 1),
            pID = pID,
            scriptID = scriptID)
          saveRDS(ret, file = pFile)
          ret
        }, error = function(e) {message(pID, ", ", scriptID, ": ", e$message); NULL})
      })))
  )
}

#' Queue parsimony search on Slurm remote
#' @inheritParams MakeSlurm
#' @importFrom ssh ssh_exec_wait
#' @seealso [`MakeSlurm()`]
#' @export
TreeSearch <- function(pID, replace = FALSE) {
  session <- SshSession()
  squeue <- SlurmQueue(session)
  jobName <- paste0("ts-", pID)
  if (jobName %in% squeue[["JobName"]]) {
    if (isTRUE(replace)) {
      ssh_exec_wait(
        session,
        paste("scancel", squeue[["JobID"]][[match(jobName, squeue[["JobName"]])]])
      )
    } else {
      return(structure(FALSE, reason = "Job already in queue"))
    }
  }
  
  template <- file.path(SlurmDir(), "TreeSearch.sh")
  newLines <- gsub("%PID%", fixed = TRUE, pID,
                   gsub("%SCRIPTBASE%", fixed = TRUE,
                        ScriptBase(pID, "by_ki"), readLines(template)))
  
  remoteFile <- paste0(RemoteDir(), "/slurm/ts-", pID, ".sh")
  .shLines <- function(lines) {
    paste(gsub("$", "\\$", fixed = TRUE, lines), collapse = "\n")
  }
  ssh_exec_wait(session, paste0("cat > ", remoteFile, " <<EOF\n",
                                     .shLines(newLines), "\nEOF"))
  command <- paste("sbatch", remoteFile)
  message(command)
  
  ssh_exec_wait(session, command)
  ssh_exec_wait(session, paste("rm", remoteFile))
  TRUE
}
