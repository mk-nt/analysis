#' Violin plot of posterior parameter medians with exp-scale y-axis
#'
#' Plots log-transformed posterior medians for rate_loss, rate_neo, and
#' tree_length as violins, with the y-axis labelled on the original
#' (non-log-transformed) scale.
#'
#' @param loss,neo,lng Summary matrices (rows = summary stats, cols = replicates)
#'   as returned by \code{sapply(..., summary)}.
#' @param true_vals Named numeric vector of true simulation values, in order
#'   \code{c(n, t, length)}.
#' @export
PlotParamViolin <- function(loss, neo, lng, true_vals) {
  vioplot::vioplot(
    log(loss["Median", ]), log(neo["Median", ]), log(lng["Median", ]),
    names = c("", "", ""),
    col = 5:3,
    axes = FALSE,
    xaxt = "n",
    yaxt = "n",
    frame.plot = FALSE
  )
  axis(1, at = 1:3,
       labels = expression("Gain-loss ratio",
                           italic("t") ~ "rate",
                           "tree length"),
       las = 1, lty = 0)
  log_ticks <- pretty(range(log(c(loss["Median", ], neo["Median", ],
                                  lng["Median", ]))))
  axis(2, at = log_ticks, labels = round(exp(log_ticks), 2), las = 2)
  abline(h = 0, lty = "dashed", col = "#888888")
  points(1:3, log(true_vals), pch = 95, cex = 2, col = 2, lwd = 3)
  points(1:3, c(median(log(loss["Median", ])),
                median(log(neo["Median", ])),
                median(log(lng["Median", ]))),
         pch = 20, cex = 0.8, col = "white")
  # segments(1:3 - 0.5, log(true_vals), 1:3 + 0.5,
  #          col = "white", lwd = 2, lty = "dashed")
}

#' @export
PackageFile <- function(...) {
  ret <- system.file(..., package = "neotrans")
  if (ret == "") {
    pkgRoot <- system.file(package = "neotrans")
    if (pkgRoot == "") {
      # Package is not installed
      file.path(getNamespaceInfo("neotrans", "path"), "inst", ...)
    } else {
      file.path(pkgRoot, ...)
    }
  } else {
    ret
  }
}

# Count number of binary characters with informative pattern
#' @export
NInformative <- function(...) {
  sum(colSums(apply(PhyDatToMatrix(ReadAsPhyDat(file.path(...))), 2, table) < 2) == 0)
}

#' Queue simulation for remote analysis
#' @param impute Logical: should we request imputation simulation, per simImpute.R?
#' @export
QueueSim <- function(seed, scriptID, cores = 16,
                     time = getOption("jobTime", (2 * 60 + 25) * 60),
                     myMem = 3000, myTmp = myMem * 8, replace = FALSE,
                     withT = TRUE, impute = FALSE) {
  session <- SshSession()
  if (isFALSE(session)) {
    return(structure(FALSE, reason = attr(session, "e")))
  }
  squeue <- SlurmQueue(session)
  simID <- sprintf("sim%03d", seed)
  
  if (sprintf("%s_%s_%s", basename(MkPath()), simID,
              scriptID) %in% squeue[["JobName"]]) {
    if (isTRUE(replace)) {
      ssh_exec_wait(
        session,
        paste("scancel", squeue[["JobID"]][[match(jobName, squeue[["JobName"]])]])
      )
    } else {
      return(structure(FALSE, reason = "Job already in queue"))
    }
  }
  
  slurmFile <- file.path(SlurmDir(), sprintf("%s_%s.sh", simID, scriptID))
  template <- file.path(SlurmDir(), "mc3sim.sh")
  newLines <- readLines(template)
  
  newLines <- gsub("%SIMNAME%", fixed = TRUE,
                   gsub("sim", "", fixed = TRUE,
                        paste0(basename(MkPath()), simID)),
                   gsub("%SIMDIR1%", fixed = TRUE,  basename(MkPath()),
                        gsub("%SIMDIR2%", fixed = TRUE, simID,
                             gsub("%SCRIPTID%", fixed = TRUE, scriptID,
                                  newLines))))
  
  if (!withT) {
    newLines <- gsub("sim-mc3.Rev", "sim-mc3n.Rev", fixed = TRUE, newLines)
  }
  if (isTRUE(impute)) {
    newLines <- c(
      gsub("sim-mc3", "imp-mc3", fixed = TRUE, newLines),
      "if ls *.states 1> /dev/null 2>&1; then",
      "git add ./%SCRIPTID%_run_*.states",
      "git commit -m \"Imputed states: %SIMDIR1% %SIMDIR2%: %SCRIPTID%\"",
      "git rebase main",
      "git push origin main",
      "else ",
      "  echo \"No .states files found. No output to commit.\"",
      "fi"
    )
  }
  
  remoteFile <- file.path(RemoteDir(), basename(slurmFile))
  
  .shLines <- function(lines) {
    paste(gsub("$", "\\$", fixed = TRUE, lines), collapse = "\n")
  }
  ssh_exec_wait(session, paste0("cat > ", remoteFile, " <<EOF\n",
                                .shLines(newLines), "\nEOF"))
  
  command <- paste0("sbatch",
                    " -n ", cores,
                    " -p ", if (time > .config$maxTime) "long" else
                      if (myMem > .config$maxSharedMem) "bigmem" else "shared",
                    " --mem=", myMem, "M",
                    " --time=", AsHMS(time),
                    " --gres=tmp:", myTmp, "M ",
                    remoteFile
  )
  message(command)
  ssh_exec_wait(session, command)
  ssh_exec_wait(session, paste("rm", remoteFile))
  TRUE
}

#' @export
SimTrees <- function(simID, scriptID, run = 1) {
  trFile <- MkPath(simID, sprintf("%s_run_%d.trees", scriptID, run))
  if (!file.exists(trFile)) {
    tarFile <- sub(".trees", ".tar.gz", trFile)
    
    if (file.exists(tarFile)) {
      # Extract the file; exdir ensures it lands in the right folder
      untar(tarFile, exdir = dirname(trFile))
    } else {
      warning("Error: Neither .trees nor .tar.gz found at ", trFile, immediate. = TRUE)
      if (isTRUE(getOption("requeue", FALSE))) {
        QueueSim(as.numeric(substr(simID, 4, 6)), scriptID,
                 myMem = if(scriptID == "sp_kv") 2888 else 3987)
      }
      return(NULL)
    }
  }
  if (!file.exists(trFile)) {
    warning("Error: tree could not be extracted for ", simID, ": ", scriptID)
    return(NULL)
  } else {
    read.tree(trFile)
  }
}

#' Fetch log file from simulation study
#' 
#' set `options("requeue" = TRUE)` to automatically re-execute any inference
#' for which a log file is missing
#' @export
FetchLogIfMissing <- function(seed, logFile) {
  simID <- sprintf("sim%03d", seed)
  localPath <- MkPath(simID, logFile)
  if (!file.exists(localPath)) {
    remotePath <- file.path(RemoteDir(), "neotrans", "inst", "simEmpirical",
                            simID, logFile)
    ok <- tryCatch({
      scp_download(SshSession(), files = remotePath, to = MkPath(simID))
      file.exists(localPath)
    }, error = function(e) FALSE)
    if (!ok) {
      remoteSimDir <- file.path(RemoteDir(), "neotrans", "inst", "simEmpirical",
                                simID)
      if (isTRUE(getOption("requeue", FALSE))) {
        ssh_exec_wait(SshSession(),
                      paste0("rm -f ", remoteSimDir, "/sp_nt_kv_run_*.tar.gz"))
        QueueSim(seed, "sp_nt_kv", myMem = 3000)
      }
      return(NULL)
    }
  }
  read.table(localPath, header = TRUE)
}
