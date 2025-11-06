#' @export
MemRequired <- function(pID, scriptID, ml, cores = 16L) {
  if (pID %in% .config$syab[4:5] && scriptID %in% c("hg2_ki", "hg_b_kv")) {
    .config$maxSharedMem * 0.9
  }
  ceiling(PredictResource(pID, scriptID, ml, cores)["mem", ] / 1024 / 1024) # Mb
}

#' @export
TmpRequired <- function(pID, scriptID, ml, cores = 16L) {
  PredictResource(pID, scriptID, ml, cores)["tmp", "upr"] # Mb
}

#' @export
EssTarget <- function(pID, scriptID) {
  convergence <- HasConverged(pID, scriptID)
  stats <- attr(convergence, "stats")
  margin <- 1.3
  if (convergence) {
    message(pID, "_", scriptID, " has converged; ESS target set to zero")
    0 # Complete
  } else if (length(stats[["ess"]]) == 0 || is.na(stats[["ess"]])) {
    margin * .config$essThreshold
  } else {
    margin * stats[["ess"]] * .config$essThreshold / min(stats$frechet, stats$med)
  }
}

#' Create a slurm job
#' Constructs and submits a slurm job to complete an analysis
#' @param pID Character specifying the project identifier to analyse
#' @param scriptID Character specifying the script to run
#' @param ml Logical; if true, calculate the marginal likelihood; if false,
#' run MCMCMC analysis
#' @param ess Numeric: Estimated sample size to attempt to obtain
#' @param replace Logical: should a corresponding job be cancelled in order to
#' run this job instead?
#' @importFrom ssh ssh_exec_wait
#' @export
MakeSlurm <- function(pID, scriptID, ml = FALSE, ess = EssTarget(pID, scriptID),
                      replace = FALSE, generous = TRUE) {
  if (grepl("\\d+\\-\\d+\\-\\d+", scriptID)) {
    return(structure(FALSE, reason = "Script ID not recognized"))
  }
  if (ml && pID %in% .config$syab[4:5] &&
      (scriptID %in% c("by_nt_kv", "hg_kv", "hg2_kv") ||
       grepl("ki$", scriptID, perl = TRUE))) {
    return(structure(FALSE, reason = "ML for _ki is futile"))
  }
  session <- SshSession()
  if (isFALSE(session)) {
    return(structure(FALSE, reason = attr(session, "e")))
  }
  squeue <- SlurmQueue(session)
  jobName <- sub(".sh", "", fixed = TRUE, basename(SlurmFile(pID, scriptID, ml)))
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
  
  if (ml) {
    cores <- Cores(pID, scriptID)
    time <- cores[["time"]]
    cores <- cores[["cores"]]
    if (generous && time == .config$maxTime && cores < .config$maxCores) {
      # Allow a little wriggle room to avoid near misses
      cores <- ceiling(128 / (ceiling(128 / cores) - 1))
    }
    if (is.na(cores)) {
      if (is.null(getOption("longRuns"))) {
        return(structure(FALSE, "reason" = "Unlikely to complete in time"))
      } else {
        cores <- 128L
        time <- .config$longTime
      }
    }
  } else {
    time <-
      if (length(list.files(AnalysisDir(pID, scriptID), "*_run_1.log"))) {
        .config$mc3ContinueTime # Standard time for MCMCMC continuation jobs
      } else {
        .config$maxTime
      }
    cores <- 16L
  }
  
  mem <- MemRequired(pID, scriptID, ml, cores)
  
  if (mem[["fit"]] > .config$maxMem) {
    return(structure(FALSE, "reason" = "Cannot allocate enough memory"))
  }
  
  myMem <- if (mem[["upr"]] > .config$bigSharedMem
               && mem[["fit"]] < .config$bigSharedMem) {
    message("Trying bigSharedMem")
    .config$bigSharedMem
  } else {
    min(mem[["upr"]], .config$maxMem)
  }
  
  myMem <- min(ceiling(myMem), .config$maxMem)
  
  if (ml && time > .config$longTime) {
    return(structure(FALSE, "reason" = "Run time too long"))
  }
  
  template <- SlurmTemplate(ml)
  slurmFile <- file.path(SlurmDir(), SlurmFile(pID, scriptID, ml))
  newLines <- gsub("%PID%", fixed = TRUE, pID,
                   gsub("%SCRIPTBASE%", fixed = TRUE,
                        ScriptBase(pID, scriptID), readLines(template)))
  
  # --args no longer supported in RevBayes 1.4.0, but serves as a placeholder
  # in the template file.
  rbLine <- grep(" --args \\d+", newLines)
  
  if (!ml && length(rbLine)) {
    newLines[[rbLine]] <- gsub(" --args \\d+",
                               sprintf(" %d", ceiling(ess)),
                               newLines[[rbLine]])
  }
  
  remoteFile <- file.path(RemoteDir(), SlurmFile(pID, scriptID, ml))
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
                    " --gres=tmp:", ceiling(TmpRequired(pID, scriptID, ml, cores)), "M ",
                    remoteFile
  )
  if (ml) {
    message(command)
  } else {
    message(command, "; ESS = ", ceiling(ess))
  }
  ssh_exec_wait(session, command)
  ssh_exec_wait(session, paste("rm", remoteFile))
  TRUE
}

#' @importFrom ssh ssh_exec_wait
#' @export
SimSlurm <- function(pID, scriptID, stat = TRUE, replace = FALSE) {
  if (grepl("\\d+\\-\\d+\\-\\d+", scriptID)) {
    return(structure(FALSE, reason = "Script ID not recognized"))
  }
  session <- SshSession()
  squeue <- SlurmQueue(session)
  mc3Name <- sub(".sh", "", fixed = TRUE, basename(SlurmFile(pID, scriptID, FALSE)))
  if (mc3Name %in% squeue[["JobName"]]) {
    return(structure(FALSE, reason = "MC3 job in queue"))
  }
  
  atStat <- if (isTRUE(stat)) "stat" else "obs"
  jobName <- paste0(atStat, "-", mc3Name)
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
  
  mem <- MemRequired(pID, scriptID, FALSE, 1)
  if (mem[["fit"]] > .config$maxMem) {
    return(structure(FALSE, "reason" = "Cannot allocate enough memory"))
  }
  
  myMem <- if (mem[["upr"]] > .config$bigSharedMem &&
               mem[["fit"]] < .config$bigSharedMem) {
    message("Trying bigSharedMem")
    .config$bigSharedMem
  } else {
    min(mem[["upr"]], .config$maxMem)
  }
  
  myMem <- min(ceiling(myMem), .config$maxMem)
  
  burninF <- ExistingResults(pID, scriptID)[["convergence"]][["burnin"]]
  
  template <- file.path(SlurmDir(), "pp-sim.sh")
  newLines <- gsub(
    "%PID%", fixed = TRUE, pID, gsub(
      "%SCRIPTID%", fixed = TRUE, scriptID, gsub(
        "%BURNIN%", fixed = TRUE, burninF, gsub(
          "%STAT%", fixed = TRUE, atStat, gsub(
            "%SCRIPTBASE%", fixed = TRUE,
            ScriptBase(pID, scriptID),
            readLines(template))))))
  
  slurmFile <- SlurmFile(pID, scriptID, FALSE)
  remoteFile <- file.path(RemoteDir(), "slurm",
                          paste0("ppsim_", basename(slurmFile)))
  .shLines <- function(lines) {
    paste(gsub("$", "\\$", fixed = TRUE, lines), collapse = "\n")
  }
  ssh_exec_wait(session, paste0("cat > ", remoteFile, " <<EOF\n",
                                .shLines(newLines), "\nEOF"))
  
  localScript <- paste0("rbScripts/ppsim_", scriptID, ".Rev")
  remoteScript <- file.path(RemoteDir(),
                            ScriptBase(pID, scriptID), basename(localScript))
  ssh_exec_wait(session, paste0("cat > ", remoteScript, " <<EOF\n",
                                .shLines(readLines(localScript)), "\nEOF"))
  command <- paste0("sbatch",
                    " -n 1",
                    " -p ", if (myMem > .config$maxSharedMem) "bigmem" else "shared",
                    " --mem=", myMem, "M",
                    " --time=", AsHMS(12 * 60 * 60),
                    " --gres=tmp:", ceiling(TmpRequired(pID, scriptID, FALSE, 1)), "M ",
                    remoteFile
  )
  message(command)
  ssh_exec_wait(session, command)
  ssh_exec_wait(session, paste("rm", remoteFile))
  TRUE
}

#' @importFrom ssh ssh_exec_wait
#' @export
TimeSlurm <- function(pID, scriptID, replace = FALSE) {
  if (grepl("\\d+\\-\\d+\\-\\d+", scriptID)) {
    return(structure(FALSE, reason = "Script ID not recognized"))
  }
  session <- SshSession()
  squeue <- SlurmQueue(session)
  jobName <- sub(".sh", "", fixed = TRUE, basename(SlurmFile(pID, scriptID, FALSE)))
  if (any(c(jobName, paste0("time-", jobName)) %in% squeue[["JobName"]])) {
    if (isTRUE(replace)) {
      ssh_exec_wait(
        session,
        paste("scancel", squeue[["JobID"]][[match(jobName, squeue[["JobName"]])]])
      )
    } else {
      return(structure(FALSE, reason = "Job already in queue"))
    }
  }
  
  time <- 60 * 60 * 12 # in seconds
  cores <- 1L
  
  mem <- MemRequired(pID, scriptID, FALSE, cores)
  
  if (mem[["fit"]] > .config$maxMem) {
    return(structure(FALSE, "reason" = "Cannot allocate enough memory"))
  }
  
  myMem <- if (mem[["upr"]] > .config$bigSharedMem && 
               mem[["fit"]] < .config$bigSharedMem) {
    message("Trying bigSharedMem")
    .config$bigSharedMem
  } else {
    min(mem[["upr"]], .config$maxMem)
  }
  
  slog <- SlurmLog()
  oom <- slog[slog$pID == pID & slog$scriptID == scriptID & slog$State == "OUT_OF_MEMORY", ]
  if (nrow(oom)) {
    lastMem <- as.numeric(
      gsub(" ", "", fixed = TRUE,
           sub("M", "000000", fixed = TRUE,
               sub("K", "000", fixed = TRUE,
                   oom[["ReqMem"]]))))
    myMem <- max(
      myMem,
      2 * lastMem / 1e6, # ~Mb, per hack above
      na.rm = TRUE
    )
  }
  myMem <- min(ceiling(myMem), .config$maxMem)
  
  # Write RevBayes script on remote server
  template <- "rbScripts/mc3time.Rev"
  newLines <- gsub(
    "%PID%", fixed = TRUE, pID, gsub(
      "%SCRIPT%", scriptID, fixed = TRUE, gsub(
        "%SCRIPTBASE%", fixed = TRUE, ScriptBase(pID, scriptID), gsub(
          "%MATRIXBASE.%", fixed = TRUE, basename(MatrixFile(pID, "")),
          readLines(template)))))
  
  remoteFile <- file.path(RemoteDir(), ScriptBase(pID, scriptID), 
                          "mc3time.Rev")
  .shLines <- function(lines) {
    paste(gsub("$", "\\$", fixed = TRUE, lines), collapse = "\n")
  }
  ssh::ssh_exec_wait(session, paste0("cat > ", remoteFile, " <<EOF\n",
                                     .shLines(newLines), "\nEOF"))
  
  # Write Slurm script to remote server
  template <- file.path(SlurmDir(), "mc3time.sh")
  newLines <- gsub(
    "%PID%", fixed = TRUE, pID, gsub(
      "%SCRIPTBASE%", fixed = TRUE, ScriptBase(pID, scriptID),
      readLines(template)))
  
  remoteFile <- file.path(RemoteDir(),
                          paste0("time-",
                                 basename(SlurmFile(pID, scriptID, FALSE))))
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
                    " --gres=tmp:", ceiling(TmpRequired(pID, scriptID, FALSE, cores)), "M ",
                    remoteFile
  )
  message(command)
  
  ssh_exec_wait(session, command)
  ssh_exec_wait(session, paste("rm", remoteFile))
  TRUE
}
