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
    margin * max(
      .config$essThreshold,
      stats[["ess"]] * .config$essThreshold / min(stats$frechet, stats$med)
    )
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
#' @param dependency Optional character giving a SLURM dependency clause to be
#' passed via `--dependency=` (e.g. `"afterany:17218353"`).  When set, the
#' "Job already in queue" guard is skipped, since a chain of jobs deliberately
#' shares a job name.  See `man sbatch` for the full grammar; this package's
#' [`MakeSlurmChain()`] supplies `afterany:<prev>` to chain checkpoint-resuming
#' runs end-to-end.
#' @param timeOverride Optional integer specifying wall-clock seconds to request
#' for this job, bypassing the prediction returned by [`Cores()`].  Used by
#' [`MakeSlurmChain()`] to slice a long analysis into uniform chunks.
#' @return On success, an integer SLURM JobID (truthy in logical contexts).
#' On failure, `FALSE` with a `"reason"` attribute.
#' @importFrom ssh ssh_exec_wait ssh_exec_internal
#' @export
MakeSlurm <- function(pID, scriptID, ml = FALSE, ess = EssTarget(pID, scriptID),
                      replace = FALSE, generous = TRUE,
                      dependency = NULL, timeOverride = NULL) {
  if (grepl("\\d+\\-\\d+\\-\\d+", scriptID)) {
    return(structure(FALSE, reason = "Script ID not recognized"))
  }
  session <- SshSession()
  if (isFALSE(session)) {
    return(structure(FALSE, reason = attr(session, "e")))
  }
  squeue <- SlurmQueue(session)
  jobName <- sub(".sh", "", fixed = TRUE, basename(SlurmFile(pID, scriptID, ml)))
  if (jobName %in% squeue[["JobName"]] && is.null(dependency)) {
    # Chained chunks legitimately share a jobName; only guard when not chaining.
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
    if (is.null(timeOverride)) {
      cores <- Cores(pID, scriptID)
      time <- cores[["time"]]
      cores <- cores[["cores"]]
      if (generous && time == .config$maxTime && cores < .config$maxCores) {
        # Allow a little wriggle room to avoid near misses
        cores <- ceiling(.config$maxCores / (ceiling(.config$maxCores / cores) - 1))
      }
      if (is.na(cores)) {
        cores <- .config$maxCores
        time <- .config$longTime
      }
    } else {
      cores <- .config$maxCores
      time <- timeOverride
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

  if (!is.null(timeOverride)) {
    time <- as.integer(timeOverride)
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

  template <- SlurmTemplate(ml)
  slurmFile <- file.path(SlurmDir(), SlurmFile(pID, scriptID, ml))
  newLines <- gsub("%PID%", fixed = TRUE, pID,
                   gsub("%SCRIPTBASE%", fixed = TRUE,
                        ScriptBase(pID, scriptID), readLines(template)))

  # RevBayes 1.4.0 takes positional args after the script name; mcmcmc.Rev
  # reads minEss as args[1] and runHours as args[2].
  rbLine <- grep("mcmcmc\\.Rev \\d+ \\d+", newLines)

  if (!ml && length(rbLine)) {
    newLines[[rbLine]] <- sub("mcmcmc\\.Rev \\d+",
                              sprintf("mcmcmc.Rev %d", ceiling(ess)),
                              newLines[[rbLine]])
  }

  remoteFile <- file.path(RemoteDir(), SlurmFile(pID, scriptID, ml))
  .shLines <- function(lines) {
    paste(gsub("$", "\\$", fixed = TRUE, lines), collapse = "\n")
  }
  ssh_exec_wait(session, paste0("cat > ", remoteFile, " <<EOF\n",
                                     .shLines(newLines), "\nEOF"))
  # --kill-on-invalid-dep=yes prevents PENDING (DependencyNeverSatisfied)
  # orphans if the predecessor never reaches a terminal state for some reason.
  dependencyClause <- if (is.null(dependency)) {
    ""
  } else {
    paste0(" --dependency=", dependency, " --kill-on-invalid-dep=yes")
  }
  command <- paste0("sbatch",
                    " -n ", cores,
                    " -p ", if (time > .config$maxTime) "long" else
                      if (myMem > .config$maxSharedMem) "bigmem" else "shared",
                    " --mem=", myMem, "M",
                    " --time=", AsHMS(time),
                    " --gres=tmp:", ceiling(TmpRequired(pID, scriptID, ml, cores)), "M",
                    dependencyClause,
                    " ",
                    remoteFile
  )
  if (ml) {
    message(command)
  } else {
    message(command, "; ESS = ", ceiling(ess))
  }
  result <- ssh_exec_internal(session, command, error = FALSE)
  ssh_exec_wait(session, paste("rm", remoteFile))
  if (result[["status"]] != 0) {
    return(structure(FALSE, reason = paste(
      "sbatch failed:", rawToChar(result[["stderr"]]))))
  }
  stdout <- rawToChar(result[["stdout"]])
  # sbatch may emit warnings ahead of the success line; extract the JobID
  # from anywhere in stdout rather than anchoring to start-of-string.
  jobIDMatch <- regmatches(stdout,
                           regexpr("Submitted batch job [0-9]+", stdout))
  if (length(jobIDMatch) == 0) {
    return(structure(FALSE, reason = paste(
      "Could not parse JobID from sbatch output:", stdout)))
  }
  # Return:
  as.integer(sub("Submitted batch job ", "", jobIDMatch, fixed = TRUE))
}

#' Submit a chain of checkpoint-resuming SLURM jobs
#'
#' Splits a long marginal-likelihood analysis into a series of shorter chunks
#' that resume from the checkpoint files written by `marginal.Rev`
#' (revbayes#991).  Each chunk depends on the previous via SLURM's
#' `--dependency=afterany:` clause, so they run strictly in series: chunk K+1
#' is held `PENDING (Dependency)` until chunk K reaches a terminal state
#' (COMPLETED, TIMEOUT, FAILED, CANCELLED, NODE_FAIL), then enters normal
#' scheduling.  At most one chunk in the chain is ever `RUNNING` at a time.
#'
#' This is useful when the projected wall-clock time exceeds the cluster's
#' single-job maximum (`.config$maxTime`) but the analysis supports
#' checkpoint resumption.  Shorter chunks also benefit from SLURM backfill,
#' typically starting sooner than a single long job.
#'
#' @inheritParams MakeSlurm
#' @param hoursPerChunk Numeric: wall-clock hours requested for each chunk.
#' Defaults to 24 — a compromise between backfill responsiveness and the
#' overhead of repeated queue waits and checkpoint reloads.
#' @param chunks Integer: number of chunks to submit.  If `NULL` (default),
#' chosen automatically from [`PredictResource()`]'s upper-bound projection
#' plus one chunk of margin, capped at 14.
#' @param margin Integer: extra chunks to append beyond the prediction.
#' Defaults to 1; ignored when `chunks` is supplied.
#' @return Invisibly returns an integer vector of submitted JobIDs.
#' On failure, `FALSE` with a `"reason"` attribute.
#' @seealso [`MakeSlurm()`]; [`SCancel()`] to terminate a chain mid-flight
#' (cancelling chunk K satisfies `afterany`, so chunk K+1 will start next).
#' @importFrom ssh ssh_exec_wait
#' @export
MakeSlurmChain <- function(pID, scriptID, ml = TRUE,
                           hoursPerChunk = 24, chunks = NULL,
                           margin = 1L, replace = FALSE) {
  if (grepl("\\d+\\-\\d+\\-\\d+", scriptID)) {
    return(structure(FALSE, reason = "Script ID not recognized"))
  }
  session <- SshSession()
  if (isFALSE(session)) {
    return(structure(FALSE, reason = attr(session, "e")))
  }

  # Cross-chain guard: a chain implies multiple jobs share a name, but no
  # two chains should ever coexist for one (pID, scriptID) — they would race
  # on checkpoint files and git pushes.
  jobName <- sub(".sh", "", fixed = TRUE, basename(SlurmFile(pID, scriptID, ml)))
  squeue <- SlurmQueue(session)
  existing <- squeue[squeue[["JobName"]] == jobName, "JobID"]
  if (length(existing) > 0) {
    if (isTRUE(replace)) {
      for (id in existing) {
        ssh_exec_wait(session, paste("scancel", id))
      }
    } else {
      return(structure(FALSE, reason = sprintf(
        "%d existing job(s) for %s; pass replace = TRUE to cancel",
        length(existing), jobName)))
    }
  }

  if (is.null(chunks)) {
    res <- PredictResource(pID, scriptID, ml = ml, cores = .config$maxCores)
    upperHours <- res["time", "upr"] / (60 * 60)
    chunks <- ceiling(upperHours / hoursPerChunk) + as.integer(margin)
    chunks <- max(2L, min(as.integer(chunks), 14L))
    message(jobName, ": estimated ", chunks, " chunk(s) of ", hoursPerChunk,
            " h (projected upper bound ", signif(upperHours, 3), " h)")
  } else {
    chunks <- as.integer(chunks)
  }

  timeSec <- as.integer(hoursPerChunk * 60 * 60)
  jobIDs <- integer(chunks)

  first <- MakeSlurm(pID, scriptID, ml = ml, timeOverride = timeSec)
  if (!is.numeric(first)) {
    return(structure(FALSE, reason = paste(
      "First chunk could not be submitted:", attr(first, "reason"))))
  }
  jobIDs[[1]] <- as.integer(first)

  for (k in seq_len(chunks - 1L) + 1L) {
    next_ <- MakeSlurm(pID, scriptID, ml = ml,
                       timeOverride = timeSec,
                       dependency = paste0("afterany:", jobIDs[[k - 1L]]))
    if (!is.numeric(next_)) {
      warning(sprintf("Chain truncated at chunk %d/%d: %s",
                      k - 1L, chunks, attr(next_, "reason")),
              immediate. = TRUE)
      break
    }
    jobIDs[[k]] <- as.integer(next_)
  }

  submitted <- jobIDs[jobIDs > 0]
  message("Submitted ", length(submitted), " chained job(s) for ", jobName,
          ": ", paste(submitted, collapse = " -> "))
  invisible(submitted)
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
  
  localScript <- file.path(dirname(SlurmDir()), "rbScripts",
                           paste0("ppsim_", scriptID, ".Rev"))
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
