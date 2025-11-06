#' Connect to remote server via SSH
#' @importFrom ssh ssh_connect
#' @export
ConnectSSH <- function() {
  tryCatch(
    ssh_connect(Sys.getenv("sshLogin"), keyfile = Sys.getenv("sshKey"),
                passwd = Sys.getenv("sshPass")),
    error = function(e) tryCatch(
      ssh_connect(Sys.getenv("sshLogin")),
      error = function(e) structure(structure(FALSE, e = e))
    )
  )
}

.SlurmTime <- function(dt) format(dt, "%Y-%m-%dT%H:%M")

ReadSacct <- function(lines) {
  widths <- nchar(strsplit(lines[[2]], " ")[[1]]) + 1L
  margins <- rbind(1 + c(0, cumsum(widths)[-length(widths)]), cumsum(widths))
  ret <- apply(margins, 2, function(lim) {
    gsub("^\\s+|\\s+$", "", substr(lines[-2], lim[[1]], lim[[2]]))
    })
  if (length(lines) < 3) {
    return(matrix(NA, 0, length(ret), dimnames = list(NULL, ret)))
  }
  dimnames(ret) <- list(NULL, ret[1, ])
  ret[-1, ]
}

Sacct <- function(slurmID = NULL, state = NULL, time = "now-24hours") {
  session <- SshSession()
  result <- ssh_exec_internal(
    session,
    paste("sacct", 
          if (!is.null(slurmID)) sprintf("-j %s", slurmID),
          "--format=\"JobID%20,JobName%24,State%16,Elapsed,End,AllocCPUS,ReqMem,MaxRSS,NodeList%12\"",
          if (!is.null(time)) {
            if (length(time) < 2) time[[2]] <- "now"
            paste("--starttime", time[[1]], "--endtime", time[[2]])
          },
          if (!is.null(state)) sprintf("--state %s", state))
  )
  
  if (result$status != 0) {
    stop("Error in sacct command: ", result)
  }
  ReadSacct(strsplit(rawToChar(result$stdout), "\n")[[1]])
}

.Filter3 <- function(lines) {
  lines <- lines[!grepl("\\d\\.\\d$", lines[, 1], ), ]
  nLines <- dim(lines)[[1]]
  ret <- lines[seq_len(nLines) %% 3 == 1, , drop = FALSE]
  ret[, "MaxRSS"] <- sub("K$", "000", # strictly we should * 1024
                          sub("M$", "000000", 
                              lines[seq_len(nLines) %% 3 == 2, "MaxRSS", 
                                    drop = FALSE]))
  ret <- ret[!startsWith(ret[, "JobName"], "ts-"), , drop = FALSE]
  data.frame(
    ret,
    "pID" = sub("^\\D*(\\d+).*$", "\\1", ret[, "JobName"], perl = TRUE),
    "scriptID" = sub("^\\D*\\d+_", "", ret[, "JobName"], perl = TRUE),
    "task" = ifelse(substr(ret[, "JobName"], 1, 3) == "ml-", "marginal", "mcmcmc")
  )
}

.HMSToS <- function(hms) {
  vapply(strsplit(hms, ":", fixed = TRUE), function(x) {
    sum(as.numeric(x) * c(60 * 60, 60, 1))
  }, numeric(1))
}

.ToSeconds <- function(times) {
  sPerDay <- 24 * 60 * 60
  ifelse(grepl("-", times, fixed = TRUE),
         as.integer(substr(times, 1, 1)) * sPerDay, 0) +
    .HMSToS(gsub("\\d+\\-", "", times))
}

#' Process analyses that have completed on remote server
#' @param days Numeric specifying how many days since last collection to
#' retrieve
#' @return `Collect()` is called for its side effects of pulling remote results
#' from the server, and updating the local cache accordingly.
#' @export
Collect <- function(days = 1) {
  sacct.csv <- file.path(OutputDir(), "sacct.csv")
  database <- if (file.exists(sacct.csv)) {
    read.csv(sacct.csv, stringsAsFactors = FALSE, 
             colClasses = c("pID" = "character", "scriptID" = "character"))
  } else {
    data.frame(
      task = character(),
      pID = character(),
      scriptID = character(),
      jobID = integer(),
      State = character(),
      endTime = numeric(),
      seconds = numeric(),
      nCPU = factor(numeric(0), levels = 2 ^ (0:7), ordered = TRUE),
      ReqMem = numeric(),
      maxRSS = numeric(),
      tmpUsed = numeric()
    )
  }
  
  startTime <- if (nrow(database) > 0) {
    as.POSIXct(database[["endTime"]][[1]], origin = "1970-01-01")
  } else {
    # Specify the first date from which to build the database
    as.POSIXct("2024-11-11T00:00", format = "%Y-%m-%dT%H:%M")
  }
  sPerDay <- 24 * 60 * 60
  window <- format(startTime + c(1, days * sPerDay),  "%Y-%m-%dT%H:%M:%S")
  
  toDB <- .Filter3(Sacct(state = "TO,OOM,CD", time = window))
  if (nrow(toDB) == 0) {
    if (window[[2]] < as.POSIXct(Sys.Date())) {
      message("No entries; extending period \U1F551")
      return(Collect(days + 1))
    } else {
      message("\U1F4A4 Nothing new to process")
      return(invisible(toDB))
    }
  }
  toDB <- cbind(toDB, "tmpUsed" = NA_real_)
  completed <- which(toDB[["State"]] == "COMPLETED")

  updated <- list()
  for (i in completed) {
    row <- toDB[i, ]
    message(row$JobID, ": ", row$JobName)
    # Check that job ran for > 60 s
    if (.ToSeconds(row$Elapsed) < 60) {
      if (row$task == "marginal") {
        if (tryCatch(
          ssh_exec_internal(
            .ssh$session,
            paste0("grep -H -E ", shQuote(".out.pp[[:space:]]+already exists"), 
                   " ", RemoteDir(), "neotrans/", row$JobName, ".out", sep = "")
          )$status, error = function(e) 99) == 0) {
          if (is.null(updated[[row$JobName]])) {
            if (!file.exists(StoneFile(row$pID, row$scriptID))) {
              UpdateRecords(row$pID, row$scriptID, searchRemote = TRUE)
            }
            if (!file.exists(StoneFile(row$pID, row$scriptID))) {
              stop("File.out.pp exists on remote but not updated locally: ", row$JobName)
            }
            nOnFile <- read.table(StoneFile(row$pID, row$scriptID))[[5]]
            if (nOnFile < 128) {
              warning(row$JobID, ": ", row$JobName, " out.pp exists; updating Rev scripts",
                      immediate. = TRUE)
              updated[[row$JobName]] <- RevBayes(row$pID, row$scriptID)
              MakeSlurm(row$pID, row$scriptID, ml = row$task == "marginal")
            } else {
              message(row$JobID, ": ", row$JobName, " already done")
            }
          }
        } else {
          warning(row$JobID, ": ", row$JobName, " succeeded in < 60 s")
        }
      } else {
        if (!UpdateRecords(row$pID, row$scriptID, searchRemote = TRUE)) {
          warning(row$JobID, ": ", row$JobName, " succeeded in < 60 s")
        }
      }
      next
    }
    toDB[i, "tmpUsed"] <- TmpUsed(row$pID, row$scriptID, row$task == "marginal")
    if (is.null(updated[[row$JobName]])) {
      updated[[row$JobName]] <- UpdateRecords(row$pID, row$scriptID,
                                              searchRemote = TRUE)
    }
  }
  
  oom <- toDB[toDB[["State"]] == "OUT_OF_MEMORY", ]
  # On exit, as we first need to update the database with the oom lines
  if (nrow(oom)) {
    on.exit(
      add = TRUE,
      apply(oom, 1, function(row) {
        message(row[["JobID"]], ": ", row[["JobName"]], " needs more memory")
        MakeSlurm(row[["pID"]], row[["scriptID"]], ml = row[["task"]] == "marginal")
      })
    )
  }
  
  timeout <- toDB[toDB[["State"]] == "TIMEOUT", ]
  # On exit, as we first need to update the database with the timeout lines
  if (nrow(timeout)) {
    on.exit(
      add = TRUE,
      apply(timeout, 1, function(row) {
        if (!UpdateRecords(row[["pID"]], row[["scriptID"]], searchRemote = TRUE,
                           makeSlurm = TRUE)) {
          message(row[["JobID"]], ": ", row[["JobName"]], " needs more time")
        }
      })
    )
  }
  
  toDB <- toDB[order(toDB[, "End"], decreasing = TRUE), ]
  newEntries <- data.frame(
    task = as.character(toDB$task),
    pID = as.character(toDB$pID),
    scriptID = as.character(toDB$scriptID),
    jobID = as.integer(toDB$JobID),
    State = as.character(toDB$State),
    endTime = as.numeric(as.POSIXct(toDB$End, format = "%Y-%m-%dT%H:%M:%S")),
    seconds = as.character(.ToSeconds(toDB$Elapsed)),
    nCPU = toDB$AllocCPUS,
    ReqMem = toDB$ReqMem,
    maxRSS = toDB$MaxRSS,
    tmpUsed = toDB$tmpUsed
  )
  if (nrow(newEntries) == 0 || 
      (nrow(newEntries) == 1 && newEntries[["endTime"]] == startTime)) {
    if (window[[2]] > as.POSIXct(Sys.Date())) {
      message("\U2713 Nothing new to process")
    } else {
      return(Collect(days + 1))
    }
  } else {
    combined <- rbind(newEntries, database)
    write.csv(combined[!duplicated(combined[c("task", "pID", "scriptID", "State")]), ],
              sacct.csv, row.names = FALSE, quote = FALSE)
  }
  
  message("Fetched ", nrow(newEntries), " new entries from ",
          as.POSIXct(window[[1]]), " to ", as.POSIXct(window[[2]]))
  if (window[[2]] > as.POSIXct(Sys.Date())) {
    on.exit(add = TRUE, message("\U2713 Up to date! \U1F600"))
  }
  invisible(newEntries)
}

#' Read log of completed SLURM jobs
#' @return `SlurmLog()` returns a data.frame with details of resource allocation
#' to previous slurm jobs.
#' @export
SlurmLog <- function() {
  df <- read.csv(file.path(OutputDir(), "sacct.csv"), stringsAsFactors = FALSE,
                 colClasses = c("task" = "factor",
                                "pID" = "factor", "scriptID" = "factor",
                                "jobID" = "numeric"))
  df[["inf"]] <- as.factor(ifelse(grepl("kv$", df[["scriptID"]]), "kv", "ki"))
  df[["model"]] <- as.factor(gsub("_k.$", "", df[["scriptID"]]))
  df
}

Timeouts <- function() {
  df <- SlurmLog()
  df <- df[!is.na(df$pID) & df$task == "marginal",
           setdiff(colnames(df), "task")]
  df$job <- paste0(df$pID, "_", df$scriptID)
  df <- df[order(df$endTime), ]
  to <- df[df$State == "TIMEOUT", ]
  cd <- df[df$State == "COMPLETED", ]
  to[!to$job %in% cd$job, ]
}

#' Read current slurm queue on remote host
#' @returns `SlurmQueue()` returns a data.frame containing details of queued
#' jobs.
#' @export
SlurmQueue <- function(session = SshSession()) {
  result <- ssh_exec_internal(
    session,
    "squeue -u $USER --format=\"%.24i %.24j %.4t %.24M\""
  )
  if (result[["status"]] != 0) {
    stop("Error in squeue command: ", result)
  }
  txt <- strsplit(rawToChar(result[["stdout"]]), "\n", fixed = TRUE)[[1]][-1]
  colNames <- c("JobID", "JobName", "status", "elapsed")
  if (length(txt) == 0) {
    as.data.frame(
      matrix(NA, 0, length(colNames), dimnames = list(NULL, colNames))
    )
  } else {
    ret <- as.data.frame(do.call(rbind, strsplit(txt, "\\s+"))[, -1, drop = FALSE])
    colnames(ret) <- colNames
    ret[["pID"]] <- sub("^\\D*(\\d+).*$", "\\1", ret[["JobName"]], perl = TRUE)
    ret[["scriptID"]] <- sub("^\\D*\\d+_", "", ret[["JobName"]], perl = TRUE)
    ret[order(ret[["JobID"]]), ]
  }
}

#' Previous attempts to run an analysis
#' @inheritParams MakeSlurm
#' @export
PastRuns <- function(pID, scriptID) {
  df <- SlurmLog()
  df[df$pID %in% pID & df$scriptID %in% scriptID, ]
}
LastRun <- PastRun <- PastRuns

#' View tail of job output
#' 
#' Displays the last `lines` lines of the output file produced by a slurm task.
#' 
#' @inheritParams MakeSlurm
#' @export
Peek <- function(pID, scriptID, ml = TRUE, lines = 20) {
  logFile <- sprintf(file.path(RemoteDir(), "neotrans", "%s%s_%s.out"),
                     if(ml) "ml-" else "", pID, scriptID)
  session <- SshSession()
  if (ssh_exec_internal(
    session,
    paste("test -f", logFile, "&& echo 1 || echo 0"))[["stdout"]][[1]] |> 
    rawToChar() == 0) {
    NA_character_
  } else {
    ssh_exec_internal(session, paste("tail -n", as.integer(lines),
                                     shQuote(logFile)))$stdout |>
      rawToChar() |> message()
  }
}

#' squeue
#' Output the current slurm queue on the remote host
#' @param k Sort criterion for job list
#' 
#' @importFrom ssh ssh_exec_wait
#' @export
SQ <- function(k = NULL) {
  session <- SshSession()
  if (is.null(k)) {
    ssh_exec_wait(session,
                  paste0("bash ", RemoteDir(), "/neotrans/squeue.sh"))
  } else if (is.numeric(k)) {
    ssh_exec_wait(
      session,
      paste("bash", file.path(RemoteDir(), "neotrans/squeue.sh"),
            "| sort -k", k)
    )
  } else {
    # i j P t M e C m R
    ssh_exec_wait(session,
                  paste("bash",
                        file.path(RemoteDir(), "neotrans/squeue.sh"),
                        paste0(k, ",")))
  }
}

#' scancel
#' 
#' Cancel a queued job on the remote host.
#' 
#' @param jobID Integer: Slurm job identifier
#' @importFrom ssh ssh_exec_wait
#' @export
SCancel <- function(jobID) {
  ssh_exec_wait(SshSession(), paste("scancel", jobID))
}

#' Predict number of cores required
#' @export
Cores <- function(pID, scriptID, cores = 8L) {
  if (pID %in% c("3448")) {
    return(c(cores = 64L, time = .config$maxTime))
  }
  if (pID %in% c("4817")) {
    # Too slow even for the 'long' service
    return(c(cores = 512L, time = .config$maxTime))
  }
  if (pID %in% c("3708", "3345", "3670")) {
    # Known to require the 'long' service
    return(c(cores = .config$maxCores, time = .config$longTime))
  }
  time <- PredictResource(pID, scriptID, ml = TRUE, cores)["time", ]
  if (time[["upr"]] < .config$maxTime * 
      if (cores >= 32) 1 else if (cores == 8L) 0.25 else 0.5) {
    request <- time[["upr"]]
    if (request > .config$maxTime * 0.6) {
      # Annoying to time out on a long run that could have completed
      request <- .config$maxTime
    }
    c(cores = cores, time = request)
  } else {
    if (cores >= .config$maxCores) {
      if (time[["fit"]] < .config$maxTime) {
        message(pID, "_", scriptID, ": > 50% chance of success with max cores")
        return(c(cores = .config$maxCores, time = .config$maxTime))
      } else if (time[["lwr"]] < .config$maxTime) { 
        message(pID, "_", scriptID, ": < 50% chance of success with max cores \U1F641")
        return(c(cores = .config$maxCores, time = .config$maxTime))
      } else {
        warning(pID, "_", scriptID, ": < 1% chance of success with max cores \U1F63F")
        return(c(cores = NA_integer_, time = Inf))
      }
    }
    proposal <- ceiling(cores * time[["lwr"]] / .config$maxTime)
    proposal <- ceiling(128 / min(ceiling(128 / proposal),
                                  floor(128 / (cores + 1))))
    Cores(pID, scriptID, min(proposal, .config$maxCores))
  }
}

TmpUsed <- function(pID, scriptID, ml = TRUE) {
  logFile <- sprintf(file.path(RemoteDir(), "%s_%s", "%s-tmpdir_usage.log"),
                     pID, scriptID, if(ml) "ml" else "mc3")
  session <- SshSession()
  if (ssh_exec_internal(
    session,
    paste("test -f", logFile, "&& echo 1 || echo 0"))[["stdout"]][[1]] |> 
    rawToChar() == 0) {
    NA_real_
  } else {
    tmpUsed <- ssh_exec_internal(session, paste("cat", shQuote(logFile)))[["stdout"]] |> 
      rawToChar() |> strsplit("\t")
    tmpUsed <- as.numeric(sub("G", "e9", fixed = TRUE,
                              sub("M", "e6", fixed = TRUE,
                                sub("K", "e3", fixed = TRUE, tmpUsed[[1]][[1]]))))
    tmpUsed / 1000 / 1000 # M
  }
}

PredictResource <- function(pID, scriptID, ml = TRUE, cores = 16,
                            verbose = FALSE) {
  df <- SlurmLog()
  df <- df[!is.na(df$pID), ]
  run <- data.frame(pID = as.factor(pID),
                    model = as.factor(sub("_k.$", "", scriptID)),
                    inf = as.factor(sub(".*(k.)$", "\\1", scriptID)),
                    nCPU = cores,
                    task = ifelse(ml, "marginal", "mcmcmc"))
  .Predict <- function(lm, mem = NA_real_) {
    predict(lm, newdata = c(run, "maxRSS" = mem),
            interval = "confidence", level = 0.99) 
  }
  
  tries <- df[which(df[["pID"]] == pID & df[["scriptID"]] == scriptID & 
                df[["task"]] == ifelse(ml, "marginal", "mcmcmc")), ]
  if (ml) {
    # `which` removes NA values
    timeouts <- tries[tries[["State"]] == "TIMEOUT", ]
    if (nrow(timeouts) > 0) {
      cpuSeconds <- mean(timeouts[["seconds"]] * timeouts[["nCPU"]])
      outFile <- sprintf(file.path(RemoteDir(), "neotrans", "%s%s_%s.out"),
                         if(ml) "ml-" else "", pID, scriptID)
      outExists <- tryCatch(
        ssh_exec_internal(SshSession(), paste("test -f", outFile))$status == 0,
        error = function(e) FALSE)
      
      if (outExists) {
        # From rbScripts/marginal.Rev
        #TODO Bonus points for reading rather than hard-coding these
        burninGen <- 12000
        runGen <- 2000 * 2 # x2 because each run includes a preBurnin of nGen
        steps <- 128
            
        # Download the file content
        outLines <- ssh_exec_internal(.ssh$session, 
                                      paste("cat", outFile))$stdout |>
          rawToChar() |>
          strsplit("\n")
        outLines <- outLines[[1]]
        progressLine <- which(outLines == "Progress:")
        if (length(progressLine) == 1 && length(outLines) >= progressLine + 2) {
          progressLine <- outLines[[progressLine + 2]]
          nStars <- nchar(progressLine)
          gensDone <- if (nStars < 68) {
            burninGen * nStars / 68
          } else {
            stepExp <- "Steps\\s+\\d+\\-\\-(\\d+)\\s+/\\s+(\\d+)\\t\\t(\\**)"
            stepsLines <- grepv(stepExp, outLines)
            if (length(stepsLines) == 0) {
              warning("Could not read ss outfile for ", pID, "_", scriptID)
              # Likely because last run was from an old marginal.Rev: update
              RevBayes(pID, scriptID)
              gensDone <- Inf
            } else {
              stepsParts <- do.call(cbind,
                                    regmatches(stepsLines,
                                               gregexec(stepExp, stepsLines)))
              starsPerChunk <- 40
              nChunks <- as.integer(stepsParts[3, 1]) / as.integer(stepsParts[2, 1])
              starsNeeded <- starsPerChunk * nChunks
              starsFound <- sum(nchar(stepsParts[4, ]))
              if (starsFound %% 40 == 0) {
                # Take an optimistic view that many runs were almost complete
                # but hadn't printed a * to the output yet
                starsFound <- starsFound + 35
              }
              burninGen + ((runGen * steps / timeouts$nCPU) * starsFound / starsNeeded)
            }
          }
          # Ceiling because in the final stage, some cores may sit idle
          stages <- ceiling(steps / cores)
          gensNeeded <- burninGen + (runGen * stages)
          sPerGen <- timeouts$seconds / gensDone
          tooShort <- sPerGen * gensNeeded
          message(pID, "_", scriptID, " timed out at ",
                  round(gensDone / gensNeeded * 100), "% after ",
                  signif(timeouts$seconds / 60 / 60, 3), " h with ", timeouts$nCPU,
                  " cores; need >", round(tooShort / 60 / 60), " h with ", cores)
        } else {
          tooShort <- cpuSeconds / cores
        }
      } else {
        tooShort <- timeouts$seconds #TODO explore why; Can we do better?
      }
    } else {
      tooShort <- 0
    }
  } else {
    # mcmcmc; can continue from last timeout
    tooShort <- 0
  }
  
  ooms <- tries[tries[["State"]] == "OUT_OF_MEMORY", ]
  tooWee <- if (nrow(ooms) > 0 && any(!is.na(ooms[["ReqMem"]]))) {
    reqMem <- as.numeric(sub("K", "e3", fixed = TRUE,
                             sub("M", "e6", fixed = TRUE,
                                 sub("G", "e9", fixed = TRUE, ooms[["ReqMem"]]))))
    max(reqMem, na.rm = TRUE)
  } else {
    0
  }
  
  
  relevant <- df[df[["State"]] == "COMPLETED" &
                   df[["seconds"]] > 180 # Probably not an aborted run
                 , ]
  if (pID %in% relevant$pID) {
    mem <- tryCatch(
      .Predict(lm(maxRSS ~ pID + model + inf + task + nCPU, data = relevant)),
      error = function(e) .Predict(lm(maxRSS ~ pID + task + nCPU, data = relevant))
    )
    # Seems to produce more sensible estimates than predicting coreSeconds and
    # multiplying by nCPU - which leads to more cores predicting MORE time...
    wallTime <- tryCatch(
      .Predict(lm(log(seconds) ~ pID + model + inf + task + nCPU, data = relevant)),
      error = function(e) .Predict(lm(log(seconds) ~ pID + task + nCPU, data = relevant))
    )
    .PredictTmp4 <- function(e) {
      .Predict(lm(tmpUsed ~ task + maxRSS, data = relevant), mem[1, "fit"])
    }
    .PredictTmp3 <- function(e) tryCatch(
      .Predict(lm(tmpUsed ~ model + task + maxRSS, data = relevant), mem[1, "fit"]),
      warning = .PredictTmp4,
      error = .PredictTmp4
    )
    .PredictTmp2 <- function(e) tryCatch(
      .Predict(lm(tmpUsed ~ pID + task + maxRSS, data = relevant), mem[1, "fit"]),
      warning = .PredictTmp3,
      error = .PredictTmp3
    )
    tmp <- tryCatch(
      .Predict(lm(tmpUsed ~ pID + model + inf + task + maxRSS, data = relevant),
               mem[1, "fit"]),
      warning = .PredictTmp2,
      error = .PredictTmp2
    )
  } else {
    mem <- .Predict(lm(maxRSS ~ model + inf + task + nCPU, data = relevant))
    wallTime <- .Predict(lm(log(seconds) ~ model + inf + task + nCPU, data = relevant))
    tmp <- tryCatch(
      .Predict(lm(tmpUsed ~ model + inf + task + maxRSS, data = relevant),
               mem[1, "fit"]),
      error = function(e) {
        rbind(setNames(quantile(relevant[["tmpUsed"]], c(0.5, 0.01, 0.99),
                                na.rm = TRUE), c("fit", "lwr", "upr")))
      })
  }
  
  RoundUpMem <- function(x) {
    proposal <- x * c(1.5, 1, 1.8) # 1.8 takes ages to get to
      # 1000 not 1024 because we might have extracted by replacing M with 000
    if (proposal[[2]] < .config$maxSharedMem * 1000 * 1000) {
      proposal <- pmin(.config$maxSharedMem * 1024 * 1024, proposal)
    }
    proposal
  }
  rbind("mem" = pmax(mem[1, ],
                     RoundUpMem(tooWee),
                     # big numbers - but jobs that request much memory sit in
                     # queues for days
                     1024 * 1024),
        "time" = pmax(exp(wallTime[1, ]),
                      tooShort * c(1, 1, 1.15), # Estimate from timeout ought to be close
                      300, na.rm = TRUE),
        "tmp" = pmax(tmp[1, ], .config$tmpMin, na.rm = TRUE))
}
