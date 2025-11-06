#' Remove burnin
#' @param values series of parameter estimates
#' @param burnin fraction or number of samples to omit from start
#' @return `BurnOff() returns `values`, without the first `burnin` samples
#' @export
BurnOff <- function(values, burnin) {
  .Keep <- function(n, burnin) {
    (if (burnin < 1) n * burnin else burnin):n
  }
  if (is.null(dim(values))) {
    values[.Keep(length(values), burnin)]
  } else {
    values[.Keep(dim(values)[[1]], burnin), ]
  }
}

#' Convert time to H:M:S format
#' @param secs time period in seconds
#' @return Time period formatted for SLURM
#' @export
AsHMS <- function(secs) {
  d <- secs %/% (24 * 3600)
  h <- secs %% (24 * 3600) %/% 3600
  m <- secs %% 3600 %/% 60
  s <- secs %% 60
  paste(if (d > 0) sprintf("%d-", d),
        sprintf("%02d:%02d:%02d", as.integer(h), as.integer(m), as.integer(s)),
        sep = "")
}

.ColourBy <- function(x, palette = "inferno") {
  n <- 512
  hcl.colors(n, palette = palette)[cut(x, n)]
}

#' @importFrom PlotTools SpectrumLegend
.LegendBy <- function(x, palette = "inferno", where = "topleft", label = NULL) {
  SpectrumLegend(
    where,
    bty = "n",
    xpd = NA,
    palette = hcl.colors(48, palette = palette),
    legend = signif(seq(max(x, na.rm = TRUE), min(x, na.rm = TRUE),
                        length.out = 5), 4),
    title = label
  )
}

.NChar <- function(path) {
  if (file.exists(path)) {
    as.integer(
      gsub(".*NCHAR\\s*=\\s*(\\d+)\\D.*", "\\1", readLines(path, 3)[[3]])
    )
  } else {
    warning("No file at ", path)
    NA_integer_
  }
}

.NTaxa <- function(path) {
  if (file.exists(path)) {
    as.integer(
      gsub(".*NTAX\\s*=\\s*(\\d+)\\D.*", "\\1", readLines(path, 3)[[3]])
    )
  } else {
    warning("No file at ", path)
    NA_integer_
  }
}

#' Number of characters with non-ambiguous state
#' @param path Path to nexus file
#' @export
.NCoded <- function(path) {
  if (file.exists(path)) {
    sum(ReadCharacters(path) %in% 0:9)
  } else {
    warning("No file at ", path)
    NA_integer_
  }
}

.ReadTable <- function(x) {
  tryCatch(
    res <- read.table(x, header = TRUE, colClasses = rep("real", 5)),
    warning = function(w) {
      res <- withCallingHandlers(
        read.table(x, header = TRUE),
        warning = function(w) invokeRestart("muffleWarning")
      )
      res[!apply(is.na(res), 1, any), ]
    }, error = function(e) {
      msg <- e[["message"]]
      
      if (msg ==  "scan() expected 'a real', got 'Iteration'") {
        read.table(x, header = TRUE)
      } else {
        nrows <- as.numeric(sub(".*?(\\d+).*", "\\1", msg, perl = TRUE)) - 2
        # Why 2, not 1? I don't know!
        #  Sometimes nrows = 10 fails with "line 11 didn't have enough elements"
        tryCatch(read.table(x, header = TRUE, nrows = nrows),
                 error = function(e) {
                   stop("Error reading ", x, ":\r\n ", e)
                 })
      }
    }
  )
}

#' Fetch results from remote server
#' @inheritParams MakeSlurm
#' @export
FetchResults <- function(pID, scriptID) {
  if (dir.exists(AnalysisDir(pID, scriptID))) {
    wd <- setwd(AnalysisDir(pID, scriptID))
    on.exit(setwd(wd))
    has_changes <- system2("git", c("status", "--porcelain", "-uno"),
                           stdout = TRUE) |> length() > 0
    if (has_changes) {
      system2("git", c("stash", "push", "-m", "temp_fetch_results"))
      on.exit(system2("git", c("stash", "pop")), add = TRUE)
    }
    
    fetchMsg <- system2("git", "fetch --depth 1", stdout = TRUE)
    if (length(fetchMsg) > 1 && 
        (grepl("forced update", fetchMsg[[2]], fixed = TRUE) ||
         grepl("-> FETCH_HEAD", fetchMsg[[2]], fixed = TRUE)
        )) {
      rebase <- system2("git", "pull --rebase", stdout = TRUE)
      if (length(rebase) && grepl("Successful", rebase[[1]])) {
        message("Fetched new results from ", pID, " ", scriptID)
        return(TRUE)
      } else {
        message("Couldn't rebase ", pID, " ", scriptID, ": \n")
        warning(paste(rebase), immediate. = TRUE)
        return(FALSE)
      }
    }
    statusMsg <- system2("git", "status", stdout = TRUE)
    if (length(statusMsg) > 1 &&
        statusMsg[[2]] == "Your branch and 'origin/main' have diverged,") {
      if (any(grepl("modified: .*_\\d\\.trees", statusMsg))) {
        system2("git", "rm --cached *_run_1.trees *_run_2.trees")
        system2("git", "commit -m \"Cache *.trees files\"")
        file.remove(list.files(pattern = c("*_run_1.trees")))
        file.remove(list.files(pattern = c("*_run_2.trees")))
        system2("git", "rebase")
        system2("git", "push")
      }
      rebase <- system2("git", "rebase", stdout = TRUE)
      if (length(rebase) && grepl("Successful", rebase[[1]])) {
        message("Fetched new results from ", pID, " ", scriptID)
        push <- system2("git", "push", stdout = TRUE)
        return(TRUE)
      } else {
        message("Couldn't rebase ", pID, " ", scriptID, ": \n")
        warning(paste(rebase), immediate. = TRUE)
        return(FALSE)
      }
    }
    setwd(wd)
  } else {
    return(.GitClone(pID, scriptID))
  }
  return(FALSE)
}

#' Read trees from cache
#' @inheritParams MakeSlurm
#' @importFrom cli cli_progress_message cli_progress_done
#' @export
ReadTrees <- function(pID, scriptID) {
  cli_progress_message(paste("Reading trees from", pID, scriptID))
  on.exit(cli_progress_done())
  gzFiles <- TreeFiles(pID, scriptID, compressed = TRUE)
  gz <- length(gzFiles)
  treeFiles <-TreeFiles(pID, scriptID, compressed = FALSE)
  
  if (gz != 0) {
    toUnzip <- if (length(treeFiles) == gz) {
      gzFiles[file.info(gzFiles)$mtime > file.info(treeFiles)$mtime]
    } else {
      gzFiles
    }
    lapply(toUnzip, untar, exdir = dirname(gzFiles[[1]]))
  }
  if (length(TreeFiles(pID, scriptID))) {
    lapply(TreeFiles(pID, scriptID), ape::read.tree)
  } else {
    warning("No tree files found for ", pID, "_", scriptID)
    list()
  }
}

#' Has an analysis converged?
#' @param pt Gelman-Rubin statistic (potential scale reduction factor) threshold;
#' analyses with PSRF > `pt` have not converged.
#' @param et Estimated sample size threshold; analyses with ESS < `et` have
#' not converged.
#' @inheritParams MakeSlurm
#' @returns `HasConverged()` returns a logical specifying whether the specified
#' analysis has converged at the specified thresholds.
#' @export
HasConverged <- function(pID, scriptID, pt = .config$psrfThreshold, et = .config$essThreshold) {
  convFile <- ConvergenceFile(pID, scriptID)
  if (!file.exists(convFile)) {
    return(structure(FALSE, reason = "No convergence file; UpdateRecords()?"))
  }
  convStats <- read.table(ConvergenceFile(pID, scriptID))
  conv <- c(psrf = convStats[["psrf"]] < pt,
            ess = convStats[["ess"]] > et,
            frechet = convStats[["frechetCorrelationESS"]] > et,
            median = convStats[["medianPseudoESS"]] > et)
  # Return:
  structure(all(conv),
            stats = convStats,
            atThreshold = conv)
}

.GitPush <- function(...) {
  std <- system2("git", paste("push", ...), stdout = NULL, stderr = TRUE)
  if (length(std) &&
      std[[1]] != "Everything up-to-date" &&
      substr(std[[1]], 0, 11) != "To https://") {
    warning(std)
  }
  std
}

#' Type of model
#' @returns Returns `TRUE` if the model is of the specified type;
#' `FALSE` otherwise.
#' @inheritParams MakeSlurm
#' @export
ModelIsHeterogeneous <- function(scriptID) {
  grepl("^hg.?_", scriptID)
}

#' @rdname ModelIsHeterogeneous
#' @inheritParams MakeSlurm
#' @export
ModelIsStationary <- function(scriptID) {
  startsWith(scriptID, "ns_")
}

.GitClone <- function(pID, scriptID) {
  clone <- suppressWarnings(
    system2("git", sprintf("clone --depth 1 https://github.com/%s/%s %s",
                           Sys.getenv("ntGithubAccount"),
                           ScriptBase(pID, scriptID),
                           dirname(ScriptFile(pID, scriptID))),
            stdout = TRUE, stderr = TRUE)
  )
  if (length(clone) == 1 && grepl("^Cloning into", clone)) {
    return(TRUE)
  } else if (length(clone) > 1 && grepl("^Cloning into", clone[[1]])
             && grepl("empty repository", clone[[2]])) {
    return(TRUE)
  } else if (length(clone) > 2 && clone[[2]] == "remote: Repository not found.") {
    return (FALSE)
  } else {
    warning(paste(clone, collapse = "\n  "), immediate. = TRUE)
  }
  return(FALSE)
}
