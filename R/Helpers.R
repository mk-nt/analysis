source("R/FilePaths.R")

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

.ColourBy <- function(x, palette = "inferno") {
  n <- 512
  hcl.colors(n, palette = palette)[cut(x, n)]
}

.LegendBy <- function(x, palette = "inferno", where = "topleft", label = NULL) {
  PlotTools::SpectrumLegend(
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
        source("R/un-double.R") # TODO REMOVE - fixes glitch
        Dedouble(x) # TODO REMOVE - fixes transient glitch
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

FetchResults <- function(pID, scriptID) {
  if (dir.exists(AnalysisDir(pID, scriptID))) {
    wd <- setwd(AnalysisDir(pID, scriptID))
    on.exit(setwd(wd))
    fetchMsg <- system2("git", "fetch --depth 1", stdout = TRUE)
    if (length(fetchMsg) > 1 && 
        grepl("forced update", fetchMsg[[2]], fixed = TRUE)) {
      rebase <- system2("git", "rebase", stdout = TRUE)
      if (length(rebase) && grepl("Successful", rebase[[1]])) {
        message("Fetched new results from ", pID, " ", scriptID)
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

HasConverged <- function(pID, scriptID, pt = psrfThreshold, et = essThreshold,
                         verbose = FALSE) {
  convFile <- ConvergenceFile(pID, scriptID)
  if (!file.exists(convFile)) {
    if (verbose) {
      message("No convergence file at ", convFile)
    }
    return(FALSE)
  }
  convStats <- read.table(ConvergenceFile(pID, scriptID))
  if (verbose) {
    print(convStats)
  }
  
  # Return:
  convStats[["psrf"]] < pt && 
    convStats[["ess"]] > et && 
    convStats[["frechetCorrelationESS"]] > et &&
    convStats[["medianPseudoESS"]] > et
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

.GitClone <- function(pID, scriptID) {
  clone <- system2("git", sprintf("clone https://github.com/%s/%s %s",
                                  githubAccount, ScriptBase(pID, scriptID),
                                  dirname(ScriptFile(pID, scriptID))),
                   stdout = TRUE, stderr = TRUE)
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
