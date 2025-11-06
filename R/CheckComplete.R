#' Effective sample size
#'
#' Calculates the effective sample size (ESS) for one or more MCMC chains,
#' discarding an initial burn-in.
#'
#' If a list of chains is supplied, `ESS()` applies the calculation to each chain
#' and summarizes results across chains using the specified summary function.
#'
#' @param values Either a numeric matrix (samples × parameters) or a list of such
#'   matrices, each representing an independent MCMC run.
#' @param burnin Integer giving the number of initial samples to discard.
#' @param Summarize Function to summarize ESS values across runs
#'   (default: \code{min}).
#'
#' @details
#' The effective sample size is estimated using \code{\link[mcmcse]{ess}} from
#' the \pkg{mcmcse} package. Burn-in samples are removed using
#' \code{BurnOff()} prior to calculation.
#'
#' @return
#' `ESS()` returns a numeric vector of effective sample sizes for each parameter.
#' When multiple runs are supplied, the summarized ESS across runs is returned.
#'
#' @seealso [PSRF()] for potential scale reduction diagnostics.
#'
#' @importFrom mcmcse ess
#' @export
ESS <- function(values, burnin, Summarize = min) {
  if (is.list(values) && !is.data.frame(values)) {
    apply(vapply(values, ESS, double(dim(values[[1]])[[2]]), burnin), 1, Summarize)
  } else {
    mcmcse::ess(BurnOff(values, burnin))
  }
}

#' Potential scale reduction factor (PSRF)
#'
#' Computes the univariate potential scale reduction factor (PSRF) for assessing
#' MCMC convergence across multiple chains.
#'
#' @param runs A list of numeric matrices, each representing an independent MCMC
#'   chain with dimensions samples × parameters.
#' @param burnin Integer giving the number of initial samples to discard.
#'
#' @details
#' The PSRF compares between-chain and within-chain variances to assess whether
#' chains have mixed adequately. Values near 1 indicate convergence.
#' This implementation follows Equation 4 of Vats & Knudson (2021),
#' *Statistical Science*, \doi{10.1214/20-STS812}.
#'
#' @return
#' `PSRF()` returns a numeric vector giving the PSRF for each parameter.
#'
#' @seealso [ESS()] for assessing effective sample sizes.
#'
#' @export
PSRF <- function(runs, burnin) {
  if (!is.list(runs)) {
    stop("`runs` must be a list in which each element corresponds to a run")
  }
  nRun <- length(runs)
  samples <- lapply(runs, BurnOff, burnin)
  if (is.null(dim(samples[[1]]))) {
    samples <- lapply(samples, matrix)
  }
  nSample <- dim(samples[[1]])[[1]]
  nPar <- dim(samples[[1]])[[2]]
  
  chainMean <- vapply(samples, colMeans, double(nPar)) # Xbari
  chainVar <- vapply(samples, apply, double(nPar), 2, var) #s2i
  allMean <- rowMeans(chainMean) # muHat
  varOfMeans <- apply(chainMean, 1, var) # B / n
  # equivalently, varOfMeans <- rowSums((chainMean - allMean) ^ 2) / (nRun - 1)
  meanChainVar <- rowMeans(chainVar) # s2
  
  sigmaHat2 <- ((nSample - 1) / nSample * meanChainVar) + varOfMeans
  # Eqn 4 in Vats & Knudson 2021, doi:10.1214/20-STS812
  psrf <- sigmaHat2 / meanChainVar
  # StratoBayes#133: Replace this calculation with Vats & Knudson's univariate PSRF
  
  psrf
}

#' Submit a continuation job
#' @param proportion How much longer should the run be?
#'   e.g. if the observed ess is 100 and we need an ESS of 300, proportion = 3
#' @param buffer Ratio to increase proportion to be sure we clear
#'   our threshold.  1.15 is a 15% buffer.
ExtendRun <- function(pID, scriptID, proportion, buffer = 1.15) {
  if (proportion < 1) {
    stop("`proportion` must be >= 1")
  }
  existing <- ExistingResults(pID, scriptID)[["convergence"]][["ess"]]
  if (is.na(existing)) {
    message("No existing ESS for ", pID, " ", scriptID)
  } else {
    MakeSlurm(pID, scriptID,
              ess = ceiling(existing * proportion * buffer / 10) * 10)
  }
}


#' Fetch ML and MCMCMC results from GitHub.
#' 
#' If new MCMCMC result are available, compute and cache convergence statistics.
#' 
#' Parameters and trees are written to cache using `UpdateResults()`, called
#' from within this function.
#' 
#' @param searchRemote Logical: if `git fetch` does not retrieve new results,
#' should we check for uncommitted changes on the remote host?
#' @param forgetCache Logical: if no new results are found, should we
#' recalculate cached statistics anyway?
#' @importFrom cli cli_progress_message cli_progress_done
#' @importFrom parallel detectCores
#' @importFrom tools md5sum
#' @importFrom treess getESSMethods treess
#' @importFrom TreeDist GetParallel RobinsonFoulds StartParallel
#' @export
UpdateRecords <- function(pID, scriptID, searchRemote = FALSE,
                          forgetCache = FALSE, makeSlurm = FALSE) {
  revDir <- AnalysisDir(pID, scriptID)
  scriptBase <- ScriptBase(pID, scriptID)
  complete <- TRUE
  updateResults <- TRUE
  
  fetched <- FetchResults(pID, scriptID)
  
  stoneOrigin <- StoneOrigin(pID, scriptID)
  .CopyStone <- function() {
    if (!file.copy(stoneOrigin, StoneFile(pID, scriptID), overwrite = TRUE)) {
      warning("Couldn't copy ", stoneOrigin)
    }
  }
  if (file.exists(stoneOrigin)) {
    if (read.table(stoneOrigin, row.names = NULL)[[5]] < 128) {
      RevBayes(pID, scriptID) # Update marginal.Rev to use more samples
      message("Attempting more accurate ml run; previous only 51 steps")
      MakeSlurm(pID, scriptID, ml = TRUE)
    }
    if (!file.exists(StoneFile(pID, scriptID))) {
      message("Copying new ML estimate for ", pID, " ", scriptID)  
      .CopyStone()
    } else if(md5sum(stoneOrigin) !=
              md5sum(StoneFile(pID, scriptID))) { 
      old <- read.table(StoneFile(pID, scriptID), row.names = NULL)
      new <- read.table(stoneOrigin, row.names = NULL)
      if (!isTRUE(all.equal(old[1:5], new[1:5]))) {
        message("Updating ML estimate for ", pID, " ", scriptID, " by ",
                signif(as.numeric(old[[1]]) - as.numeric(new[[1]]), 4),
                " from ", old[[1]], " to ", new[[1]])
        .CopyStone()
      }
    }
  
  } else {
    MakeSlurm(pID, scriptID, ml = TRUE)
    complete <- FALSE
  }
  
  if (file.exists(ConvergenceFile(pID, scriptID))) {
    if (HasConverged(pID, scriptID)) {
      if (fetched || searchRemote) {
        UpdateResults(pID, scriptID, fetch = FALSE)
        updateResults <- FALSE
      }
      
      if (!complete) {
        message("MCMCMC complete; still awaiting ML estimate for ", pID, " ",
                scriptID)
      }
      return(complete)
    } else {
      # continue: our convergence table may be outdated
    }
  }
  
  if (!dir.exists(revDir)) {
    # Check out possible result files into revDir
    clone <- suppressWarnings(system2("git", sprintf(
      "clone --depth 1 https://%s@github.com/%s/%s %s",
      Sys.getenv("MKNT_READ"),
      Sys.getenv("ntGithubAccount"),
      ScriptBase(pID, scriptID),
      revDir), stdout = TRUE))
  }
  
  if (!dir.exists(revDir)) {
    message("Couldn't clone to ", revDir)
    message(clone)
    return(FALSE)
  }

  if (!fetched) {
    message("No new results for ", pID, " ", scriptID,
            if (!searchRemote && !forgetCache) "; exiting")
    if (!searchRemote && !forgetCache) {
      return(FALSE)
    }
  }
  
  wd <- setwd(revDir)
  on.exit(setwd(wd))
  checkout <- system2("git", "checkout main", stdout = TRUE, stderr = TRUE)
  pull <- system2("git", "pull --rebase", stdout = TRUE)
  
  if (pull[[1]] != "Already up to date.") {
    message("Pull error 157 in ", pID, "_", scriptID, ": ",
            paste(pull, collapse = "\r\n"))
  }
  
  setwd(wd)

  log1 <- sprintf("%s_run_1.log", scriptBase)
  if (!file.exists(file.path(revDir, log1)) || isTRUE(searchRemote)) {
    # Check on the remote
    remoteFile <- file.path(RemoteDir(), scriptBase, log1)
    session <- SshSession()
    if (ssh_exec_wait(session, paste("test -f", remoteFile)) == 0) {
      push <- ssh_exec_wait(
        .ssh$session,
        paste("(cd", dirname(remoteFile), "&&",
              "shopt -s nullglob;",
              "if ls *.trees 1> /dev/null 2>&1; then",
              "  for file in *.trees;",
              "    do tar -czf \"${file%.trees}.tar.gz\" \"$file\";",
              "  done;",
              sprintf("  git add *.ckp %s*.log *.tar.gz &&", pID),
              "  git commit -m \"UpdateRecords(): MCMCMC output files\" &&",
              "  git push origin main;",
              "fi)"))
      if (push == 0) {
        return(UpdateRecords(pID, scriptID,
                             searchRemote = FALSE, # Already done that!
                             forgetCache = forgetCache,
                             makeSlurm = makeSlurm))
      }
    }
    
    if (makeSlurm) {
      message("Creating ", SlurmFile(pID, scriptID))
      MakeSlurm(pID, scriptID)
    }
    return(structure(FALSE, reason = "No .log files found"))
  }
  
  if (fetched || forgetCache) {
    
    pFiles <- list.files(
      path = revDir,
      pattern = sprintf("%s\\.p_run_.*\\.log$", scriptID),
      full.names = TRUE, ignore.case = TRUE
    )
    # Each row contains parameter values for a given iteration
    fileContents <- lapply(pFiles, .ReadTable)
    # Remove iteration column
    fileContents <- lapply(fileContents, function(x) {x[["Iteration"]] <- NULL; x})
    
    burnins <- seq(0.0, 0.95, by = 0.05)
    psrf <- sapply(burnins, function(b) PSRF(fileContents, b))
    converged <- apply(psrf < .config$psrfThreshold, 2, all)
    converged[is.na(converged)] <- FALSE
    if (!any(converged)) {
      message("Not converged: PSRF ", signif(min(apply(psrf, 2, max)), 6), " > ",
              signif(.config$psrfThreshold, 6))
      if (makeSlurm) {
        ExtendRun(pID, scriptID, 2)
      }
      return(structure(FALSE, reason = "Not converged"))
    }
    hasConverged <- TRUE
    
    # Calculate ESS from POOLED samples (not individual runs)
    ess <- vapply(burnins[converged],
                  function(b) ESS(fileContents, b, sum),
                  psrf[, 1])
    bestESS <- which.max(apply(ess, 2, min))
    bestBurn <- burnins[converged][[bestESS]]
    lowestESS <- min(ess[, bestESS])
    if (lowestESS < .config$essThreshold) {
      message("Undersampled: ESS ", signif(lowestESS, 3), " < ",
              round(.config$essThreshold))
      
      # RevBayes considers individual run ESS
      runESS <- min(ESS(fileContents, burnins[converged][[bestESS]], min))
      if (makeSlurm) {
        ExtendRun(pID, scriptID, .config$essThreshold / signif(runESS, 3))
      }
      hasConverged <- FALSE
    }
    
    cli_progress_message("Reading trees from file")
    trees <- lapply(ReadTrees(pID, scriptID), BurnOff, bestBurn)
    cli_progress_message("Rarifying trees")
    nTrees <- lengths(trees)
    # Very small values will crash treess; don't bother computing
    if (min(nTrees) < .config$essThreshold / 16) {
      treeEss <- c("frechetCorrelationESS" = NA_real_,
                   "medianPseudoESS" = NA_real_)
      message("Undersampled; only ", min(nTrees), " trees per run")
      # Assume run already extended during overall ESS check earlier
      hasConverged <- FALSE
    } else {
      t1000 <- lapply(trees, function(x) x[seq(1, length(x),
                                               length.out = min(1000, nTrees))])
      cli_progress_message("Calculating tree ESS")
      
      if (length(GetParallel()) == 0) {
        StartParallel(ceiling(detectCores() * 0.75))
      }
      treeESS <- do.call(
        rbind,
        treess(t1000, RobinsonFoulds, methods = getESSMethods(TRUE))
        )[, c("frechetCorrelationESS", "medianPseudoESS")]
      cli_progress_done()
      
      # Mean ensures sum > 2 × threshold
      treeEss <- colSums(treeESS)
      
      if (min(treeEss) < .config$essThreshold) {
        message("Undersampled?: Tree ESS = ",
                paste(signif(treeEss, 4), collapse = ", "))
        if (makeSlurm) {
          ExtendRun(pID, scriptID, .config$essThreshold / min(treeEss))
        }
        hasConverged <- FALSE
      }
    }
    
    write.table(data.frame(burnin = bestBurn,
                           psrf = max(psrf[, match(bestBurn, burnins)]),
                           ess = lowestESS,
                           t(treeEss)),
                file = ConvergenceFile(pID, scriptID))
  } else {
    hasConverged <- HasConverged(pID, scriptID)
  }
  
  if (hasConverged) {
    if (isTRUE(updateResults)) {
      UpdateResults(pID, scriptID)
    }
  }
  
  return(hasConverged)
}

#' Existing results
#' Read results that have already been cached using `UpdateRecords()`.
#' @param `checkRemote` Logical. If results do not exist locally, should we
#' attempt to retrieve them from the remote server?
#' @inheritParams MakeSlurm
#' @export
ExistingResults <- function(pID, scriptID, checkRemote = FALSE) {
  if (checkRemote) {
    message(paste(pID, scriptID))
  }
  revDir <- AnalysisDir(pID, scriptID)
  scriptBase <- ScriptBase(pID, scriptID)
  
  
  convergence <- if (file.exists(ConvergenceFile(pID, scriptID))) {
    read.table(ConvergenceFile(pID, scriptID))
  } else {
    c(burnin = NA, psrf = NA, ess = NA, NA, NA)
  }
  
  if (checkRemote && !file.exists(ParameterFile(pID, scriptID))) {
    UpdateResults(pID, scriptID, fetch = TRUE)
  }
  p <- if (file.exists(ParameterFile(pID, scriptID))) {
    readRDS(ParameterFile(pID, scriptID))
  } else {
    NULL
  }
  par <- if (!is.null(p)) {
    ignore <- c("Iteration", "Posterior", "Likelihood", "Prior")
    params <- p[, setdiff(colnames(p), ignore)]
    rbind(
      apply(params, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975)),
      mad = apply(params, 2, mad)
    )
  } else {
    NULL
  }
  
  ppOut <- StoneFile(pID, scriptID)
  
  list(
    marginal = if (file.exists(ppOut)) {
      suppressWarnings(read.table(ppOut))[[1]]
    } else {
      NA_real_
    },
    convergence = convergence,
    parameters = par
  )
}

#' @export
ValidateStone <- function(pID, scriptID) {
  stone <- StoneFile(pID, scriptID)
  if (file.exists(stone)) {
    stoneData <- strsplit(readLines(stone), "\t")[[1]]
    if (length(stoneData) == 5 || stoneData[[9]] == "NA") {
      if (Sys.getenv("rb.exe") == "") {
        warning(stone, "outdated, but Sys.getenv('rb.exe') not set.")
        return(stone)
      }
      
      # Initial analyses didn't calculate ML SE.
      # If it's not been calculated, calculate it now.
      ppFile <- gsub("-stone.txt", ".pp", basename(stone), fixed = TRUE)
      tmpFile <- gsub("\\\\", "/", tempfile())
      
      revDir <- dirname(DummyFile(pID, scriptID))
      
      ppPath <- file.path(revDir, ppFile)
      if (!file.exists(ppPath)) {
        .GitClone(pID, scriptID)
      }
      if (!file.exists(ppPath)) {
        warning("Could not clone repo for \"", pID, "\", \"", scriptID, "\"")
        return(NULL)
      }
      
      rbArgs <- paste("-c", file.path(getwd(), "rbScripts/marginalSE.Rev"),
                      file.path(revDir, ppFile), tmpFile)
      system2(Sys.getenv("rb.exe"), rbArgs, stdout = "rb.stdout")
      
      if (file.exists(tmpFile)) {
        revOut <- strsplit(readLines(tmpFile), "\t")[[1]]
        if (isTRUE(all.equal(revOut[1:2], stoneData[1:2]))) {
          write(paste0(c(stoneData[1:5], revOut[-(1:2)]), collapse = "\t"), stone)
          stone
        } else {
          UpdateRecords(pID, scriptID, searchRemote = FALSE,
                        forgetCache = FALSE, makeSlurm = FALSE)
          warning("Marginal likelihood calculation mismatch for ",
                  pID, " ", scriptID, ": ",
                  paste(revOut[1:2], collapse = ", "), " != ",
                  paste(stone[1:2], collapse = ", "), .immediately = TRUE)
          NULL
        }
      } else {
        message("Couldn't compute Std Err for ", pID, " ", scriptID)
        NULL
      }
    } else if (length(stoneData) != 9) {
      warning("Stone file ", stone, " has ", length(output),
              " columns; expected 9")
      NULL
    } else {
      stone
    }
  } else {
    NULL
  }
}

#' @export
GetMarginal <- function(pID, scriptID) {
  stone <- ValidateStone(pID, scriptID)
  if (is.null(stone)) {
    NA_real_
  } else {
    unlist(read.table(stone, colClasses = c("numeric", rep("NULL", 8))),
           use.names = FALSE)
  }
}

#' @export
GetMarginalSE <- function(pID, scriptID) {
  stone <- ValidateStone(pID, scriptID)
  if (is.null(stone)) {
    NA_real_
  } else {
    stdE <- unlist(read.table(stone, colClasses = c(rep("NULL", 5),
                                                    "numeric", "NULL",
                                                    "numeric", "NULL")),
                   use.names = FALSE)
    if (is.na(stdE[[2]])) {
      warning("Bootstrap estimate failed for ", pID, ", \"", scriptID, "\"")
      stdE[[1]]
    } else {
      stdE[[2]]
    }
  }
}

#' @export
GetConvStat <- function(pID, scriptID) {
  conv <- HasConverged(pID, scriptID)
  unlist(attr(conv, "stats")) %||% rep(NA_real_, 5)
}

#' @export
MarginalDiff <- function(pID, scriptID) {
  stone <- ValidateStone(pID, scriptID)
  if (is.null(stone)) {
    NA_real_
  } else {
    diff(unlist(read.table(stone,
                           colClasses = c("numeric", "numeric",
                                          rep("NULL", 7))), use.names = FALSE))
  }
}

#' @export
GetMarginals <- function(projects, models) {
  structure(
    t(`colnames<-`(
      vapply(models, function(scriptID) vapply(projects, GetMarginal, numeric(1), scriptID),
             setNames(double(length(projects)), projects)), models)),
    stdErr = t(`colnames<-`(
      vapply(models, function(scriptID) vapply(projects, GetMarginalSE, numeric(1), scriptID),
             setNames(double(length(projects)), projects)), models)))
}

#' @export
GetConvergence <- function(projects, models) {
  vapply(models, function(scriptID) vapply(projects, GetConvStat, numeric(5), scriptID),
         matrix(NA_real_, 5, length(projects), dimnames = list(
           c("burnin", "psrf", "ess", "frESS", "mpsESS"), projects)))
}

#' Difference between marginal likelihood estimates obtained by stepping stone
#' and path sampling analyses.
#' @param projects Character vector listing project identifiers to evaluate.
#' @param models Character vector specifying models to evaluate.
#' @returns `MarginalDiffs()` returns a matrix specifying the difference
#' (in log units) between the stepping stone and path sampling estimates of the
#' marginal likelihood.
#' @importFrom stats setNames
#' @export
MarginalDiffs <- function(projects, models) {
  t(`colnames<-`(
    vapply(models, function(scriptID) vapply(projects, MarginalDiff, numeric(1), scriptID),
           setNames(double(length(projects)), projects)), models))
}
