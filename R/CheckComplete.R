source("R/Helpers.R")
source("R/MakeSlurm.R")
source("R/UpdateResults.R")

if (length(TreeDist::GetParallel()) == 0) {
  TreeDist::StartParallel(ceiling(parallel::detectCores() * 0.75))
}

ESS <- function(values, burnin, Summarize = min) {
  if (is.list(values) && !is.data.frame(values)) {
    apply(vapply(values, ESS, double(dim(values[[1]])[[2]]), burnin), 1, Summarize)
  } else {
    mcmcse::ess(BurnOff(values, burnin))
  }
}

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

ExtendRun <- function(pID, scriptID) {
  revFile <- ScriptFile(pID, scriptID)
  revLines <- readLines(revFile)
  # TODO this is temporary; increasing by a fixed proportion will cause this
  # to balloon if the run has stopped for some other reason.
  newLines <- gsub("srMinESS\\(\\d+,", "srMinESS(625,", revLines)
  if (any(newLines != revLines)) {
    writeLines(newLines, revFile)
    wd <- setwd(dirname(revFile))
    on.exit(setwd(wd))
    add <- system2("git", paste("add", basename(revFile)), stdout = TRUE)
    if (length(add)) {
      warning(add, immediate. = TRUE)
    }
    commit <- system2("git", "commit -m \"srMinESS = 625\"", stdout = TRUE)
    .GitPush()
    setwd(wd)
    message("Increased ESS stopping rule for ", pID, "_", scriptID)
  } else {
    message("We should continue ", pID, "_", scriptID)
  }
  
  MakeSlurm(pID, scriptID)
}

# Cache ML estimate and convergence diagnostics to `results`.
# if `onlyIfNew`, skip calculations unless new results have been fetched.
# If convergence has not been attained, prepare a longer run for submission
# via `qsub.sh`.
UpdateRecords <- function(pID, scriptID, onlyIfNew = TRUE) {
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
    RemoveSlurm(pID, scriptID, ml = TRUE)
    if (!file.exists(StoneFile(pID, scriptID))) {
      message("Copying new ML estimate for ", pID, " ", scriptID)  
      .CopyStone()
    } else if(tools::md5sum(stoneOrigin) !=
              tools::md5sum(StoneFile(pID, scriptID))) { 
      message("Updating ML estimate for ", pID, " ", scriptID)  
      .CopyStone()
    }
  
  } else {
    MakeSlurm(pID, scriptID, ml = TRUE)
    complete <- FALSE
  }
  
  if (file.exists(ConvergenceFile(pID, scriptID))) {
    if (HasConverged(pID, scriptID)) {
      if (fetched) {
        UpdateResults(pID, scriptID, fetch = FALSE)
        updateResults <- FALSE
      }
      RemoveSlurm(pID, scriptID)
      
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
      githubAccount,
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
            if (onlyIfNew) "; exiting")
    if (onlyIfNew) {
      return(FALSE)
    }
  }
  
  wd <- setwd(revDir)
  on.exit(setwd(wd))
  checkout <- system2("git", "checkout main", stdout = TRUE, stderr = TRUE)
  pull <- system2("git", "pull --rebase", stdout = TRUE)
  
  if (pull[[1]] != "Already up to date.") {
    message("Pull 79 in ", pID, ": ", paste(pull, collapse = "\r\n"))
  }
  
  setwd(wd)


  if (!file.exists(sprintf("%s/%s_run_1.log", revDir, scriptBase))) {
    message("Creating ", SlurmFile(pID, scriptID))
    MakeSlurm(pID, scriptID)
    return(FALSE)
  }
  
  
  pFiles <- list.files(
    path = revDir,
    pattern = sprintf("%s\\.p_run_.*\\.log$", scriptID),
    full.names = TRUE, ignore.case = TRUE
  )
  fileContents <- lapply(pFiles, .ReadTable)
  fileContents <- lapply(fileContents, function(x) {x[["Iteration"]] <- NULL; x})
  
  burnins <- seq(0.0, 0.95, by = 0.05)
  psrf <- sapply(burnins, function(b) PSRF(fileContents, b))
  converged <- apply(psrf < psrfThreshold, 2, all)
  converged[is.na(converged)] <- FALSE
  if (!any(converged)) {
    message("Not converged: PSRF ", signif(min(apply(psrf, 2, max)), 6), " > ",
            signif(psrfThreshold, 6))
    ExtendRun(pID, scriptID)
    return(FALSE)
  }
  hasConverged <- TRUE
  
  ess <- vapply(burnins[converged], function(b) ESS(fileContents, b, sum), psrf[, 1])
  bestESS <- which.max(apply(ess, 2, min))
  bestBurn <- burnins[converged][[bestESS]]
  lowestESS <- min(ess[, bestESS])
  if (lowestESS < essThreshold) {
    message("Undersampled: ESS ", signif(lowestESS, 3), " < ",
            round(essThreshold))
    ExtendRun(pID, scriptID)
    hasConverged <- FALSE
  }
  
  cli::cli_progress_message("Reading trees from file")
  trees <- lapply(TreeFiles(pID, scriptID),
                  function(x) BurnOff(ape::read.tree(x), bestBurn))
  cli::cli_progress_message("Rarifying trees")
  nTrees <- lengths(trees)
  t1000 <- lapply(trees, function(x) x[seq(1, length(x),
                                           length.out = min(1000, nTrees))])
  cli::cli_progress_message("Calculating tree ESS")
  treeESS <- do.call(
    rbind,
    treess::treess(t1000, TreeDist::RobinsonFoulds,
                   methods = treess::getESSMethods(TRUE))
    )[, c("frechetCorrelationESS", "medianPseudoESS")]
  cli::cli_progress_done()
  
  # Mean ensures sum > 2 × threshold
  treeEss <- colSums(treeESS)
  if (min(treeEss) < essThreshold) {
    message("Undersampled?: Tree ESS = ",
            paste(signif(treeEss, 4), collapse = ", "))
    ExtendRun(pID, scriptID)
    hasConverged <- FALSE
  }
  
  write.table(data.frame(burnin = bestBurn,
                         psrf = max(psrf[, match(bestBurn, burnins)]),
                         ess = lowestESS,
                         t(treeEss)),
              file = ConvergenceFile(pID, scriptID))
  
  if (hasConverged) {
    RemoveSlurm(pID, scriptID)
    
    if (isTRUE(updateResults)) {
      UpdateResults(pID, scriptID)
    }
  }
  
  return(hasConverged)
}

# Read existing results from cache.
ExistingResults <- function(pID, scriptID, approx = FALSE) {
  message(paste(pID, scriptID))
  revDir <- AnalysisDir(pID, scriptID)
  scriptBase <- ScriptBase(pID, scriptID)
  
  
  convergence <- if (file.exists(ConvergenceFile(pID, scriptID))) {
    read.table(ConvergenceFile(pID, scriptID))
  } else {
    c(burnin = NA, psrf = NA, ess = NA, NA, NA)
  }
  
  pFiles <- PFiles(pID, scriptID)
  par <- if (length(pFiles)) {
    p <- do.call(rbind,
                 BurnOff(
                   lapply(pFiles, read.table, header = TRUE),
                   if (is.na(convergence[["burnin"]])) 0.5 else
                     convergence[["burnin"]])
    )
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
  if (approx && !file.exists(ppOut)) {
    # Previous stepping stone results were computed with fewer generations
    # and steps; new results have been computed where possible with the
    # computational resources available.
    # The less precise original results provide similar estimates of the Bayes
    # Factors between each method.
    ppOut <- OldStoneFile(pID, scriptID)
  }
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
