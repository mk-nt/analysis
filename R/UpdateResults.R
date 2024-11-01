source("R/Helpers.R")

UpdateResults <- function(pID, scriptID, keepTrees = 256, fetch = TRUE) {
  convFile <- ConvergenceFile(pID, scriptID)
  stoneOrigin <- StoneOrigin(pID, scriptID)
  
  if (isTRUE(fetch)) {
    FetchResults(pID, scriptID)
  }
  
  if (file.exists(stoneOrigin)) {
    if (!file.copy(stoneOrigin, StoneFile(pID, scriptID), overwrite = TRUE)) {
      warning("Couldn't copy ", stoneOrigin)
    }
  } else {
    message("Awaiting stone results for ", pID, "_", scriptID)
  }
  
  converged <- HasConverged(pID, scriptID)
  if (!converged) {
    if (!UpdateRecords(pID, scriptID)) {
      message(pID, "_", scriptID, " is not complete.")
    }
    if (!file.exists(ConvergenceFile(pID, scriptID))) {
      return(FALSE)
    }
    converged <- HasConverged(pID, scriptID)
  }
  
  convStats <- read.table(ConvergenceFile(pID, scriptID))
  burninF <- convStats[["burnin"]]
  
  pFiles <- PFiles(pID, scriptID)
  if (!length(pFiles)) {
    message("No parameter files found for ", pID, "_", scriptID)
    return(FALSE)
  }
  nFiles <- length(pFiles)
  fileContents <- lapply(pFiles, .ReadTable)
  fileContents <- lapply(fileContents, function(x) {x[["Iteration"]] <- NULL; x})
  
  posterior <- do.call(rbind, lapply(fileContents, BurnOff, burninF))
  ess <- convStats[["ess"]]
  posterior <- posterior[seq.int(1, dim(posterior)[[1]], length.out = ess), ]
  saveRDS(posterior, ParameterFile(pID, scriptID))
  
  cli::cli_progress_message(paste("Reading trees from", pID, scriptID))
  allTrees <- lapply(TreeFiles(pID, scriptID), ape::read.tree)
  cli::cli_progress_done()
  
  trees <- unlist(lapply(allTrees, BurnOff, burninF), recursive = FALSE)
  nTrees <- length(trees)
  keptTrees <- trees[seq.int(1, nTrees, length.out = min(keepTrees, nTrees))]
  ape::write.tree(keptTrees, TreeSampleFile(pID, scriptID))
  
  return(TRUE)
}
