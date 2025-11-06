# Update cached results summaries - parameter values and trees, but not
# convergence statistics - based on existing local files.
# If run hasn't converged, or fetch = TRUE, try fetching new results.
# Primarily intended to be called from within UpdateRecords().
# Note that calling with fetch = TRUE will mean that UpdateRecords won't
# recognize results as new.
#' @inheritParams MakeSlurm
#' @importFrom ape write.tree
#' @importFrom TreeDist Entropy
#' @export
UpdateResults <- function(pID, scriptID, keepTrees = 256, fetch = FALSE) {
  convFile <- ConvergenceFile(pID, scriptID)
  stoneOrigin <- StoneOrigin(pID, scriptID)
  
  if (!is.numeric(keepTrees)) {
    stop("keepTrees must be numeric")
  }
  
  if (isTRUE(fetch)) {
    FetchResults(pID, scriptID)
  }
  
  if (file.exists(stoneOrigin)) {
    if (!file.copy(stoneOrigin, StoneFile(pID, scriptID), overwrite = TRUE)) {
      warning("Couldn't copy ", stoneOrigin)
    }
  } else {
    message("Awaiting marginal likelihood estimate for ", pID, "_", scriptID)
  }
  
  converged <- HasConverged(pID, scriptID)
  if (!converged) {
    if (!UpdateRecords(pID, scriptID)) {
      stats <- attr(converged, "stats")
      if (is.null(stats)) {
        message(pID, "_", scriptID, " has not converged")
      } else {
        message(pID, "_", scriptID,
                " is not complete: PSRF = ", signif(stats$psrf, 6),
                "; ESS = ", signif(stats$ess, 5),
                "; TreeESS = ", signif(min(stats$frechet, stats$med), 4))
      }
      MakeSlurm(pID, scriptID, ml = FALSE)
    }
    if (!file.exists(ConvergenceFile(pID, scriptID))) {
      return(structure(FALSE, reason = "UpdateRecords() reports incomplete run"))
    }
    converged <- HasConverged(pID, scriptID)
  }
  
  convStats <- read.table(ConvergenceFile(pID, scriptID))
  burninF <- convStats[["burnin"]]
  
  pFiles <- PFiles(pID, scriptID)
  if (!length(pFiles)) {
    message("No parameter files found for ", pID, "_", scriptID)
    return(structure(FALSE, reason = "No parameter files found"))
  }
  nFiles <- length(pFiles)
  fileContents <- lapply(pFiles, .ReadTable) |>
    lapply(function(x) {x[["Iteration"]] <- NULL; x})

  posterior <- do.call(rbind, lapply(fileContents, BurnOff, burninF))
  ess <- convStats[["ess"]]
  posterior <- posterior[seq.int(1, dim(posterior)[[1]], length.out = ess), ]
  
  if (ModelIsHeterogeneous(scriptID)) {
    probs <- posterior[, grep("matrix_probs.", fixed = TRUE, colnames(posterior))]
    scale <- posterior[, ifelse("beta_scale" %in% colnames(posterior),
                                "beta_scale", "neo_beta")]
    meanN <- vapply(seq_along(scale), function(i) {
      betas <- fnDiscretizeBeta(scale[[i]], scale[[i]], 4)
      # n is equivalent to 1 / rate_loss in by_n_ki
      # So with n = 2.5, we get fnFreeK = [[-1.75, 1.75], [0.7, -0.7]]
      # (as 1.75 / 0.7 = 2.5)
      # This is equivalent to fnF81(Simplex(n, 1)).
      # hg_ki is specified with
      # fnF81( Simplex(abs(1 - beta_categories[i]), beta_categories[i]) )
      # n = 1 / ((1/b) - 1)
      # ns <- 1 / betas - 1
       
      ns <- if (betas[[1]] < sqrt(.Machine[["double.eps"]])) {
        # Then betas[[4]] ~ 1, which could mean that we start
        # propagating rounding errors as we divide by (1 - ~1) = ~0
        nsTmp <- betas[1:2] / (1 - betas[1:2])
        c(nsTmp, 1 / rev(nsTmp))
      } else {
        betas / (1 - betas)
      }
      exp(sum((log(ns) * probs[i, ])))
    }, numeric(1))
    h <- apply(probs, 1, Entropy)
    meanCat <- apply(probs, 1, function(x) sum(x * seq_along(x)))
    posterior <- cbind(posterior, mean_n = meanN, mean_cat = meanCat, cat_h = h)
  }
  if (ModelIsStationary(scriptID)) {
    freqs <- posterior[, startsWith(colnames(posterior), "root_freqs.")]
    posterior <- cbind(
      posterior,
      lg_root_01 = freqs[, 1] / freqs[, 2],
      root_h = apply(freqs, 1, Entropy)
    )
  }
    
  saveRDS(posterior, ParameterFile(pID, scriptID))
  
  allTrees <- ReadTrees(pID, scriptID)
  trees <- unlist(lapply(allTrees, BurnOff, burninF), recursive = FALSE)
  nTrees <- length(trees)
  keptTrees <- trees[seq.int(1, nTrees, length.out = min(keepTrees, nTrees))]
  write.tree(keptTrees, TreeSampleFile(pID, scriptID))
  
  return(TRUE)
}
