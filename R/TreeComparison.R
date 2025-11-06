#' Sample of 128 posterior trees for comparison with well-corroborated tree
#' 
#' `ComparisonTrees` samples trees from the posterior; then drops and relabel
#' leaves until the taxon set matches that in the corresponding well
#' corroborated trees
#' 
#' @inheritParams MakeSlurm
#' @importFrom ape read.tree
#' @export
ComparisonTrees <- function(scriptID, pID) {
  tsFile <- TreeSampleFile(pID, scriptID)
  if (!file.exists(tsFile)) {
    ExistingResults(pID, scriptID)
    UpdateRecords(pID, scriptID)
  }
  trees <- if (file.exists(tsFile)) {
    read.tree(tsFile)
  } else {
    message(pID, "_", scriptID, " has not yet converged; sampling")
    do.call(c,
            lapply(TreeFiles(pID, scriptID), function(x)
              read.tree(x) |>
                BurnOff(0.25) |>
                sample(128)
            )
    )
  }
  lapply(trees |> lapply(DeZZ), KeepTip,
         intersect(DeZZ(trees[[1]])[["tip.label"]],
                   wcTrees[[pID]][["tip.label"]]))
}

#' Compare trees with well corroborated tree
#' @param useCache Logical specifing whether to use cached results. Set to
#' `FALSE` to regenerate the cache, e.g. if new trees are available.
#' @importFrom cli cli_progress_done cli_progress_message
#' @importFrom digest digest
#' @importFrom Quartet QuartetStatus SimilarityToReference
#' @importFrom TreeDist ClusteringEntropy MutualClusteringInfo
#' @importFrom TreeTools KeepTip RootTree TipLabels
#' @export
TreeSimilarities <- function(pID, models, method = "quartet", sampleSize = 256,
                             useCache = TRUE) {
  if (trimws(tolower(method)) %in% c("qt", "quartet", "q")) {
    method <- "qt"
  }
  if (trimws(tolower(method)) %in% c("cid", "clustering", "cluster",
                                     "mci", "mc")) {
    method <- "cid"
  }
  hash <- digest(c(pID, models, method, sampleSize), algo = "md5")
  cacheFile <- file.path(OutputDir(), "data", paste0(hash, ".Rds"))
  if (useCache) {
    if (file.exists(cacheFile)) {
      return(setNames(readRDS(cacheFile), models))
    }
  }
  cli_progress_message(paste0(pID, ": Reading trees"))
  on.exit(cli_progress_done())
  treeSets <- lapply(models, ComparisonTrees, pID)
  models <- models[lengths(treeSets) >= sampleSize]
  treeSets <- lapply(treeSets[lengths(treeSets) >= sampleSize], sample, sampleSize)
  wcTree <- wcTrees[[pID]] |>
    KeepTip(TipLabels(treeSets[[1]][[1]])) |>
    RootTree(outgroup[[pID]])
  
  cli_progress_message("Comparing to reference")
  similarities <- if (substr(tolower(method), 1, 1) == "q") {
    lapply(treeSets, function(trees) {
      SimilarityToReference(QuartetStatus(unname(trees), wcTree),
                            normalize = TRUE)
    })
  } else {
    lapply(treeSets, function(trees) {
      MutualClusteringInfo(
        trees, wcTree,
        normalize = ClusteringEntropy(wcTree)
      )
    })
  }
  

  saveRDS(similarities, cacheFile)
  # Return:
  setNames(similarities, models)
}
