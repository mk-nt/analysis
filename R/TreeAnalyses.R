#' Compute pairwise tree distances between two model samples
#'
#' This function compares two sets of posterior trees for a given project ID and
#' models, computing a normalized Clustering Information distance using the
#' \pkg{TreeDist} package. If distance results already exist, they are loaded
#' from disk. Otherwise, the function attempts to generate missing tree samples
#' using `MakeSlurm()`, and if both are available, it computes and stores the
#' distance matrix.
#'
#' @inheritParams MakeSlurm
#' @param model1,model2 Character scalars giving the model identifiers to
#'   compare.
#'
#' @return
#' 'TreeDistances()' returns a matrix of pairwise tree distances (as computed by
#' \code{\link[TreeDist]{ClusteringInfoDistance}}) if both model samples are
#' available. If any tree sample is missing and a job is submitted via
#' `MakeSlurm()`, the function returns \code{NULL}.
#'
#' @family tree-comparisons
#' @importFrom parallel detectCores
#' @importFrom TreeTools RootTree
#' @importFrom ape read.tree
#' @importFrom TreeDist ClusteringInfoDistance GetParallel StartParallel
#' @export
TreeDistances <- function(pID, model1, model2) {
  distFile <- DistanceFile(pID, model1, model2)
  if (file.exists(distFile)) {
    readRDS(distFile)
  } else {
    if (!file.exists(TreeSampleFile(pID, model1))) {
      MakeSlurm(pID, model1, ml = FALSE)
      if (!file.exists(TreeSampleFile(pID, model2))) {
        MakeSlurm(pID, model2, ml = FALSE)
      }
      NULL
    } else if (!file.exists(TreeSampleFile(pID, model2))) {
      MakeSlurm(pID, model2, ml = FALSE)
      NULL
    } else {
      sample1 <- read.tree(TreeSampleFile(pID, model1))
      sample2 <- read.tree(TreeSampleFile(pID, model2))
      tip1 <- sample1[[1]][["tip.label"]][[1]]
      rooted1 <- RootTree(sample1, tip1)
      rooted2 <- RootTree(sample2, tip1)
      dTrees <- c(rooted1, rooted2)
      stopifnot(length(rooted1) == length(rooted2))
      
      if (length(GetParallel()) == 0) {
        StartParallel(ceiling(detectCores() * 0.75))
      }
      d <- ClusteringInfoDistance(dTrees, normalize = TRUE)
      saveRDS(d, file = distFile)
      d
    }
  }
}

#' Compute total tree lengths for all trees in a sample
#'
#' Reads the tree sample corresponding to a given project and script identifier,
#' and returns the total branch length for each tree. If the sample file does
#' not exist, a \code{MakeSlurm()} job is initiated to generate it.
#'
#' @return
#' 'TreeLengths()' returns a numeric vector giving the total branch length of
#' each tree in the sample. If the tree file is missing, it returns
#' \code{NULL} after initiating a \code{MakeSlurm()} job.
#'
#' @inheritParams MakeSlurm
#' @family tree-comparisons
#' @importFrom ape read.tree
#' @export
TreeLengths <- function(pID, scriptID) {
  if (file.exists(TreeSampleFile(pID, scriptID))) {
    trees <- read.tree(TreeSampleFile(pID, scriptID))
    colSums(sapply(trees, "[[", "edge.length"))
  } else {
    MakeSlurm(pID, scriptID, ml = FALSE)
    NULL
  }
}

#' Calculate proportional improvement in precision
#'
#' Computes the proportional change in precision between two values, depending
#' on whether improvement is expected to correspond to an increase or decrease
#' relative to a baseline.
#'
#' @param x Numeric vector of length two, giving the values to compare.
#' @param baseline Numeric scalar giving the reference baseline (default
#'   \code{-Inf}).
#' @param improvement Numeric scalar giving the direction of expected
#'   improvement (default \code{Inf}).
#'
#' @return
#' 'PrecisionIncrease()' returns a single numeric value giving the proportional
#' improvement in precision, or \code{NA_real_} if \code{x} is \code{NULL}.
#' 
#' @family tree-comparisons
#' @export
PrecisionIncrease <- function(x, baseline = -Inf, improvement = Inf) {
  if (is.null(x)) {
    NA_real_
  } else {
    switch(order(c(baseline, improvement))[[1]],
           (x[[1]] - x[[2]]) / x[[1]],
           (x[[2]] - x[[1]]) / x[[2]])
  }
}

#' Summarize dispersion of tree distances within and between samples
#'
#' Computes within- and between-sample tree distance summaries, including
#' nearest-neighbour distances, mean minimum spanning tree edge length, median
#' distances, and silhouette width, for a distance matrix computed by
#' \code{\link{TreeDistances}}.
#'
#' @param d A distance matrix, typically the output of
#'   \code{\link{TreeDistances}}.
#'
#' @return
#' 'Dispersion()' returns a list with elements:
#' \itemize{
#'   \item \code{treePairs}: a data frame of pairwise distances and comparison type
#'     ("1 vs 1", "1 vs 2", or "2 vs 2");
#'   \item \code{spread}: a matrix summarising mean nearest-neighbour and MST edge
#'     lengths, and median distances per run;
#'   \item \code{mdmd}: the median-to-median distance between the two runs;
#'   \item \code{sil}: the average silhouette width.
#' }
#'
#' @family tree-comparisons
#' @importFrom TreeDist DistanceFromMedian MeanMSTEdge MeanNN
#' @importFrom cluster silhouette
#' @export
Dispersion <- function(d) {
  if (is.null(d)) {
    return(NULL)
  }
  
  dMat <- as.matrix(d)
  # We assert in TreeDistance that n1 == n2
  n <- dim(dMat)[[1]] / 2
  runID <- rep(1:2, each = n)
  mat11 <- dMat[runID == 1, runID == 1]
  kk <- mat11[lower.tri(mat11)]
  
  mat22 <- dMat[runID == 2, runID == 2]
  nn <- mat22[lower.tri(mat22)]
  
  nk <- dMat[runID == 1, runID == 2]
  df <- data.frame(dist = c(kk, nk, nn),
                   comp = rep(c("1 vs 1", "1 vs 2", "2 vs 2"), 
                              c(length(kk), length(nk), length(nn))))
  
  medianIndex <- c(which.min(colSums(unname(mat11))),
                   which.min(colSums(unname(mat22))))
  spread <- cbind(
    mdI = medianIndex,
    mst = MeanMSTEdge(d, cluster = runID),
    nn = MeanNN(d, cluster = runID, Average = median),
    mad = DistanceFromMedian(d, cluster = runID, Average = median)
  )
  
  list(treePairs = df, spread = spread,
       mdmd = dMat[medianIndex[[1]], n + medianIndex[[2]]],
       sil = mean(silhouette(dist = d, runID)[, 3])
  )
}
