epsBF <- log(10^2) # "Decisive" support, per Jeffries 1939
epsBF <- log(10^1.5) # "Very strong" support, per Jeffries
epsBF <- log(10^0.5) # "Substantial" support, per Jeffries
epsBF <- log(10) # "Strong" support, per Jeffries 1939

#' Bayes factor support for each model, relative to best fitting
#' @param marginals Matrix specifying marginal likelihood of each model (rows)
#' for each dataset (columns)
#' @param models vector identifying the rows of `marginals` to be considered
#' 
#' @export
ModelBF <- function(marginals, models) {
  cf <- marginals[models, , drop = FALSE]
  apply(cf[, colSums(is.na(cf)) == 0, drop = FALSE], 2, function(x) x - max(x))
}

#' Combine two error terms
#' 
#' The error associated with a difference can be computed from the error
#' associated with each input variable
#' @param err1,err2 Numeric vectors giving errors to be combined.
#' @returns `.DiffErr()` returns a vector detailing the combined error for
#' each element in `err1`.
#' @export
.DiffErr <- function(err1, err2) {
  2 * sqrt((err1 * err1) + (err2 * err2))
}

#' Outwith error
#' 
#' Do two estimates differ by a noteworthy amount?
#' 
#' @param x Value to compare; typically a Bayes Factor giving difference
#' in estimated marginal likelihoods of two models
#' @param err Standard error of comparison, perhaps computed using `.DiffErr()`
#' @param eps Numeric; errors smaller than `eps` are not deemed noteworthy.
#' @returns `.OutwithError()` returns a logical vector stating whether each
#' element in `x` is greather than the corresponding `eps` and `err`
#' @export
.OutwithError <- function(x, err, eps = epsBF) {
  x > eps & x > err
}

#' Evaluate t parameter
#' 
#' Does setting a free \eqn{t} parameter improve model fit relative to a fixed
#' \eqn{t = 1}?
#' 
#' @export
EvaluatePartitioning <- function(marginals, suffix = "_ki") {
  models <- paste0(c("by", "by_t", "rm_by_t"), suffix) |>
    setNames(c("base", "t", "rm_t"))
  
  stdErr <- attr(marginals, "stdErr")
  errT <- stdErr[models[["t"]], ]
  errRmT <- stdErr[models[["rm_t"]], ]
  errBase <- stdErr[models[["base"]], ]
  
  partition <- ModelBF(marginals, models)
  
  t_best <- partition[models[["t"]], ] == 0
  
  t_really_wins <- t_best &
    .OutwithError(-partition[models[["base"]], ], .DiffErr(errBase, errT)) &
    .OutwithError(-partition[models[["rm_t"]], ], .DiffErr(errRmT, errT))
  
  t_beats_nothing <- .OutwithError(partition[models[["t"]], ] - 
                                     partition[models[["base"]], ],
                                   .DiffErr(errBase, errT))
  t_beats_rt_and_nothing <- .OutwithError(
    partition[models[["t"]], t_beats_nothing] -
      partition[models[["rm_t"]], t_beats_nothing],
    .DiffErr(errRmT, errT)[t_beats_nothing])
    
  
  rt_best <- partition[models[["rm_t"]], ] == 0
  rt_really_wins <- rt_best & 
    .OutwithError(-partition[models[["base"]], ], .DiffErr(errRmT, errBase)) &
    .OutwithError(-partition[models[["t"]], ], .DiffErr(errRmT, errT))
    
  partitioning_worse <- partition[models[["base"]], ] == 0
  
  n_projects <- length(t_really_wins)
  
  message("by_t is >eps better than unpartitioned in ", sum(t_beats_nothing),
          " / ", n_projects, " projects.")
  message("by_t is >eps better than rm_by_t in ", sum(t_beats_rt_and_nothing),
          " / ", sum(t_beats_nothing), " of these.")
  message("rm_by_t is best model by >eps in ", sum(rt_really_wins), " / ",
          n_projects)
  
  # Sanity check: Have we caught all cases?
  stopifnot(all(partitioning_worse | t_best | rt_best))
  invisible(partition)
}

#' Informative neomorphic cells
#'
#' Reads the neomorphic character matrix for a given project and extracts only
#' those characters that are informative (i.e., showing more than one repeated
#' state across taxa).  Returns per-taxon counts of neomorphic cell states,
#' together with proportions and overall counts.
#'
#' @inheritParams MakeSlurm
#'
#' @details
#' The function reads the neomorphic matrix file \code{"neo.nex"} associated with
#' the given project ID, filters out uninformative characters, and calculates for
#' each taxon:
#' \itemize{
#'   \item Number of cells with state \code{0} (\code{"0"})
#'   \item Number of cells with state \code{1} (\code{"1"})
#'   \item Number of ambiguous cells (\code{"?"})
#'   \item Proportion of \code{1} cells among those scored as \code{0} or \code{1} (\code{"p1"})
#' }
#'
#' @return
#' \code{InfNeo()} invisibly returns a list with elements:
#' \describe{
#'   \item{\code{"0"}}{Vector of per-taxon counts of cells scored as 0.}
#'   \item{\code{"1"}}{Vector of per-taxon counts of cells scored as 1.}
#'   \item{\code{"?"}}{Vector of per-taxon counts of ambiguous cells.}
#'   \item{\code{"p1"}}{Vector of proportions of 1-cells among scored cells.}
#'   \item{\code{"count"}}{Named vector giving total counts of 0-, 1-, and 0/1-cells across all taxa.}
#'   \item{\code{"p"}}{Named vector giving overall proportions of 0 and 1 states among scored cells.}
#' }
#'
#' @importFrom TreeTools ReadCharacters
#' @export
InfNeo <- function(pID) {
  neo <- ReadCharacters(MatrixFile(pID, "neo.nex"))
  neo <- neo[, apply(neo, 2, function(x) {
    sum(table(x)[as.character(0:9)] > 1, na.rm = TRUE) > 1
  })]
  ones <- rowSums(neo == 1)
  zeros <- rowSums(neo == 0)
  zeroone <- ones + zeros
  ambig <- dim(neo)[[2]] - ones - zeros
  count <- c("0" = sum(zeros), "1" = sum(ones), "01" = sum(zeros) + sum(ones))
  list("0" = zeros, "1" = ones, "?" = ambig, "p1" = ones / zeroone,
       "count" = count, "p" = count[c("0", "1")] / count["01"])
}
