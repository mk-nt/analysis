#' @export
.AmbDist <- function(path) {
  if (file.exists(path)) {
    chars <- ReadCharacters(path)
    chars[] <- chars %in% 0:9
    mode(chars) <- "logical"
    list(
      nChar = length(chars),
      nAmb = sum(!chars),
      pChar = colSums(!chars) / dim(chars)[[1]],
      pTax = rowSums(!chars) / dim(chars)[[2]]
    )
  } else {
    warning("No file at ", path)
    vector("list", 4)
  }
}

#' @export
.Ambiguate <- function(x.nex) {
  amb <- vapply(vapply(AllProjects(), MatrixFile, character(1), x.nex),
                .AmbDist, vector("list", 4))
  pAmb <- unlist(amb["nAmb", ]) / unlist(amb["nChar", ])
  mat <- PhyDatToMatrix(ReadAsPhyDat(MkPath(simID, x.nex)))
  nTip <- dim(mat)[[1]]
  nChar <- dim(mat)[[2]]
  charToCensor <- nChar * median(pAmb) * nTip
  censorTax <- sample(unlist(amb["pTax", ], use.names = FALSE), nTip)
  censorChar <- sample(unlist(amb["pChar", ], use.names = FALSE), nChar)
  # pCensor <- outer(censorChar / sum(censorChar), censorTax / sum(censorTax))
  # censored <- sample.int(length(pCensor), charToCensor, prob = as.vector(pCensor))
  censorCol <- sample.int(nChar, charToCensor,
                          p = censorChar / sum(censorChar),
                          replace = TRUE) |> 
    table()
  censorRow <- sapply(seq_along(censorCol), function(i) {
    sample.int(nTip, censorCol[[i]], censorTax, replace = FALSE)
  }) |> unlist()
  censored <- cbind(censorRow, as.numeric(rep(names(censorCol), censorCol)))
  mat[censored] <- "?"
  # hist(unlist(amb["pTax", ]))
  # hist(censorTax, col = 3, add = TRUE)
  # hist(colSums(mat == "?") / ncol(mat), col = 2, add = TRUE)
  # hist(unlist(amb["pChar", ]))
  # hist(censorChar, col = 3, add = TRUE)
  # hist(rowSums(mat == "?") / nrow(mat), add = TRUE, col = 2)
  write.nexus.data(mat, format = "standard", missing = "?",
                   file = MkPath(simID, paste0("imp_", x.nex)))
}


# Equivalently:
# .DecodeTips <- function(n, row = 2) {
#   tr0 <- readLines(MkPath(sprintf("sim%03d", n), "imp_sp_nt_kv_run_1.trees"), row)[[row]]
#   tipTop <- strsplit(tr0, "tip_")[[1]][-1]
#   tip <- gsub("^(\\d+).*", "\\1", tipTop)
#   idx <- gsub("^\\d+\\[&index=(\\d+)\\].*", "\\1", tipTop)
#   # ti <- as.numeric(tip)[order(as.numeric(idx))] # ti[2] = 10 because end_2 is tip 10
#   it <- as.numeric(idx)[order(as.numeric(tip))] # it[10] = 2 because tip 10 is end_2
# }
#' @export
.DecodeTips <- function(nTip) {
  order(order(as.character(seq_len(nTip))))
}

#' @export
.Invariant <- function(sim, nex) {
  impDat <- do.call(rbind, read.nexus.data(
    MkPath(sprintf("sim%03d", sim), nex)))
  impDat[impDat == "?"] <- NA
  impStates <- apply(`mode<-`(impDat, "integer") + 1, 2, tabulate, 2)
  structure(which(colSums(impStates == 0) > 0), nChar = ncol(impStates))
}

#' Pool runs within each replicate, then summarise per replicate
#' Each cell of `accuracy` is a list of nTip vectors (one per tip);
#' c() concatenates the two runs' values for the same tip.
#' @export
.PoolRuns <- function(r1, r2) {
  Map(c, r1, r2)
}

#' Mean accuracy across all tips and characters for one replicate
#' @export
.MeanAcc <- function(pooled) {
  mean(unlist(pooled), na.rm = TRUE)
}

