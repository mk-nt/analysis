# Before running this script:
# Set rb.exe=path/to/rb.exe using file.edit(".Renviron")
# (or usethis::edit_r_environ())

library("neotrans")

MakeInv <- function(pID) {
  neoLines <- readLines(MatrixFile(pID, "neo.nex"))
  dimLine <- startsWith(neoLines, "  DIMENSIONS ")
  if (!any(dimLine)) {
    stop("Unexpected format in neo.nex")
  }
  taxonLines <- startsWith(neoLines, "    ")
  nTaxa <- sum(taxonLines)
  stopifnot(nTaxa == as.numeric(
    gsub("^.*NTAX=(\\d+).*$", "\\1", neoLines[dimLine], perl = TRUE)))
  neoLines[dimLine] <- sub("NCHAR=(\\d+)", paste0("NCHAR=", 2 * (nTaxa + 1)),
                           neoLines[dimLine])
  neoLines[taxonLines] <- apply(cbind(
    sub("\\S+$", "", neoLines[taxonLines], perl = TRUE),
    "01",
    diag(nTaxa), 1 - diag(nTaxa)), 1, paste0, collapse = "")
  
  writeLines(neoLines, MatrixFile(pID, "01.nex"))
}

CalcPInf <- function(pID, scriptID) {
  FetchResults(pID, scriptID)
  
  if (!HasConverged(pID, scriptID)) {
    MakeSlurm(pID, scriptID, ml = FALSE)
    return(structure(FALSE, reason = "Not converged"))
  }
  revDir <- dirname(DummyFile(pID, scriptID))
  file.copy(MatrixFile(pID, "01.nex"), DummyFile(pID, scriptID), overwrite = TRUE)
  
  for (run in 1:2) {
    treeFile <- file.path(revDir,
                          paste0(pID, "_", scriptID, "_run_", run, ".trees"))
    pFile <- file.path(revDir,
                       paste0(pID, "_", scriptID,  ".p_run_", run, ".log"))
    if (!file.exists(treeFile)) {
      UpdateRecords(pID, scriptID, searchRemote = TRUE)
    }
    treepFile <- file.path(revDir,
                           paste0(pID, "_", scriptID,  "_run_", run, ".treep"))
    treeLines <- read.table(treeFile, header = TRUE)
    pLines <- read.table(pFile, header = TRUE)
    write.table(pLines[pLines[, 1] %in% treeLines[, 1], ], file = treepFile,
                row.names = FALSE)
    on.exit(unlink(treepFile), add = TRUE)
    
    
    revLines <- gsub(
      "%PID%", pID, fixed = TRUE,
      gsub(
        "%SCRIPTBASE%", paste0(pID, "_", scriptID), fixed = TRUE,
        gsub("%BURNIN%", attr(HasConverged(pID, scriptID), "stats")$burnin,
             fixed = TRUE,
             gsub("%RUN%", run, fixed = TRUE,
                  readLines(sprintf("rbScripts/pinv_%s.Rev", scriptID)))
        )
      )
    )
    writeLines(revLines, file.path(revDir, "pInv.tmp.Rev"))
    on.exit(unlink(file.path(revDir, "pInv.tmp.Rev")), add = TRUE)
    withr::with_dir(revDir,
                    system2(Sys.getenv("rb.exe"), "pInv.tmp.Rev"))
  }
}

PInf <- function(pID, scriptID) {
  pInv <- tryCatch(rbind(
    read.table(file.path(dirname(DummyFile(pID, scriptID)),
                         paste0(pID, "_", scriptID, "_run_1.pinv.log")),
               header = TRUE),
    read.table(file.path(dirname(DummyFile(pID, scriptID)),
                         paste0(pID, "_", scriptID, "_run_2.pinv.log")),
               header = TRUE)),
    warning = function(e) NULL,
    error = function(e) {
      NULL
    })
  if (is.null(pInv)) {
    matrix(NA_real_, 0, 3, dimnames = list(NULL, c("aut", "inv", "inf"))) 
  } else {
    cbind(
      aut = pInv[, 3] - rowSums(pInv[, 1:2]),
      inv = rowSums(pInv[, 1:2]),
      inf = 1 - pInv[, 3]
    )
  }
}

# 
# Compute number of neomorphic characters
# 
PresencePlot <- function(scriptID) {
  fieldSize <- vapply(KiProjects(), function(pID) {
    # Load parameters and determine modelled p(given leaf is present)
    rateLoss <- 1 / ExistingResults(pID, scriptID)[["parameters"]][
      c("2.5%", "25%", "50%", "75%", "97.5%"), "rate_loss"]
    if (length(rateLoss) == 0) {
      MakeSlurm(pID, scriptID, ml = FALSE)
      return(rep(NA_real_, 12))
    }
    statP1 <- rateLoss / (1 + rateLoss)
    rawMad <- ExistingResults(pID, "by_n_ki")[["parameters"]]["mad", "rate_loss"]
    approxMad <- rawMad / (1 + (1 / rateLoss[["50%"]])) ^ 2
    
    #pInf <- mean(PInf(pID, scriptID)[, "inf"])
    
    neo <- TreeTools::ReadCharacters(MatrixFile(pID, "neo.nex"))
    inf <- apply(neo, 2, function(x) {
      sum(table(x)[as.character(0:9)] > 1, na.rm = TRUE) > 1
    })
    ones <- rowSums(neo[, inf] == 1)
    zeros <- rowSums(neo[, inf] == 0)
    ambig <- sum(inf) - ones - zeros
    p1 <- ones / (ones + zeros)
    c(stationary = statP1, stationary.mad = rawMad,
      obs = quantile(p1, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE),
      obs.mad = mad(p1))
  }, double(12))
  
  nBin <- 84
  kiNeoChar <- .meta$nChar["neo", KiProjects()]
  ord <- order(kiNeoChar)
  col <- cut(kiNeoChar, nBin)
  cols <- hcl.colors(nBin, rev = TRUE)[col][ord]
  # Extract quantiles
  x_md <- fieldSize["stationary.50%", ord]
  # x_hi <- x_md + fieldSize["stationary.mad", ord]
  # x_lo <- x_md - fieldSize["stationary.mad", ord]
  x_lo <- fieldSize["stationary.25%", ord]
  x_hi <- fieldSize["stationary.75%", ord]
  
  y_md <- fieldSize["obs.50%", ord]
  # y_hi <- y_md + fieldSize["obs.mad", ord]
  # y_lo <- y_md - fieldSize["obs.mad", ord]
  y_lo <- fieldSize["obs.25%", ord]
  y_hi <- fieldSize["obs.75%", ord]
  
  plot(range(x_lo, x_hi, na.rm = TRUE),
       range(y_lo, y_hi, na.rm = TRUE),
       type = "n", xlab = "Stationary p(present)", ylab = "Observed p(present)",
       asp = 1, frame.plot = FALSE,
       xlim=c(0,1), ylim=c(0,1))
  abline(0, 1, col = "grey70", lty = 2)  # 1:1 line
  
  segments(x_lo, y_md, x_hi, y_md, col = cols)  # horizontal
  segments(x_md, y_lo, x_md, y_hi, col = cols)  # vertical
  points(x_md, y_md, pch = 21, bg = "white", col = cols)
  
  # text(x_md, y_md, names(x_md), pos=4, cex=0.8)
  PlotTools::SpectrumLegend("topleft", palette = hcl.colors(nBin, rev = TRUE),
                            title = "Characters",
                            legend = rev(floor(quantile(kiNeoChar))), bty = "n")
}
# for (pID in KiProjects()) CalcPInf(pID, "ns_n_ki")
# {
#   pdf("figures/presences.pdf", 6, 6)
#   PresencePlot("by_n_ki")
#   PresencePlot("ns_n_ki")
#   dev.off()
# }

if (FALSE) {
  # Summary: Autapomorphic aren't particularly variable, so 2D is enough.
  library("Ternary")
  source("R/metadata.R")
  het <- sapply(KiProjects(), PInf, "hg_b_ki")
  byn <- sapply(KiProjects(), PInf, "by_n_ki")
  par(mfrow = c(1, 1))
  h3 <- sapply(het[lengths(het) > 0], colMeans)
  n3 <- sapply(byn[lengths(byn) > 0], colMeans)
  TernaryPlot("aut", "inv", "inf", "Autapomorphic", "Invariant", "Informative",
             main = "Inferred neomorphic characters")
  TernaryPoints(h3, col = durham$purple, pch = 3)
  TernaryPoints(n3, col = durham$cyan, pch = 4)
  legend("topleft", pch = 3:4, col = c(durham$purple, durham$cyan),
         legend = modelLabel[c("hg_b_ki", "by_n_ki")], bty = "n")
}


#' @importFrom ape read.nexus.data
NInf <- function(pID, scriptID) {
  pInv <- tryCatch(rbind(
    read.table(file.path(dirname(DummyFile(pID, scriptID)),
                         paste0(pID, "_", scriptID, "_run_1.pinv.log")),
               header = TRUE),
    read.table(file.path(dirname(DummyFile(pID, scriptID)),
                         paste0(pID, "_", scriptID, "_run_2.pinv.log")),
               header = TRUE)),
    warning = function(e) NULL,
    error = function(e) {
      NULL
    })
  if (is.null(pInv)) {
    summary(1:2) * NA
  } else {
    nV <- length(read.nexus.data(MatrixFile(pID, "neo.nex"))[[1]])
    nInv <- nV / (1 - pInv[, 3])
    summary(nInv)
  }
}

if (FALSE) {
  het <- sapply(KiProjects(), NInf, "hg_b_ki")
  byn <- sapply(KiProjects(), NInf, "by_n_ki")
  par(mfrow = c(2, 1))
  boxplot(het["Median", ], byn["Median", ], notch = TRUE,
          names = modelLabel[c("hg_b_ki", "by_n_ki")],
          col = c(durham$heather, durham$cyan),
          frame.plot = FALSE,
          ylab = "Number of available characters")
  plot(het["Median", ] ~ .meta$nChar["neo", colnames(het)], frame.plot = FALSE,
       col = durham$heather, pch = 3, asp = 1,
       xlab = "Informative neomorphic characters",
       ylab = "Inferred neomorphic characters")
  points(byn["Median", ] ~ nChar["neo", colnames(het)], pch = 4,
         col = durham$cyan)
  abline(0, 1, lty = "dotted")
  plot(het["Median", ] ~ byn["Median", ], frame.plot = FALSE,
       pch = 3, asp = 1,
       main = "Inferred neomorphic characters",
       xlab = modelLabel["by_n_ki"],
       ylab = modelLabel["hg_b_ki"])
  abline(0, 1, lty = "dotted")
}
