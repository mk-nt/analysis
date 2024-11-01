library("TreeTools", quietly = TRUE)
library("colorspace")

# Plot a consensus tree for summary PDFs
PlotCons <- function(forest, title, colour) {
  cons <- SortTree(Consensus(forest, p = 0.5))
  splitFreqs <- SplitFrequency(cons, forest)
  splitP <- splitFreqs / length(forest)
  
  edgeSupport <- rep(1, nrow(cons$edge)) # Initialize trivial splits to 1
  childNode <- cons[["edge"]][, 2]
  edgeSupport[match(names(splitFreqs), childNode)] <- splitP
  hcl <- attr(as(hex2RGB(figPalette[[4]]), "polarLUV"), "coords")[1, ]
  scale <- rev(diverge_hcl(101, h = c(hcl[["H"]], 0),
                           c = c(hcl[["C"]], 0), l = c(hcl[["L"]], 90), power = 1))
  
  plot(cons, main = title, edge.col = SupportColour(edgeSupport, scale = scale),
       edge.width = 3, col.main = colour)
  LabelSplits(cons, round(100 * splitP, 1), unit = "%",
              col = SupportColour(splitP, scale = scale),
              frame = "none", pos = 3L)
  invisible(cons)
}
