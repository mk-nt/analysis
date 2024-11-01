source("R/_setup.R")
source("R/Helpers.R")
source("R/PlotCons.R")
if (length(TreeDist::GetParallel()) == 0) {
  TreeDist::StartParallel(ceiling(parallel::detectCores() *.75))
}

modelNames <- c(
  "by_ki" = "Mk: Equal rates",
  "by_t_ki" = "Mk-t: Neo/Trans",
  "by_n_ki" = "Mk-n: Loss/Gain",
  "by_nt_ki" = "Mk-nt: Neo/Trans, Loss/Gain"
)
ntkModels <- c("by_nt_ki")

AnalysisPDF <- function(pID, model1, model2, nDist = 256, overwrite = FALSE) {
  pdfFile <- PDFFile(pID, model1, model2)
  converged <- HasConverged(pID, model1) && HasConverged(pID, model2)
  distFile <- DistanceFile(pID, model1, model2, final = converged)
  if (!overwrite && file.exists(pdfFile) && file.exists(distFile)) {
    return(structure(FALSE, reason = "Already exists"))
  }
  
  if (!file.exists(ParameterFile(pID, model1)) ||
      !file.exists(ParameterFile(pID, model2)) ||
      !file.exists(TreeSampleFile(pID, model1)) ||
      !file.exists(TreeSampleFile(pID, model2))) {
    message("Results not available for ", pID, " ", model1, " vs ", model2)
    return(structure(FALSE, reason = "Results not available"))
  }
  
  res1 <- readRDS(ParameterFile(pID, model1))
  res2 <- readRDS(ParameterFile(pID, model2))
  
  tree1 <- ape::read.tree(TreeSampleFile(pID, model1))
  tree2 <- ape::read.tree(TreeSampleFile(pID, model2))
  tip1 <- tree1[[1]][["tip.label"]][[1]]
  rooted1 <- TreeTools::RootTree(tree1, tip1)
  rooted2 <- TreeTools::RootTree(tree2, tip1)
  
  pdf(pdfFile, 6.5, 8,
      title = sprintf("Model %s vs %s, project %s%s", model1, model2, pID,
                      if (converged) "" else " (preliminary)"))
  
  oPar <- par(mfrow = 1:2, mar = c(0.1, 0.1, 2, 0), cex = 0.75, xpd = NA)
  PlotCons(rooted1, sprintf("Consensus, %s %s", pID, model1), col = figPalette[[2]])
  PlotCons(rooted2, sprintf("Consensus, %s %s", pID, model2), col = figPalette[[3]])
  par(oPar)
  
  oPar <- par(mfrow = 2:1, mar = c(0, 0, 1, 0))
  d1 <- rooted1[seq(1, length(rooted1), length.out = min(length(rooted1), nDist))]
  d2 <- rooted2[seq(1, length(rooted2), length.out = min(length(rooted2), nDist))]
  dTrees <- c(d1, d2)
  n1 <- length(d1)
  n2 <- length(d2)
  d <- TreeDist::ClusteringInfoDistance(dTrees, normalize = TRUE)
  saveRDS(d, distFile)
  runID <- rep(1:2, c(n1, n2))
  plot(cmdscale(d), pch = runID + 2,
       col = figPalette[runID + 1],
       font.main = 1,
       cex = 1.5,
       asp = 1, # Preserve aspect ratio - do not distort distances
       ann = FALSE, axes = FALSE # Don't label axes: dimensions are meaningless
  )
  
  title("Mapping of tree-tree distances", line = 0)
  legend("bottomleft", modelNames[c(model1, model2)],
         bty = "n", pch = 3:4, col = figPalette[2:3], xpd = NA)
  
  
  dMat <- as.matrix(d)
  kTrees <- seq_len(n1)
  kk <- dMat[kTrees, kTrees]
  kk <- kk[lower.tri(kk)]
  
  nn <- dMat[-kTrees, -kTrees]
  nn <- nn[lower.tri(nn)]
  
  nk <- dMat[kTrees, -kTrees]
  df <- data.frame(dist = c(kk, nk, nn),
                   comp = rep(c("+ vs +", "+ vs x", "x vs x"), 
                              c(length(kk), length(nk), length(nn))))
  par(mar = c(2, 4, 0, 0), new = FALSE)
  boxplot(dist ~ comp, df, notch = TRUE, col = figPalette[c(2, 4, 3)],
          frame.plot = FALSE,
          ylab = "Clustering Info Distance", ylim = c(0, max(df$dist)))
  
  
  spread <- cbind(
    mst = TreeDist::MeanMSTEdge(d, cluster = runID),
    nn = TreeDist::MeanNN(d, cluster = runID, Average = median),
    mad = TreeDist::DistanceFromMedian(d, cluster = runID, Average = median)
  )
  
  par(fig = c(0.05, .25, 0.8, 0.95), mar = c(1, 1, 0, 0),
              new = TRUE, cex = 0.75)
  
  barplot(spread[2:1, ], # Reversed as plotted bottom-to-top
          col = figPalette[3:2], # Also necessary to reverse
          beside = TRUE, horiz = TRUE, font.main = 1, cex.main = 1,
          axisnames = TRUE, main = "Dispersion\n", xpd = NA, las = 1)
  
  par(oPar)
  plot(cluster::silhouette(dist = d, runID),
       main = "Silhouette plot",
       col = figPalette[2:3])
  
  for(ntk in list(res1, res2)[c(model1, model2) %in% ntkModels]) {
    par(mfrow = 2:1, mar = c(2, 4, 2, 0))
    ntkLoss <- ntk[[if ("rate_dollo" %in% names(ntk)) "rate_dollo" else "rate_loss"]] 
    hist(log(ntkLoss), col = 6, main = "Log(Neomorphic: Rate of loss:gain)")
    abline(v = mean(log(ntkLoss)), col = "white")
    abline(v = quantile(log(ntkLoss), c(0.025, 0.5, 0.975)), col = "darkgrey")
    abline(v = 0, lty = "solid")
    abline(v = 0, lty = "dashed", col = "white")
    legend("topright", bty = "n",
           c("Median; 95% range", "Equal rates"),
           lty = c("solid", "dashed"),
           col = c("darkgrey", "black"))
    hist(log(ntk$rate_neo), col = 5, main = "Log(Neomorphic:Transformational rate)")
    abline(v = quantile(log(ntk$rate_neo),  c(0.025, 0.5, 0.975)), col = "darkgrey")
    abline(v = mean(log(ntk$rate_neo)), col = "white")
    abline(v = 0, lty = "solid")
    abline(v = 0, lty = "dashed", col = "white")
  }
  
  par(oPar)
  .TableRow <- function(model) {
    stoneFile <- StoneFile(pID, model)
    cbind(signif(read.table(ConvergenceFile(pID, model), col.names = 
                       c("Burnin", "PSRF", "ESS", "FrechESS", "mdPsESS")), 6),
          if (file.exists(stoneFile)) {
            suppressWarnings(tryCatch(
              signif(
                read.table(stoneFile,
                           col.names = c("steppingStn", "pathSmpl",
                                         "burnin", "stepGen", "nStep"))[, 1:2],
                6),
              error = function(e) {
                x <- signif(
                  read.table(stoneFile,
                             col.names = c("steppingStn", "pathSmpl")),
                  3)
                x[] <- paste("~", x)
                x
              }))
          } else {
            cbind("steppingStn" = NA, "pathSmpl" = NA)
          }
    )
  }
  conv1 <- .TableRow(model1)
  conv2 <- .TableRow(model2)
  plot.new()
  gridExtra::grid.table(
    `rownames<-`(rbind(conv1, conv2),
                 modelNames[c(model1, model2)]),
    theme = gridExtra::ttheme_default(base_size = 8)
  )
  dev.off()
  return(TRUE)
}
