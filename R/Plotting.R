#' @importFrom stats quantile
#' @importFrom graphics points segments
#' @importFrom PlotTools SpectrumLegend
#' @export
ColErrPlot <- function(
    xLab, x_md, x_lo, x_hi,
    yLab, y_md, y_lo, y_hi,
    lgdPos = "topleft", lgdTitle = "",
    colBy, pal = hcl.colors(nBin, rev = TRUE), pch = 21,
    text = FALSE, nBin = 84,
    add = FALSE,
    asp = 1, xlim = c(0, 1), ylim = xlim, line = "1:1") {
  
  ord <- order(colBy)
  col <- if (is.numeric(colBy)) {
    lgd <- rev(floor(quantile(colBy)))
    cut(colBy, nBin)
  } else if (is.factor(colBy)) {
    nBin <- length(levels(colBy))
    lgd <- levels(colBy)[seq(nBin, 1, length.out = 4)]
    colBy
  } else {
    stop("`colBy` type not supported")
  }
  cols <- pal[col][ord]
  
  x_md <- x_md[ord]
  x_hi <- x_hi[ord]
  x_lo <- x_lo[ord]
  
  y_md <- y_md[ord]
  y_hi <- y_hi[ord]
  y_lo <- y_lo[ord]
  
  if (!isTRUE(add)) {
    plot(range(x_lo, x_hi, na.rm = TRUE),
         range(y_lo, y_hi, na.rm = TRUE),
         type = "n", xlab = xLab, ylab = yLab,
         asp = asp, frame.plot = FALSE,
         xlim = xlim, ylim = ylim)
  }
  if ("1:1" %in% line) {
    abline(0, 1, col = "grey70", lty = 2)
  }
  if ("h" %in% line) {
    abline(h = 0, col = "grey70", lty = 2)
  }
  if ("v" %in% line) {
    abline(v = 0, col = "grey70", lty = 2)
  }
  
  segments(x_lo, y_md, x_hi, y_md, col = cols, lty = "dotted")  # horizontal
  segments(x_md, y_lo, x_md, y_hi, col = cols, lty = "dotted")  # vertical
  points(x_md, y_md, pch = pch, bg = "white", col = cols)
  
  if (text) text(x_md, y_md, names(x_md), pos = 4, cex = 0.67, col = cols)
  if (!isTRUE(add)) {
    SpectrumLegend(lgdPos, palette = pal, title = lgdTitle,
                              legend = lgd, bty = "n")
  }
}

#' Spindle plot
#' 
#' `SpindlePlot()` creates a histogram arranged analogously to a violin plot
#' 
#' @param x data.frame or matrix in which columns correspond to distinct values
#' to plot.
#' @param \dots Additional arguments to `plot()`
#' @param cols Vector giving colour for each spindle
#' @param nBin Integer giving number of bins
#' @param width Numeric giving relative width of each spindle
#' @param log Logical specifying whether to log-transform the y axis
#' @param weights Numeric specifying weight of each observation, for a
#' weighted mean.
#' @param Behind Function to be executed before drawing the spindles
#' @export
SpindlePlot <- function(x, nBin = 20, width = 1, xlab = "",
                        log = FALSE, Behind = function() {},
                        clip, weights = NULL, ...) {
  centres <- seq_len(ncol(x))
  logY <- isTRUE(log) || (is.character(log) && "y" %in% unlist(strsplit(log, "")))
  Unlog <- if (logY) exp else function(x) x
  scripts <- colnames(x)
  clipped <- x
  if (!missing(clip)) {
    if (length(clip) == 1) {
      clip <- c(1 - clip, clip)
    }
    clipped[x > quantile(x, clip[[2]], na.rm = TRUE) | 
              x < quantile(x, clip[[1]], na.rm = TRUE)] <- NA
  }
  
  plot(NULL, axes = FALSE, xlab = xlab,
       xlim = c(0.5, ncol(x) + 0.5), ylim = range(clipped, na.rm = TRUE),
       log = if (logY) "y" else "", ...)
  axis(1, at = centres, labels = ModelLabel(scripts), las = 2, lwd = 0, line = -1)
  axis(2, las = 2)
  logged <- if (logY) log(clipped) else clipped
  
  Behind()
  
  histAll <- hist(as.matrix(logged), breaks = nBin * max(x, na.rm = TRUE) /
                    max(clipped, na.rm = TRUE),
                  plot = FALSE)
  for (i in seq_along(centres)) {
    hist_vals <- hist(logged[, i], breaks = histAll$breaks, plot = FALSE)
    mids <- hist_vals$mids
    dens <- hist_vals$counts / nrow(x) * width
    brks <- Unlog(hist_vals$breaks)
    
    # draw symmetrical rectangles
    for (j in seq_along(mids)) {
      if (dens[j] > 0) {
        rect(centres[i] - dens[j], brks[j],
             centres[i] + dens[j], brks[j + 1],
             col = ModelCol(scripts[[i]]),
             border = NA)
      }
    }
    
    if (!is.null(weights)) {
      # make sure weights are aligned with rows of x
      xBar <- Unlog(weighted.mean(logged[, i], weights, na.rm = TRUE))
    } else {
      xBar <- Unlog(mean(logged[, i], na.rm = TRUE))
    }
    xMed <- median(x[, i], na.rm = TRUE)
    xMad <- mad(x[, i], na.rm = TRUE)
    xIQR <- quantile(x[, i], na.rm = TRUE, c(0.25, 0.75))
    xSD <- sd(x[, i], na.rm = TRUE)
    segments(centres[i], xIQR[[1]],
             centres[i], xIQR[[2]],
             col = "#000000aa", lwd = 2)
    points(centres[i], xMed, pch = 21, col = 1, bg = "white", lwd = 2, cex = 1.4)
    points(centres[i] + 0.1, xBar, pch = 3, lwd = 2, col = 1, cex = 1.6)
  }
}

#' Plot a row of tree similarity spindles
#' 
#' A useful figure component
#' 
#' @inheritParams MakeSlurm
#' @param models Which script IDs should be included in the plot?
#' @param silID Phylopic uuid identifier for silhouette to add
#' @param Annotate1 Function to be called after first panel is drawn,
#' perhaps to add annotation or panel label
#' @importFrom graphics arrows
#' @importFrom rphylopic add_phylopic_base
#' @export
TreeSimSpindleRow <- function(pID, marginals, models, silID,
                              Annotate1 = function() {}) {
  rowCex <- 0.7
  par(cex = rowCex)
  
  clip <- switch(pID, "07206" = 0.9999, 1)
  width <- switch(pID, "07203" = 1.4, "07206" = 0.6, 0.9)
  
  TreeSimSpindle(pID, models, method = "MCI", width = width, clip = clip,
                 topMar = 2)
  Annotate1()
  lims <- par("usr")
  xWid <- diff(lims[1:2])
  xOffset <- xWid * 0.02
  yOffset <- diff(lims[3:4]) * 0.1
  silScale <- 2.84
  par(xpd = NA)
  add_phylopic_base(uuid = silID, x = lims[1] - xOffset, y = lims[3] - yOffset,
                    alpha = if (pID == "07204") 0.5 else 1,
                    width = silScale * xWid / length(models),
                    hjust = 1, vjust = 1)
  
  oPar <- par(mar = c(3.6, 3.2, 2.1, 0.2))
  heights <- setNames(marginals[models, pID] -
                        min(marginals[models, pID], na.rm = TRUE),
                      ModelLabel(models))
  stdErr <- attr(marginals, "stdErr")
  barErr <- .DiffErr(stdErr[models, pID],
                     stdErr[models, which.min(marginals[, pID])])
  bp <- barplot(heights, las = 2, col = ModelCol(models), border = NA,
                space = 0)
  
  belowZero <- heights - barErr < 0
  arrows(bp[!belowZero], heights[!belowZero] - barErr[!belowZero],
         bp[!belowZero], heights[!belowZero] + barErr[!belowZero],
         angle = 90, code = 3, length = 0.02, lwd = 1, xpd = NA)
  
  arrows(bp[belowZero], 0, bp[belowZero], heights[belowZero] + barErr[belowZero],
         angle = 90, code = 2, length = 0.02, lwd = 1, xpd = NA)
  
  noBF <- is.na(marginals[models, pID])
  text(seq_along(models)[noBF] - 0.5, par("usr")[[4]] * 0.2, "?",
       col = ModelCol(models[noBF]), font = 2)
  mtext("log Bayes Factor", 2, line = 2, cex = rowCex)
  par(oPar)
  
  # Quartet results are not materially different from MCI.
  # TreeSimSpindle(pID, models, method = "Quartet", width = width, clip = clip)
}

#' Comparative spindle plot of model Bayes factors
#'
#' `SpindleCompare()` draws side-by-side spindles comparing the marginal log
#' likelihoods of paired models. Each spindle shows the distribution of log Bayes
#' factors between the two models across datasets, with colours indicating which
#' model is favoured. Areas of non-material difference (where absolute log Bayes
#' factors are smaller than the propagated standard error) are shown in
#' semitransparent “ghost” shading.
#'
#' @param modelA Character vector naming the first set of models to compare.
#' @param modelB Character vector naming the second set of models to compare;
#'   must be the same length as `modelA`.
#' @param marginals Matrix of marginal log-likelihoods, with models as rows and
#'   datasets as columns.
#' @param stdErr Matrix of standard errors corresponding to `marginals`.
#' @param nBin Integer giving the number of histogram bins to use for each
#'   spindle (default `20`).
#' @param width Numeric scaling factor for the maximum spindle width (default `1`).
#' @param xlab Character label for the x-axis (default `"Model comparison"`).
#' @param ... Additional graphical parameters passed to [plot()].
#'
#' @details
#' Each spindle is constructed from the histogram of Bayes factors
#' \eqn{bf = marginals[modelA[i], ] - marginals[modelB[i], ]}.  
#' Bars where \eqn{|bf| < .DiffErr(stdErr[modelA[i], ], stdErr[modelB[i], ])}
#' are shown in semitransparent colour; other bars are fully opaque. Positive
#' values (favouring `modelA`) are coloured by `ModelCol(modelA)`, and negative
#' values (favouring `modelB`) by `ModelCol(modelB)`.
#'
#' @return
#' `SpindleCompare()` returns a matrix whose rows correspond to the entries
#' in `modelA` and whose columns count:
#' - `better`:  the number of times `modelA` outperforms `modelB`; 
#' - `worse`:  the number of times `modelA` is outperformed by `modelB`; 
#' - `n`:  the number of datasets compared, after removing `NA`s.
#' The function is called for its side-effect of producing a comparative spindle
#' plot.
#'
#' @export
SpindleCompare <- function(modelA, modelB, marginals, stdErr,
                           nBin = 14, width = 1, xlab = "",
                           ...) {
  stopifnot(length(modelA) == length(modelB))
  
  n <- length(modelA)
  centres <- seq_len(n)
  allBFs <- marginals[modelA, , drop = FALSE] - marginals[modelB, , drop = FALSE]
  
  # Compute global histogram breaks for consistent scaling
  global_range <- range(allBFs, na.rm = TRUE)
  raw_breaks <- seq(global_range[1], global_range[2], length.out = nBin)
  closest_to_zero <- raw_breaks[which.min(abs(raw_breaks))]
  shifted_breaks <- raw_breaks - closest_to_zero
  bin_width <- diff(shifted_breaks)[1]
  
  if (max(allBFs, na.rm = TRUE) > max(shifted_breaks)) {
    shifted_breaks <- c(shifted_breaks, max(shifted_breaks) + bin_width)
  } else if (min(allBFs, na.rm = TRUE) < min(shifted_breaks)) {
    shifted_breaks <- c(min(shifted_breaks) - bin_width, shifted_breaks)
  }
  
  plot(NULL, axes = FALSE, xlab = xlab, ylab = "log Bayes factor",
       xlim = c(0.5, n + 0.5), ylim = global_range, ...)
  text(seq_along(modelA), global_range[2], ModelLabel(modelA), pos = 2, xpd = NA)
  text(seq_along(modelB), global_range[1], ModelLabel(modelB), pos = 2, xpd = NA)
  axis(2, las = 2)
  abline(h = 0, col = "grey70", lty = "dashed")
  
  for (i in seq_along(modelA)) {
    scriptA <- modelA[[i]]
    scriptB <- modelB[[i]]
    bf <- allBFs[i, ]
    diffErr <- .DiffErr(stdErr[scriptA, ], stdErr[scriptB, ])
    
    hist_vals <- hist(bf, breaks = shifted_breaks, plot = FALSE)
    mids <- hist_vals$mids
    dens <- hist_vals$counts / length(bf) * width
    brks <- hist_vals$breaks
    
    # classify bins by sign and materiality
    material <- .OutwithError(abs(bf), diffErr)
    colFull <- ifelse(mids > 0, ModelCol(scriptA), ModelCol(scriptB))
    colGhost <- adjustcolor(colFull, alpha.f = 0.3)
    
    # first: ghost layer (all bins)
    rect(centres[i] - dens, brks[-length(brks)],
         centres[i] + dens, brks[-1],
         col = colGhost, border = NA)
    
    # second: opaque layer (only material bins)
    if (any(material, na.rm = TRUE)) {
      material_vals <- hist(bf[material], breaks = shifted_breaks, plot = FALSE)
      dens_mat <- material_vals$counts / length(bf) * width
      rect(centres[i] - dens_mat, brks[-length(brks)],
           centres[i] + dens_mat, brks[-1],
           col = colFull, border = NA)
    }
    
    # summary marks
    xMed <- median(bf, na.rm = TRUE)
    xIQR <- quantile(bf, c(0.25, 0.75), na.rm = TRUE)
    segments(centres[i], xIQR[[1]], centres[i], xIQR[[2]],
             col = "#000000aa", lwd = 2)
    points(centres[i], xMed, pch = 21, col = 1, bg = "white", lwd = 2, cex = 1.3)
  }
  
  nWin <- rowSums(.OutwithError(allBFs, .DiffErr(stdErr[modelA, ], stdErr[modelB, ])), na.rm = TRUE)
  nLose <- rowSums(.OutwithError(-allBFs, .DiffErr(stdErr[modelA, ], stdErr[modelB, ])), na.rm = TRUE)
  nA <- rowSums(!is.na(.OutwithError(allBFs, .DiffErr(stdErr[modelA, ], stdErr[modelB, ]))))
  cbind(better = nWin, worse = nLose, n = nA)
}


#' Mark meaningful differences
#' 
#' Adds lines to plots at position of `.config$eps`, i.e. the smallest
#' likelihood difference deemed to support one model over another.
#' 
#' @export
EpsLine <- function(horiz = TRUE, vert = FALSE, diag = FALSE) {
  if (horiz) {
    abline(h = 0)
    abline(h = c(.config$eps, -.config$eps), lty = "dotted")
  }
  if (vert) {
    abline(v = 0)
    abline(v = c(.config$eps, -.config$eps), lty = "dotted")
  }
  if (diag) {
    abline(a = 0, b = 1)
    abline(a = .config$eps, b = 1, lty = "dotted")
    abline(a = -.config$eps, b = 1, lty = "dotted")
  }
}

#' Plot values of _n_ from heterogeneous models
#' @param hgModelBF Numeric giving Bayes factor associated with each model,
#' compared to Mk model
#' @param hgN Numeric matrix with named rows giving 2.5, 25, 50, 75 and 97.5% 
#' quantiles for mean N value for each dataset
#' @param hgModel character specifying `scriptID` corresponding to the chosen
#' model.
#' 
#' @importFrom graphics segments
#' @export
HetNPlot <- function(hgN, hgModelBF, hgModel) {
  
  plot(hgN["50%", ], hgModelBF,
       type = "n", # just whilst texting
       xlab = paste0("Inferred mean 'n': ", modelLabel[hgModel]),
       ylab = "Model BF", xpd = NA,
       frame.plot = FALSE, pch = 4, log = "x", xlim = xLim)
  iqr <- log(hgN["75%", ]) - log(hgN["25%", ])
  segments(hgN["2.5%", ], hgModelBF, hgN["97.5%", ], xpd = NA,
           col = confidenceCol[cut(iqr, breaks = 256)])
  text(hgN["50%", ], hgModelBF, colnames(hgN), cex = 0.7, font = 2)
  abline(v = 1, lty = 2)
  abline(h = c(-1, 1) * epsBF, lty = 2)
  
  plot(hgModelBF / hgModelBF, hgModelBF, type = "n",
       frame.plot = FALSE,  log = "x", xlim = xLim,
       xlab = paste0("Characters per category: ", modelLabel[hgModel]),
       ylab = "Model BF")
  abline(v = 1,  lty = "dotted")
  for (pID in as.character(projects)) {
    param <- ExistingResults(pID, hgModel)[["parameters"]]
    if (is.null(param)) {
      message("No parameters found for ", pID)
      if (UpdateRecords(pID, hgModel)) {
        message(" > Found missing parameters.")
        param <- ExistingResults(pID, hgModel)[["parameters"]]
      } else {
        MakeSlurm(pID, hgModel, ml = FALSE)
        next
      }
    }
    beta25 <- fnDiscretizeBeta(param["25%", "beta_scale"], param["25%", "beta_scale"], 4)
    beta75 <- fnDiscretizeBeta(param["75%", "beta_scale"], param["75%", "beta_scale"], 4)
    probs <- param["50%", grep("matrix_probs.", fixed = TRUE, colnames(param))]
    segments((1 / beta25) - 1, hgModelBF[[pID]], (1 / beta75) - 1, hgModelBF[pID],
             #xpd = NA,
             lwd = probs * 10, lend = 1,
             col = gray.colors(256, rev = TRUE)[cut(probs, breaks = 256)])
  }
}

#' Heatmap of model comparisons
#' @param marginals,stdErr Matrix specifying marginal likelihood (`marginals`)
#' and standard error (`stdErr`) of each analysis.
#' Rows correspond to models; columns to datasets
#' @returns `ModelHeatmap()` is called for its side-effect of plotting a model
#' heatmap.
#' @importFrom grDevices colorRamp rgb
#' @export
ModelHeatmap <- function(marginals, stdErr) {
  models <- rownames(marginals)
  nModels  <- length(models)
  
  # initialise count matrix
  countMat <- matrix(0, nModels, nModels,
                     dimnames = list(models, models))
  
  # populate pairwise counts
  for (i in seq_len(nModels)) {
    for (j in seq_len(nModels)) {
      if (i != j) {
        diffOK <- .OutwithError(
          marginals[models[i], ] - marginals[models[j], ],
          .DiffErr(stdErr[models[i], ], stdErr[models[j], ])
        )
        countMat[i, j] <- sum(diffOK, na.rm = TRUE)
      }
    }
  }
  
  # scale each row by its maximum for colour intensity
  maxCount <- max(countMat)
  plotCols <- sapply(models, ModelCol)  # base colour for each row
  
  # create an empty plot
  par(mar = c(5, 5, 0.4, 0.4), cex = 0.8)
  image(1:nModels, 1:nModels, t(countMat[nModels:1, ]),
        col = "white", axes = FALSE, xlab = "", ylab = "",
        xlim = c(0.5, nModels + 0.5), ylim = c(0.5, nModels + 0.5))
  
  # manually draw coloured rectangles row by row
  for (i in seq_len(nModels)) {
    rowCounts <- countMat[i, ]
    rowCol <- plotCols[i]
    for (j in seq_len(nModels)) {
      intensity <- if (maxCount == 0) 0 else rowCounts[j] / maxCount
      fill <- if (intensity == 0) "white" else
        rgb(colorRamp(c("white", rowCol))(intensity), maxColorValue = 255)
      rect(j - 0.5, nModels - i + 0.5, j + 0.5, nModels - i + 1.5,
           col = fill, border = NA)
    }
  }
  
  # add labels and axes
  axis(1, at = 1:nModels, labels = ModelLabel(models), las = 2, tick = FALSE)
  axis(2, at = nModels:1, labels = ModelLabel(models), las = 2, tick = FALSE)
  
  # overlay count text with automatic contrast
  for (i in seq_len(nModels)) {
    for (j in seq_len(nModels)) {
      if (i == j) next
      val <- countMat[i, j]
      intensity <- if (maxCount == 0) 0 else val / maxCount
      txtCol <- if (intensity > 0.5) "white" else "black"
      text(j, nModels - i + 1, labels = val, col = txtCol, font = 2)
    }
  }
  
  #box()
  mtext("Model B", side = 1, line = 4)
  mtext("Model A", side = 2, line = 4)
}

#' Box plot of tree similarity
#' @inheritParams MakeSlurm
#' @importFrom cli cli_progress_done cli_progress_message
#' @export
TreeSimBox <- function(pID, models, output = NULL, sampleSize = 256,
                       method = "quartet") {
  
  outType <- substr(output, nchar(output) - 3, nchar(output))
  if (!dir.exists(dirname(output))) {
    dir.create(dirname(output)) 
  }
  switch(outType[1], # NA for length 0 if is.null(output)
         .pdf = pdf(sprintf(output, pID, method), width = 8, height = 10),
         .png = png(sprintf(output, pID, method), width = 800, height = 1000)
  )
  
  similarities <- TreeSimilarities(pID, models, method, sampleSize)
  
  cli_progress_message("Generating plot")
  oPar <- par(mar = c(6, 3.8, 0, 0), xpd = NA)
  boxplot(similarities, notch = TRUE, frame.plot = FALSE, las = 2,
          names = modelLabel[models],
          col = ModelCol(models),
          outline = FALSE, # Hide outliers
          ylab = sprintf("Similarity to well-corroborated tree (%s)", method))
  
  mtext(round(marginals7[models, pID] - 
                max(marginals7[models, pID], na.rm = TRUE), 2), 1, line = -4,
        at = seq_along(models), las = 2, cex = 0.9)
  points(vapply(similarities, mean, double(1)),
         pch = 3, col = "white", cex = 1.8)
  points(vapply(similarities, mean, double(1)),
         pch = 3, col = durham$ink, cex = 1.3)
  points(vapply(similarities, mean, double(1)),
         pch = 3, col = "white", cex = 0.05)
  
  cli_progress_done()
  print(setNames(vapply(similarities, median, 1), modelLabel[models]))
  print(setNames(vapply(similarities, mean, 1), modelLabel[models]))
  par(oPar)
  switch(outType, .pdf = dev.off(), .png = dev.off())
}

#' Spindle plot of tree similarity
#' @inheritParams MakeSlurm
#' @inheritParams SpindlePlot
#' @importFrom cli cli_progress_done cli_progress_message
#' @export
TreeSimSpindle <- function(pID, models, output = NULL, sampleSize = 256,
                           method = "quartet", nBin = 12, width = 1,
                           clip, topMar = 1) {
  similarities <- TreeSimilarities(pID, models, method, sampleSize)
  if (is.list(similarities)) {
    similarities <- do.call(cbind, similarities) |>
      `colnames<-`(models)
  }
  
  oPar <- par(mar = c(3.6, 4, topMar, 0.4), xpd = NA)
  SpindlePlot(similarities, width = width, nBin = nBin, clip = clip,
              ylab = sprintf("Similarity / normalized %s", method)
              #ylab = sprintf("%s similarity to well-corroborated tree", method)
              )
  par(oPar)
}

#' @importFrom Quartet ManyToManyQuartetAgreement QuartetDivergence
#' @export
PairwiseQuartets <- function(trees, Measure = QuartetDivergence) {
  Measure(ManyToManyQuartetAgreement(trees))
}

#' @rdname TreeSimBox
#' @importFrom cli cli_progress_message cli_progress_done
#' @importFrom graphics contour points
#' @importFrom MASS kde2d
#' @importFrom Quartet QuartetStatus SimilarityToReference
#' @importFrom Rogue QuickRogue
#' @importFrom TreeTools Consensus DropTip KeepTip LabelSplits RootTree SortTree 
#' SplitFrequency SupportColour TipLabels
#' @export
TreeSimPlot <- function(pID, models, pdf = FALSE, sampleSize = 128) {
  cli_progress_message(paste0(pID, ": Reading trees"))
  treeSets <- lapply(models, ComparisonTrees, pID)
  models <- models[lengths(treeSets) >= sampleSize]
  treeSets <- lapply(treeSets[lengths(treeSets) >= sampleSize], sample, sampleSize)
  wcTree <- wcTrees[[pID]] |>
    KeepTip(TipLabels(treeSets[[1]][[1]])) |>
    RootTree(outgroup[[pID]])
  
  cli_progress_message(paste0(pID, ": Computing distances"))
  d <- PairwiseQuartets(c(do.call(c, treeSets), list(wcTree)))
  pchType <- c(lengths(treeSets), 1)
  
  cli_progress_message("Finding rogues")
  rogues <- lapply(treeSets, QuickRogue)
  if (any(lengths(rogues)) > 5) {
    message("Found rogue taxa")
  }
  wcOrder <- TipLabels(wcTree)
  cli_progress_message("Computing consensi")
  cons <- lapply(treeSets, function(trees) {
    trees |> DropTip(rogues[[1]]$taxon[-1]) |>
      RootTree(outgroup[[pID]]) |>
      Consensus(p = 0.5) |>
      SortTree("TipLabels", order = wcOrder)
  })
  
  if (!isFALSE(pdf)) {
    if (!dir.exists(dirname(pdf))) {
      dir.create(dirname(pdf))
    }
    pdf(pdf, width = 8, height = 10)
  }
  
  nModels <- length(models)
  layout(if (nModels > 4) {
    if (nModels %% 2 == 0) {
      rbind(t(matrix(seq_len(nModels), ncol = 2)),
            c(rep(nModels + 1, nModels / 2 - 2), nModels + 2:3))
    } else {
      rbind(
        matrix(c(seq_len(nModels), nModels + 2), byrow = TRUE, nrow = 2),
        c(rep(nModels + 1, (nModels - 1) / 2), nModels + 3)
      )
    }
  } else if (nModels == 4) {
    rbind(1:3, c(4, 6, 6), c(5, 5, 7))
  } else if (nModels == 3) {
    rbind(1:3, 4:6)
  } else if (nModels == 2) {
    rbind(c(1, 2, 5), c(3, 3, 4))
  } else {
    rbind(c(1, 3), c(2, 4))
  })
  oPar <- par(mar = rep(0.5, 4))
  lapply(seq_along(models), function(i) {
    tr <- cons[[i]]
    splitFreqs <- SplitFrequency(tr, treeSets[[i]]) / sampleSize
    edgeSupport <- rep(1, nrow(tr$edge))
    childNode <- tr$edge[, 2]
    edgeSupport[match(names(splitFreqs), childNode)] <- splitFreqs
    plot(tr, tip.col = tipCol[tr[["tip.label"]]], main = models[[i]], xpd = NA,
         edge.col = SupportColour(edgeSupport), edge.width = 3)
    LabelSplits(tr, round(splitFreqs * 100, 1), unit = "%", col = SupportColour(splitFreqs),
                frame = "none", pos = 3)
  })
  
  .Median <- function(d) {
    which.min(unname(colSums(as.matrix(d))))
  }
  
  cum <- cumsum(lengths(treeSets))
  treeIndices <- rbind(cum - cum[[1]] + 1, cum + 0)
  medTree <- apply(treeIndices, 2, function(i) .Median(d[i[1]:i[2], i[1]:i[2]]))
  
  cli_progress_message("Comparing to reference")
  similarities <- lapply(treeSets, function(trees) {
    SimilarityToReference(QuartetStatus(trees, wcTree), normalize = TRUE)
  })
  medIndex <- treeIndices[1, ] + medTree - 1
  medSim <- d[medIndex, sum(1, lengths(treeSets))]
  
  cli_progress_message("Mapping distances")
  mapping <- cmdscale(d, add = TRUE)[["points"]][, 1:2]
  
  cli_progress_message("Generating plots")
  plot(mapping, asp = 1,
       col = rep(c(paste0(figPalette[1 + 1:length(models)], "88"), durham[["ink"]]),
                 pchType),
       pch = rep(c(seq_len(length(models) + 1)[-3], 3), pchType),
       cex = rep(c(rep(1, length(models)), 10) *  0.8, pchType),
       frame.plot = FALSE, axes = FALSE, ylab = "", xlab = "")
  for (i in seq_along(models)) {
    # add contour lines according to density of points of this model
    xy <- mapping[treeIndices[1, i]:treeIndices[2, i], , drop = FALSE]
    density <- kde2d(xy[, 1], xy[, 2],
                     lims = c(range(mapping[, 1]), range(mapping[, 2])))
    contour(density, add = TRUE, drawlabels = FALSE, lty = "dashed",
            col = paste0(figPalette[1 + i], "88"))
  }
  points(mapping[medIndex, 1], mapping[medIndex, 2],
         col = figPalette[1 + 1:length(models)],
         pch = "M")
  points(mapping[1 + sum(lengths(treeSets)), 1],
         mapping[1 + sum(lengths(treeSets)), 2], pch = 3, cex = 8)
  
  par(mar = c(6, 2.5, 0, 0), xpd = NA)
  boxplot(similarities, notch = TRUE, frame.plot = FALSE, las = 2,
          names = NULL,
          col = figPalette[1 + 1:length(models)],
          outline = FALSE, # Hide outliers
          ylab = "Similarity to well-corroborated tree")
  
  mtext(paste(modelLabel[models], round(marginals7[models, pID], 2)), 1,
        line = -2, at = seq_along(models), las = 2, cex = 0.7)
  points(vapply(similarities, mean, double(1)),
         pch = 3, col = "white", cex = 1.8)
  points(vapply(similarities, mean, double(1)),
         pch = 3, col = durham$ink, cex = 1.3)
  
  par(mar = rep(0.5, 4))
  plot(wcTree, main = "Well-corroborated tree",
       tip.col = tipCol[wcTree[["tip.label"]]])
  cli_progress_done()
  par(oPar)
  if (!isFALSE(pdf)) {
    dev.off()
  }
}


#' Explore whether differences in tree topology are deep or shallow
#' @importFrom cli cli_progress_done cli_progress_message
#' @importFrom TreeDist ClusteringInfoDistance
#' @importFrom TreeTools KeepTip RandomTree RootTree
#' @export
DiffDepth <- function(pID, models, pdf = FALSE, sampleSize = 128,
                      resamplings = 1024 * 4) {
  cli_progress_message(paste0(pID, ": Reading trees"))
  treeSets <- lapply(models, ComparisonTrees, pID)
  models <- models[lengths(treeSets) >= sampleSize]
  treeSets <- lapply(treeSets[lengths(treeSets) >= sampleSize], sample, sampleSize)
  wcTree <- wcTrees[[pID]] |>
    KeepTip(TipLabels(treeSets[[1]][[1]])) |>
    RootTree(outgroup[[pID]])
  
  if (!isFALSE(pdf)) {
    if (!dir.exists(dirname(pdf))) {
      dir.create(dirname(pdf))
    }
    pdf(pdf, width = 6, height = 6)
  }
  
  cli_progress_message("Resampling to benchmark")
  oPar <- par(mar = c(4.2, 4.2, 0.4, 0.4))
  wcScores <- attr(ClusteringInfoDistance(wcTree, wcTree, reportMatching = TRUE),
                   "matchedScores")
  rTrees <- replicate(resamplings, RandomTree(wcTree), simplify = FALSE)
  randomScores <- vapply(rTrees, function(tr) {
    attr(ClusteringInfoDistance(wcTree, tr, reportMatching = TRUE),
         "matchedScores")
  }, wcScores)
  
  cli_progress_message("Comparing trees")
  uqScores <- sort(unique(wcScores))
  quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  scoreRank <- vapply(treeSets, function(trees) {
    matchedScores <- vapply(trees, function(tr) {
      attr(ClusteringInfoDistance(wcTree, tr, reportMatching = TRUE),
           "matchedScores")
    }, wcScores)
    
    ranks <- t(vapply(1:nrow(matchedScores), function(i) {
      rowMeans(outer(matchedScores[i, ], randomScores[i, ], ">="))
    }, double(ncol(matchedScores))))
    
    vapply(uqScores, function(score) quantile(ranks[wcScores == score, ],
                                              quants), quants)
  }, matrix(0, length(quants), length(uqScores)))
  
  plot(scoreRank["50%", , 1] ~ uqScores,
       xlim = range(uqScores), ylim = c(min(scoreRank["50%", , ]), 1),
       xlab = "Node entropy", ylab = "Similarity to WC tree / Quantile",
       type = "n", frame.plot = FALSE)
  
  
  for (i in seq_along(models)) {
    x <- uqScores + ((i - 2) / 400)
    lines(x, scoreRank["50%", , i], lty = "dashed",
          type = "l", col = ModelCol(models[[i]]))
    points(x, scoreRank["50%", , i],
           pch = i, col = ModelCol(models[[i]]))
    # for (iqr in c("25%", "75%")) {
    #   lines(scoreRank[iqr, , i] ~ uqScores,
    #        col = figPalette[[i]], lty = "dotted")
    # }
  }
  
  legend("bottomleft", legend = modelLabel[models], bty = "n",
         lty = "dashed", pch = seq_along(models),
         col = ModelCol(models))
  
  cli_progress_done()
  par(oPar)
  if (!isFALSE(pdf)) {
    dev.off()
  }
}

#' Add label for plot region
#' @param i Index of label (1 = "(a)", or string to print
#' @export
Panel <- function(i, xOffset = 3, yOffset = 0) {
  if (is.numeric(i)) {
    i <- paste0("(", letters[i], ")")
  }
  usr <- par("usr")
  xOffset <- xOffset * strwidth("M")
  yOffset <- yOffset * strheight("M")
  
  x <- usr[[1]] - xOffset
  if (par("xlog")) {
    x <- 10 ^ x
  }
  
  y <- usr[[4]] - yOffset
  if (par("ylog")) {
    y <- 10 ^ y
  }
  
  text(x, y, i, xpd = NA, adj = c(1, 1)) # adj: right-align | top-align
}

#' Output a plot to file or graphics device
#' @param figName Figure filename (without suffix)
#' @param width,height Figure dimensions in inches
#' @param Plot function to execute, which should produce a plot.
#' 
#' @return `OutputPlot()` is called for its side-effecits: it launches a new
#' device of the type specified by `getOption("ntOutput")`
#' (supported: "png", "pdf", "device", where the latter
#' is the default and plots to the active graphics device).
#' @export
OutputPlot <- function(figName, width, height, Plot) {
  outputType <- trimws(tolower(getOption("ntOutput") %||% "device"))
  switch(outputType,
         pdf = cairo_pdf(file.path(OutputDir(), "figures",
                                   sprintf("%s.pdf", figName)),
                         width, height),
         png = png(file.path(OutputDir(), "figures",
                             sprintf("%s.png", figName)),
                   width, height, units = "in", res = 300))
  if (outputType != "device") {
    message("Creating ", sprintf("%s.%s", figName, getOption("ntOutput")))
    on.exit(dev.off())
  }
  Plot()
}
