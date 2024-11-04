source("R/_setup.R")
source("R/CheckComplete.R")

# Define analytical IDs
projects <- AllProjects()
models <- c("by_ki", "by_n_ki", "by_t_ki", "by_nt_ki")

# Set plotting variables
modelPalette <- setNames(figPalette[c("yellow", "cyan", "red", "purple")],
                         models)
modelPalette["NA"] <- figPalette[["stone"]]
modelPch <- setNames(c(1, 22, 18, 24), models)

# Decipher
modelName <- c("by_ki" = "Mk-fixed",
               "by_n_ki" = "Mk-n",
               "by_t_ki" = "Mk-t",
               "by_nt_ki" = "Mk-nt")
eqModel <- c(" " = "by_ki",
             "n " = "by_n_ki",
             " t" = "by_t_ki",
             "n t" = "by_nt_ki",
             "NA" = "NA")

# Load results from cache
results <- setNames(
  lapply(projects, function(pID) setNames(lapply(models, function(model)
    ExistingResults(pID, model, approx = TRUE)), models)), projects)
marginals <- `colnames<-`(
  vapply(results, function(x) vapply(x, `[[`, double(1), "marginal"),
         setNames(double(4), models)), projects)

# Filter the 75 original datasets to remove those for which not all results are
# available.
# Some analyses did not complete due to computational limits or errors.
allIn <- colSums(is.na(marginals)) == 0
message(sum(allIn), " marginal likelihoods available")
message("Results not available for ",
        paste(projects[!allIn[projects]], collapse = ", "))

# Filter analyses by convergence and ESS
cnv <- results[[1]][[1]][["convergence"]]
convergence <- vapply(results, function(x) {
  y <- sapply(x, `[[`, "convergence")
  matrix(as.numeric(y), nrow = dim(y)[[1]], dimnames = dimnames(y))
}, matrix(0, 5, 4))

convergence["psrf", , ]
converged <- colSums(is.na(convergence["psrf", , ])) == 0
sum(converged)
hasML <- names(allIn)[allIn]
hasConv <- names(converged)[converged]
hasESS <- names(converged)[colSums(convergence["ess", , ] >= essThreshold, na.rm = TRUE) == 4]

intersect(hasML, hasConv)
mlEss <- intersect(hasML, hasESS)
length(hasESS)
length(mlEss)
# Use the subset of results that have converged, and for which marginal
# likelihood estimates are available
projOK <- names(results) %in% mlEss

# Use a lower threshold for tree ESS
ess2percent <- 156.25 # = 2% precision; Fabreti and Höhna 2022
hasFCE <- names(converged)[colSums(convergence[4, , ] >= ess2percent,
                                   na.rm = TRUE) == 4]
hasMPE <- names(converged)[colSums(convergence[5, , ] >= ess2percent,
                                   na.rm = TRUE) == 4]

hasTrees <- union(hasFCE, hasMPE) # generous: at least one measure indicates ESS
treesOK <- names(results) %in% hasTrees
length(intersect(mlEss, hasTrees))


# Question: How often is each model preferred?
bestModel <- rep(NA, length(allIn))
bestModel[allIn] <- models[apply(marginals[, allIn], 2, which.max)]
table(bestModel[projOK])[models]

taxon <- {
  meta <- read.csv("matrices/metadata.csv")
  factor(meta[["taxon"]],
         levels = c("Euarthropoda", "Invertebrates", "Vertebrates",
                    "Non-animals"))[
                      match(projects, meta[["project"]])]
}

.ModelBy <- function(bins, title = "") {
  tab <- table(data.frame(bestModel[projOK], bins[projOK]))[models, levels(bins)]
  atX <- seq(0, 1, length.out = length(models))
  atY <- seq(0, 1, length.out = length(levels(bins)))
  image(apply(tab, 2, function(x) x / sum(x)),
        asp = 1, frame.plot = FALSE, axes = FALSE,
        col = hcl.colors(255, "viridis", rev = FALSE)[-seq_len(40)])
  text(rep(atX, length(atY)), rep(atY, each = length(atX)), tab)
  t(tab)
  tabRows <- rownames(tab)
  tabRows[tabRows == "by_ki"] <- "Mk-fixed"
  modLabel <- gsub("_", "-", gsub("by", "Mk", gsub("_ki", "", tabRows)))
  axis(1, labels = modLabel, at = atX, las = 1, lty = 0, line = -0.5)
  axis(2, labels = gsub(", Inf]", "+", fixed = TRUE,
                        gsub("(-Inf,", "<", fixed = TRUE, colnames(tab))),
       at = atY, las = 2, lty = 0, line = -3)
  mtext(title, 2, cex = 0.8, line = 2)
}

.NTBy <- function(bins, param = "n", title = character(0)) {
  
  with <- switch(
    param,
    "n" = marginals[c("by_n_ki", "by_nt_ki"), projOK],
    "t" = marginals[c("by_t_ki", "by_nt_ki"), projOK]
  )
  without <- switch(
    param,
    "n" = marginals[c("by_ki", "by_t_ki"), projOK],
    "t" = marginals[c("by_ki", "by_n_ki"), projOK]
  )
  margin <- apply(with - without, 2, max)
  
  tab <- table(data.frame(cut(margin, c(-Inf, -1, 1, Inf)),
                          bins[projOK]))[, rev(levels(bins))]
  atX <- seq(0, 1, length.out = dim(tab)[[1]])
  atY <- seq(0, 1, length.out = dim(tab)[[2]])
  
  image(apply(tab, 2, function(x) x / sum(x)),
        asp = 1, frame.plot = FALSE, axes = FALSE,
        main = paste(title, "-", param, "parameter"),
        col = hcl.colors(255, "viridis", rev = FALSE))
  text(rep(atX, length(atY)), rep(atY, each = length(atX)), tab)
  t(tab)
  axis(1, labels = c("Rejected", "Indecisive", "Preferred"),
       at = atX, las = 2, lty = 0, line = -0.5)
  axis(2, labels = gsub(", Inf]", "+", fixed = TRUE,
                        gsub("(-Inf,", "<", fixed = TRUE, colnames(tab))),
       at = atY, las = 2, lty = 0, line = -2)
}

# Question: Does the preferred model depend on the taxon?
par(mfrow = c(3, 1))
.ModelBy(taxon, "Taxon")
.NTBy(taxon, "n", "Taxon")
.NTBy(taxon, "t", "Taxon")
## Result: Vertebrate datasets slightly prefer a t parameter; equal in others


.EvenBins <- function(x, n = 4) {
  cut(x, c(-Inf, quantile(x[allIn], seq_len(n - 1) / n), Inf))
}

# Question: Does the preferred model depend on the number of taxa?
nTaxa <- vapply(vapply(AllProjects(), MatrixFile, character(1)),
                .NTaxa, integer(1))
nTaxaBins <- .EvenBins(nTaxa)
.ModelBy(nTaxaBins, "n Taxa")
.NTBy(nTaxaBins, "n", "n Taxa")
.NTBy(nTaxaBins, "t", "n Taxa")
## Result: Not particularly

# Question: Does the preferred model depend on the number of characters?
nChar <- vapply(AllProjects(), function(pID) {
  c(
    trans = .NChar(MatrixFile(pID, "trans.nex")),
    neo = .NChar(MatrixFile(pID, "neo.nex"))
  )},
  integer(2))
  
nCharBins <- .EvenBins(colSums(nChar))
.ModelBy(nCharBins, "n Characters")
.NTBy(nCharBins, "n", "n Characters")
.NTBy(nCharBins, "t", "n Characters")
## Result: Subtly stronger support for n or t with more characters




# Question: Does the preferred model reflect the ratio of t:n characters?
tnRatio <- nChar["trans", ] / nChar["neo", ]
tnBins <- .EvenBins(tnRatio)
par(mfrow = c(2, 2))
.ModelBy(tnBins, "trans:neo characters")
.ModelBy(cut(tnRatio, c(0, .5, .8, 1/.8, 1/.5, 3, Inf)), "trans:neo characters")
## Result: Mk-nt preferred more often when lower t:n ratio

.NTBy(cut(tnRatio, c(0, .5, .8, 1/.8, 1/.5, 3, Inf)), "n", "trans:neo characters")
## Result: n preferred when proportionally more neomorphic characters
.NTBy(cut(tnRatio, c(0, .5, .8, 1/.8, 1/.5, 3, Inf)), "t", "trans:neo characters")
## Result: t preferred when neomorphic characters predominate (>2:1) 
## This is probably not statistically significant though.
## Could try a chi-squared test - probably not worthwhile.


GetPar <- function(model, parameter, projects = names(results),
                   statistic = "50%") {
  FUN.VALUE <- matrix(NA_real_, length(statistic), length(parameter),
                      dimnames = list(statistic, parameter))
  vapply(results[projects], function(x) {
    par <- x[[c(model, "parameters")]]
    if (is.null(par)) {
      FUN.VALUE
    } else {
      par[statistic, parameter]
    }
  }, FUN.VALUE)
}

mlEss <- intersect(hasML, hasESS)
lossN <- GetPar("by_n_ki", "rate_loss", proj = mlEss)
lossNT <- GetPar("by_nt_ki", "rate_loss", proj = mlEss)
neoT <- GetPar("by_t_ki", "rate_neo", proj = mlEss)
neoNT <- GetPar("by_nt_ki", "rate_neo", proj = mlEss)
# Question: What's the mean value of n and t?
.Ave <- function(x) {
  message(sprintf("%s: %.2f (%.2f); %.2f (%.2f)",
                  deparse(substitute(x)),
                  mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
                  median(x, na.rm = TRUE), mad(x, na.rm = TRUE)))
}
.Ave(lossNT)
.Ave(neoNT)
.Ave(lossN)
.Ave(neoT)

FigPlot <- function(..., pch = 3) {
  plot(..., frame.plot = FALSE, pch = pch)
}

# Safety check:
# Does the estimated values of one parameters differ when the other may vary?
par(mfrow = c(2, 1))
FigPlot(lossN, lossNT, asp = 1, log = "xy",
        col = c(durham$sky, durham$cyan, durham$ink)[
          1 + colSums(convergence["ess", c("by_n_ki", "by_nt_ki"), mlEss] >
                        2 * essThreshold)])
abline(0, 1, lty = "dotted")
## Result: the value of n is not strongly affected when t is a free parameter.

FigPlot(neoT, neoNT, asp = 1, log = "xy",
        col = c(durham$sky, durham$cyan, durham$ink)[
          1 + colSums(convergence["ess", c("by_t_ki", "by_nt_ki"), mlEss] >
                        2 * essThreshold)])
abline(0, 1, lty = "dotted")
## Result: The value of t is slightly lower when n is a free parameter.



# Safety check:
# Does this differ when one parameter is fixed to 1?
FigPlot(neoT, lossN, asp = 1, log = "xy",
        col = c(durham$sky, durham$cyan, durham$ink)[
          1 + colSums(convergence["ess", c("by_t_ki", "by_n_ki"), mlEss] >
                        2 * essThreshold)])
abline(h = 1, lty = "dashed", col = durham$yellow)
abline(v = 1, lty = "dashed", col = durham$yellow)
## Result: Not really; the previous figure is a suitable summary


# Question:
# How often does unity fall outside the likely range of values of n or t?
NotOne <- function(model, parameter) {
  ranges <- GetPar(model, parameter, mlEss, c("2.5%", "97.5%"))
  isOne <- ranges["2.5%", , ] < 1 & ranges["97.5%", , ] > 1 
  message(sum(!isOne, na.rm = TRUE), " / ", sum(!is.na(isOne)), " ",
          parameter, " values unlikely to be 1")
  !isOne
}

lossNot1NT <- NotOne("by_nt_ki", "rate_loss")
neoNot1NT <- NotOne("by_nt_ki", "rate_neo")
lossNot1 <- NotOne("by_n_ki", "rate_loss")
neoNot1 <- NotOne("by_t_ki", "rate_neo")

neither1 <- lossNot1 & neoNot1
message(sum(neither1, na.rm = TRUE), " / ", sum(!is.na(neither1)),
        " neither value likely to be 1")
notOnes <- apply(rbind(ifelse(is.na(lossNot1), "NA", ifelse(lossNot1, "n", "")),
                       ifelse(is.na(neoNot1), "NA", ifelse(neoNot1, "t", ""))
                       ),
                 2, paste, collapse = " ")
resOK <- names(results) %in% mlEss
bestModel1 <- bestModel[resOK][!is.na(lossNot1) & !is.na(neoNot1)]
notOnes[grepl("NA", notOnes)] <- "NA"
notOnes1 <- notOnes[!is.na(lossNot1) & !is.na(neoNot1)]
table(notOnes1)









lossNot1NT <- NotOne("by_nt_ki", "rate_loss")
neoNot1NT <- NotOne("by_nt_ki", "rate_neo")
notOnesNT <- apply(rbind(ifelse(is.na(lossNot1NT), "NA", ifelse(lossNot1NT, "n", "")),
                       ifelse(is.na(neoNot1NT), "NA", ifelse(neoNot1NT, "t", ""))
),
2, paste, collapse = " ")
notOnesNT[grepl("NA", notOnesNT)] <- "NA"

################################################################################
# Key figure:
# What values of n and t are inferred?
################################################################################
{
  
  oPar <- par(mfrow = c(4, 1), mar = c(6, 4, 0, 0))
  FigPlot(neoNT, lossNT, log = "xy",
          #        col = c(durham$sky, durham$ink)[1 + converged["by_nt_ki", ]],
          col = ifelse(modelPch[eqModel[notOnesNT]] > 18,
                       "#ffffff88", modelPalette[eqModel[notOnesNT]]),
          bg = modelPalette[eqModel[notOnesNT]],
          xlab = "t (Neomorphic : transformational rate ratio)",
          xlim = c(1 / max(c(neoNT, 1/neoNT)), max(c(neoNT, 1/neoNT))),
          ylab = "n (Loss : gain rate ratio)",
          ylim = c(1 / max(c(lossNT, 1/lossNT)), max(c(lossNT, 1/lossNT))),
          pch = modelPch[eqModel[notOnesNT]], cex = 1.2
  )
  #text(neoNT, lossNT, AllProjects(), cex = 0.5, col = modelPalette[eqModel[notOnesNT]])
  abline(h = 1, lty = "dashed", col = durham$concrete)
  abline(v = 1, lty = "dashed", col = durham$concrete)
  legend(2.5, 0.5, c("n \u2248 1, t \u2248 1",
                     "n \u2260 1, t \u2248 1",
                     "n \u2248 1, t \u2260 1",
                     "n \u2260 1, t \u2260 1"),
         pch = modelPch[models], col = ifelse(modelPch[models] > 18,
                                              "white", modelPalette[models]),
         pt.bg = modelPalette[models],
         y.intersp = 1.8,
         bty = "n", cex = 0.8, pt.cex = 1.2)
  ## Result: neo slower than transf; loss faster than gain



################################################################################
# Key figure
################################################################################
# Question: is n != 1 likely to favour a model with an n parameter?

  par(mar = c(4, 5, 0, 2))
  tab <- table(data.frame(notOnes1, bestModel1))
  grid <- t(tab[rev(c(" ", "n ", " t", "n t")), models])
  rownames(grid)[rownames(grid) == "by_ki"] <- "Mk-fixed"
  at <- seq(0, 1, length.out = length(models))
  image(grid, asp = 1, frame.plot = FALSE, axes = FALSE,
        col = hcl.colors(255, "viridis", rev = FALSE)[-seq_len(40)])
  axis(1, labels = gsub("_", "-", gsub("by", "Mk", gsub("_ki", "", rownames(grid)))),
       at = at, las = 1, line = -0.5, lty = 0)
  mtext("Parameters differing from 1", 2, cex = 0.8, line = 3)
  axis(2, labels = gsub(" ", "", colnames(grid)),
       at = at, las = 2, lty = 0)
  
  labs <- grid
  text(rep(at, 4), rep(at, each = 4), labs)
  
  .ModelBy(cut(nChar["neo", ] / nChar["trans", ], c(-Inf, 0.5, 0.8, 1.25, Inf)),
           "Neomorphic : transformational ratio")
  mtext("Preferred model", 1, cex = 0.8, line = 2)
  .ModelBy(cut(colSums(nChar), c(-Inf, 40, 80, 140, Inf)), "Number of characters")
  
  par(oPar)
}

table(neoNot1, bestModel[resOK])
table(lossNot1, bestModel[resOK])
## Result: n and t each differ from unity about 50% of the time, independently








# Safety check: when n beats k and t beats k, does nt beat k?
nPreferred <- marginals["by_n_ki", projOK] > marginals["by_ki", projOK]
tPreferred <- marginals["by_t_ki", projOK] > marginals["by_ki", projOK]
ntPreferred <- marginals["by_nt_ki", projOK] > marginals["by_ki", projOK]
oddities <- which(nPreferred & tPreferred & !ntPreferred)
## Result: Not always, but in the exceptions the difference between model
## log-likelihoods is ~2, which might be considered "within error",
## particularly given the variability between stepping stone estimates
## (e.g. https://www.biorxiv.org/content/10.1101/2024.09.25.614964v1.full.pdf)
 
for (pID in names(oddities)) {
  ml <- vapply(models, function(m) {
    stoneFile <- StoneFile(pID, m)
    if (file.exists(stoneFile)) {
      # Warning: Incomplete final line
      suppressWarnings(read.table(stoneFile)[[1]])
    } else {
      NA_real_
    }
  }, double(1))
  
  diffs <- outer(ml, ml, "-")
  grid <- t(diffs[4:1, ])
  at <- seq(0, 1, length.out = length(models))
  image(grid, asp = 1, frame.plot = FALSE, axes = FALSE,
        main = pID,
        col = hcl.colors(255, "Green-Orange", rev = TRUE))
  axis(1, labels = models, at = at, las = 2)
  axis(2, labels = rev(models), at = at, las = 2)
  
  labs <- round(grid, 2)
  labs[labs == 0] <- ""
  text(rep(at, 4), rep(at, each = 4), labs,
       font = (grid == max(grid)) + 1)
}


table(bestModel[projOK])
modelFreeN <- c("by_n_ki", "by_nt_ki")
modelFreeT <- c("by_t_ki", "by_nt_ki")
message("Preferred model has free _n_ in ",
        sum(nPreferred), " / ", sum(projOK), " projects")
message("Preferred model has free _t_ in ",
        sum(bestModel[projOK] %in% modelFreeT), " / ", sum(projOK), " projects")

apply(t(marginals[modelFreeN, names(nPreferred)[nPreferred]]) - marginals["by_ki", names(nPreferred)[nPreferred]], 1, max)

better <- "by_nt_ki"
worse <- "by_t_ki"
Support <- function(better, worse) {
  examine <- !is.na(bestModel) & bestModel == better & projOK
  table(cut(marginals[better, examine] - marginals[worse, examine],
            c(-Inf, 0, 1, 3, 5, Inf)))
}
nSupport <- Support("by_nt_ki", "by_t_ki") + Support("by_n_ki", "by_ki")
tSupport <- Support("by_nt_ki", "by_n_ki") + Support("by_t_ki", "by_ki")

sum(nSupport)
sum(tSupport)



# Question: Does using a preferred model change tree estimation?
.PrecisionIncrease <- function(x) {
  (x[[1]] - x[[2]]) / x[[1]]
}

## Distances
.TreeStats <- function(x, cf = rep_len("by_ki", length(x)), better = "best") {
  `colnames<-`(vapply(seq_along(x), function(i) {
    index <- x[[i]]
    pID <- AllProjects(index)
    baseModel <- cf[[i]]
    if (better == "best") {
      betterModel <- bestModel[[index]]
    } else {
      betterModel <- better
    }
    dFile <- DistanceFile(
      pID, baseModel, betterModel,
      final =  file.exists(DistanceFile(pID, baseModel, betterModel))
    )
    mdmd <- if (file.exists(TreeSampleFile(pID, baseModel)) && 
        file.exists(TreeSampleFile(pID, betterModel))) {
      TreeDist::ClusteringInfoDistance(
        median(ape::read.tree(TreeSampleFile(pID, baseModel))),
        median(ape::read.tree(TreeSampleFile(pID, betterModel))),
        normalize = TRUE
      )
    } else {
      NA_real_
    }
    
    if (file.exists(dFile)) {
      d <- readRDS(dFile)
      runID <- rep(1:2, rep(attr(d, "Size") / 2, 2))
      
      c(
        mdmd = mdmd,
        mad = .PrecisionIncrease(
          TreeDist::DistanceFromMedian(d, runID, Average = median)),
        mst = .PrecisionIncrease(TreeDist::MeanMSTEdge(d, runID)),
        nn = .PrecisionIncrease(TreeDist::MeanNN(d, runID, Average = median)),
        sil = mean(cluster::silhouette(dist = d, runID)[, 3])
      )
    } else {
      # If statistics have not been computed, running the output code will
      # create them.
      # AnalysisPDF() is defined in R/AnalysisPDF.R
      message("    AnalysisPDF(", pID, ", \"", baseModel, "\", \"",
              betterModel, "\");")
      rep(NA_real_, 5)
    }
  }, c(mdmd = NA_real_, mad = NA_real_, mst = 0, nn = 0, sil = 0)),
  AllProjects(x))
}

# Subquestion: Does precision correlate to the value of n?
nHelps <- bestModel %in% c("by_n_ki", "by_nt_ki") & treesOK
nDists <- .TreeStats(which(nHelps),
                     gsub("_n_", "_", gsub("_nt_", "_t_", bestModel[nHelps])))
{
nValue <- ifelse(bestModel[nHelps] == "by_nt_ki", lossNT, lossN)
FigPlot(nValue, 100 * nDists["mad", ], type = "n",
        log = "x",
        main = "Effect on tree topology: n is favoured",
        xlab = "n", xlim = range(c(1, nValue), na.rm = TRUE),
        ylab = "Precision increase (%)",
        ylim = range(c(0, 100 * nDists[c("nn", "mst", "mad"), ]), na.rm = TRUE))
abline(h = 0, lty = "dotted", col = durham$yellow)
abline(v = 1, lty = "dotted", col = durham$yellow)
legend("left", bty = "n", cex = 0.8,
       legend = c("NN", "MST", "MAD"), pch = c(1, 3, 4),
       col = figPalette[c(2, 4, 5)])
#points(nValue, 100 * nDists["nn", ], pch = 1, col = figPalette[[2]])
#points(nValue, 100 * nDists["mst", ], pch = 3, col = figPalette[[4]])
points(nValue, 100 * nDists["mad", ], pch = 4, col = modelPalette["by_n_ki"])
summary(lm(nDists["mad", ] ~ nValue))
}
## Result: No clear pattern


# Subquestion: Does silhouette score correlate to the value of n?
FigPlot(nValue, nDists["sil", ], type = "n", log = "x",
        main = "Effect on tree topology: n is favoured",
        xlab = "n", xlim = range(c(1, nValue), na.rm = TRUE),
        ylab = "Silhouette score")
abline(h = 0, lty = "dotted", col = durham$yellow)
abline(v = 1, lty = "dotted", col = durham$yellow)
points(nValue, nDists["sil", ], pch = 2, col = modelPalette["by_n_ki"])
summary(lm(nDists["sil", ] ~ nValue))
## Result: No clear pattern

# Subquestion: Does precision correlate to the value of t?
tHelps <- bestModel %in% c("by_t_ki", "by_nt_ki") & treesOK
tDists <- .TreeStats(which(tHelps),
                     gsub("_t_", "_", gsub("_nt_", "_n_", bestModel[tHelps])))
{
tValue <- ifelse(bestModel[tHelps] == "by_nt_ki", neoNT, neoT)
FigPlot(tValue, 100 * tDists["nn", ], type = "n",
        log = "x", xlab = "t", ylab = "Precision increase (%)",
        ylim = range(c(0, 100 * tDists[c("nn", "mst", "mad"), ]), na.rm = TRUE))
abline(h = 0, lty = "dotted", col = durham$yellow)
abline(v = 1, lty = "dotted", col = durham$yellow)
legend("topright", bty = "n", cex = 0.8,
       legend = c("NN", "MST", "MAD"), pch = c(1, 3, 4),
       col = figPalette[c(2, 4, 5)])
points(tValue, 100 * tDists["nn", ], pch = 1, col = figPalette[[2]])
points(tValue, 100 * tDists["mst", ], pch = 3, col = figPalette[[4]])
points(tValue, 100 * tDists["mad", ], pch = 4, col = figPalette[[5]])
summary(lm(tDists["mad", ] ~ tValue))
}
## Result: No clear pattern

# Subquestion: Does silhouette score correlate to the value of t?
{FigPlot(tValue, tDists["sil", ], type = "n", log = "x",
        main = "Effect on tree topology: t is favoured",
        xlab = "t", xlim = range(c(1, tValue), na.rm = TRUE),
        ylab = "Silhouette score")
abline(h = 0, lty = "dotted", col = durham$yellow)
abline(v = 1, lty = "dotted", col = durham$yellow)
points(tValue, tDists["sil", ], pch = 2, col = figPalette[[1]])
}## Result: No clear pattern


# Subquestion: Does precision correlate to silhouette score?
kiBeaten <- which(!is.na(bestModel) & bestModel != "by_ki" & treesOK)
treeStats <- .TreeStats(kiBeaten)
#FigPlot(treeStats[c("nn", "mst", "mad"), ] ~ treeStats[rep("sil", 3), ],
#        pch = c(1, 3, 4), col = figPalette[c(2, 4, 5)],

{
FigPlot(100 * treeStats["mad", ] ~ treeStats["sil", ],
        pch = 3, col = modelPalette[bestModel[kiBeaten]],
        xlab = "Silhouette score", ylab = "Precision increase (%)")
abline(h = 0, lty = "dotted", col = durham$yellow)
abline(v = 0, lty = "dotted", col = durham$yellow)
legend("topright", bty = "n", cex = 0.8, legend = models[-1], pch = 3,
       col = modelPalette[models[-1]])
} ## Result: silhouettes and precision uncorrelated; neither dependent on best model


# Subquestion: Does a more preferred model result in better trees?
modelDiff <- apply(marginals[, kiBeaten], 2, max) - marginals["by_ki", kiBeaten]
{
FigPlot(modelDiff, treeStats["sil", ],# type = "n",
        pch = 3, col = modelPalette[bestModel[kiBeaten]],
        xlab = "BF(best model, Mk)", ylab = "Silhouette score")
abline(h = 0, lty = "dotted", col = durham$yellow)
abline(v = 0, lty = "dotted", col = durham$yellow)
# text(modelDiff, treeStats["sil", ], AllProjects(kiBeaten),
     # col = modelPalette[bestModel[kiBeaten]], cex = 0.7, xpd = NA)
legend("bottomright", bty = "n", cex = 0.8, legend = models[-1], pch = 3,
       col = modelPalette[models[-1]])
## Result 1: Trees differ when model is preferred

FigPlot(modelDiff, 100 * treeStats["mad", ],
        pch = 3, col = modelPalette[bestModel[kiBeaten]],
        xlab = "BF(best model, Mk)", ylab = "Precision increase (%)")
abline(h = 0, lty = "dotted", col = durham$yellow)
abline(v = 0, lty = "dotted", col = durham$yellow)
legend("topright", bty = "n", cex = 0.8, legend = models[-1], pch = 3,
       col = modelPalette[models[-1]])
}
## Result 2: Precision shows no consistent pattern

.TRUE <- function(x) ifelse(is.na(x), FALSE, x)

withN <- .TRUE(marginals["by_n_ki", ] > marginals["by_ki", ])
nDiff <- marginals["by_n_ki", withN] - marginals["by_ki", withN]
# NOTE - comparing BEST model to Mk, not n to Mk
nStats <- .TreeStats(which(withN & treesOK),
                     cf = rep_len("by_ki", sum(withN)),
                     better = "by_n_ki")
withT <- .TRUE(marginals["by_t_ki", ] > marginals["by_ki", ])
tDiff <- marginals["by_t_ki", withT] - marginals["by_ki", withT]
tStats <- .TreeStats(which(withT & treesOK), better = "by_t_ki")

.Length <- function(pID, scriptID) {
  mean(vapply(ape::read.tree(TreeSampleFile(pID, scriptID)), function(tr) {
    sum(tr[["edge.length"]])
  }, numeric(1)))
}
nLength <- vapply(AllProjects()[withN], function(pID) {
  tryCatch(.Length(pID, "by_n_ki") / .Length(pID, "by_ki"),
           error = function(e) {warning(e); NA_real_})
}, double(1))
tLength <- vapply(AllProjects()[withT], function(pID) {
  tryCatch(.Length(pID, "by_t_ki") / .Length(pID, "by_ki"),
           error = function(e) {warning(e); NA_real_})
}, double(1))

################################################################################
{
oPar <- par(mfrow = c(4, 1))
supportLines <- c(3) # strong
# Broken down by n and t
FigPlot(nDiff[colnames(nStats)], nStats["sil", ],# type = "n",
        xlab = "Model support",
        pch = 0, col = modelPalette["by_n_ki"],
        ylab = "Silhouette score")
points(tDiff[colnames(tStats)], tStats["sil", ], pch = 5,
       col = modelPalette["by_t_ki"])
abline(h = 0, lty = "dotted", col = durham$concrete)
abline(v = supportLines, lty = "dotted", col = durham$concrete)
## Result 1: Trees differ when n is preferred
## Result 1: Trees don't really differ more when t is more strongly preferred


FigPlot(nDiff[colnames(nStats)], nStats["mdmd", ],# type = "n",
        pch = 0, col = modelPalette["by_n_ki"],
        xlab = "Model support",
        ylab = "Normalized CI Distance between median trees")
points(tDiff[colnames(tStats)], tStats["mdmd", ], pch = 5,
       col = modelPalette["by_t_ki"])
abline(h = 0, lty = "dotted", col = durham$concrete)
abline(v = supportLines, lty = "dotted", col = durham$concrete)


FigPlot(nDiff[colnames(nStats)], 100 * nStats["mad", ],
        pch = 0, col = modelPalette["by_n_ki"],
        xlab = "Model support",
        ylab = "Change in precision (MAD, %)",
        ylim = c(-1, 1) * max(abs((100 * nStats["mad", ])), na.rm = TRUE))
points(tDiff[colnames(tStats)], 100 * tStats["mad", ], pch = 5,
       col = modelPalette["by_t_ki"])
abline(h = 0, lty = "dotted", col = durham$concrete)
abline(v = supportLines, lty = "dotted", col = durham$concrete)
## Result 2: Precision shows no consistent pattern

FigPlot(nDiff, nLength,# type = "n",
        pch = 0, col = modelPalette["by_n_ki"],
        xlab = "Model support (Bayes Factor, vs. Mk-fixed)",
        ylab = "Relative tree length",
        ylim = c(1/2, 2),
        log = "y", xpd = NA)
points(tDiff, tLength, pch = 5, xpd = NA,
       col = modelPalette["by_t_ki"])
abline(h = 1, lty = "dotted", col = durham$concrete)
abline(v = supportLines, lty = "dotted", col = durham$concrete)
## Tree lengths differ by factor of < 2


legend("topright", bty = "n", legend = c("n", "t"),
       title = "Free parameter",
       pch = c(0, 5), col = modelPalette[c("by_n_ki", "by_t_ki")])
par(oPar)
}


## SUPPLEMENTARY TABLE ##
library("DT")
library("htmltools")
library("htmlwidgets")

LinkMB <- function(pID) {
  mbRoot <- "https://morphobank.org/index.php/Projects/Matrices/project_id"
  sprintf("<a href=\"%s/%s\" target=\"blank\">%s</a>", mbRoot, pID, pID)
}

LinkData <- function(pID) {
  ghRoot <- "http://github.com/mk-nt"
  paste(sprintf("<a href=\"%s/%s_%s\" target=\"blank\">%s</a>", ghRoot, pID, models,
                modelName[models]), collapse = " | ")
}

Sig3 <- function(n) signif(n, 3)
Round1 <- function(n) ifelse (n == 0, "*", sprintf("%.1f", n))
BF <- function(n) Round1(n - apply(marginals[, resOK], 2, max))
BFLink <- function(model) {
  ghRoot <- "http://github.com/mk-nt"
  sprintf("<a href=\"%s/%s_%s\" target=\"blank\">%s</a>", 
          ghRoot, namesOK, model, BF(marginals[model, namesOK]))
}

intro_html <- markdown::markdownToHTML(text = "
## Summary of data and analyses

- [neo-trans/matrices](https://github.com/neo-trans/matrices) contains
  original matrices downloaded
  from MorphoBank (`projectID.nex`).
  Each character is marked as neomorphic or
  transformational in the accompanying `ProjectID_name.xlsx` file.  Characters
  that require reformulation are labelled with a proposed reformulation.
  
- These raw files were used to procedurally generate separate matrices
  for neomorphic (`*.neo.nex`) and transformational (`*.trans.nex`) characters,
  and analytical scripts, each housed in a separate GitHub repository
  and analysed on the Hamilton high performance computing cluster.
  
  Each repository contains two RevBayes scripts:
  - `marginal.Rev`, used to compute the marginal likelihood of the model, and
  - `mcmcmc.Rev`, used to conduct MCMCMC analysis.
  
  The analytical model is defined in `by_XX.Rev`.
  The repository also contains logs of parameter estimates (`.log`),
  tree samples (`.trees`), and marginal likelihood estimates (`.pp`).

## Table explanation

- The \"Project ID\" column contains a link to the original study at
  MorphoBank.
  Links to the analytical files and results are provided in the \"BF\" columns
  of the table below.
  
- Inferred values of _n_ and _t_ are marked with an asterisk where the 95%
  posterior distribution does not include 1.
  
- Columns can be sorted by clicking the arrows, or filtered using the boxes
  at the column tops.

", fragment.only = TRUE)


namesOK <- names(results)[resOK]
{# Supplementary reporting
  overview <- data.frame(
  "ID" = LinkMB(names(results[resOK])),
  Taxon = taxon[resOK],
  Taxa = nTaxa[resOK],
  Characters = colSums(nChar)[resOK],
  "Transf." = nChar["trans", resOK],
  "Neom." = nChar["neo", resOK],
  "T:N ratio" = Sig3(tnRatio[resOK]),
  "Inferred n" = paste(Sig3(ifelse(bestModel[resOK] == "by_nt_ki",
                                    lossNT[namesOK], lossN[namesOK])),
                       ifelse(lossNot1[namesOK], "*", "")),
  "Inferred t" = paste(Sig3(ifelse(bestModel[resOK] == "by_nt_ki",
                                   neoNT[namesOK], neoT[namesOK])),
                       ifelse(neoNot1[namesOK], "*", "")),
  "Optimal model" = as.factor(modelName[bestModel[resOK]]),
  "Mk-fixed" = BFLink("by_ki"),
  "Mk-n" = BFLink("by_n_ki"),
  "Mk-t" = BFLink("by_t_ki"),
  "Mk-nt" = BFLink("by_nt_ki")
)

cols <- colnames(overview)
cols <- gsub(fixed = TRUE, "Transf", "Transf.",
        gsub(fixed = TRUE, "Characters", "Chars",
        gsub(fixed = TRUE, "Neom", "Neom.",
        gsub(fixed = TRUE, "T N ratio", "Ratio",
        gsub("Inferred ([nt])", "Inferred <i>\\1</i>", perl = TRUE,
             gsub("Mk ", "BF, Mk-", fixed = TRUE,
                  gsub(".", " ", cols, fixed = TRUE)))))))

widget <- 
  datatable(overview, escape = FALSE,
            filter = "top",
            options = list(pageLength = sum(resOK), dom = "t"),
            #caption = "Project overview",
            rownames = FALSE,
            colnames = cols) |>
    formatStyle("Taxa",
                background = styleColorBar(c(0, max(overview$Taxa)), 
                                           durham$concrete, angle = 270),
                backgroundSize = "98% 88%",
                backgroundRepeat = "no-repeat",
                backgroundPosition = "center"
                ) |>
  formatStyle("Transf.",
              background = styleColorBar(c(0, max(overview$Characters)), 
                                         durham$concrete, angle = 90),
              backgroundSize = "98% 88%",
              backgroundRepeat = "no-repeat",
              backgroundPosition = "center"
  ) |>
  formatStyle("Neom.",
              background = styleColorBar(c(0, max(overview$Characters)), 
                                         durham$concrete, angle = 270),
              backgroundSize = "98% 88%",
              backgroundRepeat = "no-repeat",
              backgroundPosition = "center"
  ) |>
  formatStyle(
      "Optimal.model",
      background = styleEqual(modelName, paste0(modelPalette[models], "44")),
      backgroundSize = "98% 88%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center"
    )

saveWidget(widget, "overview.html", selfcontained = TRUE,
           title = "Overview of results")
lines <- readLines("overview.html")
writeLines(
  append(lines, intro_html, which(grepl("id=\"htmlwidget_container", lines)) - 1),
  "overview.html")
}
