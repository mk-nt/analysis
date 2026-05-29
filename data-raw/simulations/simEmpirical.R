library("TreeTools")
#library("RevBayes")
pkgload::load_all("../revbayesr")
pkgload::load_all()


### Set up simulation parameters ###

MkPath <- function(...) file.path("data-raw", "simulations", ...)
simIdxs <- 300:499
nRep <- length(simIdxs)

# Set to median values across empirical datasets
nTip <- ceiling(median(.meta$nTaxa))
# Multiply becuase datasets contain only informative characters -
# our simulated datasets will also contain variable but uninformative characters.
nNeo <- ceiling(median(.meta$nChar["neo", ]) * 1.56) # 1.48 better
nTrans <- ceiling(median(.meta$nChar["trans", ]) * 1.3) # 1.202 better
n <- 0.497
t <- 2.47
length <- 1.43 # median of medians inferred under Mk


### Simulate datasets ###

if (simulateDatasets <- FALSE) {
  for (seed in simIdxs) {
    simID <- sprintf("sim%03d", seed)
    outDir <- MkPath(simID)
    
    rb <- rbSession(c(PackageFile("rbScripts", "sim-by_nt_kv.Rev"), outDir,
                      nTip, nNeo, nTrans, n, t, length, seed))
    
    # NInformative(outDir, "neo.nex")
    # NInformative(outDir, "trans.nex")
  }
}

nInf <- sapply(simIdxs, function(seed) {
  simID <- sprintf("sim%03d", seed)
  outDir <- MkPath(simID)
  c(NInformative(outDir, "neo.nex"), NInformative(outDir, "trans.nex"))
})
apply(nInf, 1, median)

if (executeLocally <- FALSE) {
  for (script in c("sp_kv", "sp_nt_kv")) {
    if (!file.exists(MkPath(simID, paste0(script, ".trees")))) {
      rb <- rbSession(c(PackageFile("rbScripts", "sim-mc3.Rev"), outDir, script))
      rb$flush()
    }
  }
}

if (queueRemotely <- FALSE) {
  SQ() # Check current status of queue on remote
  for (i in simIdxs) {
    QueueSim(i, "sp_kv", myMem = 3000) # 2048 too small for sp_kv
    QueueSim(420, "sp_nt_kv", myMem = 2048*3, time = (2 * 60 + 25) * 60) # 4096 too small for half a dozen
  }
}

### Wait for analysis to complete before proceeding ###
### Wait for analysis to complete before proceeding ###
### Wait for analysis to complete before proceeding ###


log1 <- sapply(simIdxs, FetchLogIfMissing, logFile = "sp_nt_kv.p_run_1.log")
options("requeue" = FALSE)
log2 <- sapply(simIdxs, FetchLogIfMissing, logFile = "sp_nt_kv.p_run_2.log")
ok <- !vapply(log1, is.null, logical(1)) & !vapply(log2, is.null, logical(1))
if (any(!ok)) {
  message("Dropping sims with missing logs: ",
          paste(simIdxs[!ok], collapse = ", "))
}
stopifnot("No complete log pairs" = any(ok))
# Object format may not be in expected form if some results were downloaded fresh.
# Try running the log1 <- lines above again.
log1 <- sapply(log1[ok], identity)
log2 <- sapply(log2[ok], identity)
loss <- sapply(seq_along(log1["rate_loss", ]), function(i)
  summary(c(BurnOff(1 / log1["rate_loss", ][[i]], 0.25),
            BurnOff(1 / log2["rate_loss", ][[i]], 0.25)))
)
neo <- sapply(seq_along(log1["rate_neo", ]), function(i)
  summary(c(BurnOff(1 / log1["rate_neo", ][[i]], 0.25),
            BurnOff(1 / log2["rate_neo", ][[i]], 0.25)))
)
lng <- sapply(seq_along(log1["tree_length", ]), function(i)
  summary(c(BurnOff(log1["tree_length", ][[i]], 0.25),
            BurnOff(log2["tree_length", ][[i]], 0.25)))
)


if (file.exists(suppressWarnings(
    FetchLogIfMissing(simIdxs[[1]], logFile = "sp_kv.p_run_1.log")))
  ) {
  mkLog1 <- sapply(simIdxs, FetchLogIfMissing, logFile = "sp_kv.p_run_1.log")
  mkLog2 <- sapply(simIdxs, FetchLogIfMissing, logFile = "sp_kv.p_run_2.log")
  ok <- !vapply(mkLog1, is.null, logical(1)) & 
    !vapply(mkLog2, is.null, logical(1))
  if (any(!ok)) {
    message("Dropping sims with missing logs: ",
            paste(simIdxs[!ok], collapse = ", "))
  }
  stopifnot("No complete log pairs" = any(ok))
  mkLng <- sapply(seq_along(log1["tree_length", ]), function(i)
    summary(c(BurnOff(log1["tree_length", ][[i]], 0.25),
              BurnOff(log2["tree_length", ][[i]], 0.25)))
  )
} else {
  # Logs don't exist, calculate from tree files
  mkLng <- sapply(simIdxs, function(idx) {
    simID <- sprintf("sim%03d", idx)
    localPath <- MkPath(simID, sprintf("sp_kv_run_%s.trees", 1:2))
    if (!all(file.exists(localPath))) {
      message("Dropping sim with missing logs: ", idx)
      rep(NA_real_, 6)
    } else {
      trees1 <- read.tree(localPath[[1]])
      trees2 <- read.tree(localPath[[2]])
      summary(
        c(BurnOff(unname(colSums(sapply(trees1, `[[`, "edge.length"))), 0.25),
          BurnOff(unname(colSums(sapply(trees2, `[[`, "edge.length"))), 0.25))
      )
    }
  })
}

dev.new()
par(mar = c(8, 3, 0.1, 0.1))
PlotParamViolin(loss, neo, lng, mkLng, true_vals = c(n, t, length))
points(jitter(rep(1:3, each = ncol(loss))),
       log(c(loss["Median", ], neo["Median", ], lng["Median", ])),
       pch = 16, col = "#00000066", cex = 0.7)

if (!exists("ntErr")) {
  ntErr <- vector("list", nRep)
}
if (!exists("spErr")) {
  spErr <- vector("list", nRep)
}

for (seed in simIdxs) {
  s <- match(seed, simIdxs)
  simID <- sprintf("sim%03d", seed)
  trueTree <- read.tree(MkPath(simID, "tree.nwk"))
  
  if (is.null(ntErr[[s]])) {
    ntTrees1 <- SimTrees(simID, "sp_nt_kv", 1)
    ntTrees2 <- SimTrees(simID, "sp_nt_kv", 2)
    if(any(NTip(ntTrees1) == 1)) {
      warning("Corrupted tree file with ", simID, " at line ",
              which(NTip(ntTrees1) == 1),
              immediate. = TRUE)
    } else {
      ntTrees <- c(BurnOff(ntTrees1, 0.25), BurnOff(ntTrees2, 0.25))
      if (!is.null(ntTrees1)) {
        ntErr[[s]] <- TreeDist::ClusteringInfoDist(trueTree, ntTrees, normalize = TRUE)
      }
    }
  }
  
  if (is.null(spErr[[s]])) {
    spTrees1 <- SimTrees(simID, "sp_kv", 1)
    spTrees2 <- SimTrees(simID, "sp_kv", 2)
    spTrees <- c(BurnOff(spTrees1, 0.25), BurnOff(spTrees2, 0.25))
    if (!is.null(spTrees1)) {
      spErr[[s]] <- TreeDist::ClusteringInfoDist(trueTree, spTrees, normalize = TRUE)
    }
  }
}; cli::cli_progress_done()

stillNull <- vapply(ntErr, is.null, logical(1)) | vapply(spErr, is.null, logical(1))
sum(stillNull)

ntMeds <- vapply(ntErr[!stillNull], median, numeric(1))
spMeds <- vapply(spErr[!stillNull], median, numeric(1))
# Negative estimate means ntMeds decreases relative to spMeds
wilcox.test(ntMeds, spMeds, paired = TRUE, conf.int = TRUE)
summary(ntMeds)
summary(spMeds)

pCol <- ifelse(ntMeds - spMeds > 0, 3, ifelse(ntMeds - spMeds == 0, 1, 2))
plot(ntMeds - spMeds ~ spMeds, frame.plot = FALSE,
     xlab = "Tree divergence under sp_kv",
     ylab = "NT improvement on median tree distance",
     col = pCol, pch = 16
     )
abline(h = 0)

vioplot::vioplot(ntMeds - spMeds, frame.plot = FALSE,
                 ylab = "sp_nt_kv improvement over sp_kv; more is better")
abline(h = 0)


ntMns <- vapply(ntErr[!stillNull], mean, numeric(1))
spMns <- vapply(spErr[!stillNull], mean, numeric(1))
wilcox.test(ntMns, spMns, paired = TRUE, conf.int = TRUE)
plot(ntMns - spMns ~ spMns, frame.plot = FALSE,
     xlab = "Tree divergence under sp_kv",
     ylab = "NT improvement on mean tree distance",
     col = pCol, pch = 16
     )
abline(h = 0)

vioplot::vioplot(ntMns - spMns, frame.plot = FALSE,
                 ylab = "sp_nt_kv improvement over sp_kv; more is better")
abline(h = 0)


### Comparison of posterior splits ###

# ntPost / spPost are the posterior probabilities of each split in turn of the
# generative tree
# This is not obviously easy to interpret: it could be affected by the same
# class of problems as Robinson-Foulds distances.
if (!exists("ntPost")) {
  ntPost <- vector("list", nRep)
  prior <- vector("list", nRep)
  brLen <- vector("list", nRep)
}
if (!exists("spPost")) {
  spPost <- vector("list", nRep)
}


for (i in cli::cli_progress_along(ntPost)) {
  seed <- simIdxs[[i]]
  simID <- sprintf("sim%03d", seed)
  trueTree <- read.tree(MkPath(simID, "tree.nwk"))
  trueSplits <- as.Splits(trueTree)
  prior[[i]] <- signif(-SplitInformation(trueSplits) / log(2), 6)
  brLen[[i]] <- trueTree$edge.length[as.numeric(names(trueSplits))]
  
  if (is.null(ntPost[[i]])) {
    ntTrees1 <- SimTrees(simID, "sp_nt_kv", 1)
    ntTrees2 <- SimTrees(simID, "sp_nt_kv", 2)
    ntTrees <- c(BurnOff(ntTrees1, 0.25), BurnOff(ntTrees2, 0.25))
    if (!is.null(ntTrees1)) {
      ntPost[[i]] <- rowSums(vapply(as.Splits(ntTrees), function(sp) {
        trueSplits %in% sp
      }, logical(NSplits(trueSplits)))) / length(ntTrees)
    }
  }
  
  if (is.null(spPost[[i]])) {
    spTrees1 <- SimTrees(simID, "sp_kv", 1)
    spTrees2 <- SimTrees(simID, "sp_kv", 2)
    spTrees <- c(BurnOff(spTrees1, 0.25), BurnOff(spTrees2, 0.25))
    if (!is.null(spTrees1)) {
      spPost[[i]] <- rowSums(vapply(as.Splits(spTrees), function(sp) {
        trueSplits %in% sp
      }, logical(NSplits(trueSplits)))) / length(spTrees)
    }
  }
}; cli::cli_progress_done()

stillNull <- vapply(ntPost, is.null, logical(1)) | vapply(spPost, is.null, logical(1))
sum(stillNull)
unique_priors <- sort(unique(unlist(prior)))
n_groups      <- length(unique_priors)

ntBenefit <- unlist(ntPost[!stillNull]) - unlist(spPost[!stillNull])

wilcox.test(unlist(ntPost[!stillNull]), unlist(spPost[!stillNull]),
            paired = TRUE, conf.int = TRUE) #NSD

# plot(ntBenefit ~ unlist(prior[!stillNull]))
# plot(ntBenefit ~ unlist(brLen[!stillNull]), pch = ".", cex = 1.6,
#      col = "#00000055",
#      frame.plot = FALSE,
#      xlab = "Length of edge",
#      ylab = "Increase in posterior probability of true splits under NT model")
# abline(h = 0, col = 4, lty = "dotted")

if (!exists("ntAll")) {
  ntAll <- vector("list", nRep)
}
if (!exists("spAll")) {
  spAll <- vector("list", nRep)
}

### Posterior split sample entropy ###

nSample <- 200
for (i in cli::cli_progress_along(ntAll)) {
  seed <- simIdxs[[i]]
  simID <- sprintf("sim%03d", seed)
  trueTree <- read.tree(MkPath(simID, "tree.nwk"))
  trueSplits <- as.Splits(trueTree)
  
  if (is.null(ntAll[[i]])) {
    ntTrees1 <- SimTrees(simID, "sp_nt_kv", 1)
    ntTrees2 <- SimTrees(simID, "sp_nt_kv", 2)
    ntTrees <- c(ntTrees1, ntTrees2)
    if (!is.null(ntTrees1)) {
      if (length(ntTrees) >= nSample) {
        allS <- SplitFrequency(NULL, sample(ntTrees, nSample))
        counts <- attr(allS, "count")
        sInT <- allS %in% trueSplits
        ntAll[[i]] <- list(counts[sInT], counts[!sInT])
      } else {
        message("only ", length(ntTrees), " NT trees from ", simID)
      }
    }
  }
  
  if (is.null(spAll[[i]])) {
    spTrees1 <- SimTrees(simID, "sp_kv", 1)
    spTrees2 <- SimTrees(simID, "sp_kv", 2)
    spTrees <- c(spTrees1, spTrees2)
    if (!is.null(spTrees1)) {
      if (length(spTrees) >= nSample) {
        allS <- SplitFrequency(NULL, sample(spTrees, nSample))
        counts <- attr(allS, "count")
        sInT <- allS %in% trueSplits
        spAll[[i]] <- list(counts[sInT], counts[!sInT])
      } else {
        message("only ", length(spTrees), " Kv trees from ", simID)
      }
    }
  }
}; cli::cli_progress_done()


stillNull <- vapply(ntAll, is.null, logical(1)) |
  vapply(spAll, is.null, logical(1))
sum(stillNull)

resNT <- ntAll[!stillNull]
resSP <- spAll[!stillNull]

ayeNT <- lapply(resNT, `[[`, 1)
ayeSP <- lapply(resSP, `[[`, 1)
nayNT <- lapply(resNT, `[[`, 2)
naySP <- lapply(resSP, `[[`, 2)
allNT <- lapply(resNT, unlist)
allSP <- lapply(resSP, unlist)

ayePostNT <- vapply(ayeNT, sum, 1)
nayPostNT <- vapply(nayNT, sum, 1)

ayePostSP <- vapply(ayeSP, sum, 1)
nayPostSP <- vapply(naySP, sum, 1)

plot(ayePostNT - ayePostSP ~ ayePostSP, frame.plot = FALSE,
     ylab = sprintf("Additional true splits in %d trees from posterior sample under NT model", nSample))
abline(h = 0)
abline(h = mean(ayePostNT - ayePostSP), lty = "dotted", col = 3)
abline(h = median(ayePostNT - ayePostSP), lty = "dashed", col = 4)
wilcox.test(ayePostNT, ayePostSP, paired = TRUE, conf.int = TRUE)
# 95 percent confidence interval:
# 5.999985 27.499967
# sample estimates:
#   (pseudo)median 
# 16.99998 


extraFalses <- lengths(nayNT) - lengths(naySP)
plot(extraFalses ~ lengths(nayNT),
     frame.plot = FALSE,
     ylab = "Additional false splits under NT")
abline(h = 0)
abline(h = mean(extraFalses), lty = "dotted", col = 3)
abline(h = median(extraFalses), lty = "dashed", col = 4)
wilcox.test(lengths(nayNT), lengths(naySP), paired = TRUE, conf.int = TRUE) # NSD

hNT <- vapply(nayNT, TreeDist::Ntropy, double(1))
hSP <- vapply(naySP, TreeDist::Ntropy, double(1))

plot(hNT - hSP ~ hNT, frame.plot = FALSE,
     ylab = "Extra entropy of false splits under NT")
abline(h = 0)
abline(h = mean(hNT - hSP), lty = "dotted", col = 3)
abline(h = median(hNT - hSP), lty = "dashed", col = 4)
wilcox.test(hNT, hSP, paired = TRUE, conf.int = TRUE) # NSD


hNT <- vapply(allNT, TreeDist::Ntropy, double(1))
hSP <- vapply(allSP, TreeDist::Ntropy, double(1))

plot(hNT - hSP ~ hNT, frame.plot = FALSE,
     ylab = "Extra entropy of all splits under NT")
abline(h = 0)
abline(h = mean(hNT - hSP), lty = "dotted", col = 3)
abline(h = median(hNT - hSP), lty = "dashed", col = 4)
wilcox.test(hNT, hSP, paired = TRUE, conf.int = TRUE) # NSD


### How many false splits are afforded high probabilities? ###

qNT <- vapply(nayNT, quantile, double(3), c(0.9, 0.95, 0.99))
qSP <- vapply(naySP, quantile, double(3), c(0.9, 0.95, 0.99))

par(mfrow = c(2, 1), mar = c(4, 4, 0, 0))
hist(unlist(nayNT) / nSample * 100, xlab = "Posterior Pr / %", ylab = "False splits with this PP", main = "")
abline(v = apply(qNT / nSample * 100, 1, median, lty = "dotted"))
hist(unlist(naySP) / nSample * 100, xlab = "Posterior Pr / %", ylab = "False splits with this PP", main = "")
abline(v = apply(qSP / nSample * 100, 1, median, lty = "dotted"))

plot(qNT[1, ] - qSP[1 ,] ~ qNT[1, ], asp = 1, frame.plot = FALSE,
     xlim = range(qNT, qSP))
points(qNT[2, ] - qSP[2 ,] ~ qNT[2, ], col = 2)
points(qNT[3, ] - qSP[3 ,] ~ qNT[3, ], col = 3, pch = 3)
abline(h = 0)
rowMeans(qNT) - rowMeans(qSP)
wilcox.test(qNT["90%", ], qSP["90%", ], paired = TRUE, conf.int = TRUE) # NSD
wilcox.test(qNT["95%", ], qSP["95%", ], paired = TRUE, conf.int = TRUE) # NSD
wilcox.test(qNT["99%", ], qSP["99%", ], paired = TRUE, conf.int = TRUE) # Lower

tenNT <- vapply(nayNT, function(x) sort(x, decreasing = TRUE)[c(5, 10, 20)], double(3))
tenSP <- vapply(naySP, function(x) sort(x, decreasing = TRUE)[c(5, 10, 20)], double(3))


plot(tenNT[1, ] - tenSP[1 ,] ~ tenNT[1, ], asp = 1, frame.plot = FALSE,
     xlim = range(tenNT, tenSP))
points(tenNT[2, ] - tenSP[2 ,] ~ tenNT[2, ], col = 2)
points(tenNT[3, ] - tenSP[3 ,] ~ tenNT[3, ], col = 3, pch = 3)
abline(h = 0)
rowMeans(tenNT) - rowMeans(tenSP)


### Edge ratio comparison ###
if (!exists("ntEr")) {
  ntEr <- vector("list", nRep)
  trueEr <- vector("list", nRep)
}
if (!exists("spEr")) {
  spEr <- vector("list", nRep)
}

for (i in cli::cli_progress_along(ntEr)) {
  seed <- simIdxs[[i]]
  simID <- sprintf("sim%03d", seed)
  trueTree <- read.tree(MkPath(simID, "tree.nwk"))
  trueEr[[i]] <- EdgeRatio(trueTree)
  
  if (is.null(ntEr[[i]])) {
    ntTrees1 <- SimTrees(simID, "sp_nt_kv", 1)
    ntTrees2 <- SimTrees(simID, "sp_nt_kv", 2)
    if(any(NTip(ntTrees1) == 1)) {
      warning("Corrupted tree file with ", simID, " at line ",
              which(NTip(ntTrees1) == 1),
              immediate. = TRUE)
    } else {
      ntTrees <- c(BurnOff(ntTrees1, 0.25), BurnOff(ntTrees2, 0.25))
      if (!is.null(ntTrees1)) {
        ntEr[[i]] <- vapply(ntTrees, EdgeRatio, numeric(1))
      }
    }
  }
  
  if (is.null(spEr[[i]])) {
    spTrees1 <- SimTrees(simID, "sp_kv", 1)
    spTrees2 <- SimTrees(simID, "sp_kv", 2)
    spTrees <- c(BurnOff(spTrees1, 0.25), BurnOff(spTrees2, 0.25))
    if (!is.null(spTrees1)) {
      spEr[[i]] <-  vapply(spTrees, EdgeRatio, numeric(1))
    }
  }
}; cli::cli_progress_done()

stillNull <- vapply(ntEr, is.null, logical(1)) | vapply(spEr, is.null, logical(1))
sum(stillNull)

ntMeds <- vapply(ntEr[!stillNull], median, numeric(1))
spMeds <- vapply(spEr[!stillNull], median, numeric(1))
trMeds <- unlist(trueEr[!stillNull])
wilcox.test(ntMeds, spMeds, paired = TRUE, conf.int = TRUE) # NSD (p = 0.06)
#pCol <- ifelse(ntMeds - spMeds > 0, 3, ifelse(ntMeds - spMeds == 0, 1, 2))
plot(ntMeds ~ trMeds, frame.plot = FALSE, asp = 1, log = "xy",
     xlab = "True edge ratio",
     ylab = "Inferred edge ratio",
     col = 3, pch = 3
)
points(spMeds ~ trMeds, col = 2, pch = 3)
abline(0, 1)
