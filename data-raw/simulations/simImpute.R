library("TreeTools")
#library("RevBayes")
pkgload::load_all()
pkgload::load_all("../revbayesr")


MkPath <- function(...) file.path("data-raw", "simulations", ...)
simIdxs <- 300:499


if (FALSE) {
  # Simulate ambiguated datasets
  for (seed in simIdxs) {
    simID <- sprintf("sim%03d", seed)
    outDir <- MkPath(simID)
    
    set.seed(seed)
    .Ambiguate("neo.nex")
    set.seed(seed)
    .Ambiguate("trans.nex")
  }
  
  
  # Execute locally
  for (script in c("sp_kv", "sp_nt_kv")) {
    if (!file.exists(MkPath(simID, paste0("imp_", script, ".trees")))) {
      rb <- rbSession(c(PackageFile("rbScripts", "imp-mc3.Rev"),
                        MkPath(simID), script))
      rb$flush()
    }
  }
  
  SQ()
  for (i in simIdxs) {
    QueueSim(i, "sp_kv", impute = TRUE, myMem = 4096) # 3000 too small
    QueueSim(i, "sp_nt_kv", impute = TRUE, myMem = 4096) # 3000 too small
  }
}
### Wait for analysis to complete before proceeding ###
### Wait for analysis to complete before proceeding ###
### Wait for analysis to complete before proceeding ###



### Imputation accuracy ###

nTip <- 54
nStates1 <- sapply(sprintf("sim%03d", simIdxs), function(simID) {
  tryCatch(read.table(MkPath(simID, "imp_sp_nt_kv.neo_run_1.states"), header = TRUE)[, -1], error = function(e) message(simID, ": ", e))
})
nStates2 <- sapply(sprintf("sim%03d", simIdxs), function(simID) {
  read.table(MkPath(simID, "imp_sp_nt_kv.neo_run_2.states"), header = TRUE)[, -1]
})
tStates1 <- sapply(sprintf("sim%03d", simIdxs), function(simID) {
  read.table(MkPath(simID, "imp_sp_nt_kv.trans_run_1.states"), header = TRUE)[, -1]
})
tStates2 <- sapply(sprintf("sim%03d", simIdxs), function(simID) {
  read.table(MkPath(simID, "imp_sp_nt_kv.trans_run_2.states"), header = TRUE)[, -1]
})
vStates1 <- sapply(sprintf("sim%03d", simIdxs), function(simID) tryCatch(
  read.table(MkPath(simID, "imp_sp_kv_run_1.states"), header = TRUE)[, -1],
  error = function(e) {
    message(simID, ": ", e)
    if (isTRUE(getOption("requeue", FALSE))) {
      QueueSim(as.integer(substr(simID, 4, 6)), "sp_kv", impute = TRUE, myMem = 4096)
    }
  }))
vStates2 <- sapply(sprintf("sim%03d", simIdxs), function(simID) {
  read.table(MkPath(simID, "imp_sp_kv_run_2.states"), header = TRUE)[, -1]
}) 

nActual <- lapply(sprintf("sim%03d", simIdxs), function(simID) {
  do.call(rbind, read.nexus.data(MkPath(simID, "neo.nex")))
})
tActual <- lapply(sprintf("sim%03d", simIdxs), function(simID) {
  do.call(rbind, read.nexus.data(MkPath(simID, "trans.nex")))
})
vActual <- Map(cbind, nActual, tActual)

nAmbig <- lapply(sprintf("sim%03d", simIdxs), function(simID) {
  do.call(rbind, read.nexus.data(MkPath(simID, "imp_neo.nex"))) == "?"
})
tAmbig <- lapply(sprintf("sim%03d", simIdxs), function(simID) {
  do.call(rbind, read.nexus.data(MkPath(simID, "imp_trans.nex"))) == "?"
})
vAmbig <- Map(cbind, nAmbig, tAmbig)

accuracy <- sapply(simIdxs, function(sim) {
  s <- match(sim, simIdxs)
  key <- .DecodeTips(nTip)
  
  simNActual <- nActual[[s]]
  simTActual <- tActual[[s]]
  simVActual <- vActual[[s]]
  
  simNAmbig <- nAmbig[[s]]
  simTAmbig <- tAmbig[[s]]
  simVAmbig <- vAmbig[[s]]
  
  n1 <- nStates1[key, s] # list; entry i = csv of char states at tip i
  n2 <- nStates2[key, s]
  t1 <- tStates1[key, s]
  t2 <- tStates2[key, s]
  v1 <- vStates1[key, s]
  v2 <- vStates2[key, s]
  
  nInv <- .Invariant(sim, "imp_neo.nex")
  tInv <- .Invariant(sim, "imp_trans.nex")
  vInv <- c(nInv, attr(nInv, "nChar") + tInv)
  
  .Accuracy <- function(s1, simAmbig, simActual, invariant) {
    lapply(seq_len(nTip), function(i) {
      # Reconstructed state of each character at tip i
      # row 1: p(char j in state 0); row 2: p(j in state 2) 
      recon <- (do.call(cbind, lapply(paste0(s1[[i]], " "), strsplit, ",") |>
        unlist(recursive = FALSE, use.names = FALSE)) |>
        `mode<-`("integer") + 1) |>
        apply(1, tabulate, 2)
      # RevBayes places all invariant characters at the end;
      # we wish to retain the original numbering
      recon[, -invariant] <- recon[, seq_len(dim(recon)[[2]] - length(invariant))]
      recon[, invariant] <- NA
      recon <- t(t(recon) / colSums(recon))
      
      ambigIJ <- simAmbig[sprintf("tip_%d", i), ]
      ones <- simActual[sprintf("tip_%d", i), ambigIJ] == "1"
      recon[, ambigIJ][cbind(1 + ones, seq_len(sum(ambigIJ)))]
    })
  }
  
  list(
    n1 = .Accuracy(n1, simNAmbig, simNActual, nInv),
    n2 = .Accuracy(n2, simNAmbig, simNActual, nInv),
    t1 = .Accuracy(t1, simTAmbig, simTActual, tInv),
    t2 = .Accuracy(t2, simTAmbig, simTActual, tInv),
    v1 = .Accuracy(v1, simVAmbig, simVActual, vInv),
    v2 = .Accuracy(v2, simVAmbig, simVActual, vInv)
  )
})

accN  <- sapply(seq_along(simIdxs), function(s)
  .MeanAcc(.PoolRuns(accuracy[["n1", s]], accuracy[["n2", s]])))
accT  <- sapply(seq_along(simIdxs), function(s)
  .MeanAcc(.PoolRuns(accuracy[["t1", s]], accuracy[["t2", s]])))
accNT <- sapply(seq_along(simIdxs), function(s)
  .MeanAcc(c(.PoolRuns(accuracy[["n1", s]], accuracy[["n2", s]]),
             .PoolRuns(accuracy[["t1", s]], accuracy[["t2", s]]))))
accV  <- sapply(seq_along(simIdxs), function(s)
  .MeanAcc(.PoolRuns(accuracy[["v1", s]], accuracy[["v2", s]])))

# Compare the four accuracy measures across replicates
accDf <- data.frame(
  replicate = rep(simIdxs, 4),
  measure   = rep(c("n", "t", "nt", "v"), each = length(simIdxs)),
  accuracy  = c(accN, accT, accNT, accV)
)

# Paired Wilcoxon signed-rank test: NT joint model vs. V joint model
# accNT and accV are paired by replicate; values are bounded [0, 1]
wilcox.test(accNT, accV, paired = TRUE, conf.int = TRUE)

table(nt_better <- accNT > accV)


plot(NA, xlim = c(0.8, 2.2), ylim = range(accNT, accV),
     xaxt = "n", xlab = NA, ylab = "Mean imputation accuracy",
     frame.plot = FALSE)
axis(1, at = 1:2, labels = modelLabel[c("by_kv", "by_nt_kv")])

segments(x0 = 1, x1 = 2, y0 = accV, y1 = accNT,
         col = ifelse(nt_better, "#2166ac", "#d73027"),
         lwd = ifelse(nt_better, 1.4, 0.7))
points(rep(1, length(accV)), accV, pch = 16, cex = 0.7, col = modelCol["by_ki"])
points(rep(2, length(accNT)),  accNT,  pch = 16, cex = 0.7, col = modelCol["by_nt_ki"])





### Accuracy of tree reconstruction ###
# [There's no significant difference] #
if (!exists("ntiErr")) {
  ntiErr <- vector("list", length(simIdxs))
}
if (!exists("spiErr")) {
  spiErr <- vector("list", length(simIdxs))
}

for (seed in setdiff(simIdxs, NULL)) {
  simID <- sprintf("sim%03d", seed)
  message(simID)
  trueTree <- read.tree(MkPath(simID, "tree.nwk"))
  i <- match(seed, simIdxs)
  
  if (is.null(ntiErr[[i]])) {
    ntTrees1 <- SimTrees(simID, "imp_sp_nt_kv", 1)
    ntTrees2 <- SimTrees(simID, "imp_sp_nt_kv", 2)
    if(any(NTip(ntTrees1) == 1)) {
      warning("Corrupted tree file with ", simID, " at line ",
              which(NTip(ntTrees1) == 1),
              immediate. = TRUE)
    } else {
      ntTrees <- c(BurnOff(ntTrees1, 0.25), BurnOff(ntTrees2, 0.25))
      if (!is.null(ntTrees1)) {
        tryCatch(
          ntiErr[[i]] <- TreeDist::ClusteringInfoDist(trueTree, ntTrees,
                                                     normalize = TRUE),
          error = function(e) message(simID, ": ", e$message)
        )
      }
    }
  }
  
  if (is.null(spiErr[[i]])) {
    spTrees1 <- SimTrees(simID, "imp_sp_kv", 1)
    spTrees2 <- SimTrees(simID, "imp_sp_kv", 2)
    spTrees <- c(BurnOff(spTrees1, 0.25), BurnOff(spTrees2, 0.25))
    if (!is.null(spTrees1)) {
      spiErr[[i]] <- TreeDist::ClusteringInfoDist(trueTree, spTrees, normalize = TRUE)
    }
  }
}; cli::cli_progress_done()

stillNull <- vapply(ntiErr, is.null, logical(1)) | vapply(spiErr, is.null, logical(1))
sum(stillNull)

ntiMeds <- vapply(ntiErr[!stillNull], median, numeric(1))
spiMeds <- vapply(spiErr[!stillNull], median, numeric(1))
wilcox.test(spiMeds, ntiMeds, paired = TRUE, conf.int = TRUE)

table(nt_better)
ntiBetter <- ntiMeds < spiMeds
table(ntiBetter)
table(nt_better, ntiBetter)
pCol <- ifelse(ntiMeds - spiMeds > 0, 3, ifelse(ntiMeds - spiMeds == 0, 1, 2))
# plot(ntiMeds - spiMeds ~ spiMeds, frame.plot = FALSE,
#      xlab = "Tree divergence under sp_kv",
#      ylab = "NT improvement on median tree distance",
#      col = pCol, pch = 16
# )
# abline(h = 0)

# Paired slope plot: NT vs V tree distance per replicate
# Blue = NT closer to true tree; red = V closer
nt_better_tree <- ntiMeds < spiMeds
plot(NA, xlim = c(0.8, 2.2), ylim = range(ntiMeds, spiMeds),
     xaxt = "n", xlab = NA, ylab = "Median distance from true tree",
     frame.plot = FALSE)
axis(1, at = 1:2, labels = modelLabel[c("by_kv", "by_nt_kv")])
segments(x0 = 1, x1 = 2, y0 = spiMeds, y1 = ntiMeds,
         col = ifelse(nt_better_tree, "#2166ac66", "#d7302766"),
         lwd = ifelse(nt_better_tree, 0.7, 1.4))
points(rep(1, length(spiMeds)), spiMeds, pch = 16, cex = 0.5)
points(rep(2, length(ntiMeds)), ntiMeds, pch = 16, cex = 0.5)
points(c(0.99, 2.01), c(median(spiMeds), median(ntiMeds)), pch = 4)
segments(1, median(spiMeds), 2, median(ntiMeds), lwd = 2, lty = "dashed")
