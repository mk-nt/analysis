library("TreeTools")
library("TreeSearch")

nexusFiles <- list.files(pattern = "\\.nex$", full.names = FALSE)
chars <- do.call(cbind, lapply(nexusFiles, ReadCharacters)) |> MatrixToPhyDat()

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  args <- "Inf"
}
k <- as.numeric(args[[1]])
outFile <- gsub("^([^\\.]+).*", sprintf("\\1_k%s.trees", k), nexusFiles[[1]],
                perl = TRUE)
cat("Starting tree search with outFile = ", outFile)

startTree <- if (file.exists(outFile)) {
  message("Starting from previous results")
  read.tree(outFile, keep.multi = TRUE)[[1]]
} else {
  treeFiles <- list.files(pattern = "run_1\\.trees$", full.names = TRUE)
  if (length(treeFiles)) {
    message("Starting from last tree in ", treeFiles[[1]])
    trees <- read.tree(treeFiles[[1]], keep.multi = TRUE)
    trees[[length(trees)]]
  } else {
    message("Seeding with addition tree")
    AdditionTree(chars, concavity = k)
  }
}

if (startTree$Nnode < NTip(startTree) - 1) {
  startTree <- RootTree(startTree, 1)
}

lastScore <- TreeLength(startTree, chars, concavity = k)
best <- startTree

# Extra steps for long runners:
best <- MaximizeParsimony(chars, c(best)[[1]], concavity = k,
                          maxHits = 10, tbrIter = 1, ratchIter = 0)
write.tree(best, file = outFile)
bestScore <- TreeLength(best[[1]], chars, concavity = k)

best <- MaximizeParsimony(chars, c(best)[[1]], concavity = k,
                          maxHits = 20, tbrIter = 2, ratchIter = 0)
write.tree(best, file = outFile)
bestScore <- TreeLength(best[[1]], chars, concavity = k)

best <- MaximizeParsimony(chars, c(best)[[1]], concavity = k,
                          maxHits = 30, tbrIter = 3, ratchIter = 0)
write.tree(best, file = outFile)
bestScore <- TreeLength(best[[1]], chars, concavity = k)

best <- MaximizeParsimony(chars, c(best)[[1]], concavity = k,
                          maxHits = 40, tbrIter = 4, ratchIter = 0)
write.tree(best, file = outFile)
bestScore <- TreeLength(best[[1]], chars, concavity = k)

# Quick start
best <- MaximizeParsimony(chars, c(best)[[1]], concavity = k,
                          maxHits = 50, tbrIter = 5, ratchIter = 0)
write.tree(best, file = outFile)
bestScore <- TreeLength(best[[1]], chars, concavity = k)


eps <- sqrt(.Machine[["double.eps"]])
while (bestScore + eps < lastScore) {
  lastScore <- bestScore
  best <- MaximizeParsimony(chars, c(best)[[1]], concavity = k,
                            maxHits = 50, tbrIter = 5, ratchIter = 0)
  write.tree(best, file = outFile)
  bestScore <- TreeLength(c(best)[[1]], chars, concavity = k)
}

# Escape local optima
lastScore <- Inf
while (bestScore + eps < lastScore) {
  lastScore <- bestScore
  best <- MaximizeParsimony(chars, c(best)[[1]], concavity = k,
                            maxHits = 512, tbrIter = 2,
                            ratchIter = 30,
                            finalIter = 4
                            )
  write.tree(best, file = outFile)
  bestScore <- TreeLength(c(best)[[1]], chars, concavity = k)
}

cat("Completed tree search under k = ", k, " with score: ", bestScore)

system2("git", paste("add", outFile))
system2("git", sprintf("commit -m \"TreeSearch results, k = %s\"", k))
system2("git", "push")
