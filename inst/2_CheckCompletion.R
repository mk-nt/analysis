# This script provides a template for checking on the progress of remote
# analyses, and identifying analyses that require continuation.

library("neotrans")

statModels <- c("by_ki", "by_n_ki", "by_nn_ki", "by_t_ki", "by_nt_ki")
rmModels <- c("rm_by_n_ki", "rm_by_t_ki", "rm_by_nt_ki")
nsModels <- c("ns_ki", "ns_n_ki", "ns_t_ki", "ns_nt_ki")
hetModels <- c("hg_ki", "hg2_ki", "hg_b_ki", "hg_m_ki", "hg_bm_ki")
homModels <- c(statModels, nsModels)
kiModels <- c(statModels, rmModels, nsModels, hetModels)
nonRm <- setdiff(kiModels, rmModels)

marginals <- GetMarginals(AllProjects(), kiModels)
conv <- GetConvergence(setdiff(AllProjects(), .config$syab[-3]), kiModels)

# Validate that all analyses meet the ESS thresholds specified in the manuscript.
twoIsh <- c(8, 16, 32, 64, 99, 128, 168, 256, 512, 768, 1024, 1536, 2048, 4096,
            8192, 16384)
HistCol <- function(breaks = twoIsh) 2 + (breaks >= .config$essThreshold)
hist(conv["ess", , ], breaks = c(twoIsh, 5096), log = "x", col = HistCol())
LowEntry <- function(entry, thresh = .config$essThreshold) {
  slice <- conv[entry, , ]
  slice[is.na(slice)] <- Inf
  belowThresh <- slice < thresh
  belowThresh[apply(belowThresh, 1, any), apply(belowThresh, 2, any), drop = FALSE]
}
lowEss <- LowEntry("ess")
sum(lowEss)
for (pID in rownames(lowEss)) for (scriptID in colnames(lowEss)) {
  if (lowEss[pID, scriptID]) {
    message(pID, "_", scriptID, ": ESS ", signif(conv["ess", pID, scriptID]))
    MakeSlurm(pID, scriptID, ml = FALSE)
  }
}

hist(conv["frESS", , ], breaks = twoIsh, log = "x", col = HistCol())
lowFr <- LowEntry("frESS")
sum(lowFr)
for (pID in rownames(lowFr)) for (scriptID in colnames(lowFr)) {
  if (lowFr[pID, scriptID]) {
    message(pID, "_", scriptID, ": Frech ", signif(conv["frESS", pID, scriptID]))
    MakeSlurm(pID, scriptID, ml = FALSE)
  }
}

# There will probably only be a few of these not already covered in lowFr
mpsESS <- conv["mpsESS", , ]
mpsBreaks <- twoIsh[twoIsh < 2 * max(mpsESS, na.rm = TRUE) &
                      twoIsh > min(mpsESS, na.rm = TRUE) / 2]
hist(mpsESS, log = "x", col = HistCol(mpsBreaks),
     breaks = mpsBreaks)
lowMps <- LowEntry("mpsESS")
sum(lowMps)
for (pID in rownames(lowMps)) for (scriptID in colnames(lowMps)) {
  if (lowMps[pID, scriptID]) {
    message(pID, "_", scriptID, ": MPS ", signif(conv["mpsESS", pID, scriptID]))
    MakeSlurm(pID, scriptID, ml = FALSE)
  }
}
