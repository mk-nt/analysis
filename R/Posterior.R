#' Compare prior and posterior distributions
#' @inheritParams MakeSlurm
#' @param parameter Character string corresponding to column name in 
#' `ParameterFile(pID, scriptID), specifying which parameter to examine.
#' @param priorMean,priorSD Mean and SD of a normal prior on `parameter`.
#' @export
PriorVsPost <- function(pID, scriptID, parameter, priorMean = 0, priorSD = 2) {
  p <- if (file.exists(ParameterFile(pID, scriptID))) {
    readRDS(ParameterFile(pID, scriptID))
  } else {
    stop("No parameter file found.")
  }
  
  quants <- quantile(p[[parameter]], c(0.01, 0.99))
  xrange <- c(min(priorMean - 3 * priorSD, quants[[1]]),
              max(priorMean + 3 * priorSD, quants[[2]]))
  grid <- seq(xrange[1], xrange[2], length.out = 2000)
  priorY <- dnorm(grid, priorMean, priorSD)
  
  postD <- density(p[[parameter]], from = xrange[[1]], to = xrange[[2]])

  ylim <- range(c(postD$y, priorY))
  
  plot(postD, frame.plot = FALSE, main = "",
       xlim = xrange, ylim = ylim,
       xlab = expression(rate[neo]), ylab = "Density",
       col = "#0072B2", lwd = 2)
  
  lines(grid, priorY, col = "#D55E00", lwd = 2, lty = 2)
  abline(v = 0, lty = "dotted")
  
  legend("topright", bty = "n", lwd = 2, lty = c(1, 2),
         col = c("#0072B2", "#D55E00"),
         legend = c("Posterior", "Prior"))
}

#' Increase in information density in prior vs posterior.
#' @inheritParams PriorVsPost
#' @param projects Which projects to use
#' @param metaEntry Metadata item to regress with `parameter`
#' @returns `InformationGain()` returns the spearman's correlation of
#' `metaEntry` with `parameter`
#' @importFrom stats cor.test
#' @returns `InformationGain()` returns the results of a Spearman correlation
#' conducted with `cor.test(method = "spearman")` between the prior/posterior
#' SD ratio and the requested `metaEntry`.
#' @export
InformationGain <- function(projects = KiProjects(), scripts, parameter,
                            metaEntry, priorMean = 0, priorSD = 2) {
  
  post <- do.call(rbind, lapply(scripts, function(scriptID) {
    scriptPost <- vapply(projects, function(pID) {
      f <- ParameterFile(pID, scriptID)
      if (!file.exists(f)) return(c(NA_real_, NA_real_))
      p <- readRDS(f)
      if (!parameter %in% colnames(p)) {
        warning(parameter, " not found in ", paste(colnames(p), collapse = ", "),
                immediate. = TRUE)
      }
      vals <- p[[parameter]]
      c(mean(vals), sd(vals))
    }, numeric(2))
    
    data.frame(
      pID = colnames(scriptPost),
      postMean = scriptPost[1, ],
      postSD = scriptPost[2, ],
      scriptID = scriptID,
      row.names = NULL
    )
  }))
  
  post$Zshift <- (post$postMean - priorMean) / priorSD
  post$Rratio <- priorSD / post$postSD
  post$nChar <- .meta$nChar[match(projects, colnames(.meta$nChar))]
  post$nTaxa <- .meta$nTaxa[match(projects, names(.meta$nTaxa))]
  post$nNeo <- .meta$nChar["neo", match(projects, colnames(.meta$nChar))]
  post$cells <- post$nChar * post$nTaxa
  post$neoCells <- post$nNeo * post$nTaxa
  post$nCoded <- .meta$nCoded[match(projects, names(.meta$nCoded))]
  post$nNeoCoded <- .meta$nNeoCoded[match(projects, names(.meta$nNeoCoded))]
  
  .Decode <- function(x) switch (x,
    neoCells = "Neomorphic cells (characters × taxa)",
    cells = "Total cells (characters × taxa)",
    "nTaxa" = "Taxa",
    "nChar" = "Characters",
    "nNeo" = "Neomorphic characters",
    "nNeoCoded" = "Non-ambiguous neomorphic cells",
    "nCoded" = "Non-ambiguous cells",
    "root_freqs.1." = "P(absent) at root",
    "rate_loss" = "n",
    "rate_neo" = "t",
    Decrypt(x)
  )
  
  
  for (i in which(is.na(post$Rratio))) {
    UpdateRecords(post[i, "pID"], post[i, "scriptID"], searchRemote = TRUE)
  }
  post <- post[!is.na(post$Rratio), ]
  
  iffy <- post$Rratio < sqrt(.Machine$double.eps)
  plot(post[[metaEntry]], post$Rratio,
       log = "xy", xpd = NA,
       cex = ifelse(iffy, 0.001, 1),
       frame.plot = FALSE, pch = 16,
       col = ModelCol(post$scriptID),
       xlab = .Decode(metaEntry),
       ylab = paste("Prior / posterior SD: ", .Decode(parameter)))
  if (any(iffy)) {
    text(post[[metaEntry]][iffy], post$Rratio[iffy], post$pID[iffy],
         cex = 0.8, xpd = NA)
  }
  
  # points(size, post$Zshift, col = 2, pch = 2)
  # mtext(expression("|Posterior shift| ("*sigma[prior]*" units)"), 2, col = 2, line = 1.7)
  
  for (i in which(iffy)) {
    RewindRepo(post[i, "pID"], post[i, "scriptID"])
  }
  
  abline(h = 1, lty = 2, col = "grey")   # uninformative line
  
  suppressWarnings(cor.test(post[[metaEntry]], post$Rratio, method = "spearman"))
}


#' Restore MCMC files in a git repo to a previous state
#' 
#' Necessary because some files became corrupted: my suspicion is that
#' continuing a run on RB<1.4 with RB>1.4 causes an issue when no `.heat.ckp`
#' files are found. 
#' 
#' Updates the git repo, deletes cache files, and runs `UpdateRecords()`.
#' 
#' @inheritParams MakeSlurm
#' @param lastGoodID Git identifier denoting (first characters of) last good
#' commit ID.  A blank entry will return.
#' "n" will delete relevant files and enqueue a fresh analysis.
#' 
#' @importFrom ssh ssh_exec_wait
#' @export
RewindRepo <- function(pID, scriptID, lastGoodID) {
  
  if (missing(lastGoodID)) {
    browseURL(
      sprintf("https://github.com/neo-trans/%s_%s/blame/main/%s_%s.p_run_1.log",
              pID, scriptID, pID, scriptID))
    lastGoodID <- readline(paste(pID, scriptID, "-", "Last good ID: "))
  }
  
  if (nchar(lastGoodID) == 0) return()
  
  if (substr(tolower(lastGoodID), 1, 1) == "n") {
    # Hard-core aggression: Discard results and start from scratch
    result <- ssh_exec_wait(
      SshSession(),
      paste(sep = " && ",
            sprintf("(cd %s/%s_%s", RemoteDir(), pID, scriptID),
            "git stash",
            "git pull -r",
            "git stash pop",
            sprintf("git rm %s_%s*.ckp %s_%s*run*.log %s_%s*.tar.gz",
                    pID, scriptID, pID, scriptID, pID, scriptID),
            "git commit -m \"Remove corrupted run results\"",
            "git push)")
    )
    EnqueueML(pID, scriptID)
  } else {
    scriptBase <- paste0(c(pID, scriptID), collapse = "_")
    filesToRestore <- paste0(c("", ".ckp", "run*.log", ".trees", ".tar.gz"),
                             collapse = paste0(" ", scriptBase, "*"))
    
    result <- ssh_exec_wait(
      SshSession(),
      paste(
        sprintf("(cd %s/%s_%s", RemoteDir(), pID, scriptID),
        
        # ---- Stage 1 ----
        "&& echo '--- Stage 1: Restoring checkpoints and logs ---'",
        "&& git stash || { echo 'Stash failed'; exit 1; }",
        sprintf("&& git pull -r || { echo 'Rebase failed'; exit 1; }"),
        sprintf("&& git restore --source=%s --staged --worktree %s_%s*.ckp %s_%s*run*.log || { echo 'Restore of ckp/log failed'; exit 1; }",
                lastGoodID, pID, scriptID, pID, scriptID),
        "&& echo 'Stage 1 complete: .ckp and .log restored.'",
        
        # ---- Stage 2 ----
        "&& echo '--- Stage 2: Restoring tree files ---'",
        # First, try to restore .tar.gz; if fails, fall back to .trees
        sprintf("&& if git restore --source=%s --staged --worktree %s_%s*run*.tar.gz 2>/dev/null; then",
                lastGoodID, pID, scriptID),
        "     echo 'Restored .tar.gz files — extracting to overwrite .trees';",
        "     for tgz in *_run*.tar.gz; do",
        "         [ -f \"$tgz\" ] && tar -zxf \"$tgz\" || true;",
        "     done;",
        "   else",
        sprintf("     echo 'No .tar.gz found — restoring .trees instead';"),
        sprintf("     git restore --source=%s --staged --worktree %s_%s*run*.trees || true;",
                lastGoodID, pID, scriptID),
        "     for tr in *_run*.trees; do",
        "         base=\"${tr%.trees}\";",
        "         [ -f \"$tr\" ] && tar -zcf \"${base}.tar.gz\" \"$tr\";",
        "         git add \"${base}.tar.gz\";",
        "         git rm --cached \"$tr\" || true;",
        "     done;",
        "   fi",
        "&& echo 'Stage 2 complete: tree files handled.'",
        
        # ---- Stage 3 ----
        "&& echo '--- Stage 3: Commit and push ---'",
        sprintf("&& git commit -am \"Revert MCMC files to last known good state (%s)\" || echo 'Nothing to commit.'",
                substr(lastGoodID, 1, 7)),
        "&& git push || { echo 'Push failed'; exit 1; }",
        "&& echo 'Stage 3 complete: repository updated successfully.'",
        ")"
      )
    )
  }
  
  unlink(ConvergenceFile(pID, scriptID))
  unlink(MarginalFile(pID, scriptID))
  unlink(ParameterFile(pID, scriptID))
  unlink(ParsimonyFile(pID, scriptID))
  peFiles <- gsub("990099", ".*", ParsEvalFile(pID, scriptID, 990099, 990099))
  unlink(list.files(dirname(peFiles), basename(peFiles)))
  unlink(PPSamplerFile(pID, scriptID))
  unlink(TreeSampleFile(pID, scriptID))
  unlink(list.files(AnalysisDir(pID, scriptID), paste0(pID, "_.*.trees")))
  UpdateRecords(pID, scriptID, forgetCache = TRUE)
}
