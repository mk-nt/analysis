#' Retrieve counts of 0s and 1s from simulation output
#'
#' Read the remote or local simulation output created for a project/script and
#' return the counts of `0` and `1` for each simulated matrix contained in the
#' output. The function will attempt to fetch results if the expected local
#' output file is missing; when necessary it will inspect the remote repository
#' (via SSH), possibly compress and download large logs, or attempt a remote
#' git refresh. If the remote results are not yet available the function may
#' queue or trigger a re-run via `SimSlurm()` and return a structured `FALSE`
#' with a reason.
#'
#' @inheritParams MakeSlurm
#' @param prefix chr File-name prefix used when looking up the log; defaults to
#'   \code{"stat"}.
#'
#' @return
#' 'Simulated01s()' returns a numeric matrix whose rows are named \code{"0"}
#' and \code{"1"} and whose columns correspond to each simulated matrix found
#' in the output log. Each column gives the total count of zeroes and ones for
#' that simulated matrix. If the results are not yet available the function
#' returns \code{structure(FALSE, reason = <string>)} describing why.
#'
#' @details
#' The function performs several side effects when the local output file is
#' missing: it calls \code{FetchResults()} to fetch output, inspects the remote
#' repository via SSH, and may compress large
#' log files, download compressed logs with \code{scp_download()} and
#' uncompress them locally. If the remote git repo has diverged the function
#' will attempt a safe pull (stash/rebase/pop) and push updates.
#'
#'
#' @importFrom ssh ssh_exec_internal ssh_exec_wait scp_download
#' @importFrom stringi stri_count_fixed
#' @importFrom R.utils gunzip
#' @export
Simulated01s <- function(pID, scriptID, prefix = "stat") {
  simOut <- file.path(AnalysisDir(pID, scriptID),
                      sprintf("%s-%s_%s.out", prefix, pID, scriptID))
  if (!file.exists(simOut)) {
    FetchResults(pID, scriptID)
  }
  if (!file.exists(simOut)) {
    remoteDir <- file.path(RemoteDir(), ScriptBase(pID, scriptID))
    remoteLog <- file.path(remoteDir, basename(simOut))
    session <- SshSession()
    if (isFALSE(session %||% FALSE)) {
      return(structure(FALSE, reason = "Could not establish SSH session"))
    }
    remoteExists <- ssh_exec_internal(session,
                                      paste0("test -f ", remoteLog),
                                      error = FALSE)$status
    
    if (remoteExists == 0) {
      fileSize <- ssh_exec_internal(session,
                                    paste0("stat -c%s ", remoteLog),
                                    error = FALSE)
      if (fileSize$status == 0) {
        mb <- as.numeric(gsub("\n", "", fixed = TRUE,
                              rawToChar(fileSize$stdout))) / 1024 / 1024
        if (mb > 100) {
          ssh_exec_wait(session, paste0("gzip -kf ", remoteLog))
          compressedFile <- paste0(remoteLog, ".gz")
          expandedPath <- ssh::ssh_exec_internal(
            session,
            paste0("readlink -f ", compressedFile),
            error = FALSE)
          remoteZip <- trimws(rawToChar(expandedPath$stdout))
          scp_download(session, files = remoteZip,
                       to = AnalysisDir(pID, scriptID))
          git <- ssh_exec_internal(
            session,
            paste("(cd", remoteDir,
                  "&& git rm --cached", basename(remoteLog), "2>/dev/null;",
                  "git rm --cached",
                  gsub("stat-", "sim-", basename(remoteLog)), " 2>/dev/null;",
                  "git add", basename(remoteZip),
                  "&& git commit -m \"Add log file\"",
                  "&& git push)"),
            error = FALSE)
          rawToChar(git$stderr)
          gunzip(file.path(AnalysisDir(pID, scriptID),
                           basename(remoteZip)))
        }
      } else {
        warning(pID, " file size unavailable ", rawToChar(fileSize$stderr))
      }
    }
    
    status <- ssh_exec_internal(
      session,
      paste("git -C", file.path(RemoteDir(), ScriptBase(pID, scriptID)),
            "status"), error = FALSE)
    if (status[["status"]] != 0) {
      return(structure(FALSE, reason = paste("Failed to read Git repo:",
                                             rawToChar(status$stderr))))
    }
    looksGood <- FALSE
    if (grepl("have diverged", fixed = TRUE, rawToChar(status$stdout))) {
      squeue <- SlurmQueue()
      runningJobs <- squeue[squeue[["status"]] == "R", "JobName"]
      if (!any(grepl(ScriptBase(pID, scriptID), runningJobs, fixed = TRUE))) {
        on.exit(cli::cli_progress_done())
        cli::cli_progress_message("{pID}: Refreshing remote git repo")
        update <- ssh::ssh_exec_internal(
          session,
          paste(collapse = " && ", c(
            paste("cd", file.path(RemoteDir(), ScriptBase(pID, scriptID))),
            "git stash",
            "git pull --rebase",
            "git stash pop",
            paste0("git add ", basename(simOut)),
            "git commit -m \"Update simulated datasets\"",
            "git push")),
          error = FALSE)
        if (update[["status"]] == 0) {
          cli::cli_progress_message("{pID}: Fetching results")
          FetchResults(pID, scriptID)
          looksGood <- TRUE
        } else {
          return(structure(FALSE, reason = paste("Git update failed:",
                                                 rawToChar(update$stderr))))
        }
      }
    }
    return(structure(FALSE, reason = paste(pID, "_", scriptID, "running now")))
  }
  lines <- if (exists("remoteLines") && length(remoteLines) > 0) {
    remoteLines
  } else {
    lines <- readLines(simOut)
  }
  newMat <- grep("Simulated matrix #", lines, fixed = TRUE)
  matBounds <- rbind(newMat[-length(newMat)] + 1,
                     newMat[-1] - 1)
  apply(matBounds, 2, function(limits) {
    mat <- lines[limits[[1]]:limits[[2]]]
    taxLines <- which(!grepl("^[\\s01]+$", perl = TRUE, mat))
    
    # Optional check to verify that all simulated characters are informative
    #
    # tokenStart <- taxLines + 1
    # tokenEnd <- c(taxLines[-1] - 1, length(mat))
    # chars <- `colnames<-`(sapply(seq_along(taxLines), function(i) {
    #   (stringi::stri_paste(mat[tokenStart[i]:tokenEnd[i]], collapse = "") |>
    #     strsplit("\\s+"))[[1]]
    # })[-1, ], mat[taxLines]) |> `mode<-`("numeric")
    # stopifnot(all((colSums(apply(chars + 1, 1, tabulate) > 1)) == 2))
    
    c("0" = sum(stri_count_fixed(mat[-taxLines], "0")),
      "1" = sum(stri_count_fixed(mat[-taxLines], "1")))
  })
}
