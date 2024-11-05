source("R/Helpers.R")

extraMemory <- c(691, 3670, 3708, 3763, 4173, 4285)
bigMemory <- c(3707, 4284)
tooMuchMemory <- c(3707) # Out of memory even on bigmem
moreTime <- c(684, 691, 3708, 3670, 3763, 3804, 4173, 4285, 5228)

MakeSlurm <- function(pID, scriptID, ml = FALSE) {
  if (pID %in% tooMuchMemory) {
    RemoveSlurm(pID, scriptID, ml, "oom")
    return(structure(FALSE, "reason" = "Cannot allocate enough memory"))
  }
  
  slurmFile <- SlurmFile(pID, scriptID, ml)
  newLines <- gsub("%PID%", fixed = TRUE, pID,
                   gsub("%SCRIPTBASE%", fixed = TRUE,
                        ScriptBase(pID, scriptID), readLines(SlurmTemplate(ml))))
  timeLine <- 6
  partitionLine <- 12
  
  if (pID %in% extraMemory) {
    newLines[[5]] <- "#SBATCH --mem=246G"
    newLines[[7]] <- "#SBATCH --gres=tmp:350G"
  } else if (pID %in% bigMemory) {
    newLines[[5]] <- "#SBATCH --mem=1111G"
    newLines[[7]] <- "#SBATCH --gres=tmp:400G"
    newLines[[partitionLine]] <- "#SBATCH -p bigmem"
  }
  
  if (pID %in% moreTime) {
    if (pID %in% bigMemory) {
      warning(pID, "can't have lots of memory and lots of time.")
    } else {
      newLines[[timeLine]] <- "#SBATCH --time=6-23:35:0"
      newLines[[partitionLine]] <- "#SBATCH -p long"
    }
  }
  
  if (file.exists(slurmFile) &&
      isTRUE(all.equal(newLines, readLines(slurmFile)))) {
    structure(FALSE, reason = "Content not changed")
  } else {
    writeLines(newLines, slurmFile)
    system2("git", sprintf("add %s", slurmFile), stdout = NULL)
    system2("git", sprintf("commit -m \"%s_%s slurm files\"",
                           pID, scriptID), stdout = NULL)
    TRUE
  }
}

RemoveSlurm <- function(pID, scriptID, ml = FALSE, reason = "OK") {
  slurmFile <- SlurmFile(pID, scriptID, ml)
  if (file.exists(slurmFile)) {
    rmMsg <- system2("git", sprintf("rm %s", slurmFile), stdout = TRUE)
    msg <- switch(
      reason,
      "OK" = sprintf(
        "commit -m \"Remove %s_%s slurm files: run complete\"", pID,
        scriptID),
      "problem" = sprintf(
        "commit -m \"Remove %s_%s slurm files: run problematic\"", pID,
          scriptID),
      "oom" = sprintf(
        "commit -m \"Remove %s_%s slurm files: insufficient memory\"", pID,
        scriptID)
    )
    commit <- system2("git", msg, stdout = TRUE)
    .GitPush()
  }
}
