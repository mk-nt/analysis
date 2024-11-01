source("R/_setup.R")
source("R/PrepareMatrix.R")
source("R/RevBayes.R")

stopifnot(dir.exists(matrixDir))
# Load excel spreadsheets, ignoring temporary files
metaFiles <- grep("^Project\\d[^~]", list.files(matrixDir, "*.xlsx"),
                  value = TRUE, perl = TRUE)

# Prepare matrices and RevBayes scripts for execution
for (metaFile in metaFiles) {
  # Extract the MorphoBank project ID from the name of the Excel spreadsheet
  pID <- gsub("^.*?(\\d+).*?$", "\\1", metaFile, perl = TRUE)
  
  # Report progress
  cli::cli_progress_message(paste0("Project ", pID,
                            " (metafile ", match(metaFile, metaFiles), ")"))
  
  # PrepareMatrix() uses the Excel spreadsheet and MorphoBank Nexus file to 
  # create separate Nexus files for neomorphic and transformational partitions.
  if (PrepareMatrix(pID, metaFile)) {
    # RevBayes() copies these Nexus files and RevBayes scripts to a GitHub repo,
    # and creates slurm files to execute each job on the HPC cluster.
    RevBayes(pID, "by_ki")
    RevBayes(pID, "by_n_ki")
    RevBayes(pID, "by_t_ki")
    RevBayes(pID, "by_nt_ki")
  }
}

cli::cli_progress_done()
