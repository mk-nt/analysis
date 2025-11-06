.config <- new.env(parent = emptyenv())

#' Configuration settings
#' 
#' View and edit package configuration variables
#' 
#' @examples
#' ntConfig <- Config()
#' 
#' # View maximum cores that will be requested via SLURM
#' .config$maxCores
#'   
#' \dontrun{
#'   # Override default matrix directory
#'   ntConfig$matrixDir <- "path/to/my/matrices" 
#' }
#' 
#' 
#' @export
Config <- function() .config

.onLoad <- function(libname, pkgname) {
  
  # Relative path to directory containing matrices and metadata
  # The github repo "zlvb35/matrixes" must be checked out here
  .config$matrixDir <- system.file("matrices", package = "neotrans")
  
  # Regular expression for excel translation files
  .config$metaPattern <- "^(Project|syab)%s[^~]"
  
  # SSH account
  # To set up credentials:
  # 1. Run file.edit(".Renviron") (or usethis::edit_r_environ())
  # 2. Add a line: sshLogin=mylogin@myserver.com
  # 3. Save the file and restart R
  
  # To log in using an SSH key, on Windows:
  # 1. In the cmd prompt, run
  #     ssh-keygen -t ed25519 -f C:\Users\your_username\.ssh\id_ed25519
  # 2. Copy the key to the remote server by typing
  #     type C:\Users\your_username\.ssh\id_ed25519.pub | ssh user@remote-server "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"
  # 3. add sshKey=C:/Users/your_username/.ssh/id_ed25519 to .Renviron
  # 4. If you created a passphrase for the key, also add
  #        sshPass=<key passphrase> to .Renviron
  
  .config$tmpMin <- 32 # M, minimum temporary memory to request
  
  # Account that hosts the project analysis repositories
  Sys.setenv("ntGithubAccount" = "neo-trans")
  # Repositories will be added to this account using the `gh` command line tool.
  # You may need to run `gh auth login` at the console to authenticate.
  # Or try `system2("gh", "auth login")` from within R.
  
  
  # GitHub CLI must be installed locally,
  # and a token with the repo scope must be stored in the MTNK_WRITE
  # environment variable. Run usethis::edit_r_environ() and add the line
  # MTNK_WRITE=<your token goes here>
  gitName <- "R script"
  gitEmail <- "mk-nt@neo.trans"
  
  system2("git", sprintf("config --global user.email \"%s\"", gitEmail))
  system2("git", sprintf("config --global user.name \"%s\"", gitName))
  
  
  # Remote server slurm settings
  sPerDay <- 60 * 60 * 24
  .config$maxSharedMem <- 241000 # Mb
  hg2Adjust <- 1.2
  .config$bigSharedMem <- .config$maxSharedMem / hg2Adjust
  .config$maxMem <- 1024 ^ 2 # Mb
  .config$maxCores <- 128L
  .config$maxTime <- 3 * sPerDay
  .config$longTime <- 7 * sPerDay
  
  .config$mc3ContinueTime <- 60 * 60 * 71
  .config$syab <- paste0("0720", c(0, 2:6))
  # Maximum value of the Gelman-Rubin statistic above which runs will not
  # be deemed to have converged
  .config$psrfThreshold <- 1.02
  # Minimum estimated sample size to employ
  .config$essThreshold <- 256
  # Minimum difference between log(Bayes Factors) to consider noteworthy
  .config$eps <- log(10)
  
  .LoadMetadata()
}

#' @importFrom ssh ssh_disconnect
.onUnload <- function(libpath) {
  if (!is.null(.ssh$session)) {
    try({
      ssh_disconnect(.ssh$session)
    }, silent = TRUE)
  }
  rm(list = "session", envir = .config)
}
