# Relative path to directory containing matrices and metadata
matrixDir <- "matrices"

psrfThreshold <- 1.05
essThreshold <- 256

# Account that hosts the project analysis repositories
githubAccount <- "mk-nt"
# Repositories will be added to this account using the `gh` command line tool.
# Run `gh auth login` at the console to authenticate.


# GitHub CLI must be installed locally,
# and a token with the repo scope must be stored in the MTNK_WRITE
# environment variable. Run usethis::edit_r_environ() and add the line
# MTNK_WRITE=<your token goes here>

gitName <- "R script"
gitEmail <- "mk-nt@neo.trans"

system2("git", sprintf("config --global user.email \"%s\"", gitEmail))
system2("git", sprintf("config --global user.name \"%s\"", gitName))


# The first time the script is run, it is also necessary to install packages:
if (!requireNamespace("TreeDist", quietly = TRUE)) {
  install.packages("TreeDist")
}

# For checking whether GitHub files exist
if (!requireNamespace("curl", quietly = TRUE)) {
  install.packages("curl")
}

# For progress reporting
if (!requireNamespace("cli", quietly = TRUE)) {
  install.packages("cli")
}

# 
if (!requireNamespace("mcmcse", quietly = TRUE)) {
  install.packages("mcmcse")
}

# For PDF table plotting
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

# For tree ess
if (!requireNamespace("treess", quietly = TRUE)) {
  devtools::install_github("afmagee/treess")
}

figPalette <- unlist(durham <- list(
    "ink" = "#002A41",
    "purple" ="#68246D",
    "yellow" ="#FFD53A",
    "cyan" ="#00AEEF",
    "red" = "#BE1E2D",
    "gold" ="#AFA961",
    "heather" = "#CBA8B1",
    "stone" = "#DACDA2",
    "sky" = "#A5C8D0",
    "cedar" = "#B6AAA7",
    "concrete" ="#B3BDB1",
    "black" = "#333132",
    "white" = "#ffffff"
  ))

par("col" = durham[["ink"]], "fg" = durham[["ink"]])
