AllProjects <- function(indices = NULL, sort = TRUE) {
  x <- gsub("^\\D*(\\d+)\\D*?$", "\\1", perl = TRUE,
       list.files(dirname(MatrixFile("*")),
                  basename(MatrixFile(".*", "neo.nex$")))
  )
  x <- x[order(as.numeric(x))]
  if (is.null(indices)) {
    x
  } else {
    x[indices]
  }
  
}

AnalysisDir <- function(pID, scriptID = NULL) {
  dirname(ScriptFile(pID, scriptID))
}

ResultsDir <- function(pID, scriptID = NULL) {
  d <- sprintf("results/project%s", pID)
  if (!dir.exists(d)) {
    dir.create(d)
  }
  d
}

ConvergenceFile <- function(pID, scriptID) {
  sprintf("%s/%s-conv.txt", ResultsDir(pID, scriptID),
          ScriptBase(pID, scriptID))
}

if (!dir.exists("rbPDF")) {
  dir.create("rbPDF")
}

DistanceFile <- function(pID, model1, model2, final = TRUE) {
  sprintf("rbPDF/project%s-%s-%s.dist%s.rds", pID, model1, model2, 
          if (final) "" else ".tmp")
}

MarginalFile <- function(pID, scriptID) {
  ScriptFile(pID, scriptID, "marginal")
}

MatrixFile <- function(pID = "", ext = "nex") {
  sprintf("projects/project%s.%s", pID, ext)
}

OldStoneFile <- function(pID, scriptID) {
  paste0(StoneFile(pID, scriptID), ".approx")
}

ParameterFile <- function(pID, scriptID) {
  sprintf("%s/%s.rds", ResultsDir(pID, scriptID), ScriptBase(pID, scriptID))
}

PDFFile <- function(pID, model1, model2) {
  sprintf("rbPDF/project%s-%s-%s.pdf", pID, model1, model2)
}

PFiles <- function(pID, scriptID) {
  list.files(
    path = AnalysisDir(pID, scriptID),
    pattern = sprintf("%s\\.p_run_.*\\.log$", scriptID),
    full.names = TRUE, ignore.case = TRUE
  )
}

ScriptBase <- function(pID, scriptID) {
  sprintf("%s_%s", pID, scriptID)
}

ScriptFile <- function(pID, scriptID, analysis = "mcmcmc") {
  sprintf("../mknt-rb/%s/%s.Rev", ScriptBase(pID, scriptID), analysis)
}


SlurmTemplate <- function(ml) {
  if (ml) "slurm/marginal.sh" else "slurm/mcmcmc.sh"
}

SlurmFile <- function(pID, scriptID, ml = FALSE) {
  sprintf("slurm/%s_%s%s.sh", pID, scriptID, if (ml) "-ml" else "")
}

StoneFile <- function(pID, scriptID) {
  sprintf("%s/%s-stone.txt", ResultsDir(pID, scriptID),
          ScriptBase(pID, scriptID))
}

StoneOrigin <- function(pID, scriptID) {
  sprintf("%s/%s.out.pp", AnalysisDir(pID, scriptID), ScriptBase(pID, scriptID))
}

TreeFiles <- function(pID, scriptID) {
  list.files(
    path = AnalysisDir(pID, scriptID),
    pattern = sprintf("%s\\_run_.*\\.trees$", scriptID),
    full.names = TRUE, ignore.case = TRUE
  )
}

TreeSampleFile <- function(pID, scriptID) {
  sprintf("%s/%s.trees", ResultsDir(pID, scriptID), ScriptBase(pID, scriptID))
}
