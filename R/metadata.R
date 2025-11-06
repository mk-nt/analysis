.meta <- new.env(parent = emptyenv())

#' Project metadata
#' @return `Metadata()` returns a list of dataset metadata.
#' @export
Metadata <- function() as.list(.meta)


#' Load metadata
#' Called on .onLoad
#' @keywords internal
.LoadMetadata <- function() {
  .meta$csv <- file.path(.config$matrixDir, "metadata.csv")
  
  .meta$nChar <- vapply(AllProjects(), function(pID) {
    c(
      trans = .NChar(MatrixFile(pID, "trans.nex")),
      neo = .NChar(MatrixFile(pID, "neo.nex"))
    )},
    integer(2))
  
  
  .meta$taxon <- {
    meta <- read.csv(.meta$csv,
                     colClasses = c("character", "integer", "integer",
                                    "character", "factor", "character"))
    setNames(
      factor(meta[["taxon"]],
             levels = c("Euarthropoda", "Invertebrates", "Vertebrates",
                        "Non-animals"))[
                          match(AllProjects(), meta[["project"]])],
      AllProjects())
  }
  
  .meta$rank <- {
    meta <- read.csv(.meta$csv,
                     colClasses = c("character", "integer", "integer",
                                    "character", "factor", "character"))
    setNames(
      factor(meta[["rank"]], ordered = TRUE,
             levels = c("Species", "Subgenus", "Genus", "Subtribe", "Tribe",
                        "Subfamily", "Family", "Superfamily",
                        "Infraclass", "Subclass", "Class", "Superclass",
                        "Infraorder", "Suborder", "Order",
                        "Subphylum", "Phylum", "Superphylum",
                        "Kingdom"))[
                          match(AllProjects(), meta[["project"]])],
      AllProjects())
  }
  
  .meta$nNeoCoded <- vapply(
    gsub(".nex", ".neo.nex", fixed = TRUE,
         vapply(AllProjects(), MatrixFile, character(1))),
    .NCoded, integer(1))
  .meta$nCoded <- vapply(vapply(AllProjects(), MatrixFile, character(1)),
                         .NCoded, integer(1))
  
  .meta$nTaxa <- vapply(vapply(AllProjects(), MatrixFile, character(1)),
                        .NTaxa, integer(1))
}

#' Human-readable label for internal metadata variable
#' @export
Decrypt <- function(x) list(
  charPerTax = "Characters per taxon",
  nTaxa = "Taxa in dataset",
  nChar = "Characters in dataset",
  inferredA.50. = expression("Inferred" ~ italic("a")[0]),
  inferredLogA0.50. = expression("Inferred log " ~ italic("a")[0]),
  inferredN.50. = expression("Inferred" ~ italic("n")),
  inferredLogN.50. = expression("Inferred log " ~ italic("n")),
  
  tValue.50. = expression("Inferred" ~ italic("t")),
  logNTRatio = "log(neom. : transf. characters)"
)[[x]]
