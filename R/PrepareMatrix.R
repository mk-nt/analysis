#' Prepare matrices for phylogenetic analysis
#'
#' `PrepareMatrix()` constructs a processed morphological matrix and its
#' derivatives for downstream phylogenetic analysis.  
#' It reads a raw character matrix and corresponding metadata spreadsheet
#' from the directory specifed in `.config$matrixDir` (default: `/matrices`),
#' applies standardization rules, re-codes characters according to the
#' annotation patterns, and outputs curated NEXUS files for use in
#' subsequent inference.
#'
#' Each character is categorized as *neomorphic* (presence/absence),
#' *transformational* (multi-state), or *ignored* according to the pattern
#' specified in the metadata file.  
#' Inapplicable tokens in neomorphic characters are converted to absences
#' (0) following \insertCite{Brazeau2011}{neotrans}.  
#' Characters may also be translated according to custom mapping patterns
#' defined in the spreadsheet.
#'
#' The function creates four main output files at the path specified by
#' `MatrixFile(pID)`:
#' \describe{
#'   \item{`*.nex`}{Full processed matrix.}
#'   \item{`*-neo.nex`}{Subset containing informative neomorphic characters.}
#'   \item{`*-trans.nex`}{Subset containing informative transformational characters.}
#'   \item{`*-neo-rand.nex`, `*-trans-rand.nex`}{Matched random subsets for
#'   comparative analysis.}
#' }
#'
#' @param pID Character giving the project identifier, e.g. `"1210"` for
#' [MorphoBank project 1210](https://www.morphobank.org/index.php/Projects/ProjectOverview/project_id/1210).  
#' Identifiers beginning with `"072"` refer to matrices from
#' \doi{10.1093/sysbio/syab072}, corresponding to files named
#' `syab<pID>.nex`.
#'
#' @param overwrite Logical; whether to overwrite existing processed files.
#' Defaults to `FALSE`.
#'
#' @details
#' The function expects metadata spreadsheets formatted according to the
#' \link[=matrix-processing]{matrix-processing vignette}, containing at least
#' the columns `"Character Pattern"` and `"Unseen States"`.  
#' Character patterns define how tokens are grouped and interpreted;
#' invalid or mismatched patterns will trigger diagnostic errors.  
#' Parsimony-uninformative characters (which define no non-trivial bipartitions)
#' are excluded from subset matrices.
#'
#' The function automatically sanitizes taxon names for NEXUS compatibility
#' and ensures that neomorphic characters are strictly binary (`0`, `1`).
#' Multi-state characters are re-encoded so that their state tokens form a
#' contiguous sequence starting from zero.
#'
#' @return
#' `PrepareMatrix()` invisibly returns `TRUE` if processing succeeds, or
#' `FALSE` if matrix preparation cannot proceed (e.g., due to missing
#' informative characters).  
#' The function is called primarily for its side-effect of producing
#' standardized NEXUS matrix files in the appropriate output directory.
#'
#' @references
#' \insertAllCited{}
#' 
#' @importFrom ape write.nexus.data
#' @importFrom Rdpack reprompt
#' @importFrom TreeTools ReadCharacters
#' @importFrom readxl read_xlsx
#' @export
PrepareMatrix <- function(pID, overwrite = FALSE) {
  projectID <- sprintf("pID%s", pID)
  nexusFile <- MatrixFile(pID, "nex")
  outDir <- dirname(nexusFile)
  if (!dir.exists(outDir)) {
    dir.create(outDir)
  }
  
  if (overwrite || !file.exists(nexusFile)) {
    
    metaFile <- grep(sprintf(.config$metaPattern, pID),
                     list.files(.config$matrixDir, "*.xlsx"),
                     value = TRUE, perl = TRUE)
    if (length(metaFile) != 1) {
      stop("Couldn't identify metafile for project ", pID)
    }
    meta <- read_xlsx(paste0(.config$matrixDir, "/", metaFile))
    patterns <- meta[["Character Pattern"]]
    if (any(is.na(patterns))) {
      stop("No pattern specified for character ",
           paste(which(is.na(patterns)), collapse = ", "))
    }
    unseen <- meta[["Unseen States"]]
    hasUnseen <- Inf
    noUnseen <- -Inf
    unseen[toupper(unseen) == "YES"] <- hasUnseen
    unseen[toupper(unseen) == "NO"] <- noUnseen
    unseen <- as.numeric(unseen)
    
    chars <- ReadCharacters(paste0(.config$matrixDir, "/", basename(nexusFile)))
    nTax <- dim(chars)[[1]]
    nChars <- dim(chars)[[2]]
    taxa <- rownames(chars)
    sanitizedTaxa <- gsub("&[gl]t;|[\"\\?\\(\\)]", "_", perl = TRUE,
                          gsub(".", "", taxa, fixed = TRUE))
    
    if (nChars != length(patterns)) {
      stop(nChars, " characters in matrix; ", length(patterns),
           " entries in spreadsheet")
    }
    
    neo <- substr(tolower(patterns), 1, 3) == "neo"
    trans <- substr(tolower(patterns), 1, 5) == "trans"
    ignored <- substr(tolower(patterns), 1, 6) == "ignore"
    special <- !neo & !trans
    translations <- lapply(patterns[special], function(pat) {
      if (substr(tolower(pat), 1, 6) == "ignore") {
        list(type = character(0), k = integer(0), n = 0L)
      } else {
        parts <- trimws(strsplit(pat, ";")[[1]])
        newTypes <- tolower(parts[[1]])
        matches <- gregexpr("n|t.", newTypes)
        newCode <- regmatches(newTypes, matches)[[1]]
        newType <- substr(newCode, 1, 1)
        newKText <- substr(newCode, 2, 2)
        newKText[newKText == "?"] <- hasUnseen
        newK <- as.numeric(newKText)
        nNew <- length(newType)
        
        translate <- strsplit(parts[-1], "=")
        if (!all(lengths(translate) == 2)) {
          stop("Invalid translation pattern: ", pat)
        }
        
        
        # Return:
        tryCatch(
          list(type = newType, k = newK, n = nNew,
               from = vapply(translate, "[[", character(1), 1),
               to = vapply(translate,
                           function(x) strsplit(trimws(x[[2]]), "")[[1]],
                           character(nNew))),
          error = function(e) {
            stop("Failed to parse translation pattern ", pat)
          }
        )
      }
    })
    nNew <- rep(1, nChars)
    nNew[special] <- vapply(translations, "[[", integer(1), "n")
    index <- cumsum(nNew)
    type <- character(sum(nNew))
    type[index[neo]] <- "n"
    type[index[trans]] <- "t"
    nowSpecial <- type == ""
    
    type[type == ""] <- unlist(lapply(translations, "[[", "type"))
    
    toRewrite <- special & !ignored
    if (any(toRewrite)) {
      nSpecial <- sum(toRewrite)
      whichSpecial <- which(toRewrite)
      newChars <- do.call(rbind, lapply(seq_len(nSpecial), function(i) {
        char <- chars[, whichSpecial[[i]]]
        translate <- translations[[i]]
        to <- translate[["to"]]
        if (is.null(dim(to))) {
          to <- t(to)
        }
        ret <- to[, match(char, translate[["from"]])]
        `[<-`(ret, is.na(ret), "?")
      }))
      ret <- matrix("", nTax, sum(nNew), dimnames = list(sanitizedTaxa, NULL))
      ret[, index[!ignored]] <- chars[, !ignored]
      ret[, nowSpecial] <- t(newChars)
    } else {
      ret <- chars[, !ignored]
    }
    
    # Replace inapplicable with absence in neomorphic characters
    ret[, type == "n"][ret[, type == "n"] == "-"] <- 0
    
    # Check neomorphic characters only contain 0, 1
    illegalNeomorphic <- grepl("[23456789]", ret[, type == "n"])
    if (any(illegalNeomorphic)) {
      # Decipher where the illegal characters are
      illegalMat <- matrix(illegalNeomorphic, nrow(ret))
      charID <- colSums(illegalMat) > 0
      taxa <- rowSums(illegalMat) > 0
      stop("Illegal neomorphic characters in pID ", pID, ", char ",
           paste(which(type == "n")[charID], collapse = ", "),
           ", from input char ", 
           paste(match(which(type == "n")[charID], index), collapse = ", "),
           ";\n   taxa ",
           paste(rownames(ret)[taxa], collapse = ", ")
      )
      ret[taxa, type == "n"][, charID]
    }
    
    # Compress tokens so only k tokens are used for a k-state character
    # e.g. a character that only uses the tokens 1 and 3 will be recoded to
    # use tokens 0 and 1.
    ret <- apply(ret, 2, function(x) {
      tokens <- 0:9
      used <- vapply(tokens, grepl, TRUE, paste0(x, collapse = ""))
      for (i in seq_len(sum(used))) {
        x <- gsub(tokens[used][[i]], tokens[[i]], x, fixed = TRUE)
      }
      x
    })
    
    informative <- apply(ret, 2, function(x) {
      tab <- table(x[nchar(x) == 1 & x %in% 0:9])
      length(tab[tab > 1]) > 1
    })
    kObs <- apply(ret, 2, function(x) {
      tokens <- unique(strsplit(paste(x, collapse = ""), "")[[1]])
      length(tokens[tokens %in% 0:9])
    })
    just01 <- apply(ret, 2, function(x) {
      tokens <- unique(strsplit(paste(x, collapse = ""), "")[[1]])
      length(tokens[tokens %in% 0:1]) == 2 &&
        length(tokens[tokens %in% 2:9]) == 0
    })
    
    
    neoFile <- MatrixFile(pID, "neo.nex")
    transFile <- MatrixFile(pID, "trans.nex")
    neoRandFile <- MatrixFile(pID, "neo-rand.nex")
    transRandFile <- MatrixFile(pID, "trans-rand.nex")
    
    # Write characters to file
    .WriteChars <- function(chars, path) {
      write.nexus.data(chars, path, interleaved = FALSE, format = "STANDARD")
      # Remove ape comment - timestamp causes unnecessary commits
      writeLines(readLines(path)[-2], path)
    }
    .WriteChars(ret, nexusFile)
    if (any(informative & type == "n") &&
        any(informative & type == "t")) {
      .WriteChars(ret[, informative & type == "n"], neoFile)
      .WriteChars(ret[, informative & type == "t"], transFile)
      nNeo <- sum(informative & type == "n")
      
      # One approach would be to select two-state characters and relabel them
      # to use the states 0 and 1.
      # But a character with states 1 and 2 would otherwise be treated with M3.
      # .RelabelStates <- function(x) {
      #   tokens <- unique(strsplit(paste(x, collapse = ""), "")[[1]])
      #   starts <- tokens[tokens %in% 0:9]
      #   stopifnot(length(starts) == 2)
      #   x <- gsub(starts[[1]], "!ZERO!", x, fixed = TRUE)
      #   x <- gsub(starts[[2]], "!ONE!", x, fixed = TRUE)
      #   x <- gsub("!ZERO!", 0, x, fixed = TRUE)
      #   x <- gsub("!ONE!", 1, x, fixed = TRUE)
      #   x
      # }
      # randomNeo <- sort(sample(which(informative & kObs == 2), nNeo))
      # .WriteChars(apply(ret[, randomNeo], 2, .RelabelStates), neoRandFile)
      
      set.seed(pID)
      randomNeo <- sort(sample(which(informative & just01), nNeo))
      .WriteChars(ret[, randomNeo], neoRandFile)
      .WriteChars(ret[, informative & !tabulate(randomNeo, length(informative))]
                  , transRandFile)
    } else {
      message("Need informative neomorphic and transformational characters")
      if (createdOutDir) {
        message(outDir)
        message(dir.exists(outDir))
        system2("git", sprintf("submodule deinit %s", outDir))
        system2("git", sprintf("rm %s", outDir))
        unlink(outDir, recursive = TRUE)
        unlink(outDir)
      }
      return(FALSE)
    }

  } else {
    message(nexusFile, " already exists")
  }
  
  return(TRUE)
}
