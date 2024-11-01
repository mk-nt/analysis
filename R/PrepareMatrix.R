#' @param pID Character giving pID identifier, e.g. "123"
#' @param overwrite Logical; should existing files be overwritten?
PrepareMatrix <- function(pID, metaFile, overwrite = FALSE) {
  projectID <- sprintf("pID%s", pID)
  nexusFile <- MatrixFile(pID, "nex")
  outDir <- dirname(nexusFile)
  if (!dir.exists(outDir)) {
    dir.create(outDir)
  }
  
  if (overwrite || !file.exists(nexusFile)) {
    
    meta <- readxl::read_xlsx(paste0(matrixDir, "/", metaFile))
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
    
    chars <- TreeTools::ReadCharacters(paste0(matrixDir, "/",
                                              basename(nexusFile)))
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
        list(type = newType, k = newK, n = nNew,
             from = vapply(translate, "[[", character(1), 1),
             to = vapply(translate,
                         function(x) strsplit(trimws(x[[2]]), "")[[1]],
                         character(nNew)))
      }
    })
    nNew <- rep(1, nChars)
    nNew[special] <- vapply(translations, "[[", integer(1), "n")
    index <- cumsum(nNew)
    type <- character(sum(nNew))
    type[index[neo]] <- "n"
    type[index[trans]] <- "t"
    nowSpecial <- type == ""
    k <- numeric(sum(nNew))
    
    k[index[trans]] <- unseen[trans]
    k[index[neo]] <- NA
    k[type == ""] <- unlist(lapply(translations, "[[", "k"))
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
      illegalMat <- matrix(illegalNeomorphic, nrow(ret))
      chars <- colSums(illegalMat) > 0
      taxa <- rowSums(illegalMat) > 0
      stop("Illegal neomorphic characters in char ",
           paste(which(type == "n")[chars], collapse = ", "),
           ", taxa ",
           paste(rownames(ret)[taxa], collapse = ", ")
      )
    }
    
    # MrBayes handles up to 10 states
    k[k == hasUnseen] <- 10
    lacksUnseen <- !is.na(k) & k == noUnseen
    if (any(lacksUnseen)) {
      suppressWarnings(
        # If a character is solely "?"s then we might encounter a warning.
        # "?" is not numeric but is less than numbers
        k[lacksUnseen] <- as.numeric(apply(ret[, lacksUnseen, drop = FALSE], 2,
                                           sort, decreasing = TRUE)[1, ]) + 1
      )
    }
    kRet <- ret
    newTrans <- type == "t"
    kRet[, newTrans] <- vapply(seq_len(sum(newTrans)), function(i) {
      n <- which(newTrans)[[i]]
      x <- kRet[, n]
      tokens <- unique(unlist(strsplit(gsub("\\D", "", paste0(x)), "")))
      if (length(tokens)) {
        for (j in seq_along(tokens)) {
          x <- gsub(tokens[[j]], letters[[j]], x, fixed = TRUE)
        }
        for (j in seq_len(length(tokens) - 1)) {
          x <- gsub(letters[[j]], j - 1, x, fixed = TRUE)
        }
        gsub(letters[[length(tokens)]], k[[n]] - 1L, x, fixed = TRUE)
      } else {
        x
      }
    }, character(nTax))
    informative <- apply(kRet, 2, function(x) {
      tab <- table(x[nchar(x) == 1 & x %in% 0:9])
      length(tab[tab > 1]) > 1
    })
    
    
    neoFile <-  MatrixFile(pID, "neo.nex")
    transFile <-  MatrixFile(pID, "trans.nex")
    
    # Write characters to file
    .WriteChars <- function(chars, path) {
      ape::write.nexus.data(chars, path, interleaved = FALSE,
                            format = "STANDARD")
      # Remove ape comment - timestamp causes unnecessary commits
      writeLines(readLines(path)[-2], path)
    }
    .WriteChars(ret, nexusFile)
    if (any(informative & type == "n") &&
        any(informative & type == "t")) {
      .WriteChars(ret[, informative & type == "n"], neoFile)
      .WriteChars(ret[, informative & type == "t"], transFile)
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
