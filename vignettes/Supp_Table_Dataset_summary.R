# Render summary of data and analyses as supplementary HTML table

library("DT")
library("htmltools")
library("htmlwidgets")
library("neotrans")

statModels <- c("by_ki", "by_n_ki", "by_nn_ki", "by_t_ki", "by_nt_ki")
rmModels <- c("rm_by_n_ki", "rm_by_t_ki", "rm_by_nt_ki")
nsModels <- c("ns_ki", "ns_n_ki", "ns_t_ki", "ns_nt_ki")
hetModels <- c("hg_ki", "hg2_ki", "hg_b_ki", "hg_m_ki", "hg_bm_ki")
homModels <- c(statModels, nsModels)
kiModels <- c(statModels, rmModels, nsModels, hetModels)
nonRm <- setdiff(kiModels, rmModels)
paramModels <- c(statModels, nsModels, "hg_ki")

projects <- KiProjects()

marginals <- vapply(kiModels, function(scriptID) 
  vapply(projects, GetMarginal, numeric(1), scriptID),
  setNames(double(length(projects)), projects)) |>
  `colnames<-`(kiModels) |>
  t()

results <- lapply(projects, function(pID)
  setNames(lapply(kiModels, function(model)
    ExistingResults(pID, model)), kiModels)) |>
  setNames(projects)


allIn <- colSums(is.na(marginals)) == 0

# Check which marginal likelihoods are unavailable.
# Diagnose and record in MakeSlurm.R
marginals[, !allIn]

bestModel <- rep(NA, length(allIn))
bestModel[allIn] <- kiModels[apply(marginals[kiModels, allIn], 2, which.max)]

hasML <- names(allIn)[allIn]
hasConv <- rowSums(!sapply(paramModels, function(scriptID) sapply(projects, HasConverged, scriptID))) == 0
hasConv <- names(hasConv)[hasConv]

mlEss <- intersect(hasML, hasConv)
length(hasML)
length(hasConv)
length(mlEss)
projOK <- names(results) %in% mlEss


GetPar <- function(model, parameter, projects = names(results),
                   statistic = "50%") {
  FUN.VALUE <- matrix(NA_real_, length(statistic), length(parameter),
                      dimnames = list(statistic, parameter))
  vapply(results[projects], function(x) {
    par <- x[[c(model, "parameters")]]
    if (is.null(par)) {
      FUN.VALUE
    } else {
      par[statistic, parameter]
    }
  }, FUN.VALUE)
}

mlEss <- intersect(hasML, hasConv)
lossN <- GetPar("by_n_ki", "rate_loss", proj = mlEss)
lossNT <- GetPar("by_nt_ki", "rate_loss", proj = mlEss)
neoT <- GetPar("by_t_ki", "rate_neo", proj = mlEss)
neoNT <- GetPar("by_nt_ki", "rate_neo", proj = mlEss)

NotOne <- function(model, parameter) {
  ranges <- GetPar(model, parameter, mlEss, c("2.5%", "97.5%"))
  isOne <- ranges["2.5%", , ] < 1 & ranges["97.5%", , ] > 1 
  message(sum(!isOne, na.rm = TRUE), " / ", sum(!is.na(isOne)), " ",
          parameter, " values unlikely to be 1")
  !isOne
}

lossNot1NT <- NotOne("by_nt_ki", "rate_loss")
neoNot1NT <- NotOne("by_nt_ki", "rate_neo")
lossNot1 <- NotOne("by_n_ki", "rate_loss")
neoNot1 <- NotOne("by_t_ki", "rate_neo")


LinkMB <- function(pID) {
  mbRoot <- "https://morphobank.org/index.php/Projects/Matrices/project_id"
  sprintf("<a href=\"%s/%s\" target=\"blank\">%s</a>", mbRoot, pID, pID)
}

LinkData <- function(pID) {
  ghRoot <- "http://github.com/mk-nt"
  paste(sprintf("<a href=\"%s/%s_%s\" target=\"blank\">%s</a>", ghRoot, pID, models,
                modelName[models]), collapse = " | ")
}

Sig3 <- function(x) {
  ifelse(x == 0, 0, round(x, digits = -floor(log10(abs(x))) + 2))
}

Round1 <- function(n) ifelse (n == 0, "\U2013", sprintf("%.1f", n))
BF <- function(scriptID, marginals, invert) {
  Round1((marginals[scriptID, ] - apply(marginals, 2, max)) *
           if (invert) -1 else 1)
}
BFLink <- function(model, invert = TRUE) {
  ghRoot <- "https://github.com/mk-nt"
  sprintf("<a href=\"%s/%s_%s\" target=\"blank\">%s</a>", 
          ghRoot, namesOK, model,
          BF(model, marginals[, namesOK], invert = invert))
}

ModelLink <- function(model) {
  ghRoot <- "https://github.com/mk-nt"
  sprintf("<a href=\"%s/%s_%s\" target=\"blank\">%s</a>", 
          ghRoot, namesOK, model, ModelLabel(model))
}

intro_html <- markdown::markdownToHTML(text = "
## Summary of data and analyses

- [neo-trans/matrices](https://github.com/neo-trans/matrices) contains
  original matrices downloaded
  from MorphoBank (`project<ID>.nex`).
  Each character is marked as neomorphic or
  transformational in the accompanying `Project<ID>_<name>.xlsx` file.  Characters
  that require reformulation are labelled with a proposed reformulation.
  
- These raw files were used to procedurally generate separate matrices
  for neomorphic (`*.neo.nex`) and transformational (`*.trans.nex`) characters,
  and analytical scripts, each housed in a separate GitHub repository
  and analysed on the Hamilton high performance computing cluster, Durham
  University.
  
  Each repository contains two RevBayes scripts:
  - `marginal.Rev`, used to compute the marginal likelihood of the model, and
  - `mcmcmc.Rev`, used to conduct MCMCMC analysis.
  
  The analytical model is defined in `<model_code>.Rev`.
  Each repository also contains logs of parameter estimates (`.log`),
  tree samples (`.trees`), and marginal likelihood estimates (`.pp`).

## Table explanation

- The \"ID\" column contains a link to the original study at
  MorphoBank.
  Links to the analytical files and results are provided in the \"BF\" columns
  of the table below.
  
- Inferred values of _n_ and _t_ are marked with an asterisk where the 95%
  posterior distribution does not include 1.
  
- Columns can be sorted by clicking the arrows, or filtered using the boxes
  at the column tops.

", fragment.only = TRUE)

resOK <- names(results) %in% mlEss
namesOK <- names(results)[resOK]
{# Supplementary reporting
  overview <- data.frame(
    "ID" = LinkMB(names(results[resOK])),
    Taxon = taxon[mlEss],
    Rank = rank[mlEss],
    RankOrder = as.numeric(rank[mlEss]),
    Taxa = nTaxa[mlEss],
    Characters = colSums(nChar)[mlEss],
    "Transf." = nChar["trans", mlEss],
    "Neom." = nChar["neo", mlEss],
    "T:N ratio" = Sig3(nChar["trans", mlEss] / nChar["neo", mlEss]),
    "Inferred n" = paste(Sig3(ifelse(bestModel[resOK] == "by_nt_ki",
                                      lossNT[namesOK], lossN[namesOK])),
                         ifelse(lossNot1[namesOK], "*", "")),
    "Inferred t" = paste(Sig3(ifelse(bestModel[resOK] == "by_nt_ki",
                                     neoNT[namesOK], neoT[namesOK])),
                         ifelse(neoNot1[namesOK], "*", "")),
    "modelRaw" = as.factor(bestModel[resOK]),
    "Optimal model" = as.factor(ModelLink(bestModel[resOK])),
    "BF vs StMk" = BFLink("by_ki", invert = TRUE)
  )
  
  cols <- colnames(overview)
  cols <- gsub(fixed = TRUE, "Transf", "Transf.",
          gsub(fixed = TRUE, "Characters", "Chars",
          gsub(fixed = TRUE, "Neom", "Neom.",
          gsub(fixed = TRUE, "T N ratio", "Ratio",
          gsub("Inferred ([nt])", "Inferred <i>\\1</i>", perl = TRUE,
               gsub("Mk ", "BF, Mk-", fixed = TRUE,
                    gsub(".", " ", cols, fixed = TRUE)))))))
  
  widget <- 
    datatable(
      overview,
      escape = FALSE,
      filter = "top",
      options = list(
        pageLength = sum(resOK),
        dom = "t",
        columnDefs = list(
          list(
            targets = which(names(overview) %in% c("RankOrder", "modelRaw")) - 1,
            visible = FALSE
          ),
          list(
            orderData = which(names(overview) == "RankOrder") - 1,
            targets = which(names(overview) == "Rank") - 1)
        )
      ),
      #caption = "Project overview",
      rownames = FALSE,
      colnames = cols) |>
    formatStyle(
      "Taxon",
      background = styleEqual(names(TaxonCol(1:4)),
                              paste0(TaxonCol(1:4), "66")),
      backgroundSize = "98% 88%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center"
    ) |>
    formatStyle(
      "Rank",
      background = styleEqual(levels(rank),
                              hcl.colors(length(levels(rank)), rev = FALSE,
                                         alpha = 88 / 256)),
      backgroundSize = "98% 88%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center"
      ) |>
    formatStyle(
      "Taxa",
      background = styleColorBar(c(0, max(overview$Taxa)),
                                 durham$concrete, angle = 270),
      backgroundSize = "98% 88%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center"
      ) |>
    formatStyle(
      "Transf.",
      background = styleColorBar(c(0, max(overview$Characters)), 
                                 durham$concrete, angle = 90),
      backgroundSize = "98% 88%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center") |>
    formatStyle(
      "Neom.",
      background = styleColorBar(c(0, max(overview$Characters)),
                                 durham$concrete, angle = 270),
      backgroundSize = "98% 88%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center") |>
    formatStyle(
      "Optimal.model",
      valueColumns = "modelRaw",
      background = styleEqual(overview$modelRaw,
                              paste0(ModelCol(overview$modelRaw), "66")),
      backgroundSize = "98% 88%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center"
      )
  
  saveWidget(widget, "overview.html", selfcontained = TRUE,
             title = "Overview of results")
  lines <- readLines("overview.html")
  writeLines(
    append(lines, intro_html, which(grepl("id=\"htmlwidget_container", lines)) - 1),
    "overview.html")
}
