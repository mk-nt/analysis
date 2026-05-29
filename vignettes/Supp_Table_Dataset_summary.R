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
if (any(!allIn)) marginals[, !allIn]

bestModel <- rep(NA, length(allIn))
bestModel[allIn] <- kiModels[apply(marginals[kiModels, allIn], 2, which.max)]

hasML <- names(allIn)[allIn]
hasConv <- rowSums(!sapply(paramModels, function(scriptID) sapply(projects, HasConverged, scriptID))) == 0
if (any(!hasConv)) {
  awaiting <- names(hasConv)[!hasConv]
  message("Still waiting for: ", paste(awaiting, collapse = ", "))
  for (pID in awaiting) for (scriptID in paramModels) {
    if (!HasConverged(pID, scriptID)) {
      message("Queuing up ", pID, " ", scriptID)
      MakeSlurm(pID, scriptID)
    }
  }
}
hasConv <- names(hasConv)#[hasConv]

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
.meta <- Metadata()
{# Supplementary reporting
  overview <- data.frame(
    "ID" = LinkMB(names(results[resOK])),
    Taxon = .meta$taxon[mlEss],
    Rank = .meta$rank[mlEss],
    RankOrder = as.numeric(.meta$rank[mlEss]),
    Taxa = .meta$nTaxa[mlEss],
    Characters = colSums(.meta$nChar)[mlEss],
    "Transf." = .meta$nChar["trans", mlEss],
    "Neom." = .meta$nChar["neo", mlEss],
    "T:N ratio" = Sig3(.meta$nChar["trans", mlEss] / .meta$nChar["neo", mlEss]),
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
  
  rankLevels <- levels(.meta$rank)
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
      background = styleEqual(rankLevels,
                              hcl.colors(length(rankLevels), rev = FALSE,
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

{# Supplementary Bayes factor matrix: one row per project, one column per model.
 # Each cell links to the GitHub repository, labelled with the Bayes factor
 # against the best-supported model for that project (en-dash = best model).
  bfMatrix <- vapply(kiModels, BFLink, character(length(namesOK)),
                     invert = TRUE)
  bfTable <- data.frame(ID = LinkMB(namesOK), bfMatrix, check.names = FALSE)

  # Per-row viridis background: 0 (best) -> bright yellow, max penalty -> dark purple.
  # Same alpha as rank column in widget (88/256 ≈ 0.34) to keep text contrast.
  bfRaw <- vapply(kiModels, function(scriptID) {
    apply(marginals[, namesOK], 2, max) - marginals[scriptID, namesOK]
  }, setNames(numeric(length(namesOK)), namesOK))
  rowMax <- apply(bfRaw, 1, function(x) max(x, na.rm = TRUE))
  bfNorm <- 1 - bfRaw / ifelse(rowMax == 0, 1, rowMax)  # 1 = best, 0 = worst
  pal <- viridisLite::viridis(256, alpha = 88 / 256)
  cellColors <- matrix(
    pal[pmax(1L, pmin(256L, round(replace(bfNorm, is.na(bfNorm), NA) * 255) + 1L))],
    nrow = nrow(bfNorm), dimnames = dimnames(bfNorm)
  )
  cellColors[is.na(bfNorm)] <- "transparent"
  colorJSON <- paste0(
    "[",
    paste(apply(cellColors, 1, function(row)
      paste0("[", paste(sprintf('"%s"', row), collapse = ","), "]")),
      collapse = ","),
    "]"
  )
  rowCallback <- DT::JS(sprintf(
    'function(row, data, index) {
      var c = %s[index];
      $("td", row).each(function(i) {
        if (i > 0) $(this).css("background-color", c[i - 1]);
      });
    }',
    colorJSON
  ))

  headerColors <- c("transparent", paste0(ModelCol(kiModels), "66"))
  headerCallback <- DT::JS(sprintf(
    'function(thead, data, start, end, display) {
      if (!$(thead).find("th:first").text().trim()) return;
      var colors = [%s];
      $(thead).find("th").each(function(i) {
        $(this).css({
          "background-color": colors[i] || "",
          "font-weight": "normal",
          "writing-mode": "vertical-lr",
          "transform": "rotate(180deg)",
          "white-space": "nowrap",
          "text-align": "center",
          "vertical-align": "middle"
        });
      });
    }',
    paste(sprintf('"%s"', headerColors), collapse = ", ")
  ))

  bfWidget <-
    datatable(
      bfTable,
      escape = FALSE,
      class = "nt-bf-table display",
      options = list(
        pageLength = length(namesOK),
        dom = "t",
        autoWidth = FALSE,
        columnDefs = list(
          list(width = "42px", targets = 0L),
          list(
            width = "30px",
            type = "num",
            targets = seq_len(length(kiModels)),
            render = DT::JS("function(data, type, row, meta) {
              if (type === 'sort' || type === 'type') {
                var txt = data.replace(/<[^>]*>/g, '').trim();
                var n = parseFloat(txt);
                return isNaN(n) ? (txt === '\\u2013' ? 0 : 999999) : n;
              }
              return data;
            }")
          )
        ),
        headerCallback = headerCallback,
        rowCallback = rowCallback,
        initComplete = DT::JS('function(settings, json) {
          $(this.api().table().container()).css("font-size", "8pt");
          $("<style>").text(
            ".nt-bf-table { table-layout: fixed !important; width: auto !important } " +
            ".nt-bf-table th, .nt-bf-table td { " +
              "box-sizing: border-box !important; padding: 4px 2px !important; " +
              "width: 30px !important; min-width: 30px !important; max-width: 30px !important " +
            "} " +
            ".nt-bf-table th:first-child, .nt-bf-table td:first-child { " +
              "width: 42px !important; min-width: 42px !important; max-width: 42px !important " +
            "} " +
            ".nt-bf-table td { text-align: center; overflow: hidden; text-overflow: ellipsis; white-space: nowrap } " +
            ".nt-bf-table a { text-decoration: none } " +
            ".nt-bf-table thead th::before, .nt-bf-table thead th::after { display: none !important }"
          ).appendTo("head");
        }'),
        drawCallback = DT::JS('function(settings) {
          var api = new $.fn.dataTable.Api(settings);
          var $tbody = $(api.table().body());
          $tbody.find("tr.nt-repeat").remove();
          var $rows = $tbody.find("tr");
          var n = $rows.length;
          function makeRepeat() {
            return $(api.table().header()).find("tr").clone()
              .addClass("nt-repeat")
              .css("cursor", "default");
          }
          if (n > 1) makeRepeat().insertBefore($rows[Math.floor(n / 2)]);
          makeRepeat().appendTo($tbody);
        }')
      ),
      rownames = FALSE,
      colnames = unname(c("ID", ModelLabel(kiModels)))
    )

  saveWidget(bfWidget, "bf_matrix.html", selfcontained = TRUE,
             title = "Bayes factors by project and model")
  bf_intro_html <- markdown::markdownToHTML(text = "
## Bayes factor comparison of all models for each MorphoBank project

- Each row corresponds to one MorphoBank project (linked via the \"ID\" column).
  Each column corresponds to one of 17 analytical models tested.

- The value in each cell is the log Bayes factor of the best-supported model
  for that project relative to the column model.
  A dash indicates the best-supported model for that row.
  Each cell links to a corresponding GitHub repository containing the analytical
  scripts and artefacts.

- Cell background colour encodes relative model support _within each row_:
  bright yellow = best-supported model; dark purple = least supported.

- Columns can be sorted by clicking the column header.

", fragment.only = TRUE)
  lines <- readLines("bf_matrix.html")
  writeLines(
    append(lines, bf_intro_html,
           which(grepl("id=\"htmlwidget_container", lines)) - 1),
    "bf_matrix.html")
}
