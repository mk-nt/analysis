#' Depict predictive power of variable for output
#' @param meta Complete metadata available for prediction
#' @param predictor Metadata column to use for predicting
#' @param value Value to be predicted
#' @param rf Random forest, produced using
#' `randomForest(meta[, <all predictors>], value, importance = TRUE)`.
#' @param ylab,col Graphical parameters
#' @importFrom ppcor pcor.test
#' @export
PredictWith <- function(meta, predictor, value, rf, bold = FALSE,
                        ylab = "", col = 1) {
  plot(value ~ meta[, predictor], las = 2, frame.plot = FALSE,
       xlab = Decrypt(predictor), ylab = ylab,
       col = switch(
         predictor,
         "taxon" = TaxonCol(levels(meta[, predictor])),
         "rank" = rev(hcl.colors(5 + length(levels(meta[, "rank"])))),
         1))
  abline(h = 0, lty = "dotted")
  
  x <- do.call(partialPlot,
               list(x = rf, pred.data = meta, x.var = predictor, plot = FALSE))
  x$x
  lines(x$x, x$y, col = 2, lwd = 1.2)
  
  suppressWarnings(ct <- cor.test(as.numeric(meta[, predictor]), value,
                                  method = "spearman"))
  rho <- if (ct$p.value > 0.05) {
    message(predictor, ": No significant correlation")
    ""
  } else {
    message(predictor, ": Correlated (p = ", signif(ct$p.value, 2), "): ",
            ct$method, " ~ ", signif(ct$estimate, 2))
    sprintf("\n\U03C1 = %.2f", ct$estimate)
  }
  
  pcMeta <- meta[, vapply(meta, is.ordered, TRUE) |
                   vapply(meta, is.numeric, logical(1))]
  ordered <- vapply(pcMeta, is.ordered, TRUE)
  pcMeta[, ordered] <- sapply(meta[, vapply(pcMeta, is.ordered, TRUE)], as.numeric)
  if (predictor %in% colnames(pcMeta)) {
    suppressWarnings(
      pc <- pcor.test(value, pcMeta[, predictor],
                      pcMeta[, setdiff(colnames(pcMeta), predictor)],
                      method = "spearman")
    )
    pcor <- if (pc$p.value > 0.05) {
      message("      ... No significant partial correlation")
      ""
    } else {
      message("      ... Correlated (p = ", signif(pc$p.value, 2),
              "): pcor ~ ", signif(pc$estimate, 2))
      sprintf("pcor = %.2f", pc$estimate)
    }
  } else {
    pcor <- ""
  }
  text(par("usr")[2], par("usr")[4], pos = 2, font = 1 + bold,
       sprintf("%.0f%% IncMSE%s%s%s", round(importance(rf)[predictor, "%IncMSE"], 0),
               rho, if(rho != "" && pcor != "") "; " else
                 if (pcor == "") "" else "\n", pcor), xpd = NA)
  
  invisible(ct)
}
