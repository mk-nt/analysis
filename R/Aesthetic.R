# Labels
modelLabel <- c(
  "bestHom" = "Best homog.",
  
  "by_ki" = "StMk",
  "by_n_ki" = "StN",
  "by_nn_ki" = "StN\U2019",
  "by_t_ki" = "StT",
  "by_nt_ki" = "StNT",
  
  "rm_by_n_ki" = "StN-shuffled",
  "rm_by_t_ki" = "StT-shuffled",
  "rm_by_nt_ki" = "StNT-shuffled",
  
  "ns_ki" = "NstMk",
  "ns_n_ki" = "NstN",
  "ns_t_ki" = "NstT",
  "ns_nt_ki" = "NstNT",
  
  "bestHet" = "Best heterog.",
  "hg_ki" = "Het",
  "hg2_ki" = "Het1P",
  "hg_b_ki" = "HetB",
  "hg_m_ki" = "HetM",
  "hg_bm_ki" = "HetBM")
modelLabel[sub("ki", "kv", names(modelLabel))] <- modelLabel[names(modelLabel)]

# Palettes
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

modelCol <- c(
  # Baseline model
  by_ki = "#595959",  # neutral dark grey (less bluish than rm_ shades)
  
  bestHom = "#4F7B7B",
  bestStat = "#33A1C9", # Placeholder
  
  # Stationary, homogeneous
  by_n_ki  = "#0072B2",
  by_nn_ki = "#33A1C9",
  by_t_ki  = "#009E73",
  by_nt_ki = "#005C66",
  
  # Stationary, shuffled (“rm_”) – cooler bluish-slate range
  rm_by_n_ki  = "#6B7A8F",
  rm_by_t_ki  = "#50606F",
  rm_by_nt_ki = "#8798AB",

  # Non-stationary, homogeneous
  bestNS   = "#E69F00", # Placeholder
  
  ns_ki    = "#E69F00",
  ns_n_ki  = "#F0B429",
  ns_t_ki  = "#D17C00",
  ns_nt_ki = "#B44D00",

  # Heterogeneous
  bestHet  = "#B05FAE",  # mid magenta–violet: “average heterogeneous”
  hg_ki    = "#CC79A7",
  hg2_ki   = "#AA4E8B",
  hg_b_ki  = "#7A6FE1",
  hg_m_ki  = "#A463C5",
  hg_bm_ki = "#68246D",
  
  # Generic colours
  white    = "#ffffff"
)

#' Model colour
#' Return the colour used to plot results of this model
#' @export
ModelCol <- function(x) {
  x <- as.character(x)
  ret <- modelCol[x]
  na_idx <- is.na(ret)
  if (any(na_idx)) {
    fixed_keys <- gsub("_kv", "_ki", x[na_idx], fixed = TRUE)
    ret[na_idx] <- modelCol[fixed_keys]
  }
  na_idx <- is.na(ret)
  if (any(na_idx)) {
    ret[na_idx] <- "grey40"
  }
  ret
}

#' Model label
#' Translates a model script identifier to a human-readable label
#' @export
ModelLabel <- function(x) {
  x <- as.character(x)
  ret <- modelLabel[x]
  na_idx <- is.na(ret)
  if (any(na_idx)) {
    fixed_keys <- gsub("_kv", "_ki", x[na_idx], fixed = TRUE)
    ret[na_idx] <- modelLabel[fixed_keys]
  }
  na_idx <- is.na(ret)
  if (any(na_idx)) {
    ret[na_idx] <- x[na_idx]
  }
  ret
}

#' Taxon colour
#' Return an appropriate colour for each higher taxon
#' @export
TaxonCol <- function(x = 1:4) c(
  "Euarthropoda"  = "#C4A484",
  "Invertebrates" = "#468C9B",
  "Vertebrates"   = "#E6A3A1",
  "Non-animals"   = "#6BAF5E"
)[x]
