wcTrees <- setNames(unclass(ape::read.tree(
  system.file("wct/wellCorroboratedTrees.nwk", package = "neotrans"))),
                    paste0("0720", 0:6))

outgroup <- list(
  "07200" = "Ichthyornis",
  "07201" = "aaCrocodylia",
  "07202" = c("Tupaia", "Dermoptera"),
  "07203" = c("Didelphis", "Macropus"),
  "07204" = c("Monodelphis", "Sarcophilus"),
  "07205" = c("Monodelphis", "Sarcophilus"),
  "07206" = "Ornithorhynchus")

tipGp <- list(
  c("Didelphis", "Macropus", "aaCrocodylia", "Monodelphis", "Sarcophilus",
    "Ornithorhynchus", "Ichthyornis"),
  c("Dasypus", "Choloepus", "Bradypus", "Tamandua",
    # Livezey Zusi: Palaeognaths
    "Dromaius", "Struthio", "Rhea", "Apteryx", "Eudromia", "Casuarius",
    # Additional Palaeognaths
    "Tinamus"
    ),
  c("Procavia", "Orycteropus", "Echinops", "Rhynchocyon",  "Elephantulus",
    "Loxodonta", "Procavia", "Trichechus", "Amblysomus", "Echinops",
    # LZ: Anseriformes
    "Anser", "Gallus", "Numida", "Chauna", "Anas", "Anseranas", "Ortalis",
    "Tadorna", "Dendrocygna", "Phasianus", "Acryllium", "Crax", "Leipoa"
    ),
  c("Sorex", "Solenodon", "Canis", "Manis", "Felis", "Mustela", "Myotis",
    "Pteropus", "Vicugna", "Sus", "Bos", "Ovis", "Lama", "Hippopotamus", "Equus",
    "Tapirus", "Ceratotherium", "Talpa", "Erinaceus", "Artibeus", "Ailuropoda",
    "Mustela", "Vicugna", "Megaptera", "Tursiops",
    # LZ Gruiformes
    "Porphyrio", "Antigone", "Cariama", "Burhinus"
    ),
  c("Tupaia", "Ptilocercus", "Ochotona", "Cavia", "Anomalurus", "Aplodontia",
    "Homo", "Mus", "Rattus", "Lemur", "Hystrix", "Oryctolagus",
    "Galeopterus", "Tarsius", "Papio", "Nomascus", "Pongo", "Gorilla", "Pan",
    "Macaca", "Chlorocebus", "Callithrix", "Otolemur", "Microcebus",
    "Myocastor", "Ictidomys")
) |> lapply(unique) |> lapply(sort)

tipCol <- setNames(character(sum(lengths(tipGp))), unlist(tipGp))
tipCol[tipGp[[1]]] <- durham$ink # Originally #050301
tipCol[tipGp[[2]]] <- "#5d1114"
tipCol[tipGp[[3]]] <- "#da2426"
tipCol[tipGp[[4]]] <- "#f1501b"
tipCol[tipGp[[5]]] <- "#deae13"

#' Unify tip label format
#' @param tr tree whose labels should be rewritten
#' @return `DeZZ()` returns a tree whose labels have been rewritten in a
#' consistent format, with underscores replaecd with spaces and leading `zz`
#' prefixes removed.
#' @export
DeZZ <- function(tr) {
  lab <- tr[["tip.label"]]
  zz <- startsWith(lab, "zz")
  lab[zz] <- substr(lab[zz], 3, nchar(lab[zz]))
  tr[["tip.label"]] <- sub("([^_]+)_.*", "\\1", lab, perl = TRUE)
  tr
}
