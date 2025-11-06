#' Construct well-corroborated tree for Field et al. 2020 dataset (syab07200)

library("TreeTools")

bkTrees <- read.nexus(file.path(neotrans::OutputDir(), "wct", "BraunKimball.tre")) |>
  lapply(function(tr) {
    lab <- tr[["tip.label"]]
    lab <- sub("_I\\d+.*$", "", lab)
    lab <- sub("Colius_", "Colius.", lab)
    lab <- sub("_.+$", "", lab)
    tr[["tip.label"]] <- lab
    tr
  })

bkNames <- names(bkTrees)
# Selected subset based on perspective of Braun & Kimball 2021,
# with reference to their Fig. S1
includedTrees <- c("AllData_Unpart", "AllData_Part", "AllData_RY",
                   "Exons_RY", "NonCoding_Unpart", "NonCoding_Part",
                   "NonCoding_RY")

tax00 <- sub("_.*$", "", rownames(ReadCharacters("matrices/syab07200.nex")))
translate <- c(
  Acryllium = "Numida", # Numididae
  # Anhima = "Chauna", # Anhimidae: Already have Chauna
  Antigone = "Grus", # former synonym
  # Coturnix = "Gallus", # Phasianidae: Non-erectile clade: Already have Gallus
  #  Alectura  = "Leipoa", # Megapodiidae: Already have Leipoa
  #  Macrocephalon = "Leipoa", # Megapodiidae: Already have Leipoa
  #  Megapodius = "Leipoa", # Megapodiidae: Already have Leipoa
  Phasianus = "Bonasa", # Phasianidae: Erectile clade
  Tadorna = "Oxyura", # Anatidae - five candidates
  Ichthyornis = "Caiman", # outgroup
  
  Null = NULL
)

bk00 <- lapply(bkTrees, function (tr) {
  tr$tip.label[tr$tip.label %in% translate] <- names(translate)[match(tr$tip.label[tr$tip.label %in% translate], translate)]
  KeepTip(tr, intersect(tr$tip.label, tax00))
}) |> setNames(bkNames)

plot(Consensus(bk00, p = 1), cex = 0.8)
write.tree(Consensus(bk00, p = 1))
