source("R/pInv.R")
for (pID in KiProjects()[-seq_len(match(563, KiProjects()) - 1)]) {
  MakeInv(pID)
  CalcPInf(pID, "hg_b_ki")
}
