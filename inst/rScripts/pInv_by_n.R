source("R/pInv.R")
for (pID in KiProjects()) {
  MakeInv(pID)
  CalcPInf(pID, "by_n_ki")
}
