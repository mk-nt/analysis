source("R/_setup.R")
source("R/CheckComplete.R")

for (pID in AllProjects()) {
  # Fetch results from GitHub and compute metadata
  UpdateRecords(pID, "by_ki")
  UpdateRecords(pID, "by_n_ki")
  UpdateRecords(pID, "by_t_ki")
  UpdateRecords(pID, "by_nt_ki")
}
