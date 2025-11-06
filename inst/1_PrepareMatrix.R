# This file gives an example of how analyses might be queued given a new set of
# models.
# More details are available in 
# vignette("matrix-processing", package = "neotrans")

library("neotrans")

# Projects that are fast to analyse
# When prototyping a model it can be useful to analyse these matrices first
fast <- c("104", "175", "450", "493", "635", "675", "706", "748", "950", 
          "1113", "1210", "1271", "2131", "2553", "2800", "3199", "3244", 
          "3351", "3405", "3408", "3445", "3603", "3646", "3655", "3710", 
          "3832", "3833", "3929", "4111", "4220", "4291", "4305", "4308", 
          "4309", "4310", "4467", "4747", "4761", "4790", "4867", "4910", 
          "5099", "5186", "5201", "5230", "5255", "5268")
med <- c("157", "563", "3200", "3392", "3705", "3711", "4230", "5327")
# Codes for matrices with well corroborated trees from Asher & Smith (2022)
syab <- .config$syab
# Remaining projects
slow <- setdiff(AllProjects(), c(fast, med, syab))

# Check matrix directory exists; if not, edit the .config variables in R/zzz.R
stopifnot(dir.exists(.config$matrixDir))

# Paths to Excel spreadsheets that annotate character types,
metaFiles <- grep(sprintf(.config$metaPattern, "\\d"),
                  list.files(.config$matrixDir, "*.xlsx"),
                  value = TRUE, perl = TRUE)

# Begin MCMC analyses on all informative-conditioned projects and models.
# In practice one may wish to work with a subset first
EnqueueMC(KiProjects(),
          c("by_ki", "by_t_ki", "by_n_ki", "by_nn_ki", "by_nt_ki",
            "rm_by_n_ki", "rm_by_t_ki", "rm_by_nt_ki",
            "ns_ki", "ns_n_ki", "ns_nt_ki",
            "hg_ki", "hg2_ki", "hg_b_ki", "hg_m_ki", "hg_bm_ki"))

# Estimate marginal likelihoods for all informative-conditioned projects and
# models.
EnqueueML(KiProjects(),
          c("by_ki", "by_t_ki", "by_n_ki", "by_nn_ki", "by_nt_ki",
            "rm_by_n_ki", "rm_by_t_ki", "rm_by_nt_ki",
            "ns_ki", "ns_n_ki", "ns_nt_ki",
            "hg_ki", "hg2_ki", "hg_b_ki", "hg_m_ki", "hg_bm_ki"))
