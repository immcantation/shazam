# Load data
load("data-raw/CDR_Nuc_Mat.RData")
load("data-raw/FWR_Nuc_Mat.RData")
load("data-raw/CODON_AA_TABLE.RData")
load("data-raw/CODON_TABLE.RData")

# Save to R/sysdata.rda
devtools::use_data(CDR_Nuc_Mat, FWR_Nuc_Mat, CODON_AA_TABLE, CODON_TABLE, 
                   internal=TRUE, overwrite=TRUE)
