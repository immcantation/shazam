CODON_TABLE <- shazam:::computeCodonTable(aminoAcidClasses = NULL)
save(CODON_TABLE,file="data-raw/CODON_TABLE.RData")
