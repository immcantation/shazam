# AMINO_ACID_PROPERTIES <- 
#   structure(
#     list(
#       Name = structure(c(1, 5, 4, 6, 14, 8, 9, 10, 12, 11, 13, 3, 15, 7, 2, 16, 17, 20, 18, 19), 
#                        .Label = c("Alanine", "Arginine", "Asparagine", "Aspartic Acid", "Cysteine", 
#                                   "Glutamic Acid", "Glutamine", "Glycine", "Histidine", "Isoleucine", 
#                                   "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline", "Serine", 
#                                   "Threonine", "Tryptophan", "Tyrosine", "Valine"), class = "factor"), 
#       
#       One.Letter.Code = structure(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20), 
#                              .Label = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", 
#                                         "R", "S", "T", "V", "W", "Y"), class = "factor"), 
#       
#       Three.Letter.Code = structure(c(1, 5, 4, 7, 14, 8, 9, 10, 12, 11, 13, 3, 15, 6, 2, 16, 17, 20, 18, 19), 
#                                     .Label = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", 
#                                                "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
#                                                "Tyr", "Val"), class = "factor"), 
#       
#       Polarity = structure(c(1,2, 2, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 2, 2, 1, 1, 2), 
#                            .Label = c("Nonpolar","Polar"), class = "factor"), 
#       
#       Hydropathy = structure(c(2,2,1,1,2,3,3,2,1,2,2,1,3,1,1,3,3,2,2,3), 
#                              .Label = c("Hyrophilic","Hydrophobic", "Neutral"), class = "factor")
# 
#     ),
#     .Names = c("Name", "One.Letter.Code", "Three.Letter.Code", "Polarity", "Hydropathy"), 
#     class = "data.frame", 
#     row.names = c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
#                   "13", "14", "15", "16", "17", "18", "19", "20")
#   )
# save(AMINO_ACID_PROPERTIES, file="data-raw/AMINO_ACID_PROPERTIES.RData")

#### Amino acid classes ####

# Load data
# http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
aa_imgt <- read.csv("data-raw/IMGT_AminoAcidClasses.csv", as.is=TRUE)

# Set vectors or hydropathy, polarity and charge classes
aa_ambig <- setNames(rep(NA, 4), c("X", "-", ".", "*"))
AMINO_ACID_HYDROPATHY <- c(setNames(aa_imgt$HYDROPATHY, aa_imgt$IUPAC), aa_ambig)
AMINO_ACID_POLARITY <- c(setNames(aa_imgt$POLARITY, aa_imgt$IUPAC), aa_ambig)
AMINO_ACID_CHARGE <- c(setNames(aa_imgt$CHARGE, aa_imgt$IUPAC), aa_ambig)

devtools::use_data(AMINO_ACID_HYDROPATHY,
                   AMINO_ACID_POLARITY,
                   AMINO_ACID_CHARGE,
                   overwrite=TRUE)