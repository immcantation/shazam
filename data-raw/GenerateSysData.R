#### Region boundaries ####
#VLENGTH <- 312

#### Nucleotide characters ####
#NUCLEOTIDES <- c("A", "C", "G", "T", "N", "-", ".")

#### Codon translations ####

AMINO_ACIDS <- c("TTT"="F", "TTC"="F",
                 "TTA"="L", "TTG"="L", "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
                 "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S", "AGT"="S", "AGC"="S",
                 "TAT"="Y", "TAC"="Y",
                 "TGT"="C", "TGC"="C",
                 "TGG"="W",
                 "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
                 "CAT"="H", "CAC"="H",
                 "CAA"="Q", "CAG"="Q",
                 "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R", "AGA"="R", "AGG"="R",
                 "ATT"="I", "ATC"="I", "ATA"="I",
                 "ATG"="M",
                 "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
                 "AAT"="N", "AAC"="N",
                 "AAA"="K", "AAG"="K",
                 "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
                 "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
                 "GAT"="D", "GAC"="D", 
                 "GAA"="E", "GAG"="E",
                 "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G",
                 "TAA"="*", "TAG"="*", "TGA"="*")

#### Amino acid classes ####

# Load data
# http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
aa_imgt <- read.csv("data-raw/IMGT_AminoAcidClasses.csv", as.is=TRUE)

# Set vectors or hydropathy, polarity and charge classes
aa_ambig <- setNames(rep(NA, 4), c("X", "-", ".", "*"))
AMINO_ACIDS_HYDROPATHY <- c(setNames(aa_imgt$HYDROPATHY, aa_imgt$IUPAC), aa_ambig)
AMINO_ACIDS_POLARITY <- c(setNames(aa_imgt$POLARITY, aa_imgt$IUPAC), aa_ambig)
AMINO_ACIDS_CHARGE <- c(setNames(aa_imgt$CHARGE, aa_imgt$IUPAC), aa_ambig)

#### Load other saved data ####

load("data-raw/CDR_Nuc_Mat.RData")
load("data-raw/FWR_Nuc_Mat.RData")
load("data-raw/CODON_TABLE.RData")
load("data-raw/BAYESIAN_FITTED.RData")
load("data-raw/CONST_I.RData")

#### Save to R/sysdata.rda ####
devtools::use_data(CDR_Nuc_Mat, 
                   FWR_Nuc_Mat, 
                   CODON_TABLE,
                   AMINO_ACIDS,
                   BAYESIAN_FITTED,
                   CONST_I,
                   AMINO_ACIDS_HYDROPATHY,
                   AMINO_ACIDS_POLARITY,
                   AMINO_ACIDS_CHARGE,
                   internal=TRUE, overwrite=TRUE)
