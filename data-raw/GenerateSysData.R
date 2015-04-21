#### Nucleotide characters ####
NUCLEOTIDES <- c("A", "C", "G", "T", "N", "-", ".")

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

# Load data
load("data-raw/CDR_Nuc_Mat.RData")
load("data-raw/FWR_Nuc_Mat.RData")
load("data-raw/CODON_TABLE.RData")
load("data-raw/BAYESIAN_FITTED.RData")
load("data-raw/CONST_I.RData")

#### Save to R/sysdata.rda ####
devtools::use_data(CDR_Nuc_Mat, 
                   FWR_Nuc_Mat, 
                   CODON_TABLE,
                   NUCLEOTIDES,
                   AMINO_ACIDS,
                   BAYESIAN_FITTED,
                   CONST_I,                   
                   internal=TRUE, overwrite=TRUE)
