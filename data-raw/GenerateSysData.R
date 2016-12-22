library(shazam)

#### Region boundaries ####

VLENGTH <- 312

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

#### Distance matrices ####

HH_S1F_Distance <- calcTargetingDistance(model=HH_S1F)
HKL_S1F_Distance <- calcTargetingDistance(model=HKL_S1F)
MK_RS1NF_Distance <- calcTargetingDistance(model=MK_RS1NF)
HH_S5F_Distance <- calcTargetingDistance(model=HH_S5F)
HKL_S5F_Distance <- calcTargetingDistance(model=HKL_S5F)
MK_RS5NF_Distance <- calcTargetingDistance(model=MK_RS5NF)

e <- new.env()
load("data-raw/HS1F_v0.1.4.rda", envir=e)
load("data-raw/M1N_v0.1.4.rda", envir=e)
HS1F_Compat <- get("HS1FDistance", envir=e)
M1N_Compat <- get("M1NDistance", envir=e)
rm(e)

#### Load other saved data ####

load("data-raw/CDR_Nuc_Mat.RData")
load("data-raw/FWR_Nuc_Mat.RData")
load("data-raw/CODON_TABLE.RData")
load("data-raw/BAYESIAN_FITTED.RData")
load("data-raw/CONST_I.RData")

#### Save to R/sysdata.rda ####
devtools::use_data(NUCLEOTIDES,
                   VLENGTH,
                   CDR_Nuc_Mat, 
                   FWR_Nuc_Mat, 
                   BAYESIAN_FITTED,
                   CONST_I,
                   CODON_TABLE,
                   AMINO_ACIDS,
                   HH_S1F_Distance,
                   HH_S5F_Distance,
                   HKL_S1F_Distance,
                   HKL_S5F_Distance,
                   MK_RS1NF_Distance,
                   MK_RS5NF_Distance,
                   HS1F_Compat,
                   M1N_Compat,
                   internal=TRUE, overwrite=TRUE)
