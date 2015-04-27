# Loads and converts model data

#### M3N #####

# Smith DS, et al. Di- and trinucleotide target preferences of somatic mutagenesis 
#    in normal and autoreactive B cells. 
#    J Immunol. 1996 156:2642–52. 

#load("data-raw/MTri_Targeting.RData")

NUCLEOTIDES <- c("A", "C", "G", "T")
triSubstitution <- matrix(c(0, 0.156222928, 0.601501588, 0.242275484, 
                            0.172506739, 0, 0.241239892, 0.586253369, 
                            0.54636291, 0.255795364, 0, 0.197841727, 
                            0.290240811, 0.467680608, 0.24207858, 0),
                          nrow=4, byrow=TRUE, dimnames=list(NUCLEOTIDES, NUCLEOTIDES))
triMutabilityH <- matrix(c(0.24, 1.2, 0.96, 0.43, 2.14, 2, 1.11, 1.9, 0.85, 1.83, 2.36, 1.31, 0.82, 0.52, 0.89, 1.33, 1.4, 0.82, 1.83, 0.73, 1.83, 1.62, 1.53, 0.57, 0.92, 0.42, 0.42, 1.47, 3.44, 2.58, 1.18, 0.47, 0.39, 1.12, 1.8, 0.68, 0.47, 2.19, 2.35, 2.19, 1.05, 1.84, 1.26, 0.28, 0.98, 2.37, 0.66, 1.58, 0.67, 0.92, 1.76, 0.83, 0.97, 0.56, 0.75, 0.62, 2.26, 0.62, 0.74, 1.11, 1.16, 0.61, 0.88, 0.67, 0.37, 0.07, 1.08, 0.46, 0.31, 0.94, 0.62, 0.57, 0.29, NA, 1.44, 0.46, 0.69, 0.57, 0.24, 0.37, 1.1, 0.99, 1.39, 0.6, 2.26, 1.24, 1.36, 0.52, 0.33, 0.26, 1.25, 0.37, 0.58, 1.03, 1.2, 0.34, 0.49, 0.33, 2.62, 0.16, 0.4, 0.16, 0.35, 0.75, 1.85, 0.94, 1.61, 0.85, 2.09, 1.39, 0.3, 0.52, 1.33, 0.29, 0.51, 0.26, 0.51, 3.83, 2.01, 0.71, 0.58, 0.62, 1.07, 0.28, 1.2, 0.74, 0.25, 0.59, 1.09, 0.91, 1.36, 0.45, 2.89, 1.27, 3.7, 0.69, 0.28, 0.41, 1.17, 0.56, 0.93, 3.41, 1, 1, NA, 5.9, 0.74, 2.51, 2.24, 2.24, 1.95, 3.32, 2.34, 1.3, 2.3, 1, 0.66, 0.73, 0.93, 0.41, 0.65, 0.89, 0.65, 0.32, NA, 0.43, 0.85, 0.43, 0.31, 0.31, 0.23, 0.29, 0.57, 0.71, 0.48, 0.44, 0.76, 0.51, 1.7, 0.85, 0.74, 2.23, 2.08, 1.16, 0.51, 0.51, 1, 0.5, NA, NA, 0.71, 2.14), nrow=64,byrow=T)
triMutabilityM  <- matrix(c(1.31, 1.35, 1.42, 1.18, 2.02, 2.02, 1.02, 1.61, 1.99, 1.42, 2.01, 1.03, 2.02, 0.97, 0.53, 0.71, 1.19, 0.83, 0.96, 0.96, 0, 1.7, 2.22, 0.59, 1.24, 1.07, 0.51, 1.68, 3.36, 3.36, 1.14, 0.29, 0.33, 0.9, 1.11, 0.63, 1.08, 2.07, 2.27, 1.74, 0.22, 1.19, 2.37, 1.15, 1.15, 1.56, 0.81, 0.34, 0.87, 0.79, 2.13, 0.49, 0.85, 0.97, 0.36, 0.82, 0.66, 0.63, 1.15, 0.94, 0.85, 0.25, 0.93, 1.19, 0.4, 0.2, 0.44, 0.44, 0.88, 1.06, 0.77, 0.39, 0, 0, 0, 0, 0, 0, 0.43, 0.43, 0.86, 0.59, 0.59, 0, 1.18, 0.86, 2.9, 1.66, 0.4, 0.2, 1.54, 0.43, 0.69, 1.71, 0.68, 0.55, 0.91, 0.7, 1.71, 0.09, 0.27, 0.63, 0.2, 0.45, 1.01, 1.63, 0.96, 1.48, 2.18, 1.2, 1.31, 0.66, 2.13, 0.49, 0, 0, 0, 2.97, 2.8, 0.79, 0.4, 0.5, 0.4, 0.11, 1.68, 0.42, 0.13, 0.44, 0.93, 0.71, 1.11, 1.19, 2.71, 1.08, 3.43, 0.4, 0.67, 0.47, 1.02, 0.14, 1.56, 1.98, 0.53, 0.33, 0.63, 2.06, 1.77, 1.46, 3.74, 2.93, 2.1, 2.18, 0.78, 0.73, 2.93, 0.63, 0.57, 0.17, 0.85, 0.52, 0.31, 0.31, 0, 0, 0.51, 0.29, 0.83, 0.54, 0.28, 0.47, 0.9, 0.99, 1.24, 2.47, 0.73, 0.23, 1.13, 0.24, 2.12, 0.24, 0.33, 0.83, 1.41, 0.62, 0.28, 0.35, 0.77, 0.17, 0.72, 0.58, 0.45, 0.41), nrow=64,byrow=T)

triMutability <- triMutabilityM[, 2]
names(triMutability) <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")

# TODO:  can these both be done without Ns and uses the extend functions instead?
# Create 5-mer substitution matrix
NUCS <- c("A", "C", "G", "T","N")
triSubstitution <- cbind(triSubstitution, rep(NA, 4))
triSubstitution <- rbind(triSubstitution, rep(NA, 5))
colnames(triSubstitution) <- NUCS
rownames(triSubstitution) <- NUCS

withN_5mer <- seqinr::words(5, alphabet=seqinr::s2c("ACGTN"))
m3n_sub <- matrix(NA, nrow=5, ncol=length(withN_5mer), dimnames=list(NUCS, withN_5mer))
for(mer in withN_5mer){
    m3n_sub[, mer]<- triSubstitution[substr(mer, 3, 3), ]
}

# Create 5-mer mutability vector
withN_5mer <- seqinr::words(5, alphabet=seqinr::s2c("ACGTN"))
withN_5mer_mut <- array(NA, length(withN_5mer))
for(mer in withN_5mer){
    withN_5mer_mut[mer]<- triMutability[substr(mer,2,4)]
}

triIndexes <- sapply(withN_5mer, function(x)substr(x,2,4))
withN_5mer_mut2 <- array(NA, length(withN_5mer))
names(withN_5mer_mut2) <- withN_5mer
for(mer in withN_5mer){
    merAsChar <- seqinr::s2c(mer)
    N_Index <- grep("[N-]",merAsChar)
    if(any(N_Index)){
        if(!any(N_Index%in%3)){
            trimerAsChar <- merAsChar[2:4]
            merAsChar[N_Index] <- "."
            merAsStr <- seqinr::c2s(merAsChar)
            merIndex <- grep(merAsStr,withN_5mer)
            withN_5mer_mut2[mer] <- mean(withN_5mer_mut[withN_5mer[merIndex]],na.rm=T)
        }
    }else{
        withN_5mer_mut2[mer] <- withN_5mer_mut[mer]
    }
}

# TODO redo normalization
m3n_mut <- withN_5mer_mut2
m3n_mut <- (m3n_mut/sum(m3n_mut,na.rm=TRUE)) * sum(!is.na(m3n_mut))

# Compute targeting
m3n_tar <- createTargetingMatrix(m3n_sub, m3n_mut)
#m3n_tar <- sweep(m3n_sub, 2, mutability,`*`)

# Build and save TargetingModel object
M3NModel <- new("TargetingModel",
                name="m3n",
                description="3-mer targeting model.",
                species="Mus musculus",
                date="2015-01-06",
                citation="Smith DS, et al. J Immunol. 1996 156:2642–52",
                substitution=m3n_sub,
                mutability=m3n_mut,
                targeting=m3n_tar )
devtools::use_data(M3NModel, overwrite=TRUE)


#### HS5F ####

# Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#   based on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#   Front Immunol. 2013 4(November):358.

# Load substitution data and convert to proper data structure
# load("data-raw/HS5F_Targeting.RData")
hs5f_sub_raw <- read.delim("data-raw/HS5F_Substitution.csv", sep=" ", as.is=TRUE)
hs5f_sub <- t(as.matrix(hs5f_sub_raw[c("A", "C" ,"G", "T")]))
colnames(hs5f_sub) <- hs5f_sub_raw$Fivemer
rownames(hs5f_sub) <- c("A", "C" ,"G", "T")

# Reassign no-change substitution entries to NA
center_nuc <- gsub("..([ACGT])..", "\\1", colnames(hs5f_sub))
for (i in 1:length(center_nuc)) {
    hs5f_sub[center_nuc[i], i] <- NA
}

# Load mutability data and convert to proper data structure
hs5f_mut_raw <- read.delim("data-raw/HS5F_Mutability.csv", sep=" ", as.is=TRUE)
hs5f_mut <- setNames(hs5f_mut_raw$Mutability, hs5f_mut_raw$Fivemer)

# Normalize mutability
hs5f_mut <- hs5f_mut / sum(hs5f_mut, na.rm=TRUE)

# Extend substitution and mutability
hs5f_sub <- extendSubstitutionMatrix(hs5f_sub)
hs5f_mut <- extendMutabilityMatrix(hs5f_mut)

# Compute targeting
hs5f_tar <- createTargetingMatrix(hs5f_sub, hs5f_mut)

# Build and save TargetingModel object
HS5FModel <- new("TargetingModel",
                 name="hs5f",
                 description="5-mer targeting model based on silent mutations in functional sequences.",
                 species="Homo sapiens",
                 date="2015-01-06",
                 citation="Yaari G, et al. Front Immunol. 2013 4(November):358",
                 substitution=hs5f_sub,
                 mutability=hs5f_mut,
                 targeting=hs5f_tar)
devtools::use_data(HS5FModel, overwrite=TRUE)

#### MRS5NF #####
# Unpublished.  Shlomchik anti-NP mouse Kappa chain data.

# load("data-raw/MRS5NF_Targeting.RData")

# MRS5NFModel <- new("TargetingModel",
#                    name="hs5f",
#                    description="5-mer targeting model based on all mutations in non-functional sequences.",
#                    species="Mus musculus",
#                    date="2015-01-06",
#                    citation="Unpublished",
#                    mutability=Targeting[["Mutability"]],
#                    substitution=Targeting[['Substitution']],
#                    targeting=Targeting[['Targeting']])
# devtools::use_data(MRS5NFModel)

#### M1N #####

# Smith DS, et al. Di- and trinucleotide target preferences of somatic mutagenesis 
#    in normal and autoreactive B cells. 
#    J Immunol. 1996 156:2642–52. 

NUCS <- c('A','C','G','T','N','.','-')
M1NDistance <- matrix(c(1, 2.86, 1, 2.14, 0, 0, 0, 2.86, 0, 2.14, 1, 0, 0, 0, 1, 
                        2.14, 0, 2.86, 0, 0, 0, 2.14, 1, 2.86, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                      7, 7, dimnames=list(NUCS,NUCS))
devtools::use_data(M1NDistance, overwrite=TRUE)