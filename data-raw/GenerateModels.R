# Loads and converts model data

# Imports
library(seqinr)

#### M1NDistance #####

# Smith DS, et al. Di- and trinucleotide target preferences of somatic mutagenesis 
#    in normal and autoreactive B cells. 
#    J Immunol. 1996 156:2642-52. 

nuc_chars <- c('A','C','G','T','N','.','-')
M1NDistance <- matrix(c(0, 2.86, 1, 2.14, 0, 0, 0, 
                        2.86, 0, 2.14, 1, 0, 0, 0, 
                        1, 2.14, 0, 2.86, 0, 0, 0, 
                        2.14, 1, 2.86, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0),
                      7, 7, dimnames=list(nuc_chars, nuc_chars))
devtools::use_data(M1NDistance, overwrite=TRUE)


#### HS1FDistance #####

# Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#   based on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#   Front Immunol. 2013 4(November):358.

nuc_chars <- c('A','C','G','T','N','.','-')
HS1FDistance <- matrix(c(0, 2.08, 1, 1.75, 0, 0, 0, 
                        2.08, 0, 1.75, 1, 0, 0, 0, 
                        1, 1.75, 0, 2.08, 0, 0, 0, 
                        1.75, 1, 2.08, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0),
                      7, 7, dimnames=list(nuc_chars, nuc_chars))
devtools::use_data(HS1FDistance, overwrite=TRUE)


#### U5N ####

# 5-mer null model

# Define substitution matrix
nuc_chars <- c("A", "C" ,"G", "T")
nuc_words <- seqinr::words(5, nuc_chars)
u5n_sub <- matrix(1/3, nrow=length(nuc_chars), ncol=length(nuc_words), 
                  dimnames=list(nuc_chars, nuc_words))

# Reassign no-change substitution entries to NA
center_nuc <- gsub("..([ACGT])..", "\\1", colnames(u5n_sub))
for (i in 1:length(center_nuc)) {
    u5n_sub[center_nuc[i], i] <- NA
}

# Define mutability matrix
u5n_mut <- array(1/length(nuc_words), dim=length(nuc_words), dimnames=list(nuc_words))

# Extend substitution and mutability with N characters
u5n_sub <- extendSubstitutionMatrix(u5n_sub)
u5n_mut <- extendMutabilityMatrix(u5n_mut)

# Compute targeting
u5n_tar <- createTargetingMatrix(u5n_sub, u5n_mut)

# Build and save TargetingModel object
U5N <- new("TargetingModel",
                 name="U5N",
                 description="uniform 5-mer null model",
                 species="",
                 date="2015-05-06",
                 citation="",
                 substitution=u5n_sub,
                 mutability=u5n_mut,
                 targeting=u5n_tar)
devtools::use_data(U5N, overwrite=TRUE)

#### HH_S1F ####

# Yaari G, et al. Models of somatic hypermutation targeting and substitution 

# HS1FDistance[1:4, 1:4] normalized by row
# HH_S1F = round(apply(HS1FDistance[1:4, 1:4], 1, function(x){x/sum(x)}), 3)

# Hard-coded in case HS1FDistance gets renamed in the future
HH_S1F = matrix(data=c(0.000, 0.431, 0.207, 0.362,
                       0.431, 0.000, 0.362, 0.207,
                       0.207, 0.362, 0.000, 0.431,
                       0.362, 0.207, 0.431, 0.000),
                nrow=4, byrow=TRUE, 
                dimnames = list(c("A","C","G","T"), 
                                c("A","C","G","T")))

devtools::use_data(HH_S1F, overwrite=TRUE)

#### HKL_S1F ####

# Table 3 from 
# Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  
# Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
# Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
# Immunology,197(9), 3566–3574. http://doi.org/10.4049/jimmunol.1502263

HKL_S1F = matrix(data=c(0.00, 0.26, 0.50, 0.24,
                        0.25, 0.00, 0.30, 0.45,
                        0.53, 0.31, 0.00, 0.16,
                        0.20, 0.57, 0.23, 0.00),
                  nrow=4, byrow=TRUE, 
                  dimnames = list(c("A","C","G","T"), 
                                  c("A","C","G","T")))

devtools::use_data(HKL_S1F, overwrite=TRUE)

#### MK_RS1NF ####

# Directly from Ang Cui
MK_RS1NF = matrix(data=c(0.00, 0.17, 0.53, 0.30,
                         0.10, 0.00, 0.13, 0.77,
                         0.77, 0.14, 0.00, 0.08,
                         0.28, 0.53, 0.19, 0.00),
                nrow=4, byrow=TRUE, 
                dimnames = list(c("A","C","G","T"), 
                                c("A","C","G","T")))

devtools::use_data(MK_RS1NF, overwrite=TRUE)

#### HH_S5F ####

# Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#   based on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#   Front Immunol. 2013 4(November):358.

# Load substitution data and convert to proper data structure
hs5f_sub_raw <- read.csv("data-raw/HS5F_Substitution.csv", as.is=TRUE)
hs5f_sub <- t(as.matrix(hs5f_sub_raw[c("A", "C" ,"G", "T")]))
colnames(hs5f_sub) <- hs5f_sub_raw$Fivemer
rownames(hs5f_sub) <- c("A", "C" ,"G", "T")

# Reassign no-change substitution entries to NA
center_nuc <- gsub("..([ACGT])..", "\\1", colnames(hs5f_sub))
for (i in 1:length(center_nuc)) {
    hs5f_sub[center_nuc[i], i] <- NA
}

# Load mutability data and convert to proper data structure
hs5f_mut_raw <- read.csv("data-raw/HS5F_Mutability.csv", as.is=TRUE)
hs5f_mut <- setNames(hs5f_mut_raw$Mutability, hs5f_mut_raw$Fivemer)

# Normalize mutability
hs5f_mut <- hs5f_mut / sum(hs5f_mut, na.rm=TRUE)

# Extend substitution and mutability with N characters
hs5f_sub <- extendSubstitutionMatrix(hs5f_sub)
hs5f_mut <- extendMutabilityMatrix(hs5f_mut)

# Compute targeting
hs5f_tar <- createTargetingMatrix(hs5f_sub, hs5f_mut)

# Build and save TargetingModel object
HH_S5F <- new("TargetingModel",
                 name="HH_S5F",
                 description="5-mer targeting model based on silent (S) mutations in human heavy-chain functional sequences.",
                 species="Homo sapiens",
                 date="2015-01-06",
                 citation="Yaari G, et al. Front Immunol. 2013 4(November):358",
                 substitution=hs5f_sub,
                 mutability=hs5f_mut,
                 targeting=hs5f_tar)
devtools::use_data(HH_S5F, overwrite=TRUE)

#### HKL_S5F ####

# Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  
# Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
# Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
# Immunology,197(9), 3566–3574. http://doi.org/10.4049/jimmunol.1502263

load("data-raw/HKL_S5F_raw.rdata")

HKL_S5F <- new("TargetingModel",
                name = "HKL_S5F",
                description = "5-mer targeting model based on silent (S) mutations in human kappa & lambda light-chain functional sequences",
                species = "Homo sapiens",
                date = "2015-12-06",
                citation = "Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of Immunology, 197(9), 3566–3574. http://doi.org/10.4049/jimmunol.1502263",
                mutability = hL.mut,
                substitution = hL.sub,
                targeting = hL.tar)

devtools::use_data(HKL_S5F, overwrite=TRUE)

rm(hL.sub, hL.tar, hL.char, hL.mut)

#### MK_RS5NF #####

# Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  
# Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
# Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
# Immunology,197(9), 3566–3574. http://doi.org/10.4049/jimmunol.1502263

load("data-raw/MK_RS5NF_raw.rdata")

MK_RS5NF <- new("TargetingModel",
               name = "MK_RS5NF",
               description = "5-mer targeting model based on replacement (R) and silent (S) mutations in kappa light-chain non-functional sequences from NP-immunized mice",
               species = "Mus musculus",
               date = "2015-01-06",
               citation = "Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of Immunology, 197(9), 3566–3574. http://doi.org/10.4049/jimmunol.1502263",
               mutability = mL.mut,
               substitution = mL.sub,
               targeting = mL.tar)

devtools::use_data(MK_RS5NF, overwrite=TRUE)

rm(mL.sub, mL.tar, mL.char, mL.mut)
