# Loads and converts model data

# Imports
library(seqinr)

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

# from G Yaari (transposed his origianl "Normalized" matrix, 
# which had columns summing up to 1)
HH_S1F = matrix(data=c(0.0000, 0.3057, 0.4786, 0.2157,
                       0.2443, 0.0000, 0.3031, 0.4526,
                       0.5075, 0.3249, 0.0000, 0.1676,
                       0.2219, 0.4945, 0.2836, 0.0000),
                nrow=4, byrow=TRUE, 
                dimnames = list(c("A","C","G","T"), 
                                c("A","C","G","T")))

devtools::use_data(HH_S1F, overwrite=TRUE)

#### HKL_S1F ####

# Table 3 from 
# Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  
# Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
# Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
# Immunology,197(9), 3566-3574. http://doi.org/10.4049/jimmunol.1502263

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
# Immunology,197(9), 3566-3574. http://doi.org/10.4049/jimmunol.1502263

load("data-raw/HKL_S5F_raw.rdata")

HKL_S5F <- new("TargetingModel",
                name = "HKL_S5F",
                description = "5-mer targeting model based on silent (S) mutations in human kappa & lambda light-chain functional sequences",
                species = "Homo sapiens",
                date = "2015-12-06",
                citation = "Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of Immunology, 197(9), 3566-3574. http://doi.org/10.4049/jimmunol.1502263",
                mutability = hL.mut,
                substitution = hL.sub,
                targeting = hL.tar)

devtools::use_data(HKL_S5F, overwrite=TRUE)

rm(hL.sub, hL.tar, hL.char, hL.mut)

#### MK_RS5NF #####

# Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  
# Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
# Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
# Immunology,197(9), 3566-3574. http://doi.org/10.4049/jimmunol.1502263

load("data-raw/MK_RS5NF_raw.rdata")

MK_RS5NF <- new("TargetingModel",
               name = "MK_RS5NF",
               description = "5-mer targeting model based on replacement (R) and silent (S) mutations in kappa light-chain non-functional sequences from NP-immunized mice",
               species = "Mus musculus",
               date = "2015-01-06",
               citation = "Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,  Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of Immunology, 197(9), 3566-3574. http://doi.org/10.4049/jimmunol.1502263",
               mutability = mL.mut,
               substitution = mL.sub,
               targeting = mL.tar)

devtools::use_data(MK_RS5NF, overwrite=TRUE)

rm(mL.sub, mL.tar, mL.char, mL.mut)
