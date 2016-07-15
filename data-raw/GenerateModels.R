# Loads and converts model data

# Imports
library(seqinr)

#### M1N #####

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


#### HS1F #####

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
U5NModel <- new("TargetingModel",
                 name="u5n",
                 description="uniform 5-mer model",
                 species="",
                 date="2015-05-06",
                 citation="",
                 substitution=u5n_sub,
                 mutability=u5n_mut,
                 targeting=u5n_tar)
devtools::use_data(U5NModel, overwrite=TRUE)


#### HS5F ####

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

load("data-raw/db_gc_igk_nf_rs_targeting.rdata")

# extend substitution matrix; get 5x3125 (expect 5th row [N] to be NAs)
mrs5nf_sub <- extendSubstitutionMatrix(db_gc_igk_nf_sub_matrix_rs)
#table(colSums(mrs5nf_sub, na.rm=1)) # 625 0/NAs (expected)

# extend mutability vector
# db_gc_igk_nf_mut_matrix_rs is a df
# extendMutabilityMatrix takes a numeric vector with names of 5mers
mrs5nf_mut_vec <- db_gc_igk_nf_mut_matrix_rs$Mutability
names(mrs5nf_mut_vec) <- rownames(db_gc_igk_nf_mut_matrix_rs)
#sum(mrs5nf_mut_vec) # 1 (expected)
mrs5nf_mut <- extendMutabilityMatrix(mrs5nf_mut_vec)
#summary(mrs5nf_mut) # 625 NAs (expected)

# compute targeting; get 5x3125 (expect 5th row [N] to be NAs)
mrs5nf_tar <- createTargetingMatrix(mrs5nf_sub, mrs5nf_mut)

MRS5NFModel <- new("TargetingModel",
                   name = "rs5nf",
                   description = "5-mer targeting model based on replacement (R) and silent (S) mutations in non-functional kappa light-chain sequences of NP-immunized mice",
                   species = "Mus musculus",
                   date = "2015-01-06",
                   citation = "Cui A, et al. A model of somatic hypermutation targeting in mice based on high-throughput immunoglobulin sequencing data. (In prepration)",
                   mutability = mrs5nf_mut,
                   substitution = mrs5nf_sub,
                   targeting = mrs5nf_tar)

devtools::use_data(MRS5NFModel, overwrite=TRUE)
