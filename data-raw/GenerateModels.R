# Loads and converts model data

# TODO:  fix normalization

# Smith DS, et al. Di- and trinucleotide target preferences of somatic mutagenesis 
#    in normal and autoreactive B cells. 
#    J Immunol. 1996 156:2642–52. 
#
load("data-raw/MTri_Targeting.RData")
M3NModel <- new("TargetingModel",
                name="m3n",
                description="3-mer targeting model.",
                species="Mus musculus",
                date="2015-01-06",
                citation="Smith DS, et al. J Immunol. 1996 156:2642–52",
                mutability=Targeting[["Mutability"]],
                substitution=Targeting[['Substitution']],
                targeting=Targeting[['Targeting']])
devtools::use_data(M3NModel)

# Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#   based on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#   Front Immunol. 2013 4(November):358.
#
# load("data-raw/HS5F_Targeting.RData")
# Load substitution data and convert to proper data structure
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

# Unpublished.  Shlomchik anti-NP mouse Kappa chain data.
#
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