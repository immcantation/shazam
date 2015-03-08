# Loads and converts model data

# Smith DS, et al. Di- and trinucleotide target preferences of somatic mutagenesis 
#    in normal and autoreactive B cells. 
#    J Immunol. 1996 156:2642–52. 
#
load(system.file("extdata", "MTri_Targeting.RData", package="shm"))
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
load(system.file("extdata", "HS5F_Targeting.RData", package="shm"))
HS5FModel <- new("TargetingModel",
                 name="hs5f",
                 description="5-mer targeting model based on silent mutations in functional sequences.",
                 species="Homo sapiens",
                 date="2015-01-06",
                 citation="Yaari G, et al. Front Immunol. 2013 4(November):358",
                 mutability=Targeting[["Mutability"]],
                 substitution=Targeting[['Substitution']],
                 targeting=Targeting[['Targeting']])
devtools::use_data(HS5FModel)

# Unpublished.  Shlomchik anti-NP mouse Kappa chain data.
#
# load(system.file("extdata", "MRS5NF_Targeting.RData", package="shm"))
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