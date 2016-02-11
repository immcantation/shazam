# Imports
library(alakazam)
library(shazam)
library(profvis)

#### Load example data ####

file <- system.file("extdata", "InfluenzaDb.gz", package="shazam")
db <- readChangeoDb(file)
# Subset data for demo purposes
db <- subset(db, CPRIMER %in% c("IGHA","IGHG") & 
                 BARCODE != "RL013")

#### createSubstitutionMatrix ####

profvis({
    sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
})

#### createMutabilityMatrix ####

sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
profvis({
    mut_model <- createMutabilityMatrix(db, sub_model, model="S", multipleMutation="ignore",
                                        minNumSeqMutations=10)
})