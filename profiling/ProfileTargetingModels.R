# Imports
library(alakazam)
library(shazam)
library(profvis)

#### Load example data ####

# Subset data for demo purposes
db <- subset(ExampleDb, ISTOYPE %in% c("IGHA","IGHG") & SAMPLE == "+7d")

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

#### Plot mutability ####

profvis({
    # Plot one nucleotide in circular style
    plotMutability(HS5FModel, "C")
    
})

profvis({
    # Plot one nucleotide in barchart style
    plotMutability(HS5FModel, "T", style="bar")
})    