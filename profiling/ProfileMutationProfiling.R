# Imports
library(alakazam)
library(shazam)
library(profvis)

#### Load example data ####

# Subset data for demo purposes
db <- subset(ExampleDb, c_call == "IgG" & sample_id == "+7d")

#### observedMutations ####

profvis({
    mutations <- observedMutations(db, regionDefinition=NULL, frequency=T)
})


profvis({
    slideWindowDb(db=ExampleDb, sequenceColumn="sequence_alignment",
                  germlineColumn="germline_alignment_d_mask",
                  mutThresh=6, windowSize=10, nproc=2)
})
