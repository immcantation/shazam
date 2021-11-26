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


profvis({
    # old 21220
    # nproc=1 17990
    # nproc=2 13580
    # after adding second loop foreach
    # nproc=1 940
    # nproc=2 1290
    slideWindowTune(ExampleDb, mutThreshRange=2:4, windowSizeRange=7:9, nproc=2)
})


