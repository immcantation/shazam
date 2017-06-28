# Imports
library(alakazam)
library(shazam)
library(profvis)

#### Load example data ####

# Subset data for demo purposes
db <- subset(ExampleDb, ISOTYPE == "IgG" & SAMPLE == "+7d")

#### observedMutations ####

profvis({
    mutations <- observedMutations(db, regionDefinition=NULL, frequency=T)
})
