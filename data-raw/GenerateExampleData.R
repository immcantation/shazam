# Generate example data files

#### Imports ####

library(alakazam)

#### Influenza database ####

InfluenzaDb <- readChangeoDb("data-raw/InfluenzaDb.gz")

# Save data
devtools::use_data(InfluenzaDb,
                   overwrite=TRUE)