# DistToNearest tests
# @author  Mohamed Uduman
# @date    2015.4.24

library("shm")

dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
db <- readChangeoDb(dbPath)

# Collapses sequences by clones
# Counts the number of observed mutations
# Calculates expected frequency of mutations
# Performs BASELIne on each sequence ( i.e. on the clones)
# Returns   baseline[[1]] = the updated db
#           baseline[[2]] = list of BASELINe objects
baseline <- getBASELINe ( db,
                          testStatistic="focused",
                          regionDefinition=IMGT_V_NO_CDR3,
                          nproc=3 )


# Group PDFs by the BARCODE column
# Returns   baseline_grp[[1]] = BASELINe summary stats
#           baseline_grp[[2]] = list of BASELINe objects (for further future groupings)
baseline_grp <- getBASELINe ( db = baseline[[1]],
                              groupBy = "BARCODE",
                              db_BASELINe = baseline[[2]],
                              nproc=3 )