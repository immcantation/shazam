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
# Returns   Baseline object
baseline <- calcBaselinePdfs( db,
                              testStatistic="focused",
                              regionDefinition=IMGT_V_NO_CDR3,
                              nproc=3 )


# Group PDFs by the BARCODE column
baseline_grp <- groupBaseline( baseline, 
                               groupBy=c("BARCODE","VPRIMER"), 
                               nproc=1)

