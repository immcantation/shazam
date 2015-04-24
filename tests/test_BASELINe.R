# BASELINe tests
# @author  Mohamed Uduman
# @date    2015.04.23

library(shm)
dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
db <- readChangeoDb(dbPath)

baseline <- getBASELINePDF( db,
                            testStatistic="local",
                            regionDefinition=IMGT_V_NO_CDR3,
                            nproc=3 )

getBASELINeStats(baseline, nproc=3)


db_new <- getBASELINeStats(db=baseline[[1]],
                           list_BASELINe=baseline[[2]],
                           nproc=3)

baseline_grp <- getBASELINePDF(db=baseline[[1]],
                               list_BASELINe=baseline[[2]],
                               group="BARCODE",
                               testStatistic="local",
                               regionDefinition=IMGT_V_NO_CDR3,
                               nproc=3)

db_new_grp <- getBASELINeStats(db=baseline[[1]],
                               list_BASELINe=baseline[[2]],
                               nproc=3)