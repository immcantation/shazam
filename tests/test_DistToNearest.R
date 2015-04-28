# DistToNearest tests
# @author  Namita Gupta
# @date    2014.9.25

#### Imports ####
library(alakazam)
library(shm)

#### Load data ####
db_file1 <- system.file("extdata", "Influenza.tab", package="shm")
data_df1 <- readChangeoDb(db_file1)

#### DistToNearest steps ####
data_df1 <- distToNearest(data_df1, genotyped=T, first=F)
hist(data_df1$DIST_NEAREST, breaks=50, xlim=c(0,0.02))

