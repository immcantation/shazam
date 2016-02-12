# Imports
library(alakazam)
library(shazam)
library(profvis)

#### Load example data ####

file <- system.file("extdata", "InfluenzaDb.gz", package="shazam")
db <- readChangeoDb(file)
# Subset data for demo purposes
db <- subset(db, CPRIMER %in% c("IGHA","IGHM") & 
                 BARCODE != "RL013")
# Extracting the first entry in the sample db to use for input and germline sequences.
inputSeq <- db[1, "SEQUENCE_IMGT"]
germlineSeq <-  db[1, "GERMLINE_IMGT_D_MASK"]

#### Calculate expected mutations ####

profvis({
    calcExpectedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3,
                          mutationDefinition=HYDROPATHY_MUTATIONS)
})

#### Calculate hydropathy expected mutations over V region ####

profvis({
    db <- calcDBExpectedMutations(db,
                                  sequenceColumn="SEQUENCE_IMGT",
                                  germlineColumn="GERMLINE_IMGT_D_MASK",
                                  regionDefinition=IMGT_V_NO_CDR3,
                                  mutationDefinition=HYDROPATHY_MUTATIONS,
                                  nproc=1)
})

#### Collapse one clone ####

profvis({
    calcClonalConsensus(inputSeq, germlineSeq, 
                        regionDefinition=IMGT_V_NO_CDR3, nonTerminalOnly=FALSE)
})

#### Collapse clones ####

profvis({
    db_new <- collapseByClone(db, cloneColumn="CLONE", 
                              sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              expandedDb=FALSE,
                              regionDefinition=IMGT_V_NO_CDR3,
                              nproc=1)
})

#### Calculate Baseline ####

profvis({
    db_baseline <- calcBaseline(db, 
                                sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK", 
                                testStatistic="focused",
                                regionDefinition=IMGT_V_NO_CDR3,
                                targetingModel = HS5FModel,
                                nproc = 1)
})

#### Group Baseline ####

profvis({
    baseline_one <- groupBaseline(db_baseline, groupBy="BARCODE")
})

#### Plot Baseline density example ####

profvis({
    db_baseline <- calcBaseline(db, 
                                sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK", 
                                testStatistic="focused",
                                regionDefinition=IMGT_V_NO_CDR3,
                                targetingModel = HS5FModel,
                                nproc = 1)
     
    # Grouping the PDFs by the BARCODE and CPRIMER columns in the db, corresponding 
    # respectively to sample barcodes and the constant region isotype primers.
    baseline <- groupBaseline(db_baseline, groupBy=c("BARCODE", "CPRIMER"))
    
    # Plot mean and confidence interval
    plotBaselineDensity(baseline, "BARCODE", "CPRIMER", style="density")
    plotBaselineDensity(baseline, "BARCODE", "CPRIMER", subsetRegions="CDR", style="density")
    plotBaselineDensity(baseline, "BARCODE", "CPRIMER", facetBy="group", style="density")
})
