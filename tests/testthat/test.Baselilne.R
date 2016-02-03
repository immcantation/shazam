dbPath <- system.file("extdata", "InfluenzaDb.gz", package="shm")
db <- readChangeoDb(dbPath)

test_that("calculateTargeting and calculateMutationalPaths with regionDefinition==NULL", {
    
    ##--
    ## With regionDefinition = IMGT_V_NO_CDR3
    ##--
    ## calculateTargeting and calculateMutationalPaths 
    ## as used by calcExpectedMutations
    
    inputSeq <- db[1, "SEQUENCE_IMGT"]
    germlineSeq <-  db[1, "GERMLINE_IMGT_D_MASK"]
    
    targeting <- shm:::calculateTargeting(germlineSeq = germlineSeq, 
                                    inputSeq = inputSeq,
                                    targetingModel = HS5FModel,
                                    regionDefinition = IMGT_V_NO_CDR3)
    expect_equal(dim(targeting),c(5,IMGT_V_NO_CDR3@seqLength))
    
    obs_targeting <- targeting[3,c(1,10,50,78,312)]
    exp_targeting <- c(NA, 0.0011048048, NA, 0.0001235885, 0.0002081728)
    names(exp_targeting) <- c("G", "C", "G", "T", "T")
    expect_equal(obs_targeting, exp_targeting, tolerance=0.001)
        
    mutationalPaths <- shm:::calculateMutationalPaths(
                            germlineSeq = seqinr::c2s(colnames(targeting)), 
                            regionDefinition = IMGT_V_NO_CDR3)
    
    obs_mpath <- mutationalPaths[3,c(1,10,50,78,312)]
    exp_mpath <- c("S", "R", "S", "S", "R" )
    names(exp_mpath) <- c("G", "C", "G", "T", "T")
    expect_equal(obs_mpath, exp_mpath)
    
    ##--
    ## With regionDefinition = NULL
    ##--
    ## calculateTargeting and calculateMutationalPaths 
    ## as used by calcExpectedMutations
    
    targeting_null <- shm:::calculateTargeting(germlineSeq = germlineSeq, 
                                          inputSeq = inputSeq,
                                          targetingModel = HS5FModel,
                                          regionDefinition = NULL)
    expect_equal(dim(targeting_null),c(5,406))
    
    obs_targeting_null <- targeting_null[3,c(1,10,50,78,312)]
    exp_targeting_null <- c(NA, 0.0011048048, NA, 0.0001235885, 0.0001512495)
    names(exp_targeting_null) <- c("G", "C", "G", "T", "T")
    expect_equal(obs_targeting_null, exp_targeting_null, tolerance=0.001)
    
    mutationalPaths_null <- shm:::calculateMutationalPaths(
        germlineSeq = seqinr::c2s(colnames(targeting_null)), 
        regionDefinition = NULL)
    
    obs_mpath_null <- mutationalPaths_null[3,c(1,10,50,78,312)]
    exp_mpath_null <- c("S", "R", "S", "S", "R" )
    names(exp_mpath_null) <- c("G", "C", "G", "T", "T")
    expect_equal(obs_mpath_null, exp_mpath_null)
})


##--
## Expected
##--

test_that("calcExpectedMutations works with regionDefinition==NULL",{
    inputSeq <- db[1, "SEQUENCE_IMGT"]
    germlineSeq <-  db[1, "GERMLINE_IMGT_D_MASK"]
    
    obs_mutations <- calcExpectedMutations(inputSeq, germlineSeq, targetingModel = HS5FModel, 
                          regionDefinition = IMGT_V_NO_CDR3)
    exp_mutations <- c(0.16184071, 0.04872069, 0.57766105, 0.21177756)
    names(exp_mutations) <- c("CDR_R", "CDR_S", "FWR_R", "FWR_S")
    expect_equal(obs_mutations, exp_mutations, tolerance=0.001)
    
    obs_mutations_null <- calcExpectedMutations(inputSeq, germlineSeq, targetingModel = HS5FModel, 
                                           regionDefinition = NULL)
    exp_mutations_null <- c(0.7379522, 0.2620478)
    names(exp_mutations_null) <- c("SEQ_R", "SEQ_S")
    expect_equal(obs_mutations_null, exp_mutations_null, tolerance=0.001)
})

test_that("calcDBExpectedMutations works with regionDefinition==NULL",{
    db_subset <- subset(db, CPRIMER %in% c("IGHA","IGHM") & 
                     BARCODE %in% c("RL016","RL018","RL019","RL021"))
    db_mutations <- calcDBExpectedMutations( db_subset,
                             sequenceColumn="SEQUENCE_IMGT",
                             germlineColumn="GERMLINE_IMGT_D_MASK",
                             regionDefinition=IMGT_V_NO_CDR3,
                             nproc=1)
    
    ## Check 5 examples for each, at different positions
    ## CDR_R, first 5
    obs_cdr_r <- db_mutations$EXPECTED_CDR_R[1:5]
    exp_cdr_r <- c(0.2072431, 0.2054472, 0.1870927, 0.1507470, 0.2200854)
    expect_equal(obs_cdr_r, exp_cdr_r, tolerance=0.001)
    
    ## CDR_S 311:315
    obs_cdr_s <- db_mutations$EXPECTED_CDR_S[311:315]
    exp_cdr_s <- c(0.03682467, 0.03682467, 0.03682467, 0.03682467, 0.03682467)
    expect_equal(obs_cdr_r, exp_cdr_r, tolerance=0.001)
    
    ## FWR_R, 120:124
    obs_fwr_r <- db_mutations$EXPECTED_FWR_R[120:124]
    exp_fwr_r<- c(0.5883452, 0.5883452, 0.5883452, 0.5883452, 0.5883452)
    expect_equal(obs_fwr_r, exp_fwr_r, tolerance=0.001)
    
    ## FWR_S, 993:997
    obs_fwr_s <- db_mutations$EXPECTED_FWR_S[207:211]
    exp_fwr_s<- c(0.1950189, 0.1950189, 0.1928329, 0.1928329, 0.1950189)
    expect_equal(obs_fwr_s, exp_fwr_s, tolerance=0.001)
    
    db_mutations_null <- calcDBExpectedMutations( db_subset,
                                             sequenceColumn="SEQUENCE_IMGT",
                                             germlineColumn="GERMLINE_IMGT_D_MASK",
                                             regionDefinition=NULL,
                                             nproc=1)
    ## SEQ_R, first 5
    obs_seq_r <- db_mutations_null$EXPECTED_SEQ_R[1:5]
    exp_seq_r <- c(0.7590282, 0.7635794, 0.7611897, 0.7585786, 0.7761334)
    expect_equal(obs_seq_r, exp_seq_r, tolerance=0.001)
    
    ## SEQ_S 311:315
    obs_seq_s <- db_mutations_null$EXPECTED_SEQ_S[311:315]
    exp_seq_s <- c(0.2383252, 0.2361719, 0.2361354, 0.2400638, 0.2389921)
    expect_equal(obs_seq_s, exp_seq_s, tolerance=0.001)
    
})

##--
## Observed
##--

test_that("calcObservedMutations works with regionDefinition==NULL",{
    inputSeq <- db[1, "SEQUENCE_IMGT"]
    germlineSeq <-  db[1, "GERMLINE_IMGT_D_MASK"]
    
    obs_mutations <- calcObservedMutations(inputSeq, germlineSeq, 
                                           regionDefinition=IMGT_V_NO_CDR3,
                                           frequency=F)
    
    ## TODO: should sum(exp_mutations) and sum(exp_mutations_null) match?
    exp_mutations <- c(2, 2, 8, 1)
    names(exp_mutations) <- c("CDR_R", "CDR_S", "FWR_R", "FWR_S")
    expect_equal(obs_mutations, exp_mutations, tolerance=0.001)
    
    obs_mutations_null <- calcObservedMutations(inputSeq, germlineSeq, 
                            regionDefinition=NULL, frequency=F)
    exp_mutations_null <- c(12, 5)
    names(exp_mutations_null) <- c("SEQ_R", "SEQ_S")
    expect_equal(obs_mutations_null, exp_mutations_null, tolerance=0.001)
})

test_that("calcDBObservedMutations works with regionDefinition==NULL",{
    db_subset <- subset(db, CPRIMER %in% c("IGHA","IGHM") & 
                            BARCODE %in% c("RL016","RL018","RL019","RL021"))
    db_mutations <- calcDBObservedMutations( db_subset,
                                             sequenceColumn="SEQUENCE_IMGT",
                                             germlineColumn="GERMLINE_IMGT_D_MASK",
                                             regionDefinition=IMGT_V_NO_CDR3,
                                             nproc=1)
    ## Check 5 examples for each
    ## CDR_R, first 5
    obs_cdr_r <- db_mutations$OBSERVED_CDR_R[1:5]
    exp_cdr_r <- c(2, 0, 0, 0, 0)
    expect_equal(obs_cdr_r, exp_cdr_r)
    
    ## CDR_S 311:315
    obs_cdr_s <- db_mutations$OBSERVED_CDR_S[311:315]
    exp_cdr_s <- c(0, 0, 0, 0, 0)
    expect_equal(obs_cdr_r, exp_cdr_r)
    
    ## FWR_R, 120:124
    obs_fwr_r <- db_mutations$OBSERVED_FWR_R[120:124]
    exp_fwr_r<- c(0, 0, 0, 0, 0)
    expect_equal(obs_fwr_r, exp_fwr_r)
    
    ## FWR_S, 40:44
    obs_fwr_s <- db_mutations$OBSERVED_FWR_S[40:44]
    exp_fwr_s<- c(0, 0, 0, 1, 0)
    expect_equal(obs_fwr_s, exp_fwr_s)
     
    db_mutations_null <- calcDBObservedMutations( db_subset,
                                                  sequenceColumn="SEQUENCE_IMGT",
                                                  germlineColumn="GERMLINE_IMGT_D_MASK",
                                                  regionDefinition=NULL,
                                                  nproc=1)
    ## SEQ_R, first 5
    obs_seq_r <- db_mutations_null$OBSERVED_SEQ_R[1:5]
    exp_seq_r <- c(12, 0, 0, 0, 0)
    expect_equal(obs_seq_r, exp_seq_r)
    
    ## SEQ_S 38:42
    obs_seq_s <- db_mutations_null$OBSERVED_SEQ_S[38:42]
    exp_seq_s <- c(0, 2, 2, 0, 0)
    expect_equal(obs_seq_s, exp_seq_s)
    
})

##--
## Baseline
##--

test_that("calcBaseline works with regionDefinition==NULL", {
    set.seed(82)
    db_baseline <- calcBaseline(db, 
                                sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK", 
                                testStatistic="focused",
                                regionDefinition=IMGT_V_NO_CDR3,
                                targetingModel = HS5FModel,
                                nproc = 1,
                                calcStats = T)
    
    ## Check if the stats slot has been filled
    expect_gt(nrow(slot(db_baseline,"stats")),0)
    
    ## Check 5 examples for each, at different positions
    ## CDR_R, first 5
    obs_cdr_r <- slot(db_baseline,"db")$OBSERVED_CDR_R[1:5]
    exp_cdr_r<- c(2, 6, 5, 0, 0)
    expect_equal(obs_cdr_r, exp_cdr_r)
    
    ## CDR_S, 673:677
    obs_cdr_s <- slot(db_baseline,"db")$OBSERVED_CDR_S[673:677]
    exp_cdr_s <- c(1,  2,  2,  1,  3)
    expect_equal(obs_cdr_s, exp_cdr_s)
    
    ## FWR_R, 937:941
    obs_fwr_r <- slot(db_baseline,"db")$OBSERVED_FWR_R[937:941]
    exp_fwr_r<- c(0, 7, 14, 7, 0)
    expect_equal(obs_fwr_r, exp_fwr_r)
    
    ## FWR_S, 993:997
    obs_fwr_s <- slot(db_baseline,"db")$OBSERVED_FWR_S[993:997]
    exp_fwr_s<- c(10, 0, 0, 0, 0)
    expect_equal(obs_fwr_s, exp_fwr_s)
    
#     db_baseline_null <- calcBaseline(db, 
#                                 sequenceColumn="SEQUENCE_IMGT",
#                                 germlineColumn="GERMLINE_IMGT_D_MASK", 
#                                 testStatistic="focused",
#                                 regionDefinition=NULL,
#                                 targetingModel = HS5FModel,
#                                 nproc = 1)
})