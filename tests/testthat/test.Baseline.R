# Load test database
e1 <- new.env()
#load(file.path("tests", "data-tests", "TestDb.rda"), envir=e1)
load(file.path("..", "data-tests", "TestDb.rda"), envir=e1)
db <- get("TestDb", envir=e1)
rm(e1)

load(file.path("..", "data-tests", "ExampleDb.rda"))

test_that("calculateTargeting and calculateMutationalPaths with regionDefinition==NULL", {
    
    ##--
    ## With regionDefinition = IMGT_V
    ##--
    ## calculateTargeting and calculateMutationalPaths 
    ## as used by calcExpectedMutations
    
    inputSeq <- db[["SEQUENCE_IMGT"]][1]
    germlineSeq <-  db[["GERMLINE_IMGT_D_MASK"]][1]
    
    targeting <- shazam:::calculateTargeting(germlineSeq = germlineSeq, 
                                    inputSeq = inputSeq,
                                    targetingModel = HH_S5F,
                                    regionDefinition = IMGT_V)
    expect_equal(dim(targeting),c(5,IMGT_V@seqLength))
    
    obs_targeting <- targeting[3,c(1,10,50,78,312)]
    exp_targeting <- c(NA, 0.0011048048, NA, 0.0001235885, 0.0002081728)
    names(exp_targeting) <- c("G", "C", "G", "T", "T")
    expect_equal(obs_targeting, exp_targeting, tolerance=0.001)
        
    mutationalPaths <- shazam:::calculateMutationalPaths(
                            germlineSeq = seqinr::c2s(colnames(targeting)), 
                            regionDefinition = IMGT_V)
    
    obs_mpath <- mutationalPaths[3,c(1,10,50,78,312)]
    exp_mpath <- c("S", "R", "S", "S", "R" )
    names(exp_mpath) <- c("G", "C", "G", "T", "T")
    expect_equal(obs_mpath, exp_mpath)
    
    ##--
    ## With regionDefinition = NULL
    ##--
    ## calculateTargeting and calculateMutationalPaths 
    ## as used by calcExpectedMutations
    
    targeting_null <- shazam:::calculateTargeting(germlineSeq = germlineSeq, 
                                          inputSeq = inputSeq,
                                          targetingModel = HH_S5F,
                                          regionDefinition = NULL)
    expect_equal(dim(targeting_null),c(5,406))
    
    obs_targeting_null <- targeting_null[3,c(1,10,50,78,312)]
    exp_targeting_null <- c(NA, 0.0011048048, NA, 0.0001235885, 0.0001512495)
    names(exp_targeting_null) <- c("G", "C", "G", "T", "T")
    expect_equal(obs_targeting_null, exp_targeting_null, tolerance=0.001)
    
    mutationalPaths_null <- shazam:::calculateMutationalPaths(
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
    inputSeq <- db[["SEQUENCE_IMGT"]][1]
    germlineSeq <-  db[["GERMLINE_IMGT_D_MASK"]][1]
    
    obs_mutations <- calcExpectedMutations(inputSeq, germlineSeq, targetingModel = HH_S5F, 
                          regionDefinition = IMGT_V)
    exp_mutations <- c(0.16184071, 0.04872069, 0.57766105, 0.21177756)
    names(exp_mutations) <- c("CDR_R", "CDR_S", "FWR_R", "FWR_S")
    expect_equal(obs_mutations, exp_mutations, tolerance=0.001)
    
    obs_mutations_null <- calcExpectedMutations(inputSeq, germlineSeq, targetingModel = HH_S5F, 
                                           regionDefinition = NULL)
    exp_mutations_null <- c(0.7379522, 0.2620478)
    names(exp_mutations_null) <- c("SEQ_R", "SEQ_S")
    expect_equal(obs_mutations_null, exp_mutations_null, tolerance=0.001)
})

test_that("expectedMutations works with regionDefinition==NULL",{
    db_subset <- subset(db, CPRIMER %in% c("IGHA","IGHM") & 
                        BARCODE %in% c("RL016","RL018","RL019","RL021"))
    db_mutations <- expectedMutations(db_subset,
                             sequenceColumn="SEQUENCE_IMGT",
                             germlineColumn="GERMLINE_IMGT_D_MASK",
                             regionDefinition=IMGT_V,
                             nproc=1)
    db_mutations_df <- expectedMutations(data.frame(db_subset),
                                      sequenceColumn="SEQUENCE_IMGT",
                                      germlineColumn="GERMLINE_IMGT_D_MASK",
                                      regionDefinition=IMGT_V,
                                      nproc=1)
    expect_identical(db_mutations, db_mutations_df)
    
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
    
    db_mutations_null <- expectedMutations(db_subset,
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
                                           regionDefinition=IMGT_V,
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

test_that("observedMutations works with regionDefinition==NULL",{
    db_subset <- subset(db, CPRIMER %in% c("IGHA","IGHM") & 
                            BARCODE %in% c("RL016","RL018","RL019","RL021"))
    db_mutations <- observedMutations(db_subset,
                                      sequenceColumn="SEQUENCE_IMGT",
                                      germlineColumn="GERMLINE_IMGT_D_MASK",
                                      regionDefinition=IMGT_V,
                                      nproc=1)
    
    db_mutations_df <- observedMutations(data.frame(db_subset),
                                         sequenceColumn="SEQUENCE_IMGT",
                                         germlineColumn="GERMLINE_IMGT_D_MASK",
                                         regionDefinition=IMGT_V,
                                         nproc=1)
    expect_identical(db_mutations, db_mutations_df)
    
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
     
    db_mutations_null <- observedMutations(db_subset,
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

test_that("calcBaseline", {
    
    # collapse baseline
    db_clonal <- collapseClones(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                method="thresholdedFreq", minimumFrequency=0.6,
                                includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
    
    
    ## With the full db, and CLONE field
    ## Mutations calculated CLONAL_SEQUENCE vs GERMLINE_IMGT_D_MASK    
    #set.seed(723)
    db_baseline <- calcBaseline(db=db_clonal, 
                                sequenceColumn="CLONAL_SEQUENCE",
                                germlineColumn="CLONAL_GERMLINE", 
                                testStatistic="focused",
                                regionDefinition=IMGT_V,
                                targetingModel = HH_S5F,
                                nproc = 1,
                                calcStats = T)
    
    ## Check if the stats slot has been filled
    expect_gt(nrow(slot(db_baseline,"stats")),0)
    
    db_baseline_rowSums <- rowSums(
        slot(db_baseline,"db")[1:5,grep("OBSERVED",colnames(slot(db_baseline,"db")))])

# Commented because working on the code, when clones collapsed.
    expect_equivalent(
        db_baseline_rowSums,
        c(13, 19, 55, 0, 0))
    
#     ## Check 5 examples for each, at different positions
#     ## CDR_R, first 5
#     obs_cdr_r <- slot(db_baseline,"db")$OBSERVED_CDR_R[1:5]
#     exp_cdr_r<- c(2, 6, 17, 18, 0)
#     expect_equal(obs_cdr_r, exp_cdr_r)
#     
#     ## CDR_S, 673:677
#     obs_cdr_s <- slot(db_baseline,"db")$OBSERVED_CDR_S[673:677]
#     exp_cdr_s <- c(3, 5, 5, 5, 0)
#     expect_equal(obs_cdr_s, exp_cdr_s)
#     
#     ## FWR_R, 937:941
#     obs_fwr_r <- slot(db_baseline,"db")$OBSERVED_FWR_R[937:941]
#     exp_fwr_r<- c(0, 7, 14, 7, 0)
#     expect_equal(obs_fwr_r, exp_fwr_r)
#     
#     ## FWR_S, 993:997
#     obs_fwr_s <- slot(db_baseline,"db")$OBSERVED_FWR_S[993:997]
#     exp_fwr_s<- c(10, 0, 0, 0, 0)
#     expect_equal(obs_fwr_s, exp_fwr_s)

    
    ## If we trim the sequences to 1:312 (this is the region IMGT_V)
    ## we expect to get the same results setting regionDefinition=NULL
    ## Setting CLONE to NULL
    ## TODO: test with CLONE
    
    ## aux trimming function
    trimToLength <- function(inputSeq, germlineSeq, rdLength) {
        # Removing IMGT gaps (they should come in threes)
        # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
        germlineSeq <- gsub("\\.\\.\\.", "XXX", germlineSeq)
        #If there is a single gap left convert it to an N
        germlineSeq <- gsub("\\.", "N", germlineSeq)
        # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
        germlineSeq <- gsub("XXX", "...", germlineSeq)
        
        # Removing IMGT gaps (they should come in threes)
        # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
        inputSeq <- gsub("\\.\\.\\.", "XXX", inputSeq)
        #If there is a single gap left convert it to an N
        inputSeq <- gsub("\\.", "N", inputSeq)
        # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
        inputSeq <- gsub("XXX", "...", inputSeq)    
        
        # Trim the input and germline sequence to the shortest
        len_inputSeq <- nchar(inputSeq)
        len_germlineSeq <- nchar(germlineSeq)
        rdLength <- 312
        
        len_shortest <- min(c(len_inputSeq, len_germlineSeq, rdLength), na.rm=TRUE)
        
        c_inputSeq <- seqinr::s2c(inputSeq)[1:len_shortest]
        c_germlineSeq <- seqinr::s2c(germlineSeq)[1:len_shortest]
        
        # If the sequence and germline (which now should be the same length) is shorter
        # than the rdLength, pad it with Ns
        if(len_shortest<rdLength){
            fillWithNs <- array("N",rdLength-len_shortest)
            c_inputSeq <- c( c_inputSeq, fillWithNs)
            c_germlineSeq <- c( c_germlineSeq, fillWithNs)
        }        
        
        return(list("inputSeq"=paste(c_inputSeq,collapse=""),
                    "germlineSeq"=paste(c_germlineSeq,collapse="")))
    }
    
    db_clonal_trim <- db_clonal[1:5, ]
    
    for (i in 1:nrow(db_clonal_trim)) {
        trim_seqs <- trimToLength(db_clonal_trim$CLONAL_SEQUENCE[i], 
                                  db_clonal_trim$CLONAL_GERMLINE[i], 312)
        db_clonal_trim$CLONAL_SEQUENCE[i] <- trim_seqs$inputSeq
        db_clonal_trim$CLONAL_GERMLINE[i] <- trim_seqs$germlineSeq
    }
    
    db_1_to_5 <- db_clonal[1:5,]
    
    #set.seed(283)
    db_baseline <- calcBaseline(db_1_to_5, 
                                sequenceColumn="CLONAL_SEQUENCE",
                                germlineColumn="CLONAL_GERMLINE", 
                                testStatistic="focused",
                                regionDefinition=IMGT_V,
                                targetingModel = HH_S5F,
                                nproc = 1,
                                calcStats = T)
    
    #set.seed(2935)
    db_baseline_trim_null <- calcBaseline(db_clonal_trim, sequenceColumn="CLONAL_SEQUENCE",
                                          germlineColumn="CLONAL_GERMLINE",
                                          testStatistic="focused",
                                          regionDefinition=NULL,
                                          targetingModel = HH_S5F,
                                          nproc = 1, 
                                          calcStats=T)
 
    total_trim_null <- rowSums(cbind(slot(db_baseline_trim_null,"db")$OBSERVED_SEQ_S,
        slot(db_baseline_trim_null,"db")$OBSERVED_SEQ_R))
    
    total_baseline <- rowSums(slot(db_baseline,"db")[,grep("OBSERVED.*", colnames(slot(db_baseline,"db")))])
    
    expect_equivalent(total_trim_null, total_baseline)
    
    ## Should match observedMutations, with the full seqs and region 
    ## IMGT_V and the trimmed seqs and region NULL
    obs_mu <- observedMutations(db_1_to_5,
                                sequenceColumn="CLONAL_SEQUENCE",
                                germlineColumn="CLONAL_GERMLINE",
                                regionDefinition=IMGT_V,
                                nproc=1)
    expect_equivalent(
        rowSums(obs_mu[,grep("OBSERVED", colnames(obs_mu))]), 
        total_baseline
        )
    
    obs_mu <- observedMutations(db_clonal_trim,
                                sequenceColumn="CLONAL_SEQUENCE",
                                germlineColumn="CLONAL_GERMLINE",
                                regionDefinition=NULL,
                                nproc=1)
    expect_equivalent(
        rowSums(obs_mu[,grep("OBSERVED", colnames(obs_mu))]), 
        total_baseline
    )
    
    # db_baseline_null <- calcBaseline(db_clonal,
    #                             sequenceColumn="CLONAL_SEQUENCE",
    #                             germlineColumn="CLONAL_GERMLINE",
    #                             testStatistic="focused",
    #                             regionDefinition=NULL,
    #                             targetingModel = HH_S5F,
    #                             nproc = 1)
    # ## Check 5 examples for each, at different positions
    # ## CDR_R, first 5
    # obs_seq_r <- slot(db_baseline_null,"db")$OBSERVED_SEQ_R[1:5]
    # exp_seq_r<- c(12, 18, 74, 73)
    # expect_equal(obs_seq_r, exp_seq_r)
    # 
    # obs_seq_s <- slot(db_baseline_null,"db")$OBSERVED_SEQ_S[1:5]
    # exp_seq_s<- c(5, 6, 34, 46, 1)
    # expect_equal(obs_seq_r, exp_seq_r)

    
})


test_that("Test groupBaseline", {
    # Subset example data from alakazam
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG"))
                      
    # collapse baseline
    db_clonal <- collapseClones(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                method="thresholdedFreq", minimumFrequency=0.6,
                                includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
    
    
    # Calculate BASELINe
    baseline <- calcBaseline(db=db_clonal, 
                             sequenceColumn="CLONAL_SEQUENCE",
                             germlineColumn="CLONAL_GERMLINE",
                             testStatistic="focused",
                             regionDefinition=IMGT_V,
                             targetingModel=HH_S5F,
                             nproc=1)
    
    baseline_df <- calcBaseline(data.frame(db_clonal), 
                             sequenceColumn="CLONAL_SEQUENCE",
                             germlineColumn="CLONAL_GERMLINE",
                             testStatistic="focused",
                             regionDefinition=IMGT_V,
                             targetingModel=HH_S5F,
                             nproc=1)
    
    ## expect same results with tbl and data.frame
    expect_identical(baseline, baseline_df)
    
    # Group PDFs by sample
    grouped1 <- groupBaseline(baseline, groupBy="SAMPLE")
    pdf1 <- slot(grouped1, "pdfs")
    sigma1 <- slot(grouped1,"stats")$BASELINE_SIGMA
    
    expect_equal(range(pdf1$CDR[1,]), c(0,5.018), tolerance=0.01)
    expect_equal(range(pdf1$CDR[2,]), c(0,7.333), tolerance=0.01)
    expect_equal(sigma1, c(-0.263, -0.693, -0.100, -0.694), tolerance=0.01)
    
    # Group PDFs by both sample (between variable) and isotype (within variable)
    grouped2 <- groupBaseline(baseline, groupBy=c("SAMPLE", "ISOTYPE"))
    pdf2 <- slot(grouped2, "pdfs")
    sigma2 <- slot(grouped2,"stats")$BASELINE_SIGMA
    
    expect_equal(range(pdf2$CDR[1,]), c(0,3.643), tolerance=0.01)
    expect_equal(range(pdf2$CDR[2,]), c(0,3.539), tolerance=0.01)    
    expect_equal(sigma2, 
                 c(-0.31, -0.59, -0.20, -0.82, -0.15, -0.78, -0.08, -0.65), tolerance=0.01)
    
    # Collapse previous isotype (within variable) grouped PDFs into sample PDFs
    grouped3 <- groupBaseline(grouped2, groupBy="SAMPLE")
    pdf3 <- slot(grouped3, "pdfs")
    sigma3 <- slot(grouped3,"stats")$BASELINE_SIGMA
    
    expect_equal(range(pdf3$CDR[1,]), c(0,4.975), tolerance=0.01)
    expect_equal(range(pdf3$CDR[2,]), c(0,7.319), tolerance=0.01)    
    expect_equal(sigma3, 
                 c( -0.25, -0.72, -0.10, -0.69), tolerance=0.01)
    
})
