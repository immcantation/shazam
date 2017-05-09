#load(file.path("tests", "data-tests", "ExampleDb.rda"))
load(file.path("..", "data-tests", "ExampleDb.rda"))

test_that("collapseClones", {
    # example data
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
    # clone with most sequences: CLONE 3128
    
    # Build clonal consensus for the full sequence
    set.seed(7)
    clones.1 <- collapseClones(db, nproc=1, expandedDb=FALSE)
    set.seed(7)
    clones.2 <- collapseClones(db, nproc=1, expandedDb=TRUE)
    
    for (clone in unique(db$CLONE)) {
        # expect CLONAL_SEQUENCE for all seqs in the same clone to be the same
        # expect CLONAL_GERMLINE for all seqs in the same clone to be the same
        # use result from expandedDb=TRUE to test
        expect_equal(length(unique(clones.2[clones.2$CLONE==clone, "CLONAL_SEQUENCE"])), 1)
        expect_equal(length(unique(clones.2[clones.2$CLONE==clone, "CLONAL_GERMLINE"])), 1)
        
        # expect result from expandedDb=TRUE to be the same as that from expandedDb=FALSE
        expect_equal(clones.1[clones.1$CLONE==clone, ], clones.2[clones.2$CLONE==clone, ][1, ])
    }
}) 


test_that("binMutationsByRegion", {
    set.seed(8)
    numbOfMutations <- sample(3:10, 1) 
    set.seed(60)
    posOfMutations <- sort(sample(330, numbOfMutations))
    set.seed(13)
    mutation_types <- sample(c("R", "S"), length(posOfMutations), replace=TRUE)
    mutations_array <- array(mutation_types, dimnames=list(posOfMutations))
    
    observed_bin <- binMutationsByRegion(mutations_array, regionDefinition=IMGT_V)
    expected_bin <- c(1, 0, 3, 2)
    names(expected_bin) <- c("CDR_R", "CDR_S", "FWR_R", "FWR_S")
    expect_equal(observed_bin, expected_bin)
    
    observed_bin <- binMutationsByRegion(mutations_array, regionDefinition=NULL)
    expected_bin <- c(4, 2)
    names(expected_bin) <- c("SEQ_R", "SEQ_S")
    expect_equal(observed_bin,expected_bin)
})


test_that("observedMutations, charge mutations", {

    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")

    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                             germlineColumn="GERMLINE_IMGT_D_MASK",
                             regionDefinition=IMGT_V,
                             mutationDefinition=CHARGE_MUTATIONS,
                             nproc=1)
    
    expect_equal(db_obs$OBSERVED_CDR_R[1:10], c(0, 0, 0, 0, 0, 0, 2, 2, 2, 0))
    expect_equal(db_obs$OBSERVED_CDR_S[1:10], c(0, 0, 0, 0, 0, 0, 3, 3, 3, 16))
    expect_equal(db_obs$OBSERVED_FWR_R[1:10], c(0, 1, 1, 0, 0, 0, 0, 0, 0, 3))
    expect_equal(db_obs$OBSERVED_FWR_S[1:10], c(0, 6, 6, 0, 0, 0, 10, 10, 10, 14))
    
})

test_that("calcObservedMutations, hydropathy", {
    in_seq <- ExampleDb[1, "SEQUENCE_IMGT"]
    germ_seq <-  ExampleDb[1, "GERMLINE_IMGT_D_MASK"]
     
    #' # Identify all mutations in the sequence
    expect_equivalent(calcObservedMutations(in_seq, germ_seq), c(0,2))
     
    # Identify only mutations in the V segment minus CDR3
    expect_equivalent(
        calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V),
        NA)
      
    # Identify mutations by change in hydropathy class
    expect_equivalent(
        calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
                           mutationDefinition=HYDROPATHY_MUTATIONS, frequency=TRUE),
        NA)
    
    ## This seq is only V
    v_seq <- ".......................................................................................CGTGTG......CTCTTCCGATCTAAGTCTTCCAAAGTACTANACGGGGGACTCCTGTGCCCCACCATGGACACACTTTGCTCCACGCTC.........CTGCTGCTGACCATCCCTTCATGGGGCTTG...TCCCAGATCACCATCTCGAGGGACACCTCCAAAAACCAGGTGATCCTTTCAATGACCAACATGGGCCCTCTGGACACAGCCACATACTTCTGTGCACACAGAC............................................................"
    germ_seq <- "CAGATCACCTTGAAGGAGTCTGGTCCT...ACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCACCTTCTCTGGGTTCTCACTCAGC......ACTAGTGGAGTGGGTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGGCCCTGGAGTGGCTTGCACTCATTTATTGGGAT.........GATGATAAGCGCTACAGCCCATCTCTGAAG...AGCAGGCTCACCATCACCAAGGACACCTCCAAAAACCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGACACAGCCACATATTACTGTGCACACAGACNNNNNNNNNNNNNNNNNTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG	"
    expect_equivalent(
        sum(calcObservedMutations(v_seq, germ_seq, regionDefinition=NULL, freq=F)),
        sum(calcObservedMutations(v_seq, germ_seq, regionDefinition=IMGT_V, freq=F)))
    
})

test_that("observedMutations, combine", {
    
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA") & SAMPLE == "+7d")[1:25,]
    
    ##
    ## With regionDefinition==NULL
    ##
    
    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=F,
                                combine=F,
                                regionDefinition=NULL,
                                mutationDefinition=NULL,
                                nproc=1)
    
    db_freq <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=T,
                                combine=F,
                                regionDefinition=NULL,
                                mutationDefinition=NULL,
                                nproc=1)
    
    db_obs_combined <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                          germlineColumn="GERMLINE_IMGT_D_MASK",
                                          frequency=F,
                                          combine=T,
                                          regionDefinition=NULL,
                                          mutationDefinition=NULL,
                                          nproc=1)    
    
    db_freq_combined <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 frequency=T,
                                 combine=T,
                                 regionDefinition=NULL,
                                 mutationDefinition=NULL,
                                 nproc=1)
    
    ## When using the whole sequence, the sum of OBSERVED_SEQ_S and OBSERVED_SEQ_R
    ## should match OBSERVED 
    expect_equal(rowSums(db_obs[,grep("OBSERVED",colnames(db_obs))]),
                 db_obs_combined$OBSERVED)

    ## When using the whole sequence, the sum of MU_FREQ_SEQ_S and MU_FREQ_SEQ_R
    ## should match MU_FREQ 
    expect_equal(rowSums(db_freq[,grep("MU_FREQ",colnames(db_freq))]),
                 db_freq_combined$MU_FREQ)
    
    ##
    ## With regionDefinition==IMGT_V
    ##
    
    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=F,
                                combine=F,
                                regionDefinition=IMGT_V,
                                mutationDefinition=NULL,
                                nproc=1)
    
    db_freq <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 frequency=T,
                                 combine=F,
                                 regionDefinition=IMGT_V,
                                 mutationDefinition=NULL,
                                 nproc=1)
    
    db_obs_combined <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                         germlineColumn="GERMLINE_IMGT_D_MASK",
                                         frequency=F,
                                         combine=T,
                                         regionDefinition=IMGT_V,
                                         mutationDefinition=NULL,
                                         nproc=1)    
    ## Get the number of nonN positions, the denomitator for the frequency
    db_obs_denominator <- sapply(1:dim(db)[1], 
                        function(i) {
                            mu_count <- calcObservedMutations(inputSeq = db[i, "SEQUENCE_IMGT"],
                                                              germlineSeq = db[i,"GERMLINE_IMGT_D_MASK"],
                                                              regionDefinition = IMGT_V,
                                                              returnRaw = T, freq=F)
                            mu_den <- sum(mu_count$nonN)
                        })
    
    db_freq_combined <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                          germlineColumn="GERMLINE_IMGT_D_MASK",
                                          frequency=T,
                                          combine=T,
                                          regionDefinition=IMGT_V,
                                          mutationDefinition=NULL,
                                          nproc=1)
    
    ## When using IMGT_V, the sum of OBSERVED_SEQ_S and OBSERVED_SEQ_R
    ## should match OBSERVED. There is only one denominator, nonN-SEQ
    expect_equal(rowSums(db_obs[,grep("OBSERVED",colnames(db_obs))]),
                 db_obs_combined$OBSERVED)
    
    ## When not using the whole sequence, the sum of the mutation frequencies
    ## may not match MU_FREQ, because CDR mutations and FWR mutations use their own
    ## denominators (nonN-CDR and nonN-FWR)
    expect_equal(db_obs_combined$OBSERVED/db_obs_denominator,
        db_freq_combined$MU_FREQ)
    
})

test_that("calcExpecteddMutations, hydropathy", {
    
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
    
    # Calculate hydropathy expected mutations over V region
    db_exp <- expectedMutations(db,
                               sequenceColumn="SEQUENCE_IMGT",
                               germlineColumn="GERMLINE_IMGT_D_MASK",
                               regionDefinition=IMGT_V,
                               mutationDefinition=HYDROPATHY_MUTATIONS,
                               nproc=1)    
    expect_equal(db_exp$EXPECTED_CDR_R[1:10],
                c(0.123, 0.114, 0.114, 0.131, 0.131, 0.131, 0.118, 0.118, 0.118, 0.139),
                tolerance=0.001
    )
    expect_equal(db_exp$EXPECTED_CDR_S[10:20],
                 c(0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.150),
                 tolerance=0.001
    )
    expect_equal(db_exp$EXPECTED_FWR_R[20:30],
                 c(0.318, 0.318, 0.318, 0.318, 0.332, 0.332, 0.323, 0.323, 0.323, 0.323, 0.315),
                 tolerance=0.001
    )
    expect_equal(db_exp$EXPECTED_FWR_S[30:40],
                 c(0.436, 0.436, 0.437, 0.436, 0.439, 0.462, 0.429, 0.484, 0.462, 0.462, 0.433),
                 tolerance=0.001
    )
})    


test_that("observedMutations overwrites with a warning pre-existing mutation counts/freqs", {
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
    
    ## Counts
    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=F,
                                combine=F,
                                regionDefinition=NULL,
                                mutationDefinition=NULL,
                                nproc=1)
    expect_warning(observedMutations(db_obs, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=F,
                                combine=F,
                                regionDefinition=NULL,
                                mutationDefinition=NULL,
                                nproc=1),
                   "Columns OBSERVED_SEQ_R, OBSERVED_SEQ_S exist and will be overwritten")
    ## Counts Combine
    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=F,
                                combine=T,
                                regionDefinition=NULL,
                                mutationDefinition=NULL,
                                nproc=1)
    expect_warning(observedMutations(db_obs, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     frequency=F,
                                     combine=T,
                                     regionDefinition=NULL,
                                     mutationDefinition=NULL,
                                     nproc=1),
                   "Columns OBSERVED exist and will be overwritten")    
    
    ## Counts  regionDefinition
    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=F,
                                combine=F,
                                regionDefinition=IMGT_V,
                                mutationDefinition=NULL,
                                nproc=1)
    expect_warning(observedMutations(db_obs, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     frequency=F,
                                     combine=F,
                                     regionDefinition=IMGT_V,
                                     mutationDefinition=NULL,
                                     nproc=1),
                   "Columns OBSERVED_CDR_R, OBSERVED_CDR_S, OBSERVED_FWR_R, OBSERVED_FWR_S exist and will be overwritten")  
    ## Counts combine regionDefinition
    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=F,
                                combine=T,
                                regionDefinition=IMGT_V,
                                mutationDefinition=NULL,
                                nproc=1)
    expect_warning(observedMutations(db_obs, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     frequency=F,
                                     combine=T,
                                     regionDefinition=IMGT_V,
                                     mutationDefinition=NULL,
                                     nproc=1),
                   "Columns OBSERVED exist and will be overwritten")    
    
    
    ## Freq
    db_freq <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                frequency=T,
                                combine=F,
                                regionDefinition=NULL,
                                mutationDefinition=NULL,
                                nproc=1)
    expect_warning(observedMutations(db_freq, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     frequency=T,
                                     combine=F,
                                     regionDefinition=NULL,
                                     mutationDefinition=NULL,
                                     nproc=1),
                   "Columns MU_FREQ_SEQ_R, MU_FREQ_SEQ_S exist and will be overwritten")
    
    ## Freq Combine
    db_freq <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 frequency=T,
                                 combine=T,
                                 regionDefinition=NULL,
                                 mutationDefinition=NULL,
                                 nproc=1)
    expect_warning(observedMutations(db_freq, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     frequency=T,
                                     combine=T,
                                     regionDefinition=NULL,
                                     mutationDefinition=NULL,
                                     nproc=1),
                   "Columns MU_FREQ exist and will be overwritten")
})



test_that("expectedMutations overwrites with a warning pre-existing values", {
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
    
    db_exp <- expectedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                regionDefinition=NULL,
                                mutationDefinition=NULL,
                                nproc=1)
    expect_warning(expectedMutations(db_exp, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     regionDefinition=NULL,
                                     mutationDefinition=NULL,
                                     nproc=1),
                   "Columns EXPECTED_SEQ_R, EXPECTED_SEQ_S exist and will be overwritten")
    
    db_exp <- expectedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                regionDefinition=IMGT_V,
                                mutationDefinition=NULL,
                                nproc=1)
    expect_warning(expectedMutations(db_exp, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     regionDefinition=IMGT_V,
                                     mutationDefinition=NULL,
                                     nproc=1),
                   "Columns EXPECTED_CDR_R, EXPECTED_CDR_S, EXPECTED_FWR_R, EXPECTED_FWR_S exist and will be overwritten")
    
})