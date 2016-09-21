#load(file.path("tests", "data-tests", "ExampleDb.rda"))
load(file.path("..", "data-tests", "ExampleDb.rda"))

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
