test_that("binMutationsByRegion", {
    set.seed(8)
    numbOfMutations <- sample(3:10, 1) 
    set.seed(60)
    posOfMutations <- sort(sample(330, numbOfMutations))
    set.seed(13)
    mutation_types <- sample(c("R", "S"), length(posOfMutations), replace=TRUE)
    mutations_array <- array(mutation_types, dimnames=list(posOfMutations))
    
    observed_bin <- binMutationsByRegion(mutations_array, regionDefinition=IMGT_V_NO_CDR3)
    expected_bin <- c(1, 0, 3, 2)
    names(expected_bin) <- c("CDR_R", "CDR_S", "FWR_R", "FWR_S")
    expect_equal(observed_bin, expected_bin)
    
    observed_bin <- binMutationsByRegion(mutations_array, regionDefinition=NULL)
    expected_bin <- c(4, 2)
    names(expected_bin) <- c("SEQ_R", "SEQ_S")
    expect_equal(observed_bin,expected_bin)
})
