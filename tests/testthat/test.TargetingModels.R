load(file.path("..", "data-tests", "ExampleDb.rda"))

test_that("createSubstitutionMatrix", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    # Create model using only silent mutations
    
    sub <- createSubstitutionMatrix(db, model="S")
    set.seed(311)
    idx <- sample(length(sub), 10)
    expected <- c(0.621, 0.314, NA, NA, 0.16, NA, 0.531, 0.152, 0.314, 0.057)
    expect_equal(sub[idx], expected, tolerance=0.001)
})


test_that("createMutabilityMatrix", {
   
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    # Create model using only silent mutations
    sub_model <- createSubstitutionMatrix(db, model="S")
    expect_warning(mut_model <- createMutabilityMatrix(db, sub_model, model="S"),
                   "Insufficient number of mutations")
    
    set.seed(61)
    idx <- sample(length(mut_model), 10)
    expected_mut <- c(0, 0, 0.003, 0.002, 0, 0.001, 0.001, 0, 0.004, 0)
    expect_equal(as.vector(mut_model[idx]), expected_mut, tolerance=0.001)
})

test_that("extendSubstitutionMatrix", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    # Create model using only silent mutations
    sub_model <- createSubstitutionMatrix(db, model="S")
    ext_model <- extendSubstitutionMatrix(sub_model)
    
    set.seed(4)
    idx <- sample(length(ext_model), 10)
    expected_ext <- c(NA, NA, NA, NA, 0.621, NA, NA, 0.547, 0.547, NA)
    expect_equal(ext_model[idx], expected_ext, tolerance=0.001)
})

test_that("extendMutabilityMatrix", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    # Create model using only silent mutations
    sub_model <- createSubstitutionMatrix(db, model="S")
    expect_warning(mut_model <- createMutabilityMatrix(db, sub_model, model="S"),
                   "Insufficient number of mutations"
    )
    ext_model <- extendMutabilityMatrix(mut_model)
    
    set.seed(31)
    idx <- sample(length(ext_model), 10)
    expected_ext <- c(0, 0, 0, 0, 0.002, 0.001, 0.001, 0.003, 0, 0)
    expect_equal(as.vector(ext_model[idx]), expected_ext, tolerance=0.001)
})

test_that("createTargetingMatrix", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    # Create 4x1024 models using only silent mutations
    sub_model <- createSubstitutionMatrix(db, model="S")
    expect_warning(mut_model <- createMutabilityMatrix(db, sub_model, model="S"),
                   "Insufficient number of mutations")
     
    # Extend substitution and mutability to including Ns (5x3125 model)
    sub_model <- extendSubstitutionMatrix(sub_model)
    mut_model <- extendMutabilityMatrix(mut_model)
     
    # Create targeting model from substitution and mutability
    tar_mat <- createTargetingMatrix(sub_model, mut_model)    
    
    set.seed(48)
    idx <- sample(length(tar_mat), 20)
    expected_tar <- c(NA, NA, 0, NA, NA, NA, NA, NA, NA, 0, 0, 0.001, 0, NA, NA, 0, NA, 0, NA, NA)
    expect_equal(tar_mat[idx], expected_tar, tolerance=0.001)
})    
    
test_that("createTargetingModel", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")

    expect_warning(model <- createTargetingModel(db, model="S", multipleMutation="ignore"),
                   "Insufficient number of mutations")
    
    set.seed(9)
    idx <- sample(length(model), 20)
    expected_model <- c(0, 0.001, 0, NA, 0, 0, 0, NA, NA, NA, 0.001, 0.003, NA, 0.001, NA, NA, NA, NA, 0, NA)
    expect_equal(model@targeting[idx], expected_model, tolerance=0.001)    
})    

test_that("calcTargetingDistance", {
    
    dist <- calcTargetingDistance(HH_S5F)
    
    set.seed(15)
    idx <- sample(length(dist), 20)
    expected_dist <- c(0.937, NA, 0.931, 1.03, NA, 1.032, 0.963, 0, 0.938, 0, NA, NA, 0, 0.959, 0.942, 0, NA, 0.995, 0.855, 0)
    expect_equal(dist[idx], expected_dist, tolerance=0.001)    
})

test_that("rescaleMutability", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    # Create model and rescale mutabilities
    model <- createTargetingModel(db, model="S", multipleMutation="ignore")
    mut <- shazam:::rescaleMutability(model)    
    
    set.seed(95)
    idx <- sample(length(mut), 20)
    expected_mut <- c(0, 0.247, 1.248, 0, 0, NA, 0, 0.555, 0.27, 0, NA, 2.035, NA, 0.921, 0.775, 2.035, NA, 0, 0, NA)
    expect_equal(as.vector(mut[idx]), expected_mut, tolerance=0.001)   
})

test_that("removeCodonGaps", {
    test_db <- data.frame("A"=c("ATC","...","ATC", "ATC...", "ATC..."), "B"=c("ATC","...","...", "ATC...","...ATC"))
    obs <- removeCodonGaps(test_db)
    exp <- matrix(c("ATC","ATC","","","","","ATC","ATC","", ""), ncol=2, byrow = T)
    expect_equal(obs, exp)
})

test_that("calcOBservedMutations and listObservedMutations reproduce same results", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    lom <- mapply(shazam:::listMutations, db[["SEQUENCE_IMGT"]], db[["GERMLINE_IMGT_D_MASK"]], 
           "independent", "RS", USE.NAMES=FALSE)
    
#     cob <- mapply(calcObservedMutations,db[["SEQUENCE_IMGT"]], db[["GERMLINE_IMGT_D_MASK"]])
#     om <- observedMutations(db, "SEQUENCE_IMGT", "GERMLINE_IMGT_D_MASK")
#     
#     expect_equal(length(lom[[1]]), sum(cob[,1]))
#     expect_equal(length(lom[[1]]), sum(om[1, c("OBSERVED_SEQ_R", "OBSERVED_SEQ_S")]))
#     
    ## Use only V
    om_v <- observedMutations(db, "SEQUENCE_IMGT", "GERMLINE_IMGT_D_MASK", 
                              regionDefinition=IMGT_V)
    ## TODO check more seqs
    expect_equal(length(lom[[1]]), sum(om_v[1, grep("OBSERVED_", colnames(om_v))]))
    
})