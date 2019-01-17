load(file.path("..", "data-tests", "ExampleDb.rda"))

#### createSubstitutionMatrix ####

test_that("createSubstitutionMatrix", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    # Create model using only silent mutations
    
    sub <- createSubstitutionMatrix(db, model="S")
    set.seed(311)
    idx <- sample(length(sub), 10)
    expected <- c(0.621, 0.314, NA, NA, 0.16, NA, 0.531, 0.152, 0.314, 0.057)
    expect_equal(sub[idx], expected, tolerance=0.001)
})

#### createMutabilityMatrix ####

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

#### extendSubstitutionMatrix ####

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

#### extendMutabilityMatrix ####

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

#### createTargetingMatrix ####

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

#### createTargetingModel ####

test_that("createTargetingModel", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")

    expect_warning(model <- createTargetingModel(db, model="S", multipleMutation="ignore"),
                   "Insufficient number of mutations")
    
    set.seed(9)
    idx <- sample(length(model@targeting), 20)
    expected_model <- c(NA, NA, NA, 0, 0.00125, NA, 0.000139, 0, 0.00073, 0.00009, NA, 0, 0, 0.00782, 0, NA, NA, 0.000989, NA, 0)
    expect_equal(model@targeting[idx], expected_model, tolerance=0.001)    
})

#### calcTargetingDistance ####

test_that("calcTargetingDistance", {
    
    dist <- calcTargetingDistance(HH_S5F)
    
    set.seed(15)
    idx <- sample(length(dist), 20)
    expected_dist <- c(0.94, 0, 0.93, 1.03, 0, 1.03, 0.96, 0, 0.94, 0, 0, 0, 0, 0.96, 0.94, 0, 0, 1.00, 0.85, 0)
    expect_equal(dist[idx], expected_dist, tolerance=0.001)
})

#### rescaleMutability ####

test_that("rescaleMutability", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    # Create model and rescale mutabilities
    expect_warning(model <- createTargetingModel(db, model="S", multipleMutation="ignore"),
                   "Insufficient number of mutations")
    mut <- shazam:::rescaleMutability(model)    
    
    set.seed(95)
    idx <- sample(length(mut), 20)
    expected_mut <- c(0, 0.247, 1.248, 0, 0, NA, 0, 0.555, 0.27, 0, NA, 2.035, NA, 0.921, 0.775, 2.035, NA, 0, 0, NA)
    expect_equal(as.vector(mut[idx]), expected_mut, tolerance=0.001)   
})

#### removeCodonGaps ####

test_that("removeCodonGaps", {
    test_db <- data.frame("A"=c("ATC","...","ATC", "ATC...", "ATC..."), "B"=c("ATC","...","...", "ATC...","...ATC"))
    obs <- shazam:::removeCodonGaps(test_db)
    exp <- matrix(c("ATC","ATC","","","","","ATC","ATC","", ""), ncol=2, byrow = T)
    expect_equal(obs, exp)
})

#### calcObservedMutations and listObservedMutations ####

test_that("calcObservedMutations and listObservedMutations reproduce same results", {
    
    db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")
    
    lom <- mapply(shazam:::listMutations, db[["SEQUENCE_IMGT"]], db[["GERMLINE_IMGT_D_MASK"]], 
           "independent", "RS", USE.NAMES=FALSE)
    
#     cob <- mapply(calcObservedMutations,db[["SEQUENCE_IMGT"]], db[["GERMLINE_IMGT_D_MASK"]])
#     om <- observedMutations(db, "SEQUENCE_IMGT", "GERMLINE_IMGT_D_MASK")
#     
#     expect_equal(length(lom[[1]]), sum(cob[,1]))
#     expect_equal(length(lom[[1]]), sum(om[1, c("MU_COUNT_SEQ_R", "MU_COUNT_SEQ_S")]))
#     
    ## Use only V
    om_v <- observedMutations(db, "SEQUENCE_IMGT", "GERMLINE_IMGT_D_MASK", 
                              regionDefinition=IMGT_V)
    ## TODO check more seqs
    expect_equal(length(lom[[1]]), sum(om_v[1, grep("MU_COUNT_", colnames(om_v))]))
    
})

#### makeDegenerate5merSub ####

test_that("makeDegenerate5merSub with extended=FALSE using HH_S1F", {
    degenerate5mer = makeDegenerate5merSub(HH_S1F, extended=FALSE)
    centers = sapply(colnames(degenerate5mer), substr, start=3, stop=3)
    
    test_HH_S1F = HH_S1F
    diag(test_HH_S1F) = NA
    
    for (nuc in c("A", "T", "G", "C")) {
        colIdx = which(centers==nuc)
        for (i in colIdx) {
            expect_equal(test_HH_S1F[nuc, ],
                         degenerate5mer[, i])
        }
    }
})

test_that("makeDegenerate5merSub with extended=TRUE using HH_S1F", {
    degenerate5mer = makeDegenerate5merSub(HH_S1F, extended=TRUE)
    centers = sapply(colnames(degenerate5mer), substr, start=3, stop=3)
    
    test_HH_S1F = HH_S1F
    diag(test_HH_S1F) = NA
    
    # 5mers with non-N central 1mer (incl. eg. ATTTN) should have same rates as 
    # corresponding central 1mer
    for (nuc in c("A", "T", "G", "C")) {
        colIdx = which(centers==nuc)
        for (i in colIdx) {
            expect_equal(test_HH_S1F[nuc, ],
                         degenerate5mer[-5, i])
        }
    }
    
    # 5mers with N central 1mer should have all NA rates
    colIdx = which(centers=="N")
    for (i in colIdx) {
        expect_equal(sum(is.na(degenerate5mer[, i])), 5)
    }
})

#### makeAverage1merSub ####

test_that("makeAverage1merSub using HH_S1F", {
    degenerate5mer = makeDegenerate5merSub(HH_S1F, extended=FALSE)
    average1mer = makeAverage1merSub(degenerate5mer)
    
    test_HH_S1F = HH_S1F
    diag(test_HH_S1F) = NA
    test_HH_S1F = test_HH_S1F[match(rownames(average1mer), rownames(test_HH_S1F)),
                              match(colnames(average1mer), colnames(test_HH_S1F))]
    
    expect_equal(test_HH_S1F, average1mer)
})

#### makeDegenerate5merMut ####

test_that("makeDegenerate5merMut with extended=FALSE using an example", {
    example1merMut <- c(A=0.2, T=0.1, C=0.4, G=0.3)
    degenerate5mer = makeDegenerate5merMut(example1merMut, extended=FALSE)
    centers = sapply(colnames(degenerate5mer), substr, start=3, stop=3)
    
    for (nuc in c("A", "T", "G", "C")) {
        idx = which(centers==nuc)
        for (i in idx) {
            expect_equal(example1merMut[nuc], degenerate5mer[idx])
        }
    }
})

test_that("makeDegenerate5merMut with extended=TRUE using an example", {
    example1merMut <- c(A=0.2, T=0.1, C=0.4, G=0.3)
    degenerate5mer = makeDegenerate5merMut(example1merMut, extended=TRUE)
    centers = sapply(colnames(degenerate5mer), substr, start=3, stop=3)
    
    # 5mers with non-N central 1mer (incl. eg. ATTTN) should have same rates as 
    # corresponding central 1mer
    for (nuc in c("A", "T", "G", "C")) {
        idx = which(centers==nuc)
        for (i in idx) {
            expect_equal(example1merMut[nuc], degenerate5mer[idx])
        }
    }
    
    # 5mers with N central 1mer should have all NA rates
    idx = which(centers=="N")
    for (i in idx) {
        expect_equal(NA, degenerate5mer[idx])
    }
})

#### makeAverage1merMut ####

test_that("makeAverage1merMut using an example", {
    example1merMut <- c(A=0.2, T=0.1, C=0.4, G=0.3)
    degenerate5mer = makeDegenerate5merMut(example1merMut, extended=FALSE)
    average1mer = makeAverage1merMut(degenerate5mer)
    
    expect_equal(example1merMut[match(names(average1mer), 
                                      names(example1merMut))], 
                 average1mer)
})
