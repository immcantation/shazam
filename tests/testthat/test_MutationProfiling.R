#load(file.path("tests", "data-tests", "ExampleDb.rda"))
load(file.path("..", "data-tests", "ExampleDb.rda"))

#### collapseClones ####
test_that("collapseClones", {
    # example data
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
    # clone with most sequences: CLONE 3128
    
    # Build clonal consensus for the full sequence
    set.seed(7)
    clones.1 <- collapseClones(db, method="thresholdedFreq", minimumFrequency=0.2, breakTiesStochastic=TRUE, nproc=1, expandedDb=FALSE)
    set.seed(7)
    clones.1.df <- collapseClones(db, method="thresholdedFreq", minimumFrequency=0.2, breakTiesStochastic=TRUE, nproc=1, expandedDb=FALSE)
    expect_identical(clones.1.df, clones.1)
    
    set.seed(7)
    clones.2 <- collapseClones(db, method="thresholdedFreq", minimumFrequency=0.2, breakTiesStochastic=TRUE, nproc=1, expandedDb=TRUE)
    
    for (clone in unique(db$CLONE)) {
        # expect CLONAL_SEQUENCE for all seqs in the same clone to be the same
        # expect CLONAL_GERMLINE for all seqs in the same clone to be the same
        # use result from expandedDb=TRUE to test
        expect_equal(length(unique(clones.2[clones.2$CLONE==clone, "CLONAL_SEQUENCE"])), 1)
        expect_equal(length(unique(clones.2[clones.2$CLONE==clone, "CLONAL_GERMLINE"])), 1)
        
        # expect result from expandedDb=TRUE to be the same as that from expandedDb=FALSE
        expect_identical(clones.1[clones.1$CLONE==clone, ], clones.2[clones.2$CLONE==clone, ][1,])
    }
}) 

#### binMutationsByRegion ####
test_that("binMutationsByRegion", {
    
    set.seed(8)
    numbOfMutations <- sample(3:10, 1) 
    set.seed(60)
    posOfMutations <- sort(sample(330, numbOfMutations))
    set.seed(13)
    mutations_array <- matrix(0, nrow=2, ncol=numbOfMutations, dimnames=list(c("R", "S"), posOfMutations))
    mutations_array["R", ] = sample(x=0:10, size=numbOfMutations, replace=TRUE)
    mutations_array["S", ] = sample(x=0:10, size=numbOfMutations, replace=TRUE)
    
    observed_bin <- shazam:::binMutationsByRegion(mutations_array, regionDefinition=IMGT_V)
    expected_bin <- c(2, 8, 22, 31)
    names(expected_bin) <- c("CDR_R", "CDR_S", "FWR_R", "FWR_S")
    expect_equal(observed_bin, expected_bin)
    
    observed_bin <- shazam:::binMutationsByRegion(mutations_array, regionDefinition=NULL)
    expected_bin <- c(24, 39)
    names(expected_bin) <- c("SEQ_R", "SEQ_S")
    expect_equal(observed_bin, expected_bin)
})

#### observedMutations, charge mutations ####
test_that("observedMutations, charge mutations", {
    
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")[1:10, ]
    
    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                regionDefinition=IMGT_V,
                                mutationDefinition=CHARGE_MUTATIONS,
                                nproc=1)
    
    expect_equal(db_obs$MU_COUNT_CDR_R, c(0, 0, 0, 0, 0, 0, 2, 2, 2, 0))
    expect_equal(db_obs$MU_COUNT_CDR_S, c(0, 0, 0, 0, 0, 0, 3, 3, 3, 16))
    expect_equal(db_obs$MU_COUNT_FWR_R, c(0, 1, 1, 0, 0, 0, 0, 0, 0, 3))
    expect_equal(db_obs$MU_COUNT_FWR_S, c(0, 6, 6, 0, 0, 0, 10, 10, 10, 14))
    
})

#### calcObservedMutations, hydropathy ####
test_that("calcObservedMutations, hydropathy", {
    in_seq <- ExampleDb[1, "SEQUENCE_IMGT"]
    germ_seq <-  ExampleDb[1, "GERMLINE_IMGT_D_MASK"]
    
    #' # Identify all mutations in the sequence
    expect_equivalent(calcObservedMutations(in_seq, germ_seq), c(0,2))
    
    # Identify only mutations in the V segment minus CDR3
    expect_equivalent(
        calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V),
        c(NA, NA, NA, NA))
    
    # Identify mutations by change in hydropathy class
    expect_equivalent(
        calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
                              mutationDefinition=HYDROPATHY_MUTATIONS, frequency=TRUE),
        c(NA, NA, NA, NA))
    
    ## This seq is only V
    v_seq <- ".......................................................................................CGTGTG......CTCTTCCGATCTAAGTCTTCCAAAGTACTANACGGGGGACTCCTGTGCCCCACCATGGACACACTTTGCTCCACGCTC.........CTGCTGCTGACCATCCCTTCATGGGGCTTG...TCCCAGATCACCATCTCGAGGGACACCTCCAAAAACCAGGTGATCCTTTCAATGACCAACATGGGCCCTCTGGACACAGCCACATACTTCTGTGCACACAGAC............................................................"
    germ_seq <- "CAGATCACCTTGAAGGAGTCTGGTCCT...ACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCACCTTCTCTGGGTTCTCACTCAGC......ACTAGTGGAGTGGGTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGGCCCTGGAGTGGCTTGCACTCATTTATTGGGAT.........GATGATAAGCGCTACAGCCCATCTCTGAAG...AGCAGGCTCACCATCACCAAGGACACCTCCAAAAACCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGACACAGCCACATATTACTGTGCACACAGACNNNNNNNNNNNNNNNNNTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG	"
    expect_equivalent(
        sum(calcObservedMutations(v_seq, germ_seq, regionDefinition=NULL, freq=F)),
        sum(calcObservedMutations(v_seq, germ_seq, regionDefinition=IMGT_V, freq=F)))
    
})

#### observedMutations, combine ####
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
    
    ## When using the whole sequence, the sum of MU_COUNT_SEQ_S and MU_COUNT_SEQ_R
    ## should match MU_COUNT 
    expect_equal(rowSums(db_obs[,grep("MU_COUNT",colnames(db_obs))]),
                 db_obs_combined$MU_COUNT)
    
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
    
    ## When using IMGT_V, the sum of MU_COUNT_SEQ_S and MU_COUNT_SEQ_R
    ## should match MU_COUNT. There is only one denominator, nonN-SEQ
    expect_equal(rowSums(db_obs[,grep("MU_COUNT",colnames(db_obs))]),
                 db_obs_combined$MU_COUNT)
    
    ## When not using the whole sequence, the sum of the mutation frequencies
    ## may not match MU_FREQ, because CDR mutations and FWR mutations use their own
    ## denominators (nonN-CDR and nonN-FWR)
    expect_equal(db_obs_combined$MU_COUNT/db_obs_denominator,
                 db_freq_combined$MU_FREQ)
    
})

#### expectedMutations, hydropathy ####
test_that("expectedMutations, hydropathy", {
    
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
    
    # Calculate hydropathy expected mutations over V region
    db_exp <- expectedMutations(db,
                                sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT_D_MASK",
                                regionDefinition=IMGT_V,
                                mutationDefinition=HYDROPATHY_MUTATIONS,
                                nproc=1)    
    expect_equal(db_exp$MU_EXPECTED_CDR_R[1:10],
                 c(0.123, 0.114, 0.114, 0.131, 0.131, 0.131, 0.118, 0.118, 0.118, 0.139),
                 tolerance=0.001
    )
    expect_equal(db_exp$MU_EXPECTED_CDR_S[10:20],
                 c(0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.116, 0.150),
                 tolerance=0.001
    )
    expect_equal(db_exp$MU_EXPECTED_FWR_R[20:30],
                 c(0.318, 0.318, 0.318, 0.318, 0.332, 0.332, 0.323, 0.323, 0.323, 0.323, 0.315),
                 tolerance=0.001
    )
    expect_equal(db_exp$MU_EXPECTED_FWR_S[30:40],
                 c(0.436, 0.436, 0.437, 0.436, 0.439, 0.462, 0.429, 0.484, 0.462, 0.462, 0.433),
                 tolerance=0.001
    )
})    

#### observedMutations, warning ####
test_that("observedMutations overwrites with a warning pre-existing mutation counts/freqs", {
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")[1:10, ]
    
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
                   "Columns MU_COUNT_SEQ_R, MU_COUNT_SEQ_S exist and will be overwritten")
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
                   "Columns MU_COUNT exist and will be overwritten")    
    
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
                   "Columns MU_COUNT_CDR_R, MU_COUNT_CDR_S, MU_COUNT_FWR_R, MU_COUNT_FWR_S exist and will be overwritten")  
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
                   "Columns MU_COUNT exist and will be overwritten")    
    
    
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


#### expectedMutations, warning ####
test_that("expectedMutations overwrites with a warning pre-existing values", {
    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")[1:10, ]
    
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
                   "Columns MU_EXPECTED_SEQ_R, MU_EXPECTED_SEQ_S exist and will be overwritten")
    
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
                   "Columns MU_EXPECTED_CDR_R, MU_EXPECTED_CDR_S, MU_EXPECTED_FWR_R, MU_EXPECTED_FWR_S exist and will be overwritten")
    
})

#### calcObservedMutations, no mutations ####
test_that("calcObservedMutations, when no mutations found", {
    
    in_seq <- ExampleDb[1, "SEQUENCE_IMGT"]
    germ_seq <- in_seq
    
    #' Should return c(NA,NA), not c(NA)
    expect_equivalent(calcObservedMutations(in_seq, germ_seq, regionDefinition = NULL), c(NA, NA))
    expect_equivalent(calcObservedMutations(in_seq, germ_seq, regionDefinition = IMGT_V), c(NA, NA, NA, NA))
    
    inputSeq <- ".......................GGGA...GGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTAGCTACGACATGCACTGGGTCCGCCAAGCTACAGGAAAAGGTCTGGAGTGGGTCTCAGCTATTGGTACTGCT.........GGTGACACATACTATCCAGGCTCCGTGAAG...GGCCGATTCACCATCTCCAGAGAAAATGCCAAGAACTCCTTGTATCTTCAAATGAACAGCCTGAGAGCCGGGGACACGGCTGTGTATTACTGTGCAAGAGATAAGGACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG"
    germlinSeq <- "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTAGCTACGACATGCACTGGGTCCGCCAAGCTACAGGAAAAGGTCTGGAGTGGGTCTCAGCTATTGGTACTGCT.........GGTGACACATACTATCCAGGCTCCGTGAAG...GGCCGATTCACCATCTCCAGAGAAAATGCCAAGAACTCCTTGTATCTTCAAATGAACAGCCTGAGAGCCGGGGACACGGCTGTGTATTACTGTGCAAGAGANNNNNACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG"
    expect_equivalent(calcObservedMutations(in_seq, germ_seq, frequency=TRUE, regionDefinition = NULL, mutationDefinition = NULL, returnRaw = F),
                      c(NA, NA))
    
})

#### observedMutations, no mutations ####
test_that("observedMutations, when no mutations found", {
    
    in_seq <- ExampleDb[["SEQUENCE_IMGT"]][1]
    germ_seq <- in_seq
    
    expect_equivalent(observedMutations(data.frame(
        "SEQUENCE_IMGT"=c(in_seq, in_seq),
        "GERMLINE_IMGT_D_MASK"=c(in_seq, in_seq)
    ), regionDefinition = NULL)[,c("MU_COUNT_SEQ_R", "MU_COUNT_SEQ_S")], data.frame(c(0, 0), c(0, 0)))
})

##### Tests 1A-1H (calcObservedMutations & observedMutations)
##### Tests 2A-2E (calcClonalConsensusHelper, calcClonalConsensus, collapseClones) 
##### are for changed made during commits pushed during June 2-June 12 2017

#### calcObservedMutations 1A ####
test_that("calcObservedMutations, 1A, without ambiguous characters, length is multiple of 3", {
    # 7 codons exactly
    #set.seed(1835)
    #seqinr::c2s(sample(x=shazam:::NUCLEOTIDES[1:4], size=21, replace=TRUE))
    
    #obsv: "TAT ATA ATC -GT CAG CTC TCG" 
    #germ: "TAT TAT ATA GGT CTT CNC AAC" 
    #region  W   Y   Y   W   Y   W   W 
    #codon   1   2   3   4   5   6   7
    
    # 1st codon: 0 mutation: (1 na because of no change)
    # 2nd codon: 3 mutations: TAT->AAT (R, 4); TAT->TTT (R, 5); TAT->TAA (Stop)
    # 3rd codon: 1 mutation: ATA->ATC (S, 9)
    # 4th codon: 0 mutation: GGT->-GT (na because of - in -GT)
    # 5th codon: 2 mutation: CTT->CAT (R, 14); CTT->CTG (S, 15) 
    # 6th codon: 0 mutation: CNC->CTC (na because of N in CNC)
    # 7th codon: 3 mutations: AAC->TAC (R, 19); AAC->ACC (R, 20); AAC->AAG (R, 21)
    
    obsv = "TATATAATC-GTCAGCTCTCG" 
    germ = "TATTATATAGGTCTTCNCAAC"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor(rep(c("W","Y","Y","W","Y","W","W"), each=3)))
    
    ##### both uppercase and lowercase input should work
    upperResult = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                        frequency=FALSE, returnRaw=TRUE, 
                                        regionDefinition=NULL)
    lowerResult = calcObservedMutations(inputSeq=tolower(obsv), germlineSeq=tolower(germ), 
                                        frequency=FALSE, returnRaw=TRUE, 
                                        regionDefinition=NULL)
    expect_equal(upperResult, lowerResult)
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=TRUE, 
                                                regionDefinition=NULL)
    # based on annotations above 
    exp.raw.noRegDef.position = c(4, 5, 9, 14, 15, 19, 20, 21)
    exp.raw.noRegDef.R = c(1,1,0,1,0,1,1,1)
    exp.raw.noRegDef.S = c(0,0,1,0,1,0,0,0)
    
    expect_equal(freqF.rawT.noRegDef$pos$position, exp.raw.noRegDef.position)
    expect_equal(freqF.rawT.noRegDef$pos$R, exp.raw.noRegDef.R)
    expect_equal(freqF.rawT.noRegDef$pos$S, exp.raw.noRegDef.S)
    expect_equal(freqF.rawT.noRegDef$pos$region, rep("SEQ", length(exp.raw.noRegDef.position)))
    # $nonN is named; use expect_equivalent instead of expect_equal
    exp.noRegDef.nonN.Dash.Dot = sum( seqinr::s2c(obsv) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equivalent(freqF.rawT.noRegDef$nonN, exp.noRegDef.nonN.Dash.Dot)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=FALSE, 
                                                regionDefinition=NULL)
    expect_equal(freqF.rawF.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R), 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)))
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=FALSE, 
                                                regionDefinition=NULL)
    freqT.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=TRUE, 
                                                regionDefinition=NULL)
    
    expect_equal(freqT.rawF.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    expect_equal(freqT.rawT.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=TRUE, 
                                              regionDefinition=regDef)
    exp.raw.regDef.region = c("Y","Y","Y","Y","Y","W","W","W")
    expect_equal(freqF.rawT.regDef$pos$region, exp.raw.regDef.region)
    # counts should be the same
    expect_identical(freqF.rawT.regDef$pos[, 1:3], freqF.rawT.noRegDef$pos[, 1:3])
    
    exp.regDef.W.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    exp.regDef.Y.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[regDef@boundaries=="Y"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[regDef@boundaries=="Y"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equal(freqF.rawT.regDef$nonN, c("W"=exp.regDef.W.nonN.Dash.Dot, "Y"=exp.regDef.Y.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    
    exp.raw.regDef.W.R = sum(subset(freqF.rawT.regDef$pos, region=="W")$R)
    exp.raw.regDef.W.S = sum(subset(freqF.rawT.regDef$pos, region=="W")$S)
    exp.raw.regDef.Y.R = sum(subset(freqF.rawT.regDef$pos, region=="Y")$R)
    exp.raw.regDef.Y.S = sum(subset(freqF.rawT.regDef$pos, region=="Y")$S)
    
    expect_equal(freqF.rawF.regDef, c("W_R"=exp.raw.regDef.W.R,
                                      "W_S"=exp.raw.regDef.W.S,
                                      "Y_R"=exp.raw.regDef.Y.R,
                                      "Y_S"=exp.raw.regDef.Y.S))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    freqT.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=TRUE, 
                                              regionDefinition=regDef)
    expect_equal(freqT.rawF.regDef, c("W_R"=exp.raw.regDef.W.R/exp.regDef.W.nonN.Dash.Dot,
                                      "W_S"=exp.raw.regDef.W.S/exp.regDef.W.nonN.Dash.Dot,
                                      "Y_R"=exp.raw.regDef.Y.R/exp.regDef.Y.nonN.Dash.Dot,
                                      "Y_S"=exp.raw.regDef.Y.S/exp.regDef.Y.nonN.Dash.Dot))
    
})

#### calcObservedMutations 1B ####
test_that("calcObservedMutations, 1B, with ambiguous characters, length is multiple of 3", {
    #obsv:   "TAT ATA WSC -GT CDG CTC TCG" 
    #germ:   "TAT TAT ATA GGT CTT CNC AAC" 
    #region    W   Y   Y   W   Y   W   W 
    #codon     1   2   3   4   5   6   7
    #R+S mut      **  ***      **     ***
    
    # 1st codon: 0 mutation: (1 na because of no change)
    # 2nd codon: 3 mutations: TAT->AAT (R, 4); TAT->TTT (R, 5); TAT->TAA (Stop)
    # 3rd codon: 4 mutations: 
    # ATA->WTA => ATA->ATA (na), ATA->TTA (R, 7)
    # ATA->ASA => ATA->ACA (R, 8), ATA->AGA (R, 8)
    # ATA->ATC => ATA->ATC (S, 9)
    # 4th codon: 0 mutation: GGT->-GT (na because of - in -GT)
    # 5th codon: 3 mutations: 
    # CTT->CDT => CTT->CAT (R, 14), CTT->CGT (R, 14), CTT->CTT (na)
    # CTT->CTG => CTT->CTG (S, 15)
    # 6th codon: 0 mutation: CNC->CTC (na because of N in CNC)
    # 7th codon: 3 mutations: AAC->TAC (R, 19); AAC->ACC (R, 20); AAC->AAG (R, 21)
    
    obsv = "TATATAWSC-GTCDGCTCTCG" 
    germ = "TATTATATAGGTCTTCNCAAC"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor(rep(c("W","Y","Y","W","Y","W","W"), each=3)))
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=TRUE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    # based on annotations above 
    exp.raw.noRegDef.position = c(4, 5, 7, 8, 9, 14, 15, 19, 20, 21)
    exp.raw.noRegDef.R = c(1,1,1,2,0,2,0,1,1,1)
    exp.raw.noRegDef.S = c(0,0,0,0,1,0,1,0,0,0)
    
    expect_equal(freqF.rawT.noRegDef$pos$position, exp.raw.noRegDef.position)
    expect_equal(freqF.rawT.noRegDef$pos$R, exp.raw.noRegDef.R)
    expect_equal(freqF.rawT.noRegDef$pos$S, exp.raw.noRegDef.S)
    expect_equal(freqF.rawT.noRegDef$pos$region, rep("SEQ", length(exp.raw.noRegDef.position)))
    # $nonN is named; use expect_equivalent instead of expect_equal
    # 8 non-N positions w/o mutation
    # 3,5,4,3 combinations after expanding ambiguous codons at codons 2, 3, 5, 7 resp.
    exp.noRegDef.nonN.Dash.Dot = 8 + sum(c(3,5,4,3))
    expect_equivalent(freqF.rawT.noRegDef$nonN, exp.noRegDef.nonN.Dash.Dot)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=FALSE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    expect_equal(freqF.rawF.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R), 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)))
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=FALSE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    freqT.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=TRUE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    
    expect_equal(freqT.rawF.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    expect_equal(freqT.rawT.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=TRUE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    exp.raw.regDef.region = c("Y","Y","Y","Y","Y","Y","Y","W","W","W")
    expect_equal(freqF.rawT.regDef$pos$region, exp.raw.regDef.region)
    # counts should be the same
    expect_identical(freqF.rawT.regDef$pos[, 1:3], freqF.rawT.noRegDef$pos[, 1:3])
    
    # 7 non-N positions w/o mutation in region "W"
    # 3 combinations after expanding ambiguous codons at codon 7 
    exp.regDef.W.nonN.Dash.Dot = 7 + 3
    # 1 non-N positions w/o mutation in region "Y"
    # 3,5,4 combinations after expanding ambiguous codons at codons 2, 3, 5 resp.
    exp.regDef.Y.nonN.Dash.Dot = 1 + sum(c(3,5,4))
    
    expect_equal(freqF.rawT.regDef$nonN, c("W"=exp.regDef.W.nonN.Dash.Dot, "Y"=exp.regDef.Y.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=FALSE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    
    exp.raw.regDef.W.R = sum(subset(freqF.rawT.regDef$pos, region=="W")$R)
    exp.raw.regDef.W.S = sum(subset(freqF.rawT.regDef$pos, region=="W")$S)
    exp.raw.regDef.Y.R = sum(subset(freqF.rawT.regDef$pos, region=="Y")$R)
    exp.raw.regDef.Y.S = sum(subset(freqF.rawT.regDef$pos, region=="Y")$S)
    
    expect_equal(freqF.rawF.regDef, c("W_R"=exp.raw.regDef.W.R,
                                      "W_S"=exp.raw.regDef.W.S,
                                      "Y_R"=exp.raw.regDef.Y.R,
                                      "Y_S"=exp.raw.regDef.Y.S))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=FALSE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    freqT.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=TRUE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    expect_equal(freqT.rawF.regDef, c("W_R"=exp.raw.regDef.W.R/exp.regDef.W.nonN.Dash.Dot,
                                      "W_S"=exp.raw.regDef.W.S/exp.regDef.W.nonN.Dash.Dot,
                                      "Y_R"=exp.raw.regDef.Y.R/exp.regDef.Y.nonN.Dash.Dot,
                                      "Y_S"=exp.raw.regDef.Y.S/exp.regDef.Y.nonN.Dash.Dot))
    
})

#### calcObservedMutations 1C ####
test_that("calcObservedMutations, 1C, without ambiguous characters, length is not multiple of 3", {
    # 6 codons + 1 two-nucleotide overhang
    # non-triplet overhang should be ignored
    
    #obsv: "TAT ATA ATC -GT CAG CTC TC" 
    #germ: "TAT TAT ATA GGT CTT CNC AA" 
    #region  W   Y   Y   W   Y   W   W 
    #codon   1   2   3   4   5   6   7
    
    # 1st codon: 0 mutation: (1 na because of no change)
    # 2nd codon: 3 mutations: TAT->AAT (R, 4); TAT->TTT (R, 5); TAT->TAA (Stop)
    # 3rd codon: 1 mutation: ATA->ATC (S, 9)
    # 4th codon: 0 mutation: GGT->-GT (na because of - in -GT)
    # 5th codon: 2 mutation: CTT->CAT (R, 14); CTT->CTG (S, 15) 
    # 6th codon: 0 mutation: CNC->CTC (na because of N in CNC)
    # 7th codon: 0 mutation: non-triplet overhang ignored
    
    obsv = "TATATAATC-GTCAGCTCTC" 
    germ = "TATTATATAGGTCTTCNCAA"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor(c( rep(c("W","Y","Y","W","Y","W"), each=3),
                                                         "W", "W" ))) # hard-coded overhang
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=TRUE, 
                                                regionDefinition=NULL)
    # based on annotations above 
    exp.raw.noRegDef.position = c(4, 5, 9, 14, 15)
    exp.raw.noRegDef.R = c(1,1,0,1,0)
    exp.raw.noRegDef.S = c(0,0,1,0,1)
    
    expect_equal(freqF.rawT.noRegDef$pos$position, exp.raw.noRegDef.position)
    expect_equal(freqF.rawT.noRegDef$pos$R, exp.raw.noRegDef.R)
    expect_equal(freqF.rawT.noRegDef$pos$S, exp.raw.noRegDef.S)
    expect_equal(freqF.rawT.noRegDef$pos$region, rep("SEQ", length(exp.raw.noRegDef.position)))
    # $nonN is named; use expect_equivalent instead of expect_equal
    # hard-coded; length is not multiple of 3; non-triplet overhang ignored
    exp.noRegDef.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[1:18] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[1:18] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equivalent(freqF.rawT.noRegDef$nonN, exp.noRegDef.nonN.Dash.Dot)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=FALSE, 
                                                regionDefinition=NULL)
    expect_equal(freqF.rawF.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R), 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)))
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=FALSE, 
                                                regionDefinition=NULL)
    freqT.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=TRUE, 
                                                regionDefinition=NULL)
    
    expect_equal(freqT.rawF.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    expect_equal(freqT.rawT.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=TRUE, 
                                              regionDefinition=regDef)
    exp.raw.regDef.region = c("Y","Y","Y","Y","Y")
    expect_equal(freqF.rawT.regDef$pos$region, exp.raw.regDef.region)
    # counts should be the same
    expect_identical(freqF.rawT.regDef$pos[, 1:3], freqF.rawT.noRegDef$pos[, 1:3])
    
    # hard-coded; length is not multiple of 3; non-triplet overhang ignored
    exp.regDef.W.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[1:18][regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[1:18][regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    exp.regDef.Y.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[1:18][regDef@boundaries=="Y"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[1:18][regDef@boundaries=="Y"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equal(freqF.rawT.regDef$nonN, c("W"=exp.regDef.W.nonN.Dash.Dot, "Y"=exp.regDef.Y.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    
    exp.raw.regDef.W.R = sum(subset(freqF.rawT.regDef$pos, region=="W")$R)
    exp.raw.regDef.W.S = sum(subset(freqF.rawT.regDef$pos, region=="W")$S)
    exp.raw.regDef.Y.R = sum(subset(freqF.rawT.regDef$pos, region=="Y")$R)
    exp.raw.regDef.Y.S = sum(subset(freqF.rawT.regDef$pos, region=="Y")$S)
    
    expect_equal(freqF.rawF.regDef, c("W_R"=exp.raw.regDef.W.R,
                                      "W_S"=exp.raw.regDef.W.S,
                                      "Y_R"=exp.raw.regDef.Y.R,
                                      "Y_S"=exp.raw.regDef.Y.S))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    freqT.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=TRUE, 
                                              regionDefinition=regDef)
    expect_equal(freqT.rawF.regDef, c("W_R"=exp.raw.regDef.W.R/exp.regDef.W.nonN.Dash.Dot,
                                      "W_S"=exp.raw.regDef.W.S/exp.regDef.W.nonN.Dash.Dot,
                                      "Y_R"=exp.raw.regDef.Y.R/exp.regDef.Y.nonN.Dash.Dot,
                                      "Y_S"=exp.raw.regDef.Y.S/exp.regDef.Y.nonN.Dash.Dot))
    
})

#### calcObservedMutations 1D ####
test_that("calcObservedMutations, 1D, without ambiguous characters, only 1 codon, no mutation", {
    # 1 codon only
    
    ##### no mutation
    
    #obsv: "TAT" 
    #germ: "TAT" 
    #region  W 
    #codon   1 
    
    # 1st codon: 0 mutation
    
    obsv = "TAT" 
    germ = "TAT"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor( rep("W", each=3) ))
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=TRUE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    
    expect_equal(freqF.rawT.noRegDef$pos, NA)
    # because there's no ambiguous char, expect "and" to produce same result as "eitherOr"
    # test with result produced using "eitherOr" (default for when there's no ambiguous char)
    exp.noRegDef.nonN.Dash.Dot = sum( seqinr::s2c(obsv) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equivalent(freqF.rawT.noRegDef$nonN, exp.noRegDef.nonN.Dash.Dot)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=FALSE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    expect_equal(freqF.rawF.noRegDef, c("SEQ_R"=NA, "SEQ_S"=NA))
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=FALSE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    freqT.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=TRUE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    
    expect_equal(freqT.rawF.noRegDef, c("SEQ_R"=NA, "SEQ_S"=NA))
    expect_equal(freqT.rawT.noRegDef, c("SEQ_R"=NA, "SEQ_S"=NA))
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=TRUE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    
    expect_equal(freqF.rawT.regDef$pos, NA)
    
    # because there's no ambiguous char, expect "and" to produce same result as "eitherOr"
    # test with result produced using "eitherOr" (default for when there's no ambiguous char)
    exp.regDef.W.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equal(freqF.rawT.regDef$nonN, c("W"=exp.regDef.W.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=FALSE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    
    expect_equal(freqF.rawF.regDef, c("W_R"=NA, "W_S"=NA))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=FALSE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    freqT.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=TRUE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    expect_equal(freqT.rawF.regDef, c("W_R"=NA, "W_S"=NA))
    
})

#### calcObservedMutations 1E ####
test_that("calcObservedMutations, 1E, with ambiguous characters, only 1 codon, with mutations", {
    # 1 codon only
    
    ##### 2 mutations
    
    #obsv: "THT" 
    #germ: "TAT" 
    #region  W 
    #codon   1 
    
    # 1st codon: 2 mutation:
    # TAT->THT => TAT->TAT (na); TAT->TCT (R, 2); TAT->TTT (R, 2)
    
    obsv = "THT" 
    germ = "TAT"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor( rep("W", each=3) ))
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=TRUE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    
    # based on annotations above 
    exp.raw.noRegDef.position = 2
    exp.raw.noRegDef.R = 2
    exp.raw.noRegDef.S = 0
    
    expect_equal(freqF.rawT.noRegDef$pos$position, exp.raw.noRegDef.position)
    expect_equal(freqF.rawT.noRegDef$pos$R, exp.raw.noRegDef.R)
    expect_equal(freqF.rawT.noRegDef$pos$S, exp.raw.noRegDef.S)
    expect_equal(freqF.rawT.noRegDef$pos$region, rep("SEQ", length(exp.raw.noRegDef.position)))
    
    # 2 non-N positions w/o mutation
    # 3 combinations after expanding ambiguous codon at codon 1
    exp.noRegDef.nonN.Dash.Dot = 2 + 3
    expect_equivalent(freqF.rawT.noRegDef$nonN, exp.noRegDef.nonN.Dash.Dot)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=FALSE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    expect_equal(freqF.rawF.noRegDef, c("SEQ_R"=exp.raw.noRegDef.R, 
                                        "SEQ_S"=exp.raw.noRegDef.S))
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=FALSE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    freqT.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=TRUE, 
                                                regionDefinition=NULL, ambiguousMode="and")
    
    expect_equal(freqT.rawF.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    expect_equal(freqT.rawT.noRegDef, c("SEQ_R"=sum(exp.raw.noRegDef.R)/exp.noRegDef.nonN.Dash.Dot, 
                                        "SEQ_S"=sum(exp.raw.noRegDef.S)/exp.noRegDef.nonN.Dash.Dot))
    
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=TRUE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    
    exp.raw.regDef.region = "W"
    expect_equal(freqF.rawT.regDef$pos$region, exp.raw.regDef.region)
    # counts should be the same
    expect_identical(freqF.rawT.regDef$pos[, 1:3], freqF.rawT.noRegDef$pos[, 1:3])
    
    # 2 non-N positions w/o mutation in region "W"
    # 3 combinations after expanding ambiguous codon at codon 1
    exp.regDef.W.nonN.Dash.Dot = 2 + 3 
    expect_equal(freqF.rawT.regDef$nonN, c("W"=exp.regDef.W.nonN.Dash.Dot))
    
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=FALSE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    
    exp.raw.regDef.W.R = sum(subset(freqF.rawT.regDef$pos, region=="W")$R)
    exp.raw.regDef.W.S = sum(subset(freqF.rawT.regDef$pos, region=="W")$S)
    
    expect_equal(freqF.rawF.regDef, c("W_R"=exp.raw.regDef.W.R,
                                      "W_S"=exp.raw.regDef.W.S))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=FALSE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    freqT.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=TRUE, 
                                              regionDefinition=regDef, ambiguousMode="and")
    
    expect_equal(freqT.rawF.regDef, c("W_R"=exp.raw.regDef.W.R/exp.regDef.W.nonN.Dash.Dot,
                                      "W_S"=exp.raw.regDef.W.S/exp.regDef.W.nonN.Dash.Dot))
    
})

#### calcObservedMutations 1F ####
test_that("calcObservedMutations, 1F, less than 1 codon", {
    # non-triplet overhang should be ignored
    # in this case, NA
    
    #obsv: "TA" 
    #germ: "TA" 
    #region  W 
    #codon   1 
    
    obsv = "TA" 
    germ = "TA"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor( rep("W", each=2) ))
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=TRUE, 
                                                regionDefinition=NULL)
    
    expect_equal(freqF.rawT.noRegDef$pos, NA)
    
    # hard-coded; non-triplet overhang ignored
    exp.noRegDef.nonN.Dash.Dot = NA
    expect_equivalent(freqF.rawT.noRegDef$nonN, exp.noRegDef.nonN.Dash.Dot)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=FALSE, 
                                                regionDefinition=NULL)
    expect_equal(freqF.rawF.noRegDef, c("SEQ_R"=NA, "SEQ_S"=NA))
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=FALSE, 
                                                regionDefinition=NULL)
    freqT.rawT.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=TRUE, returnRaw=TRUE, 
                                                regionDefinition=NULL)
    
    expect_equal(freqT.rawF.noRegDef, c("SEQ_R"=NA, "SEQ_S"=NA))
    expect_equal(freqT.rawT.noRegDef, c("SEQ_R"=NA, "SEQ_S"=NA))
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=TRUE, 
                                              regionDefinition=regDef)
    
    expect_equal(freqF.rawT.regDef$pos, NA)
    
    # hard-coded; non-triplet overhang ignored
    exp.regDef.W.nonN.Dash.Dot = NA
    expect_equal(freqF.rawT.regDef$nonN, c("W"=exp.regDef.W.nonN.Dash.Dot))
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    
    expect_equal(freqF.rawF.regDef, c("W_R"=NA, "W_S"=NA))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    freqT.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=TRUE, 
                                              regionDefinition=regDef)
    expect_equal(freqT.rawF.regDef, c("W_R"=NA, "W_S"=NA))
    
})

#### calcObservedMutations 1G ####
test_that("calcObservedMutations, 1G, with ambiguous characters in germline", {
    
    #obsv:   THS AGM SAA C  C  K  A  A  G
    #germ:   TMC RGC SAT C  C  K  A  T  G
    #pos     123 456 789 10 11 12 13 14 15
    #R+Smut   ** * *   *             *
    #reg     WWW YYY PPP Y  Y  Y  P  P  P
    
    # non-N, non-dash, non-dot positions w/o mutations: 9
    
    # pos 2
    # "and": R=4; S=0; nonN=6
    # "eitherOr": na
    # TMC->THC =>
    # (M)  (H)  
    # TAC->TAC: na
    # TAC->TCC: R
    # TAC->TTC: R
    # TCC->TAC: R
    # TCC->TCC: na
    # TCC->TTC: R
    
    # pos 3
    # "and": R=3; S=1; nonN=8
    # "eitherOr": na
    # TMC->TMS =>
    # (M)  (MS)
    # TAC->TAC: na
    # TAC->TAG: Stop
    # TAC->TCC: R
    # TAC->TCG: R
    # TCC->TAC: R
    # TCC->TAG: stop
    # TCC->TCC: na
    # TCC->TCG: S
    
    # pos 4
    # "and": R=1; S=0; nonN=2
    # "eitherOr": na
    # RGC->AGC =>
    #(R)
    # AGC->AGC: na
    # GGC->AGC: R
    
    # pos 6
    # "and": R=5; S=1; nonN=8
    # "eitherOr": na
    # RGC->RGM =>
    #(R)  (R M)
    # AGC->AGA: R
    # AGC->AGC: na
    # AGC->GGA: R
    # AGC->GGC: R
    # GGC->AGA: R
    # GGC->AGC: R
    # GGC->GGA: S
    # GGC->GGC: na
    
    # pos 9
    # "and": R=4; S=0; nonN=4
    # "eitherOr": R
    # SAT->SAA =>
    #(S)  (S)
    # CAT->CAA: R
    # CAT->GAA: R
    # GAT->CAA: R
    # GAT->GAA: R
    
    # pos 14
    # "and": R=1; S=0; nonN=1
    # "either": R=1; S=0
    # ATG->AAG: R
    
    obsv = "THSAGMSAACCKAAG" 
    germ = "TMCRGCSATCCKATG"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor( rep(c("W", "Y", "P", "Y", "P"), each=3) ))
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    freqF.rawT.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    # based on annotations above 
    exp.raw.noRegDef.position.and = c(2,3,4,6,9,14)
    exp.raw.noRegDef.R.and = c(4,3,1,5,4,1)
    exp.raw.noRegDef.S.and = c(0,1,0,1,0,0)
    
    exp.raw.noRegDef.position.eor = c(9,14)
    exp.raw.noRegDef.R.eor = c(1,1)
    exp.raw.noRegDef.S.eor = c(0,0)
    
    expect_equal(freqF.rawT.noRegDef.and$pos$position, exp.raw.noRegDef.position.and)
    expect_equal(freqF.rawT.noRegDef.and$pos$R, exp.raw.noRegDef.R.and)
    expect_equal(freqF.rawT.noRegDef.and$pos$S, exp.raw.noRegDef.S.and)
    expect_equal(freqF.rawT.noRegDef.and$pos$region, rep("SEQ", length(exp.raw.noRegDef.position.and)))
    
    expect_equal(freqF.rawT.noRegDef.eor$pos$position, exp.raw.noRegDef.position.eor)
    expect_equal(freqF.rawT.noRegDef.eor$pos$R, exp.raw.noRegDef.R.eor)
    expect_equal(freqF.rawT.noRegDef.eor$pos$S, exp.raw.noRegDef.S.eor)
    expect_equal(freqF.rawT.noRegDef.eor$pos$region, rep("SEQ", length(exp.raw.noRegDef.position.eor)))
    
    # $nonN is named; use expect_equivalent instead of expect_equal
    exp.noRegDef.nonN.Dash.Dot.and = 9 + sum(c(6,8,2,8,4,1)) # 9 non-N positions w/o mutations
    expect_equivalent(freqF.rawT.noRegDef.and$nonN, exp.noRegDef.nonN.Dash.Dot.and)
    
    exp.noRegDef.nonN.Dash.Dot.eor = sum( seqinr::s2c(obsv) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                              seqinr::s2c(germ) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equivalent(freqF.rawT.noRegDef.eor$nonN, exp.noRegDef.nonN.Dash.Dot.eor)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    expect_equal(freqF.rawF.noRegDef.and, c("SEQ_R"=sum(exp.raw.noRegDef.R.and), 
                                            "SEQ_S"=sum(exp.raw.noRegDef.S.and)))
    
    freqF.rawF.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    expect_equal(freqF.rawF.noRegDef.eor, c("SEQ_R"=sum(exp.raw.noRegDef.R.eor), 
                                            "SEQ_S"=sum(exp.raw.noRegDef.S.eor)))
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    freqT.rawT.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    
    freqT.rawF.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    freqT.rawT.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    
    
    expect_equal(freqT.rawF.noRegDef.and, c("SEQ_R"=sum(exp.raw.noRegDef.R.and)/exp.noRegDef.nonN.Dash.Dot.and, 
                                            "SEQ_S"=sum(exp.raw.noRegDef.S.and)/exp.noRegDef.nonN.Dash.Dot.and))
    expect_equal(freqT.rawT.noRegDef.and, c("SEQ_R"=sum(exp.raw.noRegDef.R.and)/exp.noRegDef.nonN.Dash.Dot.and, 
                                            "SEQ_S"=sum(exp.raw.noRegDef.S.and)/exp.noRegDef.nonN.Dash.Dot.and))
    
    expect_equal(freqT.rawF.noRegDef.eor, c("SEQ_R"=sum(exp.raw.noRegDef.R.eor)/exp.noRegDef.nonN.Dash.Dot.eor, 
                                            "SEQ_S"=sum(exp.raw.noRegDef.S.eor)/exp.noRegDef.nonN.Dash.Dot.eor))
    expect_equal(freqT.rawT.noRegDef.eor, c("SEQ_R"=sum(exp.raw.noRegDef.R.eor)/exp.noRegDef.nonN.Dash.Dot.eor, 
                                            "SEQ_S"=sum(exp.raw.noRegDef.S.eor)/exp.noRegDef.nonN.Dash.Dot.eor))
    
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    exp.raw.regDef.region.and = c("W","W","Y","Y","P","P")
    expect_equal(freqF.rawT.regDef.and$pos$region, exp.raw.regDef.region.and)
    
    freqF.rawT.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    exp.raw.regDef.region.eor = c("P","P")
    expect_equal(freqF.rawT.regDef.eor$pos$region, exp.raw.regDef.region.eor)
    
    # counts should be the same
    expect_identical(freqF.rawT.regDef.and$pos[, 1:3], freqF.rawT.noRegDef.and$pos[, 1:3])
    expect_identical(freqF.rawT.regDef.eor$pos[, 1:3], freqF.rawT.noRegDef.eor$pos[, 1:3])
    
    expect_equal(freqF.rawT.regDef.and$nonN, c("P"=sum(c(4,1,4)), # 4 non-N P positions w/o mutation
                                               "W"=sum(c(1,6,8)), # 1 non-N P positions w/o mutation
                                               "Y"=sum(c(2,8,4)))) # 2 non-N P positions w/o mutation
    
    exp.regDef.W.nonN.Dash.Dot.eor = sum( seqinr::s2c(obsv)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                              seqinr::s2c(germ)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    exp.regDef.Y.nonN.Dash.Dot.eor = sum( seqinr::s2c(obsv)[regDef@boundaries=="Y"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                              seqinr::s2c(germ)[regDef@boundaries=="Y"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    exp.regDef.P.nonN.Dash.Dot.eor = sum( seqinr::s2c(obsv)[regDef@boundaries=="P"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                              seqinr::s2c(germ)[regDef@boundaries=="P"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equal(freqF.rawT.regDef.eor$nonN, c("P"=exp.regDef.P.nonN.Dash.Dot.eor, 
                                               "W"=exp.regDef.W.nonN.Dash.Dot.eor, 
                                               "Y"=exp.regDef.Y.nonN.Dash.Dot.eor))
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    freqF.rawF.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    
    exp.raw.regDef.W.R.and = sum(subset(freqF.rawT.regDef.and$pos, region=="W")$R)
    exp.raw.regDef.W.S.and = sum(subset(freqF.rawT.regDef.and$pos, region=="W")$S)
    exp.raw.regDef.Y.R.and = sum(subset(freqF.rawT.regDef.and$pos, region=="Y")$R)
    exp.raw.regDef.Y.S.and = sum(subset(freqF.rawT.regDef.and$pos, region=="Y")$S)
    exp.raw.regDef.P.R.and = sum(subset(freqF.rawT.regDef.and$pos, region=="P")$R)
    exp.raw.regDef.P.S.and = sum(subset(freqF.rawT.regDef.and$pos, region=="P")$S)
    
    expect_equal(freqF.rawF.regDef.and, c("P_R"=exp.raw.regDef.P.R.and,
                                          "P_S"=exp.raw.regDef.P.S.and,
                                          "W_R"=exp.raw.regDef.W.R.and,
                                          "W_S"=exp.raw.regDef.W.S.and,
                                          "Y_R"=exp.raw.regDef.Y.R.and,
                                          "Y_S"=exp.raw.regDef.Y.S.and))
    
    exp.raw.regDef.W.R.eor = sum(subset(freqF.rawT.regDef.eor$pos, region=="W")$R)
    exp.raw.regDef.W.S.eor = sum(subset(freqF.rawT.regDef.eor$pos, region=="W")$S)
    exp.raw.regDef.Y.R.eor = sum(subset(freqF.rawT.regDef.eor$pos, region=="Y")$R)
    exp.raw.regDef.Y.S.eor = sum(subset(freqF.rawT.regDef.eor$pos, region=="Y")$S)
    exp.raw.regDef.P.R.eor = sum(subset(freqF.rawT.regDef.eor$pos, region=="P")$R)
    exp.raw.regDef.P.S.eor = sum(subset(freqF.rawT.regDef.eor$pos, region=="P")$S)
    
    expect_equal(freqF.rawF.regDef.eor, c("P_R"=exp.raw.regDef.P.R.eor,
                                          "P_S"=exp.raw.regDef.P.S.eor,
                                          "W_R"=exp.raw.regDef.W.R.eor,
                                          "W_S"=exp.raw.regDef.W.S.eor,
                                          "Y_R"=exp.raw.regDef.Y.R.eor,
                                          "Y_S"=exp.raw.regDef.Y.S.eor))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    freqT.rawT.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    
    freqT.rawF.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    freqT.rawT.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    
    expect_equal(freqT.rawF.regDef.eor, freqT.rawT.regDef.eor)
    expect_equal(freqT.rawF.regDef.and, freqT.rawT.regDef.and)
    
    expect_equivalent(freqT.rawF.regDef.and, c("P_R"=exp.raw.regDef.P.R.and/freqF.rawT.regDef.and$nonN["P"],
                                               "P_S"=exp.raw.regDef.P.S.and/freqF.rawT.regDef.and$nonN["P"],
                                               "W_R"=exp.raw.regDef.W.R.and/freqF.rawT.regDef.and$nonN["W"],
                                               "W_S"=exp.raw.regDef.W.S.and/freqF.rawT.regDef.and$nonN["W"],
                                               "Y_R"=exp.raw.regDef.Y.R.and/freqF.rawT.regDef.and$nonN["Y"],
                                               "Y_S"=exp.raw.regDef.Y.S.and/freqF.rawT.regDef.and$nonN["Y"]))
    
})

#### calcObservedMutations 1H ####
test_that("calcObservedMutations, 1H, eitherOr vs. and when there is no ambiguous chars", {
    # expect result from "eitherOr" and that from "and" to be the same
    
    obsv = "TATATAATC-GTCAGCTCTCG" 
    germ = "TATTATATAGGTCTTCNCAAC"
    
    regDef = createRegionDefinition(name="mock", 
                                    boundaries=factor(rep(c("W","Y","Y","W","Y","W","W"), each=3)))
    
    ##### freq=F, returnRaw=T, regDef=NULL
    # recall that only R+S mutations are recorded; na or Stop are dropped
    freqF.rawT.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    freqF.rawT.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    
    expect_equal(freqF.rawT.noRegDef.and, freqF.rawT.noRegDef.eor)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    freqF.rawF.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    
    freqF.rawF.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=FALSE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    expect_equal(freqF.rawF.noRegDef.and, freqF.rawF.noRegDef.eor)
    
    ##### freq=T, returnRaw=F/T, regDef=NULL
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    freqT.rawT.noRegDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="and")
    
    freqT.rawF.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=FALSE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    freqT.rawT.noRegDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                    frequency=TRUE, returnRaw=TRUE, 
                                                    regionDefinition=NULL, ambiguousMode="eitherOr")
    
    expect_equal(freqT.rawF.noRegDef.and, freqT.rawF.noRegDef.eor)
    expect_equal(freqT.rawT.noRegDef.and, freqT.rawT.noRegDef.eor)
    expect_equal(freqT.rawF.noRegDef.and, freqT.rawT.noRegDef.eor)
    
    ##### freq=F, returnRaw=T, with regDef
    freqF.rawT.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    
    freqF.rawT.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    
    expect_equal(freqF.rawT.regDef.and, freqF.rawT.regDef.eor)
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    freqF.rawF.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=FALSE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    
    expect_equal(freqF.rawF.regDef.and, freqF.rawF.regDef.eor)
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    freqT.rawT.regDef.eor = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="eitherOr")
    
    freqT.rawF.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=FALSE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    freqT.rawT.regDef.and = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                  frequency=TRUE, returnRaw=TRUE, 
                                                  regionDefinition=regDef, ambiguousMode="and")
    
    expect_equal(freqT.rawF.regDef.and, freqT.rawF.regDef.eor)
    expect_equal(freqT.rawT.regDef.and, freqT.rawT.regDef.eor)
    expect_equal(freqT.rawF.regDef.and, freqT.rawT.regDef.eor)
    
})

#### observedMutations 1I ####
test_that("observedMutations, 1I, using mock data from 1A through 1G, ambiguousMode=and", {
    # pull sequences from 1A-1F
    testDb = data.frame(obsv=c("TATATAATC-GTCAGCTCTCG", # 1A
                               "TATATAWSC-GTCDGCTCTCG", # 1B
                               "TATATAATC-GTCAGCTCTC",  # 1C
                               "TAT", # 1D
                               "THT", # 1E
                               "TA",  # 1F
                               "THSAGMSAACCKAAG"), # 1G
                        germ=c("TATTATATAGGTCTTCNCAAC", # 1A
                               "TATTATATAGGTCTTCNCAAC", # 1B
                               "TATTATATAGGTCTTCNCAA",  # 1C
                               "TAT", # 1D
                               "TAT", # 1E
                               "TA",  # 1F
                               "TMCRGCSATCCKATG"), # 1G
                        stringsAsFactors=FALSE)
 
    # pull results from 1A-1F
    # these are results from calcObservedMutations which have been tested
    exp.noRegDef.R = c(6, 10, 3, NA, 2, NA, 18)
    exp.noRegDef.S = c(2, 2, 2, NA, 0, NA, 2)
    exp.noRegDef.nonN = c(19, 23, 16, 3, 5, NA, 38)
    
    # convert NA to 0
    exp.noRegDef.R[is.na(exp.noRegDef.R)] = 0
    exp.noRegDef.S[is.na(exp.noRegDef.S)] = 0
    
    # frequency=F, combine=F, no regDef
    freqF.combF.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=F, combine=F, ambiguousMode="and")
    
    expect_equal(freqF.combF.noRegDef$MU_COUNT_SEQ_R, exp.noRegDef.R)
    expect_equal(freqF.combF.noRegDef$MU_COUNT_SEQ_S, exp.noRegDef.S)
    
    # frequency=F, combine=T, no regDef
    freqF.combT.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=F, combine=T, ambiguousMode="and")
    
    expect_equal(freqF.combT.noRegDef$MU_COUNT, 
                 exp.noRegDef.R+exp.noRegDef.S)
    
    # frequency=T, combine=F, no regDef
    freqT.combF.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=T, combine=F, ambiguousMode="and")
    
    # $nonN has NA; division gets NA; but observedMutations converts all NA to 0
    exp.noRegDef.R.freq = exp.noRegDef.R/exp.noRegDef.nonN
    exp.noRegDef.R.freq[is.na(exp.noRegDef.R.freq)] = 0
    
    exp.noRegDef.S.freq = exp.noRegDef.S/exp.noRegDef.nonN
    exp.noRegDef.S.freq[is.na(exp.noRegDef.S.freq)] = 0
    
    expect_equal(freqT.combF.noRegDef$MU_FREQ_SEQ_R, 
                 exp.noRegDef.R.freq)
    expect_equal(freqT.combF.noRegDef$MU_FREQ_SEQ_S, 
                 exp.noRegDef.S.freq)
    
    # frequency=T, combine=T, no regDef
    freqT.combT.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=T, combine=T, ambiguousMode="and")
    
    # $nonN has NA; division gets NA; but observedMutations converts all NA to 0
    exp.noRegDef.comb.freq = (exp.noRegDef.R+exp.noRegDef.S)/exp.noRegDef.nonN
    exp.noRegDef.comb.freq[is.na(exp.noRegDef.comb.freq)] = 0
    
    expect_equal(freqT.combT.noRegDef$MU_FREQ, 
                 exp.noRegDef.comb.freq)
    
})

#### observedMutations 1J ####
test_that("observedMutations, 1J, using mock data from 1A through 1G, ambiguousMode=eitherOr", {
    # pull sequences from 1A-1F
    testDb = data.frame(obsv=c("TATATAATC-GTCAGCTCTCG", # 1A
                               "TATATAWSC-GTCDGCTCTCG", # 1B ambiguous char
                               "TATATAATC-GTCAGCTCTC",  # 1C
                               "TAT", # 1D
                               "THT", # 1E ambiguous char
                               "TA",  # 1F
                               "THSAGMSAACCKAAG"), # 1G ambiguous char
                        germ=c("TATTATATAGGTCTTCNCAAC", # 1A
                               "TATTATATAGGTCTTCNCAAC", # 1B
                               "TATTATATAGGTCTTCNCAA",  # 1C
                               "TAT", # 1D
                               "TAT", # 1E ambiguous char
                               "TA",  # 1F
                               "TMCRGCSATCCKATG"), # 1G ambiguous char
                        stringsAsFactors=FALSE)
    
    # pull results from 1A-1F
    # these are results from calcObservedMutations which have been tested
    exp.noRegDef.R = c(6, 6, 3, NA, NA, NA, 2)
    exp.noRegDef.S = c(2, 2, 2, NA, NA, NA, 0)
    exp.noRegDef.nonN = c(19, 19, 16, 3, 3, NA, 15)
    
    # convert NA to 0
    exp.noRegDef.R[is.na(exp.noRegDef.R)] = 0
    exp.noRegDef.S[is.na(exp.noRegDef.S)] = 0
    
    # frequency=F, combine=F, no regDef
    freqF.combF.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=F, combine=F, ambiguousMode="eitherOr")
    
    expect_equal(freqF.combF.noRegDef$MU_COUNT_SEQ_R, exp.noRegDef.R)
    expect_equal(freqF.combF.noRegDef$MU_COUNT_SEQ_S, exp.noRegDef.S)
    
    # frequency=F, combine=T, no regDef
    freqF.combT.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=F, combine=T, ambiguousMode="eitherOr")
    
    expect_equal(freqF.combT.noRegDef$MU_COUNT, 
                 exp.noRegDef.R+exp.noRegDef.S)
    
    # frequency=T, combine=F, no regDef
    freqT.combF.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=T, combine=F, ambiguousMode="eitherOr")
    
    # $nonN has NA; division gets NA; but observedMutations converts all NA to 0
    exp.noRegDef.R.freq = exp.noRegDef.R/exp.noRegDef.nonN
    exp.noRegDef.R.freq[is.na(exp.noRegDef.R.freq)] = 0
    
    exp.noRegDef.S.freq = exp.noRegDef.S/exp.noRegDef.nonN
    exp.noRegDef.S.freq[is.na(exp.noRegDef.S.freq)] = 0
    
    expect_equal(freqT.combF.noRegDef$MU_FREQ_SEQ_R, 
                 exp.noRegDef.R.freq)
    expect_equal(freqT.combF.noRegDef$MU_FREQ_SEQ_S, 
                 exp.noRegDef.S.freq)
    
    # frequency=T, combine=T, no regDef
    freqT.combT.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=T, combine=T, ambiguousMode="eitherOr")
    
    # $nonN has NA; division gets NA; but observedMutations converts all NA to 0
    exp.noRegDef.comb.freq = (exp.noRegDef.R+exp.noRegDef.S)/exp.noRegDef.nonN
    exp.noRegDef.comb.freq[is.na(exp.noRegDef.comb.freq)] = 0
    
    expect_equal(freqT.combT.noRegDef$MU_FREQ, 
                 exp.noRegDef.comb.freq)
    
})

#### calcClonalConsensus 2A ####
test_that("calcClonalConsensusHelper, 2A, miscellaneous", {
    ##### only 1 seq
    seq1 = "ATGCATGCATGCA"
    # region def spanning nucleotides 1 through 12
    regDef1 = createRegionDefinition(boundaries=factor(rep(c("W", "Y"), each=6)))
    # no region def
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seq1),
                 list(cons=seq1, muFreq=NULL))
    # with region def
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seq1, 
                                                    lenLimit=regDef1@seqLength),
                 list(cons=substr(seq1, 1, regDef1@seqLength), muFreq=NULL))
    
    ##### multiple identical seqs
    # no region def
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=rep(seq1, 7)),
                 list(cons=seq1, muFreq=NULL))
    # with region def
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=rep(seq1, 7), 
                                                    lenLimit=regDef1@seqLength),
                 list(cons=substr(seq1, 1, regDef1@seqLength), muFreq=NULL))
})

#### calcClonalConsensus 2B ####
test_that("calcClonalConsensusHelper, 2B, methods = thresholdedFreq, mostCommon, catchAll", {
    # seq1: A T G C A T G C A T  -  G  .  N  T  C  
    # seq2: A T G G A T C G N T  -  A  .  G  N  C  G  C
    # seq3: A C T G A C T . A T  -  T  .  T  A  .  N 
    # seq4: T C G A A C C T A G  .  C  .  G  -  G  A  A  A
    # seq5: C G T A A T G - A A  -  A  .  C  G  A  T  T  A
    # pos   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
    
    # frequencies (by pos)
    # 1: A 0.6 T 0.2 C 0.2 | thresh=0.6: A | thresh=0.4: A | mostCommon: A | catchAll: H
    # 2: T 0.4 G 0.2 C 0.4 | thresh=0.6: N | thresh=0.4: TC=Y | mostCommon: TC=Y | catchAll: B
    # 3: T 0.4 G 0.6 | thresh=0.6: G | thresh=0.4: TG=K | mostCommon: G | catchAll: K
    # 4: A 0.4 G 0.4 C 0.2 | thresh=0.6: N | thresh=0.4: AG=R | mostCommon: AG=R | catchAll: V
    # 5: A 1 | thresh=0.6: A | thresh=0.4: A | mostCommon: A | catchAll: A
    # 6: T 0.6 C 0.4 | thresh=0.6: T | thresh=0.4: TC=Y | mostCommon: T | catchAll: Y
    # 7: T 0.2 G 0.4 C 0.4 | thresh=0.6: N | thresh=0.4: GC=S | mostCommon: GC=S | catchAll: B
    # 8: T 0.2 G 0.2 C 0.2 - 0.2 . 0.2 | thresh=0.6: N | thresh=0.4:N | mostCommon: TGC-.=TGC=B | catchAll: TGC-.=TGC=B
    # 9: A 0.8 N 0.2 | thresh=0.6: A | thresh=0.4: AN=A | mostCommon: A | catchAll: AN=A
    #10: A 0.2 T 0.6 G 0.2 | thresh=0.6: T | thresh=0.4: T | mostCommon: T | catchAll: ATG=D
    #11: - 0.8 . 0.2 | thresh=0.6: - | thresh=0.4: - | mostCommon: - | catchAll: -.=-
    #12: A 0.4 T 0.2 G 0.2 C 0.2 | thresh=0.6: N | thresh=0.4: A | mostCommon: A | catchAll: ATGC=N
    #13: . 1 | thresh=0.6: . | thresh=0.4: . | mostCommon: . | catchAll: .
    #14: T 0.2 G 0.4 C 0.2 N 0.2 | thresh=0.6: N | thresh=0.4: G | mostCommon: G | catchAll: TGCN=TGC=B
    #15: A 0.2 T 0.2 G 0.2 N 0.2 - 0.2 | thresh=0.6: N | thresh=0.4: N | mostCommon: ATGN-=ATG=D | catchAll: ATGN-=ATG=D
    #16: A 0.2 G 0.2 C 0.4 . 0.2 | thresh=0.6: N | thresh=0.4: C | mostCommon: C | catchAll: AGC.=AGC=V
    #17: (majority=floor(5/2)=2; 4>2 seqs have info at pos 17)
    #    A 0.25 T 0.25 G 0.25 N 0.25 | thresh=0.6: N | thresh=0.4: N | mostCommon: ATGN=ATG=D | catchAll: ATGN=ATG=D
    #18: (majority=floor(5/2)=2; 3>2 seqs have info at pos 18)
    #    A 1/3 T 1/3 C 1/3 | thresh=0.6: N | thresh=0.4: N | mostCommon: ATC=H | catchAll: ATC=H
    #19: freq doesn't matter b/c longest possible length is 17 (2<=2 seqs have info at pos 18)
    
    seqs1 = c("ATGCATGCAT-G.NTC",
              "ATGGATCGNT-A.GNCGC",
              "ACTGACT.AT-T.TA.N",
              "TCGAACCTAG.C.G-GAAA",
              "CGTAATG-AA-A.CGATTA")
    
    # region def spanning nucleotides 1 through 14
    regDef1 = createRegionDefinition(boundaries=factor(rep(c("W", "Y"), each=7)))
    
    ##### staggering lengths of sequences in seqs1 also allow for 
    # testing on correct output length
    
    ##### mostCommon
    ### ties resolved deterministically by representing ties using ambiguous chars 
    ## no region definition
    mostCommon.ambi.noRegDef = "AYGRATSBAT-A.GDCDH"
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="mostCommon", includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=mostCommon.ambi.noRegDef, muFreq=NULL))
    # when both includeAmbiguous and breakTiesStochastic are TRUE, includeAmbiguous takes precedence
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="mostCommon", includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=TRUE),
                 list(cons=mostCommon.ambi.noRegDef, muFreq=NULL))
    
    ## with region definition
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="mostCommon", includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=substr(mostCommon.ambi.noRegDef, 1, regDef1@seqLength), muFreq=NULL))
    # when both includeAmbiguous and breakTiesStochastic are TRUE, includeAmbiguous takes precedence
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="mostCommon", includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=TRUE),
                 list(cons=substr(mostCommon.ambi.noRegDef, 1, regDef1@seqLength), muFreq=NULL))
    
    ### ties resolved stochastically
    # all possible combinations
    mostCommon.sto.noRegDef = expand.grid("A", c("T", "C"), "G", c("A","G"), "A",
                                          "T", c("G","C"), c("T","G","C", "-", "."), "A", "T",
                                          "-", "A", ".", "G", c("A", "T", "G", "N", "-"),
                                          "C", c("A", "T", "G", "N"), c("A", "T", "C"))
    mostCommon.sto.noRegDef = apply(mostCommon.sto.noRegDef, 1, paste, collapse="")
    ## run stochastically 100 times without region definition
    test.mostCommon.sto.noRegDef = replicate(100, 
                                             shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                                                mtd="mostCommon", includeAmbiguous=FALSE, 
                                                                                breakTiesStochastic=TRUE)$cons)
    # for debugging: manually look at chars at each position
    #table(sapply(test.mostCommon.sto.noRegDef, function(x){substr(x,18,18)}))
    
    # result from each of 100 runs should be one of the possible combinations
    expect_true(all(test.mostCommon.sto.noRegDef %in% mostCommon.sto.noRegDef))
    
    ## run stochastically 100 times with region definition
    test.mostCommon.sto.regDef = replicate(100, 
                                           shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                                              mtd="mostCommon", includeAmbiguous=FALSE, 
                                                                              breakTiesStochastic=TRUE)$cons)
    # for debugging: manually look at chars at each position
    #table(sapply(test.mostCommon.sto.regDef, function(x){substr(x,12,12)}))
    
    # result from each of 100 runs should be one of the possible combinations
    expect_true(all(test.mostCommon.sto.regDef %in% substr(mostCommon.sto.noRegDef, 1, regDef1@seqLength)))
    
    ### resolve ties deterministcally by taking first char in the order of ATGCN-.
    mostCommon.det.1st.noRegDef = "ATGAATGTAT-A.GACAA"
    # no region definition
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="mostCommon", includeAmbiguous=FALSE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=mostCommon.det.1st.noRegDef, muFreq=NULL))
    
    # with region definitioin
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="mostCommon", includeAmbiguous=FALSE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=substr(mostCommon.det.1st.noRegDef,1,regDef1@seqLength), muFreq=NULL))
    
    ##### thresholdedFreq
    ## no region definition
    thresh0.6.ambi.noRegDef = "ANGNATNNAT-N.NNNNN"
    thresh0.4.ambi.noRegDef = "AYKRAYSNAT-A.GNCNN"
    # thresh 0.6
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="thresholdedFreq", minFreq=0.6,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=thresh0.6.ambi.noRegDef, muFreq=NULL))
    # thresh 0.4
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="thresholdedFreq", minFreq=0.4,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=thresh0.4.ambi.noRegDef, muFreq=NULL))    
    # when both includeAmbiguous and breakTiesStochastic are TRUE, includeAmbiguous takes precedence
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="thresholdedFreq", minFreq=0.6,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=TRUE),
                 list(cons=thresh0.6.ambi.noRegDef, muFreq=NULL))
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="thresholdedFreq", minFreq=0.4,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=TRUE),
                 list(cons=thresh0.4.ambi.noRegDef, muFreq=NULL))
    
    ## with region definition
    # thresh 0.6
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="thresholdedFreq", minFreq=0.6,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=substr(thresh0.6.ambi.noRegDef, 1, regDef1@seqLength), muFreq=NULL))
    # thresh 0.4
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="thresholdedFreq", minFreq=0.4,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=substr(thresh0.4.ambi.noRegDef, 1, regDef1@seqLength), muFreq=NULL))
    # when both includeAmbiguous and breakTiesStochastic are TRUE, includeAmbiguous takes precedence
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="thresholdedFreq", minFreq=0.6,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=TRUE),
                 list(cons=substr(thresh0.6.ambi.noRegDef, 1, regDef1@seqLength), muFreq=NULL))
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="thresholdedFreq", minFreq=0.4,
                                                    includeAmbiguous=TRUE, 
                                                    breakTiesStochastic=TRUE),
                 list(cons=substr(thresh0.4.ambi.noRegDef, 1, regDef1@seqLength), muFreq=NULL))
    
    
    ### ties resolved stochastically
    # all possible combinations
    thresh0.6.sto.noRegDef = expand.grid("A", "N", "G", "N", "A",
                                         "T", "N", "N", "A", "T",
                                         "-", "N", '.', "N", "N",
                                         "N", "N", "N")
    thresh0.6.sto.noRegDef = apply(thresh0.6.sto.noRegDef, 1, paste, collapse="")
    
    thresh0.4.sto.noRegDef = expand.grid("A", c("T", "C"), c("T", "G"), c("A", "G"), "A",
                                         c("T", "C"), c("G","C"), "N", c("A", "N"), "T",
                                         "-", "A", ".", "G", "N", 
                                         "C", "N", "N")
    thresh0.4.sto.noRegDef = apply(thresh0.4.sto.noRegDef, 1, paste, collapse="")
    
    ## run stochastically 100 times without region definition
    test.thresh0.6.sto.noRegDef = replicate(100, 
                                            shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                                               mtd="thresholdedFreq", minFreq=0.6,
                                                                               includeAmbiguous=FALSE, 
                                                                               breakTiesStochastic=TRUE)$cons)
    test.thresh0.4.sto.noRegDef = replicate(100, 
                                            shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                                               mtd="thresholdedFreq", minFreq=0.4,
                                                                               includeAmbiguous=FALSE, 
                                                                               breakTiesStochastic=TRUE)$cons)
    # for debugging: manually look at chars at each position
    #table(sapply(test.thresh0.4.sto.noRegDef, function(x){substr(x,18,18)}))
    
    # result from each of 100 runs should be one of the possible combinations
    expect_true(all(test.thresh0.6.sto.noRegDef %in% thresh0.6.sto.noRegDef))
    expect_true(all(test.thresh0.4.sto.noRegDef %in% thresh0.4.sto.noRegDef))
    
    ## run stochastically 100 times with region definition
    test.thresh0.6.sto.regDef = replicate(100, 
                                          shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                                             mtd="thresholdedFreq", minFreq=0.6,
                                                                             includeAmbiguous=FALSE, 
                                                                             breakTiesStochastic=TRUE)$cons)
    test.thresh0.4.sto.regDef = replicate(100, 
                                          shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                                             mtd="thresholdedFreq", minFreq=0.4,
                                                                             includeAmbiguous=FALSE, 
                                                                             breakTiesStochastic=TRUE)$cons)
    
    # for debugging: manually look at chars at each position
    #table(sapply(test.thresh0.4.sto.regDef, function(x){substr(x,12,12)}))
    
    # result from each of 100 runs should be one of the possible combinations
    expect_true(all(test.thresh0.6.sto.regDef %in% substr(thresh0.6.sto.noRegDef, 1, regDef1@seqLength)))
    expect_true(all(test.thresh0.4.sto.regDef %in% substr(thresh0.4.sto.noRegDef, 1, regDef1@seqLength)))
    
    ### resolve ties deterministcally by taking first char in the order of ATGCN-.
    thresh0.6.det.1st.noRegDef = "ANGNATNNAT-N.NNNNN"
    thresh0.4.det.1st.noRegDef = "ATTAATGNAT-A.GNCNN"
    ## no region definition
    # thresh 0.6
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="thresholdedFreq", minFreq=0.6,
                                                    includeAmbiguous=FALSE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=thresh0.6.det.1st.noRegDef, muFreq=NULL))
    # thresh 0.4
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, 
                                                    mtd="thresholdedFreq", minFreq=0.4,
                                                    includeAmbiguous=FALSE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=thresh0.4.det.1st.noRegDef, muFreq=NULL))
    
    ## with region definitioin
    # thresh 0.6
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="thresholdedFreq", minFreq=0.6,
                                                    includeAmbiguous=FALSE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=substr(thresh0.6.det.1st.noRegDef,1,regDef1@seqLength), muFreq=NULL))
    # thresh 0.4
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, 
                                                    mtd="thresholdedFreq", minFreq=0.4,
                                                    includeAmbiguous=FALSE, 
                                                    breakTiesStochastic=FALSE),
                 list(cons=substr(thresh0.4.det.1st.noRegDef,1,regDef1@seqLength), muFreq=NULL))
    
    ##### catchAll
    catchAll.noRegDef = "HBKVAYBBAD-N.BDVDH"
    # no region definition
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=NULL, mtd="catchAll"),
                 list(cons=catchAll.noRegDef, muFreq=NULL))
    # with region definition
    expect_equal(shazam:::calcClonalConsensusHelper(seqs=seqs1, lenLimit=regDef1@seqLength, mtd="catchAll"),
                 list(cons=substr(catchAll.noRegDef,1,regDef1@seqLength), muFreq=NULL))
    
    ##### longest possible length 
    ## note that this had also been tested by the tests above (by seqs1 having staggering lengths)
    # total number of seqs = 4
    # floor(n/2) = floor(4/2) = floor(2) = 2
    # longest possible length = number of positions at which >2 seqs have information
    seqs2 = c("ATGC",
              "AAG",
              "AT",
              "A")
    # expect longest possible length = 2 = length of consensus
    expect_equal(nchar(shazam:::calcClonalConsensusHelper(seqs=seqs2, mtd="catchAll")$cons), 
                 2)
})

#### calcClonalConsensus 2C ####
test_that("calcClonalConsensusHelper, 2C, methods = mostMutated, leastMutated", {
    # seq1: DUPCOUNT=37; CONSCONT=25; ERR=0.3
    # obsv: [full length=15] 1 R; nonN=15; muFreq = 1/15
    # obsv: ATG CAT GCA TGC ATA 
    # germ: ATG CAT GCA TGC ATG 
    
    # seq2: DUPCOUNT=9; CONSCONT=8; ERR=0.076
    # obsv: [full length=18] 1 R; 1 S; nonN=18; muFreq = 2/18
    # obsv: ATG CAT GCG TGC ATA CGT
    # germ: ATG CAT GCA TGC ATG CGT
    
    # seq3: DUPCOUNT=37; CONSCONT=25; ERR=0.23
    # obsv: [full length=17] 1 S; nonN=15; muFreq = 1/15
    # obsv: ATG CAC GCG TGC ATG CC
    # germ: ATG CAC GCA TGC ATG CC
    # germ: 2nd codon CAC instead of CAT
    
    # seq4: DUPCOUNT=34; CONSCONT=25; ERR=0.44
    # obsv: [full length=15] 1 S; nonN=15; muFreq = 1/15
    # obsv: ATG CAC GCG TGC ATG 
    # germ: ATG CAC GCA TGC ATG 
    # germ: 2nd codon CAC instead of CAT
    
    # seq5: DUPCOUNT=11; CONSCONT=20; ERR=0.14 
    # obsv: [full length=19] 1 R; 2 S; nonN=18 (ignoring non-triplet overhang); muFreq = 3/18
    # obsv: ATG CAT GCG TGT ATA CGC G
    # germ: ATG CAT GCA TGC ATG CGT G
    
    # seq6: DUPCOUNT=11; CONSCONT=20; ERR=0.12
    # obsv: [full length=20] 2 R; 1 S; nonN=18; muFreq = 3/18
    # obsv: ATG CAT GCA TGT ATA CGT GA
    # germ: ATT CAT GCA TGC ATG CGT GA
    # germ: 1st codon ATT instead of ATG
    
    # seq7: DUPCOUNT=11; CONSCONT=16; ERR=0.17
    # obsv: [full length=18] 2 R; 1 S; nonN=18; muFreq = 3/18
    # obsv: ATG CAT GCA TGT ATA CGT 
    # germ: ATT CAT GCA TGC ATG CGT
    # germ: 1st codon ATT instead of ATG
    
    testDb = data.frame(obsv=c("ATGCATGCATGCATA",      # seq1
                               "ATGCATGCGTGCATACGT",   # seq2
                               "ATGCACGCGTGCATGCC",      # seq3
                               "ATGCACGCGTGCATG",      # seq4
                               "ATGCATGCGTGTATACGTG",  # seq5
                               "ATGCATGCATGTATACGTGA", # seq6
                               "ATGCATGCATGTATACGT"    # seq7
    ),
    germ=c("ATGCATGCATGCATG",      # seq1
           "ATGCATGCATGCATGCGT",   # seq2
           "ATGCACGCATGCATGCC",      # seq3
           "ATGCACGCATGCATG",      # seq4
           "ATGCATGCATGCATGCGTG",  # seq5
           "ATTCATGCATGCATGCGTGA", # seq6
           "ATTCATGCATGCATGCGT"    # seq7
    ),
    DUPCOUNT=c(37,9,37,34,11,11,11),
    CONSCOUNT=c(25,8,25,25,20,20,16),
    ERR=c(0.3, 0.076, 0.23, 0.44, 0.14, 0.12, 0.17),
    MUTFREQ=c(1/15,2/18,1/15,1/15,3/18,3/18,3/18),
    stringsAsFactors=FALSE)
    
    #observedMutations(testDb, "obsv", "germ", frequency=F, combine=F)
    #observedMutations(testDb, "obsv", "germ", frequency=T, combine=T)
    #calcObservedMutations(inputSeq = testDb$obsv[4], germlineSeq = testDb$germ[4], frequency=F, returnRaw=T)
    
    # check mutation frequency
    expect_true(all(observedMutations(testDb, "obsv", "germ", frequency=T, combine=T)$MU_FREQ ==
                        testDb$MUTFREQ))
    
    
    ##### mostMutated
    ### resolve ties stochastically
    most.sto.possible = testDb$obsv[5:7]
    most.sto.1 = replicate(100, 
                           shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                              mtd="mostMutated", 
                                                              breakTiesStochastic=TRUE,
                                                              breakTiesByColumns=NULL,
                                                              lenLimit=NULL)$cons)
    # when both breakTiesStochastic and breakTiesByColumns TRUE, breakTiesStochastic takes precedence
    most.sto.2 = replicate(100, 
                           shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                              mtd="mostMutated", 
                                                              breakTiesStochastic=TRUE,
                                                              breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT", "ERR"), 
                                                                                      c(max,max,min)),
                                                              lenLimit=NULL)$cons)
    
    expect_true(all(most.sto.1 %in% most.sto.possible))
    expect_true(all(most.sto.2 %in% most.sto.possible))
    # check length
    expect_true(all( nchar(most.sto.1) %in% nchar(most.sto.possible) ))
    expect_true(all( nchar(most.sto.2) %in% nchar(most.sto.possible) ))
    
    ### resolve ties by columns
    # 3 columns; able to resolve
    # DUPCOUNT gives 5&6&7; CONSCOUNT gives 5&6; ERR gives 6
    most.byCol.1 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                      mtd="mostMutated", 
                                                      breakTiesStochastic=FALSE,
                                                      breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT", "ERR"), 
                                                                              c(max,max,min)),
                                                      lenLimit=NULL)$cons
    expect_equal(most.byCol.1, testDb$obsv[6])
    
    # 1 column; able to resolve
    # ERR gives 6
    most.byCol.2 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                      mtd="mostMutated", 
                                                      breakTiesStochastic=FALSE,
                                                      breakTiesByColumns=list(c("ERR"), 
                                                                              c(min)),
                                                      lenLimit=NULL)$cons
    expect_equal(most.byCol.2, testDb$obsv[6])
    
    # 2 columns; unable to resolve; returns sequence that appears first
    # DUPCOUNT gives 5&6&7; CONSCOUNT gives 5&6; 5 appears before 6, returns 5
    most.byCol.3 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                      mtd="mostMutated", 
                                                      breakTiesStochastic=FALSE,
                                                      breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT"), 
                                                                              c(max,max)),
                                                      lenLimit=NULL)$cons
    expect_equal(most.byCol.3, testDb$obsv[5])
    
    ### resolve ties deterministically by returning sequence that appears first (index 5)
    most.det.1 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                    mtd="mostMutated", 
                                                    breakTiesStochastic=FALSE,
                                                    breakTiesByColumns=NULL,
                                                    lenLimit=NULL)$cons
    expect_equal(most.det.1, testDb$obsv[5])
    
    ### check length when lenLimit is supplied
    most.det.2 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                    mtd="mostMutated", 
                                                    breakTiesStochastic=FALSE,
                                                    breakTiesByColumns=NULL,
                                                    lenLimit=7)$cons
    expect_equal(most.det.2, substr(testDb$obsv[5], 1, 7))
    
    
    ##### leastMutated
    ### resolve ties stochastically
    least.sto.possible = testDb$obsv[c(1,3,4)]
    least.sto.1 = replicate(100, 
                            shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                               mtd="leastMutated", 
                                                               breakTiesStochastic=TRUE,
                                                               breakTiesByColumns=NULL,
                                                               lenLimit=NULL)$cons)
    # when both breakTiesStochastic and breakTiesByColumns TRUE, breakTiesStochastic takes precedence
    least.sto.2 = replicate(100, 
                            shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                               mtd="leastMutated", 
                                                               breakTiesStochastic=TRUE,
                                                               breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT", "ERR"), 
                                                                                       c(max,max,min)),
                                                               lenLimit=NULL)$cons)
    
    expect_true(all(least.sto.1 %in% least.sto.possible))
    expect_true(all(least.sto.2 %in% least.sto.possible))
    # check length
    expect_true(all( nchar(least.sto.1) %in% nchar(least.sto.possible) ))
    expect_true(all( nchar(least.sto.2) %in% nchar(least.sto.possible) ))
    
    ### resolve ties by columns
    # 3 columns; able to resolve
    # DUPCOUNT gives 1&3; CONSCOUNT gives 1&3; ERR gives 3
    least.byCol.1 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                       mtd="leastMutated", 
                                                       breakTiesStochastic=FALSE,
                                                       breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT", "ERR"), 
                                                                               c(max,max,min)),
                                                       lenLimit=NULL)$cons
    expect_equal(least.byCol.1, testDb$obsv[3])
    
    # 1 column; able to resolve
    # ERR gives 3
    least.byCol.2 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                       mtd="leastMutated", 
                                                       breakTiesStochastic=FALSE,
                                                       breakTiesByColumns=list(c("ERR"), 
                                                                               c(min)),
                                                       lenLimit=NULL)$cons
    expect_equal(least.byCol.2, testDb$obsv[3])
    
    # 2 columns; unable to resolve; returns sequence that appears first
    # DUPCOUNT gives 1&3; CONSCOUNT gives 1&3; 1 appears before 3, returns 1
    least.byCol.3 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                       mtd="leastMutated", 
                                                       breakTiesStochastic=FALSE,
                                                       breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT"), 
                                                                               c(max,max)),
                                                       lenLimit=NULL)$cons
    expect_equal(least.byCol.3, testDb$obsv[1])
    
    ### resolve ties deterministically by returning sequence that appears first (index 1)
    least.det.1 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                     mtd="leastMutated", 
                                                     breakTiesStochastic=FALSE,
                                                     breakTiesByColumns=NULL,
                                                     lenLimit=NULL)$cons
    expect_equal(least.det.1, testDb$obsv[1])
    
    ### check length when lenLimit is supplied
    least.det.2 = shazam:::calcClonalConsensusHelper(seqs=testDb$obsv, muFreqColumn="MUTFREQ", db=testDb,
                                                     mtd="leastMutated", 
                                                     breakTiesStochastic=FALSE,
                                                     breakTiesByColumns=NULL,
                                                     lenLimit=7)$cons
    expect_equal(least.det.2, substr(testDb$obsv[1], 1, 7))
})

#### calcClonalConsensus 2D ####
test_that("calcClonalConsensus, 2D", {
    ##### same testDb from test 2C for calcClonalConsensusHelper
    testDb = data.frame(obsv=c("ATGCATGCATGCATA",      # seq1
                               "ATGCATGCGTGCATACGT",   # seq2
                               "ATGCACGCGTGCATGCC",      # seq3
                               "ATGCACGCGTGCATG",      # seq4
                               "ATGCATGCGTGTATACGTG",  # seq5
                               "ATGCATGCATGTATACGTGA", # seq6
                               "ATGCATGCATGTATACGT"    # seq7
    ),
    germ=c("ATGCATGCATGCATG",      # seq1
           "ATGCATGCATGCATGCGT",   # seq2
           "ATGCACGCATGCATGCC",      # seq3
           "ATGCACGCATGCATG",      # seq4
           "ATGCATGCATGCATGCGTG",  # seq5
           "ATTCATGCATGCATGCGTGA", # seq6
           "ATTCATGCATGCATGCGT"    # seq7
    ),
    DUPCOUNT=c(37,9,37,34,11,11,11),
    CONSCOUNT=c(25,8,25,25,20,20,16),
    ERR=c(0.3, 0.076, 0.23, 0.44, 0.14, 0.12, 0.17),
    MUTFREQ=c(1/15,2/18,1/15,1/15,3/18,3/18,3/18),
    stringsAsFactors=FALSE)
    
    ##### run calcClonalConsensus
    # no region definition
    test.result.noRegDef = shazam:::calcClonalConsensus(db=testDb, 
                                                        sequenceColumn="obsv", 
                                                        germlineColumn="germ", 
                                                        muFreqColumn="MUTFREQ",
                                                        regionDefinition=NULL, 
                                                        method="mostMutated", 
                                                        minimumFrequency=NULL, includeAmbiguous=FALSE,
                                                        breakTiesStochastic=FALSE, 
                                                        breakTiesByColumns=list(c("ERR"), c(min)))
    
    # with region definition
    test.regDef = createRegionDefinition(boundaries=factor(rep(c("W","Y"), each=6)))
    
    test.result.regDef = shazam:::calcClonalConsensus(db=testDb, 
                                                      sequenceColumn="obsv", 
                                                      germlineColumn="germ", 
                                                      muFreqColumn="MUTFREQ",
                                                      regionDefinition=test.regDef, 
                                                      method="mostMutated", 
                                                      minimumFrequency=NULL, includeAmbiguous=FALSE,
                                                      breakTiesStochastic=FALSE, 
                                                      breakTiesByColumns=list(c("ERR"), c(min)))
    
    ##### returned object should be a list 
    expect_true(is.list(test.result.noRegDef))
    expect_true(is.list(test.result.regDef))
    # of length 3
    expect_equal(length(test.result.noRegDef), 3)
    expect_equal(length(test.result.regDef), 3)
    # with entries named inputCons, germlineCons, and inputMuFreq
    expect_equal(names(test.result.noRegDef), c("inputCons", "germlineCons", "inputMuFreq"))
    expect_equal(names(test.result.regDef), c("inputCons", "germlineCons", "inputMuFreq"))
    # entries should be character, character, and numeric
    expect_equal(unlist(lapply(test.result.noRegDef, class)), 
                 c(inputCons="character", germlineCons="character", inputMuFreq="numeric"))
    expect_equal(unlist(lapply(test.result.regDef, class)), 
                 c(inputCons="character", germlineCons="character", inputMuFreq="numeric"))
    
    ##### germlineCons should be generated from mostCommon method using the same additional parameter setting
    # generate mostCommon germline using calcClonalConsensusHelper (tested in test 2B)
    exp.germ.noRegDef = shazam:::calcClonalConsensusHelper(seqs=testDb$germ, mtd="mostCommon",
                                                           includeAmbiguous=FALSE,
                                                           breakTiesStochastic=FALSE,
                                                           lenLimit=NULL)
    exp.germ.regDef = shazam:::calcClonalConsensusHelper(seqs=testDb$germ, mtd="mostCommon",
                                                         includeAmbiguous=FALSE,
                                                         breakTiesStochastic=FALSE,
                                                         lenLimit=test.regDef@seqLength)
    
    expect_equal(test.result.noRegDef$germlineCons, exp.germ.noRegDef$cons)
    expect_equal(test.result.regDef$germlineCons, exp.germ.regDef$cons)
    
    ##### inputCons and germlineCons lengths should be the same
    expect_equal(nchar(test.result.noRegDef$germlineCons), nchar(test.result.noRegDef$inputCons))
    expect_equal(nchar(test.result.regDef$germlineCons), nchar(test.result.regDef$inputCons))
    
    ##### inputCons and germlineCons lengths should both be <= regDef@seqLength if supplied
    expect_true(nchar(test.result.regDef$germlineCons) <= test.regDef@seqLength)
    expect_true(nchar(test.result.regDef$inputCons) <= test.regDef@seqLength)
    
})

#### collapseClones ####
test_that("collapseClones, 2E", {
    ##### same testDb from test 2C for calcClonalConsensusHelper
    testDb = data.frame(obsv=c("ATGCATGCATGCATA",      # seq1
                               "ATGCATGCGTGCATACGT",   # seq2
                               "ATGCACGCGTGCATGCC",      # seq3
                               "ATGCACGCGTGCATG",      # seq4
                               "ATGCATGCGTGTATACGTG",  # seq5
                               "ATGCATGCATGTATACGTGA", # seq6
                               "ATGCATGCATGTATACGT"    # seq7
    ),
    germ=c("ATGCATGCATGCATG",      # seq1
           "ATGCATGCATGCATGCGT",   # seq2
           "ATGCACGCATGCATGCC",      # seq3
           "ATGCACGCATGCATG",      # seq4
           "ATGCATGCATGCATGCGTG",  # seq5
           "ATTCATGCATGCATGCGTGA", # seq6
           "ATTCATGCATGCATGCGT"    # seq7
    ),
    DUPCOUNT=c(37,9,37,34,11,11,11),
    CONSCOUNT=c(25,8,25,25,20,20,16),
    ERR=c(0.3, 0.076, 0.23, 0.44, 0.14, 0.12, 0.17),
    MUTFREQ=c(1/15,2/18,1/15,1/15,3/18,3/18,3/18),
    stringsAsFactors=FALSE)
    
    ##### create 3 identical clones from testDb
    testDb.clone = rbind(testDb, testDb, testDb)
    testDb.clone$CLONE = rep(c("124", "39", "5"), each=nrow(testDb))
    
    ##### check input checks/warnings
    # method
    expect_error(collapseClones(db=testDb.clone, method="nonExistingMethod"))
    
    # minimumFrequency for thresholdedFreq
    expect_error(collapseClones(db=testDb.clone, method="thresholdedFreq", minimumFrequency=NULL),
                 "minimumFrequency must be a numeric value.")
    expect_error(collapseClones(db=testDb.clone, method="thresholdedFreq", minimumFrequency=1.3),
                 "minimumFrequency must be between 0 and 1.")
    expect_error(collapseClones(db=testDb.clone, method="thresholdedFreq", minimumFrequency=(-1.3)),
                 "minimumFrequency must be between 0 and 1.")
    
    # includeAmbiguous & breakTiesStochastic for methods other than catchAll
    expect_error(collapseClones(db=testDb.clone, method="mostCommon", 
                                includeAmbiguous=NULL, breakTiesStochastic=FALSE),
                 "includeAmbiguous must be TRUE or FALSE.")
    expect_error(collapseClones(db=testDb.clone, method="mostCommon", 
                                includeAmbiguous=FALSE, breakTiesStochastic=NULL),
                 "breakTiesStochastic must be TRUE or FALSE.")
    
    # breakTiesByColumns and muFreqColumn for methods most/leastMutated
    expect_error(collapseClones(db=testDb.clone, method="mostMutated", 
                                breakTiesByColumns="notList"),
                 "breakTiesByColumns must be a list.")
    expect_error(collapseClones(db=testDb.clone, method="mostMutated", 
                                breakTiesByColumns=list(1,2,3)),
                 "breakTiesByColumns must be a nested list of length 2.")
    expect_error(collapseClones(db=testDb.clone, method="mostMutated", 
                                breakTiesByColumns=list(c("col1", "col2"), c(max,max,min))),
                 "Nested vectors in breakTiesByColumns must have the same lengths.")
    expect_error(collapseClones(db=testDb.clone, method="mostMutated", 
                                breakTiesByColumns=list(c(12,22,3), c(max,max,min))),
                 "The first vector in breakTiesByColumns must contain column names.")
    expect_error(collapseClones(db=testDb.clone, method="mostMutated", 
                                breakTiesByColumns=list(c("col1", "col2"), c("max","max"))),
                 "The second vector in breakTiesByColumns must contain functions.")
    expect_error(collapseClones(db=testDb.clone, method="mostMutated", 
                                breakTiesByColumns=list(c("col1", "col2"), c(max,max))),
                 "All column named included in breakTiesByColumns must be present in db.")
    expect_error(collapseClones(db=testDb.clone, method="mostMutated", muFreqColumn="notInDb",
                                breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT"), c(max,max))),
                 "If specified, muFreqColumn must be a column present in db.")
    
    # check mutual exclusivitiy
    expect_message(collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                                  method="mostCommon", includeAmbiguous=TRUE, breakTiesStochastic=TRUE),
                   "includeAmbiguous and breakTiesStochastic are mutually exclusive. When both TRUE, includeAmbiguous will take precedence.")
    #expect_message(collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
    #                              method="mostCommon", includeAmbiguous=FALSE, breakTiesStochastic=FALSE),
    #               "When both includeAmbiguous and breakTiesStochastic are FALSE, ties are broken in the order of 'A', 'T', 'G', 'C', 'N', '.', and '-'.")
    expect_message(collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                                  method="mostCommon", includeAmbiguous=TRUE, breakTiesStochastic=FALSE,
                                  breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT"), c(max,max))),
                   "breakTiesByColumns is ignored when method is thresholdedFreq or mostCommon.")
    
    expect_message(collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                                  method="mostMutated", breakTiesStochastic=TRUE,
                                  breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT"), c(max,max))),
                   "breakTiesStochastic and breakTiesByColumns are mutually exclusive. When both set, breakTiesStochastic will take precedence.")
    #expect_message(collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
    #                              method="mostMutated", breakTiesStochastic=FALSE,
    #                              breakTiesByColumns=NULL),
    #               "When breakTiesStochastic is FALSE and breakTiesByColumns is NULL, ties are broken by taking the sequence that appears earlier in the data.frame.")
    expect_message(collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                                  method="mostMutated", breakTiesStochastic=FALSE,
                                  breakTiesByColumns=list(c("DUPCOUNT", "CONSCOUNT"), c(max,max)),
                                  includeAmbiguous=TRUE),
                   "includeAmbiguous is ignored when method is mostMutated or leastMutated.")
    
    
    ##### check calling observedMutations when muFreqColumn is not specifiied
    # resolve ties determinisically by taking seq that appears first 
    test.mut.most = collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                                   muFreqColumn=NULL,
                                   method="mostMutated", breakTiesStochastic=FALSE,
                                   breakTiesByColumns=NULL, expandedDb=FALSE)
    expect_equal(test.mut.most[["CLONAL_SEQUENCE_MUFREQ"]],
                 rep(testDb$MUTFREQ[5], 3)) # from test 2D
    
    test.mut.least = collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                                    muFreqColumn=NULL,
                                    method="leastMutated", breakTiesStochastic=FALSE,
                                    breakTiesByColumns=NULL, expandedDb=FALSE)
    expect_equal(test.mut.least[["CLONAL_SEQUENCE_MUFREQ"]],
                 rep(testDb$MUTFREQ[1], 3)) # from test 2D
    
    ##### check expandedDb
    # resolve ties determinisically by taking seq that appears first 
    test.expF = collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                               muFreqColumn="MUTFREQ",
                               method="mostMutated", breakTiesStochastic=FALSE,
                               breakTiesByColumns=NULL, expandedDb=FALSE)
    test.expT = collapseClones(db=testDb.clone, cloneColumn="CLONE", sequenceColumn="obsv", germlineColumn="germ",
                               muFreqColumn="MUTFREQ",
                               method="mostMutated", breakTiesStochastic=FALSE,
                               breakTiesByColumns=NULL, expandedDb=TRUE)
    
    ### dimension
    # nrow of unexpanded should correspond to number of clones
    expect_equal(nrow(test.expF), 3) # 3 clones
    # there should be 3 unique clones in unexpanded
    expect_equal(length(unique(test.expF[["CLONE"]])), 3) # 3 clones
    
    # nrow of expanded should correspond to nrow of input db
    expect_equal(nrow(test.expT), nrow(testDb.clone)) 
    # number of seqs in each unique clone in expanded should correspond to that in input db
    expect_equal(as.vector(table(test.expT[["CLONE"]])), rep(nrow(testDb), 3)) 
    
    
    ### muFreqColumn supplied
    # when expanded, MUTFREQ for all sequences in input db should remain unchanged
    expect_equal(length(unique(test.expT[["MUTFREQ"]])), length(unique(testDb$MUTFREQ)))
    # when unexpended, sequenceColumn and corresponding MUTFREQ retained for each clone may
    # not necessarily be the most mutated; so nothing to test here
    
    ### CLONAL_SEQUENCE_MUFREQ
    # $CLONAL_SEQUENCE_MUFREQ is the mut freq of most mutated consensus for a given clone
    # since all 3 clones here are artifically made the same, this number should be the same for all clones
    expect_equal(length(unique(test.expF[["CLONAL_SEQUENCE_MUFREQ"]])), 1)
    
    # $CLONAL_SEQUENCE_MUFREQ is the mut freq of most mutated consensus for a given clone
    # this number should be the same for all seqs in the same clone
    # since all 3 clones here are artifically made the same, this number should be the same for all seqs
    expect_equal(length(unique(test.expT[["CLONAL_SEQUENCE_MUFREQ"]])), 1)
    
})

#### mutationType ####
test_that("mutationType", {
    # R, S, Stop, na
    
    ##### no ambiguous char in input
    # no need to specify ambiguousMode
    # if specified, there should be no effect
    
    # TGG (Trp) -> TGG (Trp); expect na 1
    expect_equivalent(shazam:::mutationType("TGG", "TGG"), c(0,0,0,1))
    expect_equivalent(shazam:::mutationType("TGG", "TGG", ambiguousMode="eitherOr"), c(0,0,0,1))
    expect_equivalent(shazam:::mutationType("TGG", "TGG", ambiguousMode="and"), c(0,0,0,1))
    
    # TGG (Trp) -> TAG (Stop); expect Stop 1
    expect_equivalent(shazam:::mutationType("TGG", "TAG"), c(0,0,1,0))
    expect_equivalent(shazam:::mutationType("TGG", "TAG", ambiguousMode="eitherOr"), c(0,0,1,0))
    expect_equivalent(shazam:::mutationType("TGG", "TAG", ambiguousMode="and"), c(0,0,1,0))
    
    
    # TGG (Trp) -> TCG (Ser); expect R 1
    expect_equivalent(shazam:::mutationType("TGG", "TCG"), c(1,0,0,0))
    expect_equivalent(shazam:::mutationType("TGG", "TCG", ambiguousMode="eitherOr"), c(1,0,0,0))
    expect_equivalent(shazam:::mutationType("TGG", "TCG", ambiguousMode="and"), c(1,0,0,0))
    
    # TGC (Cys) -> TGT (Cys); expect S 1
    expect_equivalent(shazam:::mutationType("TGC", "TGT"), c(0,1,0,0))
    expect_equivalent(shazam:::mutationType("TGC", "TGT", ambiguousMode="eitherOr"), c(0,1,0,0))
    expect_equivalent(shazam:::mutationType("TGC", "TGT", ambiguousMode="and"), c(0,1,0,0))
    
    ##### ambiguous char in input
    # ambiguousMode matters
    # if not specified, default is eitherOr
    
    # TGG (Trp) -> TAG (Stop), TGG (Trp) -> TTG (Leu)
    # "and": expects R 1 + Stop 1
    # "eitherOr": expects R 1
    expect_equivalent(shazam:::mutationType("TGG", "TWG"), c(1,0,0,0))
    expect_equivalent(shazam:::mutationType("TGG", "TWG", ambiguousMode="eitherOr"), c(1,0,0,0))
    expect_equivalent(shazam:::mutationType("TGG", "TWG", ambiguousMode="and"), c(1,0,1,0))
    
    # TGG (Trp) -> TCG (Ser), TGG (Trp) -> TGG (Trp)
    # "and": expects R 1 + na 1
    # "eitherOr": expects na 1
    expect_equivalent(shazam:::mutationType("TGG", "TSG"), c(0,0,0,1))
    expect_equivalent(shazam:::mutationType("TGG", "TSG", ambiguousMode="eitherOr"), c(0,0,0,1))
    expect_equivalent(shazam:::mutationType("TGG", "TSG", ambiguousMode="and"), c(1,0,0,1))
    
    # TGG (Trp) -> TCA (Ser)
    # TGG (Trp) -> TCC (SER)
    # TGG (Trp) -> TGA (Stop)
    # TGG (Trp) -> TGC (Cys)
    # "and": expects R 3 + stop 1
    # "eitherOr": expects R 1
    expect_equivalent(shazam:::mutationType("TGG", "TSM"), c(1,0,0,0))
    expect_equivalent(shazam:::mutationType("TGG", "TSM", ambiguousMode="eitherOr"), c(1,0,0,0))
    expect_equivalent(shazam:::mutationType("TGG", "TSM", ambiguousMode="and"), c(3,0,1,0))
    
    ##### With classes
    classes <- HYDROPATHY_MUTATIONS@classes
    #HYDROPATHY_MUTATIONS@codonTable[, "TTT"]
    expect_equivalent(shazam:::mutationType("TTT", "TTC", aminoAcidClasses=classes),
                      c(0,1,0,0))
    expect_equivalent(shazam:::mutationType("TTT", "TTA", aminoAcidClasses=classes),
                      c(0,1,0,0))
    expect_equivalent(shazam:::mutationType("TTT", "TCT", aminoAcidClasses=classes),
                      c(1,0,0,0))
    expect_equivalent(shazam:::mutationType("TTT", "TGA", aminoAcidClasses=classes),
                      c(0,0,1,0)) # TGA=Stop
    
})

#### nucs2IUPAC ####
test_that("nucs2IUPAC", {
    expect_equivalent(shazam:::nucs2IUPAC(c("A", "T")), "W")
    expect_equivalent(shazam:::nucs2IUPAC(c("A", "T", "G", "C")), "N")
    expect_equivalent(shazam:::nucs2IUPAC(c("C", "T", "G")), "B")
    expect_equivalent(shazam:::nucs2IUPAC(c("C", "T", "G", "G")), "B")
})

#### chars2Ambiguous ####
test_that("chars2Ambiguous", {
    expect_equivalent(shazam:::chars2Ambiguous(c("A", "T")), "W")
    expect_equivalent(shazam:::chars2Ambiguous(c("A", "T", "N")), "W")
    expect_equivalent(shazam:::chars2Ambiguous(c("A", "T", "G", "C")), "N")
    expect_equivalent(shazam:::chars2Ambiguous(c("A", "T", "G", "C", "N")), "N")
    expect_equivalent(shazam:::chars2Ambiguous(c("A", "T", "G", "C", "-")), "N")
    expect_equivalent(shazam:::chars2Ambiguous(c(".", "-")), "-")
    expect_equivalent(shazam:::chars2Ambiguous(c(".", "N")), "N")
    expect_equivalent(shazam:::chars2Ambiguous(c(".", "A", "T")), "W")
})

#### IUPAC2nucs ####
test_that("IUPAC2nucs", {
    expect_equivalent(shazam:::IUPAC2nucs(code="N", excludeN=T), "N")
    expect_equivalent(shazam:::IUPAC2nucs(code="N", excludeN=F), c("A", "C", "G", "T"))
    expect_equivalent(shazam:::IUPAC2nucs(code="S", excludeN=T), c("C", "G"))
    expect_equivalent(shazam:::IUPAC2nucs(code="S", excludeN=F), c("C", "G"))
})
