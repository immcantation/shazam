#load(file.path("tests", "data-tests", "ExampleDb.rda"))
load(file.path("..", "data-tests", "ExampleDb.rda"))

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


test_that("observedMutations, charge mutations", {

    db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")[1:10, ]

    db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                             germlineColumn="GERMLINE_IMGT_D_MASK",
                             regionDefinition=IMGT_V,
                             mutationDefinition=CHARGE_MUTATIONS,
                             nproc=1)
    
    expect_equal(db_obs$OBSERVED_CDR_R, c(0, 0, 0, 0, 0, 0, 2, 2, 2, 0))
    expect_equal(db_obs$OBSERVED_CDR_S, c(0, 0, 0, 0, 0, 0, 3, 3, 3, 16))
    expect_equal(db_obs$OBSERVED_FWR_R, c(0, 1, 1, 0, 0, 0, 0, 0, 0, 3))
    expect_equal(db_obs$OBSERVED_FWR_S, c(0, 6, 6, 0, 0, 0, 10, 10, 10, 14))
    
})

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

test_that("expecteddMutations, hydropathy", {
    
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

test_that("observedMutations, when no mutations found", {
    
    in_seq <- ExampleDb[["SEQUENCE_IMGT"]][1]
    germ_seq <- in_seq
    
    expect_equivalent(observedMutations(data.frame(
        "SEQUENCE_IMGT"=c(in_seq, in_seq),
        "GERMLINE_IMGT_D_MASK"=c(in_seq, in_seq)
    ), regionDefinition = NULL)[,c("OBSERVED_SEQ_R", "OBSERVED_SEQ_S")], data.frame(c(0, 0), c(0, 0)))
})

##### Tests 1A-1H added on June 9 2017

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

test_that("calcObservedMutations, 1B, with ambiguous characters, length is multiple of 3", {
    #obsv: "TAT ATA WSC -GT CDG CTC TCG" 
    #germ: "TAT TAT ATA GGT CTT CNC AAC" 
    #region  W   Y   Y   W   Y   W   W 
    #codon   1   2   3   4   5   6   7
    
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
                                                regionDefinition=NULL)
    # based on annotations above 
    exp.raw.noRegDef.position = c(4, 5, 7, 8, 9, 14, 15, 19, 20, 21)
    exp.raw.noRegDef.R = c(1,1,1,2,0,2,0,1,1,1)
    exp.raw.noRegDef.S = c(0,0,0,0,1,0,1,0,0,0)
    
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
    exp.raw.regDef.region = c("Y","Y","Y","Y","Y","Y","Y","W","W","W")
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
    exp.raw.regDef.region = c("Y","Y","Y","Y","Y")
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
                                                regionDefinition=NULL)
    
    expect_equal(freqF.rawT.noRegDef$pos, NA)
    exp.noRegDef.nonN.Dash.Dot = sum( seqinr::s2c(obsv) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
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
    
    exp.regDef.W.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
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
                                                regionDefinition=NULL)
    
    # based on annotations above 
    exp.raw.noRegDef.position = 2
    exp.raw.noRegDef.R = 2
    exp.raw.noRegDef.S = 0
    
    expect_equal(freqF.rawT.noRegDef$pos$position, exp.raw.noRegDef.position)
    expect_equal(freqF.rawT.noRegDef$pos$R, exp.raw.noRegDef.R)
    expect_equal(freqF.rawT.noRegDef$pos$S, exp.raw.noRegDef.S)
    expect_equal(freqF.rawT.noRegDef$pos$region, rep("SEQ", length(exp.raw.noRegDef.position)))
    
    exp.noRegDef.nonN.Dash.Dot = sum( seqinr::s2c(obsv) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equivalent(freqF.rawT.noRegDef$nonN, exp.noRegDef.nonN.Dash.Dot)
    
    ##### freq=F, returnRaw=F, regDef=NULL
    
    freqF.rawF.noRegDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                                frequency=FALSE, returnRaw=FALSE, 
                                                regionDefinition=NULL)
    expect_equal(freqF.rawF.noRegDef, c("SEQ_R"=exp.raw.noRegDef.R, 
                                        "SEQ_S"=exp.raw.noRegDef.S))
    
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
    
    exp.raw.regDef.region = "W"
    expect_equal(freqF.rawT.regDef$pos$region, exp.raw.regDef.region)
    # counts should be the same
    expect_identical(freqF.rawT.regDef$pos[, 1:3], freqF.rawT.noRegDef$pos[, 1:3])
    
    exp.regDef.W.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
    expect_equal(freqF.rawT.regDef$nonN, c("W"=exp.regDef.W.nonN.Dash.Dot))
    
    
    ##### freq=F, returnRaw=F, with regDef
    freqF.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=FALSE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    
    exp.raw.regDef.W.R = sum(subset(freqF.rawT.regDef$pos, region=="W")$R)
    exp.raw.regDef.W.S = sum(subset(freqF.rawT.regDef$pos, region=="W")$S)
    
    expect_equal(freqF.rawF.regDef, c("W_R"=exp.raw.regDef.W.R,
                                      "W_S"=exp.raw.regDef.W.S))
    
    ##### freq=T, returnRaw=F/T, with regDef
    # returnRaw value doesn't matter if frequency=T
    freqT.rawF.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=FALSE, 
                                              regionDefinition=regDef)
    freqT.rawT.regDef = calcObservedMutations(inputSeq=obsv, germlineSeq=germ, 
                                              frequency=TRUE, returnRaw=TRUE, 
                                              regionDefinition=regDef)
    
    expect_equal(freqT.rawF.regDef, c("W_R"=exp.raw.regDef.W.R/exp.regDef.W.nonN.Dash.Dot,
                                      "W_S"=exp.raw.regDef.W.S/exp.regDef.W.nonN.Dash.Dot))
    
})

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
    exp.noRegDef.nonN.Dash.Dot = sum( seqinr::s2c(obsv) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ) %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
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
    
    exp.regDef.W.nonN.Dash.Dot = sum( seqinr::s2c(obsv)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] &
                                          seqinr::s2c(germ)[regDef@boundaries=="W"] %in% shazam:::NUCLEOTIDES_AMBIGUOUS[1:14] )
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

test_that("calcObservedMutations, 1G, with ambiguous characters in germline", {
    # expect fail
    obsv = "THT" 
    germ = "TWT"
    
    expect_error(calcObservedMutations(inputSeq=obsv, germlineSeq=germ),
                 "germlineSeq cannot contain ambiguous characters.")
    
    # test with lowercase sequences too
    # should fail too
    obsv = "ggg" 
    germ = "gsg"
    
    expect_error(calcObservedMutations(inputSeq=obsv, germlineSeq=germ),
                 "germlineSeq cannot contain ambiguous characters.")
    
})

test_that("observedMutations, 1H, using mock data from 1A through 1F", {
    # pull sequences from 1A-1F
    testDb = data.frame(obsv=c("TATATAATC-GTCAGCTCTCG", # 1A
                               "TATATAWSC-GTCDGCTCTCG", # 1B
                               "TATATAATC-GTCAGCTCTC",  # 1C
                               "TAT", # 1D
                               "THT", # 1E
                               "TA"), # 1F
                        germ=c("TATTATATAGGTCTTCNCAAC", # 1A
                               "TATTATATAGGTCTTCNCAAC", # 1B
                               "TATTATATAGGTCTTCNCAA",  # 1C
                               "TAT", # 1D
                               "TAT", # 1E
                               "TA"), # 1F
                        stringsAsFactors=FALSE)
    
    # pull results from 1A-1F
    # these are results from calcObservedMutations which have been tested
    exp.noRegDef.R = c(6, 10, 3, NA, 2, NA)
    exp.noRegDef.S = c(2, 2, 2, NA, 0, NA)
    exp.noRegDef.nonN = c(19, 19, 18, 3, 3, 2)
    
    # convert NA to 0
    exp.noRegDef.R[is.na(exp.noRegDef.R)] = 0
    exp.noRegDef.S[is.na(exp.noRegDef.S)] = 0
    
    # frequency=F, combine=F, no regDef
    freqF.combF.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=F, combine=F)
    
    expect_equal(freqF.combF.noRegDef$OBSERVED_SEQ_R, exp.noRegDef.R)
    expect_equal(freqF.combF.noRegDef$OBSERVED_SEQ_S, exp.noRegDef.S)
    
    # frequency=F, combine=T, no regDef
    freqF.combT.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=F, combine=T)
    
    expect_equal(freqF.combT.noRegDef$OBSERVED, 
                 exp.noRegDef.R+exp.noRegDef.S)
    
    # frequency=T, combine=F, no regDef
    freqT.combF.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=T, combine=F)
    
    expect_equal(freqT.combF.noRegDef$MU_FREQ_SEQ_R, 
                 exp.noRegDef.R/exp.noRegDef.nonN)
    expect_equal(freqT.combF.noRegDef$MU_FREQ_SEQ_S, 
                 exp.noRegDef.S/exp.noRegDef.nonN)
    
    # frequency=T, combine=T, no regDef
    freqT.combT.noRegDef = observedMutations(db=testDb, sequenceColumn="obsv", 
                                             germlineColumn="germ",
                                             regionDefinition=NULL, 
                                             frequency=T, combine=T)
    
    expect_equal(freqT.combT.noRegDef$MU_FREQ, 
                 (exp.noRegDef.R+exp.noRegDef.S)/exp.noRegDef.nonN)
    
})
