# Load test database
e1 <- new.env()
#load(file.path("tests", "data-tests", "TestDb.rda"), envir=e1)
load(file.path("..", "data-tests", "TestDb.rda"), envir=e1)
db <- get("TestDb", envir=e1)
rm(e1)

#### distToNearest - cross with hh_s1f ####

test_that("Test cross distToNearest with model hh_s1f", {
    ## Reproduce example
    db <- subset(db, CPRIMER %in% c("IGHA", "IGHM") & 
                     BARCODE %in% c("RL016", "RL018", "RL019", "RL021"))
    db_nrow <- nrow(db)
    db2 <- dplyr::bind_rows(db,db)
    db2$SAMPLE<- "S2"
    db2$SAMPLE[1:db_nrow] <- "S1"
    db2$DONOR <- "D1"
    
    test_idx <- c(1:5,10:15,300:305)
    
    ## Test hh_s1f 
    
    dist_hs1f <- distToNearest(db, sequenceColumn="JUNCTION",
                               vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                               model="hh_s1f", first=FALSE, normalize="len",
                               locusColumn="LOCUS")
    ## Test if the updated function reproduces results
    expect_equal(dist_hs1f$dist_nearest[test_idx],
                 c(NA, NA, NA, NA, 
                   0.4060, 0.4465, 0.3976, 0.3482, 0.3064, 0.3064, 0.4295, 
                   0.0436, 0.1216, 0.3786, 0.3875, 0.1216, 0.3700),
                 tolerance=0.005)
    
    ## There's only one donor, hence cross-donor will return all NA
    dist_hs1f_cross_donor <- distToNearest(db2, sequenceColumn="JUNCTION",
                                           vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                           model="hh_s1f", first=FALSE, normalize="none", nproc=1, cross="DONOR",
                                           locusColumn="LOCUS")

    expect_true(all(is.na(dist_hs1f_cross_donor$cross_dist_nearest)))
    
    ## fields=NULL and fields=DONOR should give same results
    cross_dist_hs1f <- distToNearest(db2, sequenceColumn="JUNCTION",
                                     vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                     model="hh_s1f", first=FALSE, normalize="len",nproc=1,
                                     cross="SAMPLE",
                                     locusColumn="LOCUS")
    
    
    cross_dist_hs1f_donor <- distToNearest(db2, sequenceColumn="JUNCTION",
                                           vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                           model="hh_s1f", first=FALSE, normalize="len",nproc=1,
                                           cross="SAMPLE", fields="DONOR",
                                           locusColumn="LOCUS")
    
    expect_equal(cross_dist_hs1f, cross_dist_hs1f_donor)
    expect_equal(cross_dist_hs1f$cross_dist_nearest[test_idx],
                 c(NA,NA,NA,NA,0.4040, 0.4447, 0.3963, 0.3469, 0.3050, 0.3050,
                   0.4284, 0.0435, 0.1212, 0.3771, 0.3862, 0.1212, 0.3687),
                 tolerance=0.005)
    ## Check cross
    ## Crossing should reproduce the same results as not crossed
    ## because both donors have the same db
    expect_equal(cross_dist_hs1f$cross_dist_nearest[1:nrow(dist_hs1f)],
                 dist_hs1f$dist_nearest, tolerance=0.005)
    
    ## Inroduce row as Sample 3, very similar to rows 1 and 316
    ## This will be the best
    db3 <- db2
    db3[nrow(db3),] <- db3[1,]
    db3[nrow(db3),"JUNCTION"] <- sub("^T","A",db3[nrow(db3),"JUNCTION"])
    db3[nrow(db3),"SAMPLE"] <- "S3"
    # t(db3[c(1,316,nrow(db3)),c("JUNCTION","V_CALL","J_CALL","DONOR","SAMPLE")] )
    
    db2_1_316_630 <- distToNearest(db2[c(1,316,630),], sequenceColumn="JUNCTION",
                                   vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                   model="hh_s1f", first=FALSE, normalize="len",cross="SAMPLE",
                                   locusColumn="LOCUS")
    ## Exactly same seq, returns NA
    expect_equal(db2_1_316_630$cross_dist_nearest,c(NA,NA,NA))
    
    ## One seq has been edited, will return distance values
    db3_1_316_630 <- distToNearest(db3[c(1,316,630),], sequenceColumn="JUNCTION",
                                   vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                   model="hh_s1f", first=FALSE, normalize="len",cross="SAMPLE",
                                   locusColumn="LOCUS")
    expect_equal(db3_1_316_630$cross_dist_nearest,c(0.0175,0.0175,0.0175), tolerance=0.005)
})

#### distToNearest - cross with hh_s5f ####

test_that("Test cross distToNearest with model hh_s5f", {
    ## Reproduce vignette
    db <- subset(db, CPRIMER %in% c("IGHA","IGHM") & 
                     BARCODE %in% c("RL016","RL018","RL019","RL021"))
    
    db_nrow <- nrow(db)
    db2 <- dplyr::bind_rows(db,db)
    db2$SAMPLE<- "S2"
    db2$SAMPLE[1:db_nrow] <- "S1"
    db2$DONOR <- "D1"
    
    test_idx <- c(1:5,10:15,300:305)
    
    ## Test hh_s5f 
    
    dist_hs5f <- distToNearest(db, sequenceColumn="JUNCTION",
                               vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                               model="hh_s5f", first=FALSE, normalize="none", nproc=1,
                               locusColumn="LOCUS")
    
    ## Test if the updated function reproduces results
    expect_equal(dist_hs5f$dist_nearest[test_idx],
                 c( NA, NA, NA, NA, 37.2436, 35.1829, 31.5680, 27.6521, 15.0128,
                    15.0128, 16.8262, 2.9351, 8.2891, 27.9191, 27.2217, 8.2891, 27.0871
                 ),
                 tolerance=0.005)
    
    ## There's only one donor, hence cross-donor will return all NA
    dist_hs5f_cross_donor <- distToNearest(db2, sequenceColumn="JUNCTION",
                                           vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                           model="hh_s5f", first=FALSE, normalize="none", nproc=1, cross="DONOR",
                                           locusColumn="LOCUS")
    expect_true(all(is.na(dist_hs5f_cross_donor$cross_dist_nearest)))
    
    ## fields=NULL and fields=DONOR should give same results
    cross_dist_hs5f <- distToNearest(db2, sequenceColumn="JUNCTION",
                                     vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                     model="hh_s5f", first=FALSE, normalize="none",nproc=1,
                                     cross="SAMPLE", locusColumn="LOCUS")
    cross_dist_hs5f_donor <- distToNearest(db2, sequenceColumn="JUNCTION",
                                           vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                           model="hh_s5f", first=FALSE, normalize="none",nproc=1,
                                           cross="SAMPLE", fields="DONOR", locusColumn="LOCUS")
    expect_equal(cross_dist_hs5f, cross_dist_hs5f_donor)
    
    expect_equal(cross_dist_hs5f$cross_dist_nearest[test_idx],
                 c( NA, NA, NA, NA, 37.2436, 35.1829, 31.5680, 27.6521, 15.0128,
                    15.0128, 16.8262, 2.9351, 8.2891, 27.9191, 27.2217, 8.2891, 27.0871
                 ),
                 tolerance=0.005)
    ## Check cross
    ## Crossing shoud reproduce the same results as not crossed
    ## because both donors have the same db
    expect_equal(cross_dist_hs5f$cross_dist_nearest[1:nrow(dist_hs5f)],
                 dist_hs5f$dist_nearest, tolerance=0.005)
    
    ## Inroduce row as Sample 3, very similar to rows 1 and 316
    ## This will be the best
    db3 <- db2
    db3[nrow(db3),] <- db3[1,]
    db3[nrow(db3),"JUNCTION"] <- sub("^T","A",db3[nrow(db3),"JUNCTION"])
    db3[nrow(db3),"SAMPLE"] <- "S3"
    # t(db3[c(1,316,nrow(db3)),c("JUNCTION","V_CALL","J_CALL","DONOR","SAMPLE")] )
    
    db2_1_316_630_hs5f <- distToNearest(db2[c(1,316,630),], sequenceColumn="JUNCTION",
                                        vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                        model="hh_s5f", first=FALSE, normalize="none",cross="SAMPLE",
                                        locusColumn="LOCUS")
    ## Exactly same seq, returns NA
    expect_equal(db2_1_316_630_hs5f$cross_dist_nearest,c(NA,NA,NA))
    
    ## One seq has been edited, will return distance values
    db3_1_316_630_hs5f <- distToNearest(db3[c(1,316,630),], sequenceColumn="JUNCTION",
                                        vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                        model="hh_s5f", first=FALSE, normalize="none",cross="SAMPLE",
                                        locusColumn="LOCUS")
    expect_equal(db3_1_316_630_hs5f$cross_dist_nearest,c(1.001,1.001,1.001), tolerance=0.005)
    
    seq1 <- c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", 
              "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
    seq2 <- c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", 
              "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
    
    targeting_distance <- calcTargetingDistance(HH_S5F)
    dist <- shazam:::dist5Mers(seq1, seq2, targeting_distance)
    expect_equal(dist, 1.0574, tolerance = 0.005)
    
    ## seq2[1] with a non valid character "S"
    ## Expect error
    seq2[1] <- "NNSCG"
    expect_error(shazam:::dist5Mers(seq1, seq2, targeting_distance), "Character not found")
    
    ## Length normalized
    dist_hs5f_len <- distToNearest(db, sequenceColumn="JUNCTION",
                                   vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                   model="hh_s5f", first=FALSE, normalize="len", nproc=1,
                                   locusColumn="LOCUS")
    
    expected_len_dist <- dist_hs5f$dist_nearest/nchar(dist_hs5f$JUNCTION)
    expect_equal(dist_hs5f_len$dist_nearest, expected_len_dist, tolerance=0.001)
    
})

#### distToNearest - unrecognized characters ####

test_that("Test distToNearest with unrecognized characters", {
    db2 <- subset(db, CPRIMER %in% c("IGHA","IGHM") &
                     BARCODE %in% c("RL016","RL018","RL019","RL021"))

    ## Create a junction with unrecognized char
    ## This seq belongs to clone of size 1
    db2$JUNCTION[3] <- gsub("T","Z",db2$JUNCTION[5] )
    expect_warning(
        dist_hs1f <- distToNearest(db2, sequenceColumn="JUNCTION",
                                   vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                   model="hh_s1f", first=FALSE, normalize="len",
                                   locusColumn="LOCUS"),
        "Invalid sequence characters"
    )
    
    ## Create a junction with unrecognized char
    ## This seq belongs to clone of size > 1
    db2$JUNCTION[5] <- gsub("T","Z",db$JUNCTION[5] )
    expect_warning(
        dist_hs1f <- distToNearest(db2, sequenceColumn="JUNCTION",
                                   vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                   model="hh_s1f", first=FALSE, normalize="len",
                                   locusColumn="LOCUS"),
        "Invalid sequence characters"
    )

})

#### distToNearest - model='aa'  ####
test_that("Test distToNearest with stop codon model='aa' ", {
    # Create a toy dataframe
    juncs <- c("ACGTACGTACGTACGTACGTACGTACGTATCGT", 
               "ACGTACGTACGTACGTACGTACGTACGTATAAT", 
               "ACGTACGTACGTACGTACGTACGTACGTATNNN", 
               "ACGTACGTACGTACGTACGTACGTACGTAT---",
               "TGGAACGTACGTACGTACGTACGTACGTACGTT",
               "---A--GTACGTACGTACGTACGTACGTAT---")
    vcall <- "Homsap IGHV3-49*03 F"
    jcall <- "Homsap IGHJ1*01 F"
    df <- data.frame(SEQUENCE_ID=c(1:length(juncs)),
                     V_CALL= rep(vcall, length(juncs)),
                     J_CALL= rep(jcall, length(juncs)),
                     JUNCTION=juncs,
                     stringsAsFactors=F)
    df$LOCUS <- getLocus(df$V_CALL)
    
    # calculate the ditances with normalization
    df <- distToNearest(df, sequenceColumn="JUNCTION", vCallColumn="V_CALL", 
                        jCallColumn="J_CALL", model="aa", locusColumn="LOCUS")
    expect_equal(df$dist_nearest, c(0.0303,0.0303,0.0606,0.0606,0.0606,NA))
    
    # calculate the ditances without normalization
    df <- distToNearest(df, sequenceColumn="JUNCTION", vCallColumn="V_CALL", 
                        jCallColumn="J_CALL", model="aa", normalize = "none",
                        locusColumn="LOCUS")
    expect_equal(df$dist_nearest, c(1,1,2,2,2,NA))
})

#### distToNearest - tibbles ####

test_that("Test distToNearest returns the same result with data.frame and tibble", {
    db2 <- subset(db, CPRIMER %in% c("IGHA","IGHM") &
                      BARCODE %in% c("RL016","RL018","RL019","RL021"))
    expect_equivalent(distToNearest(data.frame(db2), sequenceColumn="JUNCTION",
                                    vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                    locusColumn="LOCUS"),
                      distToNearest(tibble::as_tibble(db2), sequenceColumn="JUNCTION",
                                    vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                                    locusColumn="LOCUS"))
})

#### distToNearest - heavy and light chain ####

test_that("distToNearest, bulk H:L", {
    
    df <- data.frame(sequence_id=1:4,
                     v_call= c('IGHV1-01', 'IGHV1-01','IGKV1-01', 'IGKV1-01'),
                     j_call= c('IGHJ1-01','IGHJ1-01','IGKJ1-01','IGKJ1-01'),
                     junction=c(
                         'TGTAAAAAATGG',
                         'TGTAAACCCTGG',
                         'TGTCAAAAATGG',
                         'TGTCCCAAATGG'
                     ),
                     stringsAsFactors=F)
    df$locus <- getLocus(df$v_call)
    # pwd<-pairwiseDist(df$junction); colnames(pwd)<-rownames(pwd)<-df$junction
    
    dtn_h <- distToNearest(df, sequenceColumn="junction",
                  vCallColumn="v_call", jCallColumn="j_call",
                  locusColumn="locus", locusValues = "IGH")
    
    dtn_h_dfsubset <- distToNearest(df %>% filter(locus == "IGH"), sequenceColumn="junction",
                           vCallColumn="v_call", jCallColumn="j_call",
                           locusColumn="locus", locusValues = "IGH")
    
    # IGK not analyzed
    expect_equal(dtn_h$dist_nearest, c(0.25, 0.25, NA, NA))
                 
    dtn_hl <- distToNearest(df, sequenceColumn="junction",
                  vCallColumn="v_call", jCallColumn="j_call",
                  locusColumn="locus", locusValues = c("IGH","IGK"))
    
    # IGH and IGK analyzed
    expect_equal(dtn_hl$dist_nearest,c(0.25, 0.25, 0.1667, 0.1667 ))
    
    expect_equal(dtn_h[1:2,], dtn_hl[1:2,])
})


#### distToNearest - light chain ####

test_that("distToNearest, single-cell mode with VH:VL paired input", {
    
    # Load test database
    # data_t1: 2nd row contains 2 heavy chains within the same cell 
    #          (not allowed by distToNearest)
    # data_t2: removed multi heavy chain in 2nd row
    load(file.path("..", "data-tests", "db_sc.rda"))
    
    # expect error with data_t1 
    # cell "B" has 2 HCs (disallowed, as indicated in doc)
    expect_warning(expect_error(distToNearest(db=data_t1, sequenceColumn="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL", 
                               model="ham", normalize="len", symmetry="avg",
                               first=FALSE, VJthenLen=FALSE, keepVJLgroup=TRUE, 
                               cellIdColumn="CELL_ID", locusColumn="LOCUS", onlyHeavy=FALSE),
                               "multiple heavy chains found. One heavy chain per cell is expected."),
                               "onlyHeavy = FALSE is deprecated. Running as if onlyHeavy = TRUE")
    
    expect_error(distToNearest(db=data_t1, sequenceColumn="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL", 
                               model="ham", normalize="len", symmetry="avg",
                               first=FALSE, VJthenLen=FALSE, keepVJLgroup=TRUE, 
                               cellIdColumn="CELL_ID", locusColumn="LOCUS", onlyHeavy=TRUE),
                               "multiple heavy chains found. One heavy chain per cell is expected.")
    
    ## with data_t2, expect smooth run
    # Request group using HC VJL & LC VJL, but LC VJL is now deprected and will use onlyHevy=TRUE
    # first=F
    expect_warning(dtn1 <- distToNearest(db=data_t2, sequenceColumn="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL", 
                         model="ham", normalize="len", symmetry="avg",
                         first=FALSE, VJthenLen=FALSE, keepVJLgroup=TRUE, 
                         cellIdColumn="CELL_ID", locusColumn="LOCUS", onlyHeavy=FALSE),
                   "onlyHeavy = FALSE is deprecated. Running as if onlyHeavy = TRUE",
                   fixed=TRUE)
    # first=T
    # onlyHeavy=FALSe is deprecated
    expect_warning(dtn2 <- distToNearest(db=data_t2, sequenceColumn="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL", 
                         model="ham", normalize="len", symmetry="avg",
                         first=TRUE, VJthenLen=FALSE, keepVJLgroup=TRUE, 
                         cellIdColumn="CELL_ID", locusColumn="LOCUS", onlyHeavy=FALSE),
                         "onlyHeavy = FALSE is deprecated. Running as if onlyHeavy = TRUE",
                         fixed=TRUE)
    
    # group using HC VJL only, without LC VJL
    # first=F
    dtn3 <- distToNearest(db=data_t2, sequenceColumn="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL", 
                          model="ham", normalize="len", symmetry="avg",
                          first=FALSE, VJthenLen=FALSE, keepVJLgroup=TRUE, 
                          cellIdColumn="CELL_ID", locusColumn="LOCUS", onlyHeavy=TRUE)

    # first=T
    dtn4 <- distToNearest(db=data_t2, sequenceColumn="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL", 
                          model="ham", normalize="len", symmetry="avg",
                          first=TRUE, VJthenLen=FALSE, keepVJLgroup=TRUE, 
                          cellIdColumn="CELL_ID", locusColumn="LOCUS", onlyHeavy=TRUE)
    
    # Due to deprecation of onlyHeavy
    expect_identical(dtn1, dtn3)
    expect_identical(dtn2, dtn4)
    
    # calcualte expected
    
    dtn1_expect <- rep(NA, nrow(dtn1))
    for (i in 1:nrow(dtn1)) {
        
        if (dtn1[["LOCUS"]][i]=="IGH") {
            curGrp <- dtn1[["vjl_group"]][i]
            curIdx <- which(dtn1[["vjl_group"]]==curGrp & dtn1[["LOCUS"]]=="IGH")
            curJunc <- dtn1[["JUNCTION"]][i]
            
            if (length(curIdx)==1) {
                dtn1_expect[i] <- NA
            } else {
                curJuncsAll <- dtn1[["JUNCTION"]][curIdx]
                stopifnot(length(unique(nchar(curJuncsAll)))==1)
                # non-self, non-identical
                bool <- curIdx!=i & curJuncsAll!=curJunc
                if (!any(bool)) {
                    dtn1_expect[i] = NA
                } else {
                    curJuncsCf = curJuncsAll[bool]
                    tmpDists <- rep(NA, length(curJuncsCf))
                    for (j in 1:length(curJuncsCf)) {
                        tmpDists[j] <- alakazam::seqDist(seq1=curJunc, seq2=curJuncsCf[j])
                    }
                    dtn1_expect[i] <- min(tmpDists / nchar(curJuncsAll)[1])
                }
            }
        }
        # if not IGH, left as NA
        
    }
    
    expect_equal(dtn1_expect, dtn1[["dist_nearest"]], tolerance=0.0001)
    
    dtn2_expect = rep(NA, nrow(dtn2))
    
    for (i in 1:nrow(dtn2)) {
        
        if (dtn2[["LOCUS"]][i]=="IGH") {
            curGrp <- dtn2[["vjl_group"]][i]
            curIdx <- which(dtn2[["vjl_group"]]==curGrp & dtn2[["LOCUS"]]=="IGH")
            curJunc <- dtn2[["JUNCTION"]][i]
            
            if (length(curIdx)==1) {
                dtn2_expect[i] <- NA
            } else {
                curJuncsAll <- dtn2[["JUNCTION"]][curIdx]
                stopifnot(length(unique(nchar(curJuncsAll)))==1)
                # non-self, non-identical
                bool <- curIdx!=i & curJuncsAll!=curJunc
                if (!any(bool)) {
                    dtn2_expect[i] <- NA
                } else {
                    curJuncsCf <- curJuncsAll[bool]
                    tmpDists <- rep(NA, length(curJuncsCf))
                    for (j in 1:length(curJuncsCf)) {
                        tmpDists[j] <- alakazam::seqDist(seq1=curJunc, seq2=curJuncsCf[j])
                    }
                    dtn2_expect[i] <- min(tmpDists / nchar(curJuncsAll)[1])
                }
            }
        }
        # if not IGH, left as NA
    }
    
    expect_equal(dtn2_expect, dtn2[["dist_nearest"]], tolerance=0.0001)
    
    dtn3_expect = rep(NA, nrow(dtn3))
    
    for (i in 1:nrow(dtn3)) {
        
        if (dtn3[["LOCUS"]][i]=="IGH") {
            curGrp <- dtn3[["vjl_group"]][i]
            curIdx <- which(dtn3[["vjl_group"]]==curGrp & dtn3[["LOCUS"]]=="IGH")
            curJunc <- dtn3[["JUNCTION"]][i]
            
            if (length(curIdx)==1) {
                dtn3_expect[i] <- NA
            } else {
                curJuncsAll <- dtn3[["JUNCTION"]][curIdx]
                stopifnot(length(unique(nchar(curJuncsAll)))==1)
                # non-self, non-identical
                bool <- curIdx!=i & curJuncsAll!=curJunc
                if (!any(bool)) {
                    dtn3_expect[i] <- NA
                } else {
                    curJuncsCf <- curJuncsAll[bool]
                    tmpDists <- rep(NA, length(curJuncsCf))
                    for (j in 1:length(curJuncsCf)) {
                        tmpDists[j] <- alakazam::seqDist(seq1=curJunc, seq2=curJuncsCf[j])
                    }
                    dtn3_expect[i] <- min(tmpDists / nchar(curJuncsAll)[1])
                }
            }
        }
        # if not IGH, left as NA
    }
    
    expect_equal(dtn3_expect, dtn3[["dist_nearest"]], tolerance=0.0001)
    
    dtn4_expect = rep(NA, nrow(dtn4))
    
    for (i in 1:nrow(dtn4)) {
        
        if (dtn4[["LOCUS"]][i]=="IGH") {
            curGrp <- dtn4[["vjl_group"]][i]
            curIdx <- which(dtn4[["vjl_group"]]==curGrp & dtn4[["LOCUS"]]=="IGH")
            curJunc <- dtn4[["JUNCTION"]][i]
            
            if (length(curIdx)==1) {
                dtn4_expect[i] <- NA
            } else {
                curJuncsAll <- dtn4[["JUNCTION"]][curIdx]
                stopifnot(length(unique(nchar(curJuncsAll)))==1)
                # non-self, non-identical
                bool <- curIdx!=i & curJuncsAll!=curJunc
                if (!any(bool)) {
                    dtn4_expect[i] <- NA
                } else {
                    curJuncsCf <- curJuncsAll[bool]
                    tmpDists <- rep(NA, length(curJuncsCf))
                    for (j in 1:length(curJuncsCf)) {
                        tmpDists[j] <- alakazam::seqDist(seq1=curJunc, seq2=curJuncsCf[j])
                    }
                    dtn4_expect[i] <- min(tmpDists / nchar(curJuncsAll)[1])
                }
            }
        }
        # if not IGH, left as NA
    }
    
    expect_equal(dtn4_expect, dtn4[["dist_nearest"]], tolerance=0.0001)
    
})

#### findThreshold ####

test_that("Test findThreshold", {
    
    db <- distToNearest(db, sequenceColumn="JUNCTION",
                        vCallColumn="V_CALL", jCallColumn="J_CALL", 
                        model="ham", first=FALSE, normalize="len", nproc=1,
                        locusColumn="LOCUS")
    
    gmm_output <- findThreshold(db$dist_nearest, method="gmm", model="gamma-gamma", edge=0.9)
    expect_equal(gmm_output@threshold, 0.12, tolerance=0.01)

    dens_output <- findThreshold(db$dist_nearest, method="dens")
    expect_equal(dens_output@threshold, 0.14, tolerance=0.01)
})

#### calcTargetingDistance ####
    
test_that("Test distance, Change-O tests", {
    
    # aux function. Given a distance matrix,
    # get the min distance value for each row
    .getMin <- function(mat) {
        apply(mat,1,function(x) {
            gt0 <- which(x>0)
            if (length(gt0>0)) {
                min(x[gt0])
            } else {
                NA
            }
        })
    }
    
    nt_seq <- c("TGTGCAAGGGGGCCA",
                "TGTGCAAGGGGGCCA",
                "TGTATTTGGGGGCCA",
                "ACACTTGCCACTGAT",
                "NNNNNNNNNNNNTGA",
                "NNNNNNNNNNNNNNN")
    nt_seq_short <- c("TGTGCAAGG",
                     "TGTGCAAGG",
                     "TGTATTTGG",
                     "NNNNNNNNN")
    aa_seq <- alakazam::translateDNA(nt_seq)

    # Amino acid distance from sequence 1 to sequence 1,..,5
    aa_len <- 5.0
    aa_ham <- c(0.0, 0.0, 2.0, 5.0, 1, 0.0)    
    
    # Nucleotide distances sequence 1 to sequence 1,..,6
    nt_len <- 15.0
    # 1 vs 3
    #   A-C = 0; A-G = 1; A-T = 2;
    #   C-G = 0; C-T = 1; G-T = 0
    # 1 vs 4
    #   A-C = 1; A-G = 2; A-T = 4;
    #   C-G = 6; C-T = 1; G-T = 1
    # 1 vs 5
    #   A-C = 0; A-G = 0; A-T = 0;
    nt_ham <- c(0.0, 0.0, 4.0, 15.0, 2.0, 0.0)
    nt_hh_s1f <- c(0.0,
                   0.0,
                   0.0 + 1*0.64 + 2*1.16 + 0.0 + 1*0.64 + 0.0,
                   1*1.21 + 2*0.64 + 4*1.16 + 6*1.16 + 1*0.64 + 1*1.21,
                   0.0 + 0.0 + 0.0 + 0.0 + 1*1.16 + 1*0.64 + 0.0,
                   0.0)
    nt_mk_rs1nf <- c(0.0,
                     0.0,
                     0.0 + 1*0.32 + 2*1.17 + 0.0 + 1*0.32 + 0.0,
                     1*1.51 + 2*0.32 + 4*1.17 + 6*1.17 + 1*0.32 + 1*1.51,
                     0.0 + 0.0 + 0.0 + 0.0 + 1*1.17 + 1*0.32 + 0.0,
                     0.0)
    
    # 5-mer models use shorter sequence
    nt_short_len <- 9.0
    
    HH_S1F_Distance <- calcTargetingDistance(model=HH_S1F)
    HH_S5F_Distance <- calcTargetingDistance(model=HH_S5F)
    MK_RS5NF_Distance <- calcTargetingDistance(model=MK_RS5NF)
    MK_RS1NF_Distance <- calcTargetingDistance(model=MK_RS1NF)
    
    #   GTGCA-GTATT = [0.97, 0.84]
    #   TGCAA-TATTT = [0.93, 0.83]
    #   GCAAG-ATTTG = [1.08, 1.16]
    #   CAAGG-TTTGG = [0.91, 1.07]
    nt_hh_s5f_avg <- c(0.0,
                       0.0,
                       mean(c(0.97, 0.84)) + 
                           mean(c(0.93, 0.83)) +
                           mean(c(1.08, 1.16)) +
                           mean(c(0.91, 1.07)),
                       0.0)
    
    nt_hh_s5f_min <- c(0.0,
                       0.0,
                           min(c(0.97, 0.84)) +
                           min(c(0.93, 0.83)) +
                           min(c(1.08, 1.16)) + 
                           min(c(0.91, 1.07)),
                       0.0)
    
    expect_equal(shazam:::dist5Mers( "GTGCA", "GTATT",HH_S5F_Distance, symmetry="raw"),
        c(0.97, 0.84))
    expect_equal(shazam:::dist5Mers( "TGCAA", "TATTT",HH_S5F_Distance, symmetry="raw"),
                 c(0.93, 0.83))
    expect_equal(shazam:::dist5Mers( "GCAAG", "ATTTG",HH_S5F_Distance, symmetry="raw"),
                 c(1.08, 1.16))
    expect_equal(shazam:::dist5Mers( "CAAGG", "TTTGG",HH_S5F_Distance, symmetry="raw"),
                 c(0.91, 1.07))
    
    seq1 <- c( "GTGCA", "TGCAA", "GCAAG", "CAAGG")
    seq2 <-  c("GTATT", "TATTT", "ATTTG", "TTTGG")
    shazam:::dist5Mers( seq1, seq2, HH_S5F_Distance, symmetry="raw")
    
    seq1_to_seq2 <- (diag(HH_S5F_Distance[substr(seq2, 3, 3), seq1]))
    seq1_to_seq2
    seq2_to_seq1 <- (diag(HH_S5F_Distance[substr(seq1, 3, 3), seq2]))
    seq2_to_seq1
    
    # mk_rs5nf
    #   GTGCA-GTATT = [0.71, 0.77]
    #   TGCAA-TATTT = [0.71, 0.93]
    #   GCAAG-ATTTG = [1.05, 1.03]
    #   CAAGG-TTTGG = [1.08, 1.13]
    nt_mk_rs5nf_avg <- c(0.0,
                         0.0,
                         mean(c(0.71, 0.77)) +
                             mean(c(0.71, 0.93)) +
                             mean(c(1.05, 1.03)) + 
                             mean(c(1.08, 1.13)),
                         0.0)
    nt_mk_rs5nf_min <- c(0.0,
                        0.0,
                        min(c(0.71, 0.77)) + min(c(0.71, 0.93)) +
                        min(c(1.05, 1.03)) + min(c(1.08, 1.13)),
                        0.0)
    
    # aa ham
    expect_equivalent(pairwiseDist(aa_seq, getAAMatrix())[1,],
                      aa_ham)
    
    # aa ham with cross
    expect_equal(
        nearestDist(nt_seq, model="aa"),
        nearestDist(nt_seq, model="aa",crossGroups = 1:length(nt_seq))
    )
    
    # nt    
    expect_equivalent(pairwiseDist(nt_seq, getDNAMatrix())[1,],
                      nt_ham)

    # nt len
    expect_equivalent(.getMin(pairwiseDist(nt_seq, getDNAMatrix())),
                 shazam:::nearestDist(nt_seq))
    expect_equal(.getMin(pairwiseDist(nt_seq, getDNAMatrix())/nt_len),
                      shazam:::nearestDist(nt_seq, model="ham", normalize = "len"),
                      tolerance=0.001, check.attributes = FALSE)
    
    # hh_s1f
    expect_equivalent(pairwiseDist(nt_seq, HH_S1F_Distance)[1,],
                      nt_hh_s1f)
    
    # hh_s1f len 
    expect_equivalent(.getMin(pairwiseDist(nt_seq, HH_S1F_Distance)),
                          shazam:::nearestDist(nt_seq, model="hh_s1f"))
    
    expect_equal(.getMin(pairwiseDist(nt_seq, HH_S1F_Distance)/nt_len),
                 shazam:::nearestDist(nt_seq, model="hh_s1f", normalize = "len"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # mk_rs1nf
    expect_equivalent(pairwiseDist(nt_seq, MK_RS1NF_Distance)[1,],
                      nt_mk_rs1nf)
    
    # mk_rs1nf len
    expect_equivalent(.getMin(pairwiseDist(nt_seq, MK_RS1NF_Distance)),
                      shazam:::nearestDist(nt_seq, model="mk_rs1nf"))
    
    expect_equal(.getMin(pairwiseDist(nt_seq, MK_RS1NF_Distance)/nt_len),
                 shazam:::nearestDist(nt_seq, model="mk_rs1nf", normalize = "len"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # hh_s5f avg
    hh_s5f_avg_dm <- shazam:::pairwise5MerDist(nt_seq_short, 
                                               HH_S5F_Distance, symmetry="avg")
    expect_equivalent(hh_s5f_avg_dm[1,],
                      nt_hh_s5f_avg)
    
    expect_equal(.getMin(hh_s5f_avg_dm),
                 shazam:::nearestDist(nt_seq_short, model="hh_s5f", 
                                      symmetry="avg", normalize = "none"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # hh_s5f avg len
    expect_equal(.getMin(hh_s5f_avg_dm/nt_short_len),
                 shazam:::nearestDist(nt_seq_short, model="hh_s5f", 
                                      symmetry="avg", normalize = "len"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # hh_s5f min
    hh_s5f_min_dm <- shazam:::pairwise5MerDist(nt_seq_short, 
                              HH_S5F_Distance, symmetry="min")
    expect_equivalent(hh_s5f_min_dm[1,],
                      nt_hh_s5f_min)  
    expect_equal(.getMin(hh_s5f_min_dm),
                 shazam:::nearestDist(nt_seq_short, model="hh_s5f", 
                                      symmetry="min", normalize = "none"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # hh_s5f min len
    expect_equal(.getMin(hh_s5f_min_dm/nt_short_len),
                 shazam:::nearestDist(nt_seq_short, model="hh_s5f", 
                                      symmetry="min", normalize = "len"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # mk_rs5nf min
    mk_rs5nf_min_dm <- shazam:::pairwise5MerDist(nt_seq_short, 
                              MK_RS5NF_Distance, symmetry="min")
    expect_equivalent(mk_rs5nf_min_dm[1,],
                      nt_mk_rs5nf_min) 
    expect_equal(.getMin(mk_rs5nf_min_dm),
                 shazam:::nearestDist(nt_seq_short, model="mk_rs5nf", 
                                      symmetry="min", normalize = "none"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # mk_rs5nf min len
    expect_equal(.getMin(mk_rs5nf_min_dm/nt_short_len),
                 shazam:::nearestDist(nt_seq_short, model="mk_rs5nf", 
                                      symmetry="min", normalize = "len"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # mk_rs5nf avg
    mk_rs5nf_avg_dm <- shazam:::pairwise5MerDist(nt_seq_short, 
                              MK_RS5NF_Distance, symmetry="avg")
    expect_equivalent(mk_rs5nf_avg_dm[1,],
                      nt_mk_rs5nf_avg)    
    expect_equal(.getMin(mk_rs5nf_avg_dm),
                 shazam:::nearestDist(nt_seq_short, model="mk_rs5nf", 
                                      symmetry="avg", normalize = "none"),
                 tolerance=0.001, check.attributes = FALSE)
    
    # mk_rs5nf avg len TODO
    expect_equal(.getMin(mk_rs5nf_avg_dm/nt_short_len),
                 shazam:::nearestDist(nt_seq_short, model="mk_rs5nf", 
                                      symmetry="avg", normalize = "len"),
                 tolerance=0.001, check.attributes = FALSE)
})

test_that("distToNearest fields applied before gene groups", {
    
    db <- data.frame(
        subject_id=c("S1","S1","S1","S2","S2"),
        v_call=c("IGHV1-1*01","IGHV1-1*01","IGHV1-2*01","IGHV1-1*01, IGHV1-2*01","IGHV1-2*01" ),
        j_call=c("IGHJ1*01","IGHJ1*01","IGHJ1*01","IGHJ1*01","IGHJ1*01"),
        junction=c("TGTAAAAAATGG","TGTAAAAAATGG","TGTAAAACCTGG","TGTAAACCCTGG","TGTAAACCCTGG"),
        locus=c("IGH","IGH","IGH","IGH","IGH"),
        cell_id=1:5
    )
    
    expect_error(dtn <- distToNearest(db, first=F, fields="subject_id", VJthenLen = F),
                   "specify the cell_id")
    
    dtn <- distToNearest(db, first=F, fields="subject_id", VJthenLen = F, cellId="cell_id")

    # S1 should have 2 groups, and S2 one.
    expect_equal(length(unique(dtn$vjl_group)), 3)
})

### distToNearest mix sc-bulk

test_that("distToNearest onlyHeavy=FALSE is deprecated", {
    expect_warning(distToNearest(db, cellIdColumn = NULL, first=F, onlyHeavy=FALSE,
                  sequenceColumn = "JUNCTION",
                  vCallColumn = "V_CALL", jCallColumn = "J_CALL", locusColumn = "LOCUS"),
                  "onlyHeavy = FALSE is deprecated. Running as if onlyHeavy = TRUE")
})

test_that("distToNearest mix single cell and bulk", {

    db <- data.frame(
         subject_id=c("S1","S1","S1","S1","S1", "S1", "S1", "S2", "S3", "S3"),
         v_call=c("IGHV1-1*01","IGHV1-1*01","IGHV1-2*01","IGHV1-1*01,IGHV1-2*01","IGHV1-2*01", "IGKV1-1*01", "IGKV1-1*01", "IGKV1-1*01", "IGHV1-1*01", "IGKV1-1*01"),
         j_call=c("IGHJ1*01","IGHJ1*01","IGHJ1*01","IGHJ1*01","IGHJ1*01","IGKJ1*01", "IGKJ1*01", "IGKJ1*01", "IGHJ1*01", "IGKJ1*01"),
         junction=c("TGTAAAAAATGG","TGTAAAAATTGG","TGTAAAACCTGG","TGTAAACCCTGG","TGTAAACCCTGG","TGTCCCCCCTGG","TGTCCCCCGTGG", "TGTCCCCAATGG", "TGTAAAAAATGG", "TGTCCCCCCTGG"),
         locus=c("IGH","IGH","IGH","IGH","IGH","IGK","IGK","IGK", "IGH","IGK"),
         cell_id=c(1,2,3,NA,NA,1,NA, NA, 4, 4),
         junction_length=12
    )
    # mixed
    # bulk light chain sequences are not analyzed because they
    # are not grouped by groupGenes.
    # alakazam::groupGenes(db, cell_id = "cell_id", first=F)
    # IGH specified
    expect_message(
        expect_warning( 
            dtn_m_h <- distToNearest(db, 
                                            cellIdColumn = "cell_id", first=F, onlyHeavy=TRUE,
                          sequenceColumn = "junction",
                          vCallColumn = "v_call", jCallColumn = "j_call", 
                          locusColumn = "locus", locusValues=c("IGH")),
            "The vj_group column contains NA values"
        ),
        "Running in single-cell mode..+Note: [0-9]+ light/short.*$", 
        fixed=FALSE
    )
    expect_equal(dtn_m_h$dist_nearest, 
                 c(0.0833, 0.0833, 0.0833, 0.0833, 0.0833, NA, NA, NA, 0.0833, NA))
    
    #IGH and IGK specified
    expect_message(
        expect_warning(dtn_m_hk <- distToNearest(db, cellIdColumn = "cell_id", first=F, onlyHeavy=TRUE,
                            sequenceColumn = "junction",
                            vCallColumn = "v_call", jCallColumn = "j_call", 
                            locusColumn = "locus", locusValues=c("IGH","IGK")),
                    "The vj_group column contains NA values"),
        "Running in single-cell mode..+Note: [0-9]+ light/short.*Note: locusValues.*$", 
        fixed=FALSE
    )
    # the single-cell light chains should have NA dist_nearest
    expect_equal(dtn_m_hk$dist_nearest, 
                 c(0.0833, 0.0833, 0.0833, 0.0833, 0.0833, NA, NA, NA, 0.0833, NA)       
    )

    # IGK specified
    expect_message(
        expect_warning(dtn_m_k <- distToNearest(db, cellIdColumn = "cell_id", first=F, onlyHeavy=TRUE,
                            sequenceColumn = "junction",
                            vCallColumn = "v_call", jCallColumn = "j_call", 
                            locusColumn = "locus", locusValues=c("IGK")),
                    "The vj_group column contains NA values"),
        "Running in single-cell mode..+Note: [0-9]+ light/short.*Note: locusValues.*$", 
        fixed=FALSE
    )
    # All NA dist_nearest
    expect_equal(dtn_m_k$dist_nearest, 
                rep(NA, nrow(dtn_m_k))      
    )
    
    # sc only
    # alakazam::groupGenes(db %>% dplyr::filter(!is.na(cell_id)), 
    #                      cell_id = "cell_id", first=F)
    # IGH specified
    expect_message(
        dtn_sc_h <- distToNearest(db %>% dplyr::filter(!is.na(cell_id)), 
                          cellIdColumn = "cell_id", first=F, onlyHeavy=TRUE,
                          sequenceColumn = "junction",
                          vCallColumn = "v_call", jCallColumn = "j_call", 
                          locusColumn = "locus", locusValues=c("IGH")),
        "Running in single-cell mode..+Note: [0-9]+ light/short.*$", 
        fixed=FALSE
    )
    expect_equal(
        dtn_sc_h$dist_nearest, c(0.0833, 0.0833, NA,NA, 0.0833, NA)
    )
    
    # IGH and IGK specified
    expect_message(
        dtn_sc_hk <- distToNearest(db %>% dplyr::filter(!is.na(cell_id)), 
                        cellIdColumn = "cell_id", first=F, onlyHeavy=TRUE,
                        sequenceColumn = "junction",
                        vCallColumn = "v_call", jCallColumn = "j_call", 
                        locusColumn = "locus", locusValues=c("IGH","IGK")),    
        "Running in single-cell mode..+Note: [0-9]+ light/short.*Note: locusValues.*$", 
        fixed=FALSE
    )
    expect_equal(
        dtn_sc_hk$dist_nearest, c(0.0833, 0.0833, NA,NA, 0.0833, NA)
    )

    # IGK specified
    expect_message(
        dtn_sc_k <- distToNearest(db %>% dplyr::filter(!is.na(cell_id)), 
                            cellIdColumn = "cell_id", first=F, onlyHeavy=TRUE,
                            sequenceColumn = "junction",
                            vCallColumn = "v_call", jCallColumn = "j_call", 
                            locusColumn = "locus", locusValues=c("IGK")),
    "Running in single-cell mode..+Note: [0-9]+ light/short.*Note: locusValues.*$", 
    fixed=FALSE
    )
    expect_equal(
        dtn_sc_k$dist_nearest, c(NA, NA, NA,NA, NA, NA)
    )    
    
    # bulk only
    # alakazam::groupGenes(db %>% dplyr::mutate(cell_id=NULL), 
    #                      cell_id = NULL, first=F)
    # IGH specified
    expect_message(
        dtn_b_h <- distToNearest(db %>% dplyr::mutate(cell_id=NULL), 
                          cellIdColumn = NULL, first=F, onlyHeavy=TRUE,
                          sequenceColumn = "junction",
                          vCallColumn = "v_call", jCallColumn = "j_call", 
                          locusColumn = "locus", locusValues=c("IGH")),
        "Running in non-single-cell mode", 
        fixed=FALSE
    )
    expect_equal(
        dtn_b_h$dist_nearest,
        c(0.0833, 0.0833, 0.0833, 0.0833, 0.0833, NA, NA, NA, 0.0833, NA)
    )
    
    # IGH and IGK specified
    expect_message(
        dtn_b_hk <- distToNearest(db %>% dplyr::mutate(cell_id=NULL), 
                        cellIdColumn = NULL, first=F, onlyHeavy=TRUE,
                        sequenceColumn = "junction",
                        vCallColumn = "v_call", jCallColumn = "j_call", 
                        locusColumn = "locus", locusValues=c("IGH", "IGK")),
        "Running in non-single-cell mode.", 
        fixed=FALSE
    )
    expect_equal(
        dtn_b_hk$dist_nearest,
        c(0.0833, 0.0833, 0.0833, 0.0833, 0.0833,  0.0833,  0.0833, 0.1667, 0.0833, 0.0833)
    )

    # IGK specified
    expect_message(
        dtn_b_k <- distToNearest(db %>% dplyr::mutate(cell_id=NULL), 
                            cellIdColumn = NULL, first=F, onlyHeavy=TRUE,
                            sequenceColumn = "junction",
                            vCallColumn = "v_call", jCallColumn = "j_call", 
                            locusColumn = "locus", locusValues=c("IGK")),
        "Running in non-single-cell mode.", 
        fixed=FALSE
    )
    expect_equal(
        dtn_b_k$dist_nearest,
        c(NA, NA, NA, NA, NA, 0.0833,  0.0833, 0.1667, NA, 0.0833)
    )    
})

#### AIRR migration tests ####

test_that("distToNearest & findThreshold, AIRR migration", {
    
    # ExampleDb
    load(file.path("..", "data-tests", "ExampleDb.rda")) 
    # ExampleDb_airr
    load(file.path("..", "data-tests", "ExampleDb_airr.rda")) 
    
    db_c <- subset(ExampleDb, SAMPLE == "-1h")
    db_a <- subset(ExampleDb_airr, sample == "-1h")
    
    # distToNearest
    
    dist_c <- distToNearest(db_c, sequenceColumn="JUNCTION",
                            vCallColumn="V_CALL_GENOTYPED", jCallColumn="J_CALL",
                            model="ham", first=FALSE, VJthenLen=TRUE, normalize="len",
                            locusColumn="LOCUS")
    
    dist_a <- distToNearest(db_a, sequenceColumn="junction", 
                            vCallColumn="v_call_genotyped", jCallColumn="j_call",
                            model="ham", first=FALSE, VJthenLen=TRUE, normalize="len",
                            locusColumn="locus")
    
    expect_equal(dist_c[["dist_nearest"]], dist_a[["dist_nearest"]])
    
    # findThreshold, gmm
    
    # gmm method seems to be stochastic 
    # using set.seed does not help
    
    #set.seed(2945)
    #o_gmm_c <- findThreshold(dist_c[["dist_nearest"]], method="gmm", model="gamma-gamma", cutoff="user", spc=0.99)
    #set.seed(2945)
    #o_gmm_a <- findThreshold(dist_a[["dist_nearest"]], method="gmm", model="gamma-gamma", cutoff="user", spc=0.99)
    
    #expect_equal(o_gmm_c@threshold, o_gmm_a@threshold)
    
    # findThreshold, density
    o_den_c <- findThreshold(dist_c[["dist_nearest"]], method="density")
    o_den_a <- findThreshold(dist_a[["dist_nearest"]], method="density")
    
    expect_equal(o_den_c@xdens, o_den_a@xdens)
    expect_equal(o_den_c@ydens, o_den_a@ydens)
    expect_equal(o_den_c@threshold, o_den_a@threshold)
})

