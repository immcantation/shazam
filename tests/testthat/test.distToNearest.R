# Load test database
e1 <- new.env()
#load(file.path("tests", "data-tests", "TestDb.rda"), envir=e1)
load(file.path("..", "data-tests", "TestDb.rda"), envir=e1)
db <- get("TestDb", envir=e1)
rm(e1)

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
    
    dist_hs1f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
                                    model="hh_s1f", first=FALSE, normalize="len")
    ## Test if the updated function reproduces results
    expect_equal(dist_hs1f$DIST_NEAREST[test_idx],
                 c(NA,NA,NA,NA,0.4040, 0.4447, 0.3963, 0.3469, 0.3050, 0.3050,
                   0.4284, 0.0435, 0.1212, 0.3771, 0.3862, 0.1212, 0.3687),
                 tolerance=0.005)
    
    ## There's only one donor, hence cross-donor will return all NA
    dist_hs1f_cross_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                           model="hh_s1f", first=FALSE, normalize="none", nproc=1, cross="DONOR")

    expect_true(all(is.na(dist_hs1f_cross_donor$CROSS_DIST_NEAREST)))
    
    ## fields=NULL and fields=DONOR should give same results
    cross_dist_hs1f <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                     model="hh_s1f", first=FALSE, normalize="len",nproc=1,
                                     cross="SAMPLE")
    
    
    cross_dist_hs1f_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                               model="hh_s1f", first=FALSE, normalize="len",nproc=1,
                               cross="SAMPLE", fields="DONOR")
    
    expect_equal(cross_dist_hs1f, cross_dist_hs1f_donor)
    expect_equal(cross_dist_hs1f$CROSS_DIST_NEAREST[test_idx],
                 c(NA,NA,NA,NA,0.4040, 0.4447, 0.3963, 0.3469, 0.3050, 0.3050,
                   0.4284, 0.0435, 0.1212, 0.3771, 0.3862, 0.1212, 0.3687),
                 tolerance=0.005)
    ## Check cross
    ## Crossing shoud reproduce the same results as not crossed
    ## because both donors have the same db
    expect_equal(cross_dist_hs1f$CROSS_DIST_NEAREST[1:nrow(dist_hs1f)],
                 dist_hs1f$DIST_NEAREST, tolerance=0.005)
    
    ## Inroduce row as Sample 3, very similar to rows 1 and 316
    ## This will be the best
    db3 <- db2
    db3[nrow(db3),] <- db3[1,]
    db3[nrow(db3),"JUNCTION"] <- sub("^T","A",db3[nrow(db3),"JUNCTION"])
    db3[nrow(db3),"SAMPLE"] <- "S3"
    # t(db3[c(1,316,nrow(db3)),c("JUNCTION","V_CALL","J_CALL","DONOR","SAMPLE")] )
    
    db2_1_316_630 <- distToNearest(db2[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                  model="hh_s1f", first=FALSE, normalize="len",cross="SAMPLE")
    ## Exactly same seq, returns NA
    expect_equal(db2_1_316_630$CROSS_DIST_NEAREST,c(NA,NA,NA))
    
    ## One seq has been edited, will return distance values
    db3_1_316_630 <- distToNearest(db3[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                             model="hh_s1f", first=FALSE, normalize="len",cross="SAMPLE")
    expect_equal(db3_1_316_630$CROSS_DIST_NEAREST,c(0.0175,0.0175,0.0175), tolerance=0.005)
})

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
    
    dist_hs5f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
                               model="hh_s5f", first=FALSE, normalize="none", nproc=1)

    ## Test if the updated function reproduces results
    expect_equal(dist_hs5f$DIST_NEAREST[test_idx],
                 c( NA, NA, NA, NA, 37.2436, 35.1829, 31.5680, 27.6521, 15.0128,
                    15.0128, 16.8262, 2.9351, 8.2891, 27.9191, 27.2217, 8.2891, 27.0871
                 ),
                 tolerance=0.005)
    
    ## There's only one donor, hence cross-donor will return all NA
    dist_hs5f_cross_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                               model="hh_s5f", first=FALSE, normalize="none", nproc=1, cross="DONOR")
    expect_true(all(is.na(dist_hs5f_cross_donor$CROSS_DIST_NEAREST)))
    
    ## fields=NULL and fields=DONOR should give same results
    cross_dist_hs5f <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                     model="hh_s5f", first=FALSE, normalize="none",nproc=1,
                                     cross="SAMPLE")
    cross_dist_hs5f_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                           model="hh_s5f", first=FALSE, normalize="none",nproc=1,
                                           cross="SAMPLE", fields="DONOR")
    expect_equal(cross_dist_hs5f, cross_dist_hs5f_donor)
    
    expect_equal(cross_dist_hs5f$CROSS_DIST_NEAREST[test_idx],
                 c( NA, NA, NA, NA, 37.2436, 35.1829, 31.5680, 27.6521, 15.0128,
                    15.0128, 16.8262, 2.9351, 8.2891, 27.9191, 27.2217, 8.2891, 27.0871
                 ),
                 tolerance=0.005)
    ## Check cross
    ## Crossing shoud reproduce the same results as not crossed
    ## because both donors have the same db
    expect_equal(cross_dist_hs5f$CROSS_DIST_NEAREST[1:nrow(dist_hs5f)],
                 dist_hs5f$DIST_NEAREST, tolerance=0.005)
    
    ## Inroduce row as Sample 3, very similar to rows 1 and 316
    ## This will be the best
    db3 <- db2
    db3[nrow(db3),] <- db3[1,]
    db3[nrow(db3),"JUNCTION"] <- sub("^T","A",db3[nrow(db3),"JUNCTION"])
    db3[nrow(db3),"SAMPLE"] <- "S3"
    # t(db3[c(1,316,nrow(db3)),c("JUNCTION","V_CALL","J_CALL","DONOR","SAMPLE")] )
    
    db2_1_316_630_hs5f <- distToNearest(db2[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                                   model="hh_s5f", first=FALSE, normalize="none",cross="SAMPLE")
    ## Exactly same seq, returns NA
    expect_equal(db2_1_316_630_hs5f$CROSS_DIST_NEAREST,c(NA,NA,NA))
    
    ## One seq has been edited, will return distance values
    db3_1_316_630_hs5f <- distToNearest(db3[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                                   model="hh_s5f", first=FALSE, normalize="none",cross="SAMPLE")
    expect_equal(db3_1_316_630_hs5f$CROSS_DIST_NEAREST,c(1.001,1.001,1.001), tolerance=0.005)
    
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
    dist_hs5f_len <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
                               model="hh_s5f", first=FALSE, normalize="len", nproc=1)
    
    expected_len_dist <- dist_hs5f$DIST_NEAREST/nchar(dist_hs5f$JUNCTION)
    expect_equal(dist_hs5f_len$DIST_NEAREST, expected_len_dist, tolerance=0.001)
    
})

test_that("Test distToNearest with unrecognized characters", {
    db2 <- subset(db, CPRIMER %in% c("IGHA","IGHM") &
                     BARCODE %in% c("RL016","RL018","RL019","RL021"))

    ## Create a junction with unrecognized char
    ## This seq belongs to clone of size 1
    db2$JUNCTION[3] <- gsub("T","Z",db2$JUNCTION[5] )
    expect_warning(
        dist_hs1f <- distToNearest(db2, model="hh_s1f", first=FALSE, normalize="len"),
        "Invalid sequence characters"
    )
    
    ## Create a junction with unrecognized char
    ## This seq belongs to clone of size > 1
    db2$JUNCTION[5] <- gsub("T","Z",db$JUNCTION[5] )
    expect_warning(
        dist_hs1f <- distToNearest(db2, model="hh_s1f", first=FALSE, normalize="len"),
        "Invalid sequence characters"
    )

})

test_that("Test findThreshold", {
    
    db <- distToNearest(db, model="ham", first=FALSE, normalize="len", nproc=1)
    
    gmm_output <- findThreshold(db$DIST_NEAREST, method="gmm", cutEdge=0.9)
    expect_equal(gmm_output@threshold, 0.0896, tolerance=0.01)

    dens_output <- findThreshold(db$DIST_NEAREST, method="dens", cutEdge=0.9)
    expect_equal(dens_output@threshold, 0.114, tolerance=0.01)
})
