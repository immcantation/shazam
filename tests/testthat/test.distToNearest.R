test_that("Test cross distToNearest with model hs1f", {
    ## Reproduce example
    db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
                     BARCODE %in% c("RL016","RL018","RL019","RL021"))
    db_nrow <- nrow(db)
    db2 <- dplyr::bind_rows(db,db)
    db2$SAMPLE<- "S2"
    db2$SAMPLE[1:db_nrow] <- "S1"
    db2$DONOR <- "D1"
    
    test_idx <- c(1:5,10:15,300:305)
    
    ## Test hs1f 
    
    dist_hs1f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
                                    model="hs1f", first=FALSE, normalize="length")
    dist_hs1f_rcpp <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
                               model="hs1f", first=FALSE, normalize="length", rcpp=T)
    expect_equal(dist_hs1f, dist_hs1f_rcpp)
    ## Test if the updated function reproduces results
    expect_equal(dist_hs1f$DIST_NEAREST[test_idx],
                 c(NA,NA,NA,NA,0.6396,0.7027,0.6479,0.5563,0.4879,
                   0.4879,0.6998,0.0700,0.1955,0.6071,0.6145,0.1955,0.5974),
                 tolerance=0.001)
    
    ## There's only one donor, hence cross-donor will return all NA
    dist_hs1f_cross_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                           model="hs1f", first=FALSE, normalize="none", nproc=1, cross="DONOR")
    dist_hs1f_cross_donor_rcpp <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                           model="hs1f", first=FALSE, normalize="none", nproc=1, cross="DONOR",
                                           rcpp=T) 
    expect_equal(dist_hs1f_cross_donor, dist_hs1f_cross_donor_rcpp)
    
    expect_true(all(is.na(dist_hs1f_cross_donor$CROSS_DIST_NEAREST)))
    
    ## fields=NULL and fields=DONOR should give same results
    cross_dist_hs1f <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                     model="hs1f", first=FALSE, normalize="length",nproc=1,
                                     cross="SAMPLE")
    
    
    cross_dist_hs1f_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                               model="hs1f", first=FALSE, normalize="length",nproc=1,
                               cross="SAMPLE", fields="DONOR")
    
    expect_equal(cross_dist_hs1f, cross_dist_hs1f_donor)
    expect_equal(cross_dist_hs1f$CROSS_DIST_NEAREST[test_idx],
                 c(NA,NA,NA,NA,0.6396,0.7027,0.6479,0.5563,0.4879,
                   0.4879,0.6998,0.0700,0.1955,0.6071,0.6145,0.1955,0.5974),
                 tolerance=0.001)
    ## Check cross
    ## Crossing shoud reproduce the same results as not crossed
    ## because both donors have the same db
    expect_equal(cross_dist_hs1f$CROSS_DIST_NEAREST[1:nrow(dist_hs1f)],
                 dist_hs1f$DIST_NEAREST, tolerance=0.001)
    
    ## Inroduce row as Sample 3, very similar to rows 1 and 316
    ## This will be the best
    db3 <- db2
    db3[nrow(db3),] <- db3[1,]
    db3[nrow(db3),"JUNCTION"] <- sub("^T","A",db3[nrow(db3),"JUNCTION"])
    db3[nrow(db3),"SAMPLE"] <- "S3"
    # t(db3[c(1,316,nrow(db3)),c("JUNCTION","V_CALL","J_CALL","DONOR","SAMPLE")] )
    
    db2_1_316_630 <- distToNearest(db2[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                  model="hs1f", first=FALSE, normalize="length",cross="SAMPLE")
    ## Exactly same seq, returns NA
    expect_equal(db2_1_316_630$CROSS_DIST_NEAREST,c(NA,NA,NA))
    
    ## One seq has been edited, will return distance values
    db3_1_316_630 <- distToNearest(db3[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                             model="hs1f", first=FALSE, normalize="length",cross="SAMPLE")
    expect_equal(db3_1_316_630$CROSS_DIST_NEAREST,c(0.0265,0.0265,0.0265), tolerance=0.001)
})

test_that("Test cross distToNearest with model hs5f", {
    ## Reproduce vignette
    db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
                     BARCODE %in% c("RL016","RL018","RL019","RL021"))
    
    db_nrow <- nrow(db)
    db2 <- dplyr::bind_rows(db,db)
    db2$SAMPLE<- "S2"
    db2$SAMPLE[1:db_nrow] <- "S1"
    db2$DONOR <- "D1"
    
    test_idx <- c(1:5,10:15,300:305)
    
    ## Test hs1f 
    
    dist_hs5f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
                               model="hs5f", first=FALSE, normalize="none", nproc=1)

    ## Test if the updated function reproduces results
    expect_equal(dist_hs5f$DIST_NEAREST[test_idx],
                 c( NA, NA, NA, NA, 37.2436, 35.1829, 31.5680, 27.6521, 15.0128,
                    15.0128, 16.8262, 2.9351, 8.2891, 27.9191, 27.2217, 8.2891, 27.0871
                 ),
                 tolerance=0.001)
    
    ## There's only one donor, hence cross-donor will return all NA
    dist_hs5f_cross_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                               model="hs5f", first=FALSE, normalize="none", nproc=1, cross="DONOR")
    expect_true(all(is.na(dist_hs5f_cross_donor$CROSS_DIST_NEAREST)))
    
    ## fields=NULL and fields=DONOR should give same results
    cross_dist_hs5f <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                     model="hs5f", first=FALSE, normalize="none",nproc=1,
                                     cross="SAMPLE")
    cross_dist_hs5f_donor <- distToNearest(db2, vCallColumn="V_CALL_GENOTYPED", 
                                           model="hs5f", first=FALSE, normalize="none",nproc=1,
                                           cross="SAMPLE", fields="DONOR")
    expect_equal(cross_dist_hs5f, cross_dist_hs5f_donor)
    
    expect_equal(cross_dist_hs5f$CROSS_DIST_NEAREST[test_idx],
                 c( NA, NA, NA, NA, 37.2436, 35.1829, 31.5680, 27.6521, 15.0128,
                    15.0128, 16.8262, 2.9351, 8.2891, 27.9191, 27.2217, 8.2891, 27.0871
                 ),
                 tolerance=0.001)
    ## Check cross
    ## Crossing shoud reproduce the same results as not crossed
    ## because both donors have the same db
    expect_equal(cross_dist_hs5f$CROSS_DIST_NEAREST[1:nrow(dist_hs5f)],
                 dist_hs5f$DIST_NEAREST, tolerance=0.001)
    
    ## Inroduce row as Sample 3, very similar to rows 1 and 316
    ## This will be the best
    db3 <- db2
    db3[nrow(db3),] <- db3[1,]
    db3[nrow(db3),"JUNCTION"] <- sub("^T","A",db3[nrow(db3),"JUNCTION"])
    db3[nrow(db3),"SAMPLE"] <- "S3"
    # t(db3[c(1,316,nrow(db3)),c("JUNCTION","V_CALL","J_CALL","DONOR","SAMPLE")] )
    
    db2_1_316_630_hs5f <- distToNearest(db2[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                                   model="hs5f", first=FALSE, normalize="none",cross="SAMPLE")
    ## Exactly same seq, returns NA
    expect_equal(db2_1_316_630_hs5f$CROSS_DIST_NEAREST,c(NA,NA,NA))
    
    ## One seq has been edited, will return distance values
    db3_1_316_630_hs5f <- distToNearest(db3[c(1,316,630),], vCallColumn="V_CALL_GENOTYPED", 
                                   model="hs1f", first=FALSE, normalize="none",cross="SAMPLE")
    expect_equal(db3_1_316_630_hs5f$CROSS_DIST_NEAREST,c(1.75,1.75,1.75), tolerance=0.001)
})
