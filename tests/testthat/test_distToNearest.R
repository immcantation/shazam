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
    
    dist_hs1f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
                               model="hh_s1f", first=FALSE, normalize="len")
    ## Test if the updated function reproduces results
    expect_equal(dist_hs1f$DIST_NEAREST[test_idx],
                 c(NA, NA, NA, NA, 
                   0.4060, 0.4465, 0.3976, 0.3482, 0.3064, 0.3064, 0.4295, 
                   0.0436, 0.1216, 0.3786, 0.3875, 0.1216, 0.3700),
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

#### distToNearest - unrecognized characters ####

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
                     JUNCTION=juncs)
    
    # calculate the ditances with normalization
    df <- distToNearest(df, model="aa")
    expect_equal(df$DIST_NEAREST, c(0.0303,0.0303,0.0606,0.0606,0.0606,NA))
    
    # calculate the ditances without normalization
    df <- distToNearest(df, model="aa", normalize = "none")
    expect_equal(df$DIST_NEAREST, c(1,1,2,2,2,NA))
})

#### distToNearest - tibbles ####

test_that("Test distToNearest returns the same result with data.frame and tibble", {
    db2 <- subset(db, CPRIMER %in% c("IGHA","IGHM") &
                      BARCODE %in% c("RL016","RL018","RL019","RL021"))
    expect_equivalent(distToNearest(data.frame(db2)),
                      distToNearest(tibble::as_tibble(db2)))
})

#### findThreshold ####

test_that("Test findThreshold", {
    
    db <- distToNearest(db, model="ham", first=FALSE, normalize="len", nproc=1)
    
    gmm_output <- findThreshold(db$DIST_NEAREST, method="gmm", model="gamma-gamma", edge=0.9)
    expect_equal(gmm_output@threshold, 0.12, tolerance=0.01)

    dens_output <- findThreshold(db$DIST_NEAREST, method="dens")
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
