test_that("makeRegion", {
    # data("ExampleDb", package="alakazam")
    load(file.path("..", "data-tests", "ExampleDb_airr.rda"))
    junction_length_10 <- ExampleDb[10,"junction_length"]
    sequence_imgt_10 <- ExampleDb[10,"sequence_alignment"]
    seq_10_reg_def <- makeRegion(juncLength = junction_length_10, 
                                 sequenceImgt = sequence_imgt_10, 
                                 regionDefinition = IMGT_VDJ_BY_REGIONS)
    junction_length_45 <- ExampleDb[45,"junction_length"]
    sequence_imgt_45 <- ExampleDb[45,"sequence_alignment"]
    seq_45_reg_def <- makeRegion(juncLength = junction_length_45, 
                                 sequenceImgt = sequence_imgt_45, 
                                 regionDefinition = IMGT_VDJ)
    junction_length_333 <- ExampleDb[333,"junction_length"]
    sequence_imgt_333 <- ExampleDb[333,"sequence_alignment"]
    seq_reg_def_IMGT_V <- makeRegion(juncLength = junction_length_333, 
                                         sequenceImgt = sequence_imgt_333, 
                                         regionDefinition = IMGT_V)
    seq_reg_def_IMGT_V_BY_CODONS <- makeRegion(juncLength = junction_length_10, 
                                                   sequenceImgt = sequence_imgt_10, 
                                                   regionDefinition = IMGT_V_BY_CODONS)
    seq_reg_def_IMGT_V_BY_REGIONS <- makeRegion(juncLength = junction_length_45, 
                                                   sequenceImgt = sequence_imgt_45, 
                                                   regionDefinition = IMGT_V_BY_REGIONS)
    seq_reg_def_IMGT_V_BY_SEGMENTS <- makeRegion(juncLength = junction_length_333, 
                                                 sequenceImgt = sequence_imgt_333, 
                                                 regionDefinition = IMGT_V_BY_SEGMENTS)
    
    expect_that(seq_10_reg_def, is_a("RegionDefinition"))
    expect_equal(seq_10_reg_def@name, "IMGT_VDJ_BY_REGIONS")
    expect_equal(seq_45_reg_def@name, "IMGT_VDJ")
    expect_equal(seq_10_reg_def@description, 
                 "IMGT numbering scheme defining the V(D)J segment by individual cdr1/2/3 and fwr1/2/3/4 regions")
    expect_equal(seq_45_reg_def@description, 
                 "IMGT numbering scheme defining the V(D)J segment only by cdr/fwr regions, including cdr3 and fwr4")
    expect_equal(as.character(seq_10_reg_def@boundaries[1:312]), 
                 as.character(IMGT_V_BY_REGIONS@boundaries[1:312]))
    expect_equal(as.character(seq_45_reg_def@boundaries[1:312]), 
                 as.character(IMGT_V@boundaries[1:312]))
    expect_equal(as.character(seq_10_reg_def@boundaries[313:387]), rep("cdr3",75))
    expect_equal(as.character(seq_45_reg_def@boundaries[313:360]), rep("cdr",48))
    expect_equal(as.character(seq_10_reg_def@boundaries[388:nchar(sequence_imgt_10)]), 
                 rep("fwr4",nchar(sequence_imgt_10)-387))
    expect_equal(as.character(seq_45_reg_def@boundaries[361:nchar(sequence_imgt_45)]), 
                 rep("fwr",nchar(sequence_imgt_45)-360))
    expect_equal(seq_10_reg_def@seqLength,unname(nchar(sequence_imgt_10)))
    expect_equal(seq_45_reg_def@seqLength,unname(nchar(sequence_imgt_45)))
    expect_equal(seq_10_reg_def@regions,
                 c("cdr1", "cdr2", "cdr3", "fwr1", "fwr2", "fwr3", "fwr4"))
    expect_equal(seq_45_reg_def@regions, c("cdr", "fwr"))
    expect_equal(seq_10_reg_def@labels, 
                 c("cdr1_r", "cdr1_s", "cdr2_r", "cdr2_s", "cdr3_r", "cdr3_s", 
                   "fwr1_r", "fwr1_s", "fwr2_r", "fwr2_s", "fwr3_r", "fwr3_s", 
                   "fwr4_r", "fwr4_s"))
    expect_equal(seq_45_reg_def@labels, c("cdr_r", "cdr_s", "fwr_r", "fwr_s"))
    expect_equal(seq_10_reg_def@citation, IMGT_V_BY_REGIONS@citation)
    expect_equal(seq_10_reg_def@citation, IMGT_V_BY_REGIONS@citation)
    
    expect_equal(seq_reg_def_IMGT_V,IMGT_V)
    expect_equal(seq_reg_def_IMGT_V_BY_CODONS,IMGT_V_BY_CODONS)
    expect_equal(seq_reg_def_IMGT_V_BY_REGIONS,IMGT_V_BY_REGIONS)
    expect_equal(seq_reg_def_IMGT_V_BY_SEGMENTS,IMGT_V_BY_SEGMENTS)
})

test_that("getCloneRegion", {
    data("ExampleDb")
    clone_2834_reg <- getCloneRegion(clone_num = 2834, db = ExampleDb, 
                                   seq_col = "sequence_alignment", clone_col = "clone_id", 
                                   regionDefinition = IMGT_VDJ_BY_REGIONS)
    seq_343 <- ExampleDb[343,]
    seq_343_reg_def <- makeRegion(juncLength = 48, 
                                 sequenceImgt = seq_343[,"germline_alignment"], 
                                 regionDefinition = IMGT_VDJ_BY_REGIONS)
    
    clone_3227_reg <- getCloneRegion(clone_num = 3227, db = ExampleDb, 
                                     seq_col = "sequence_alignment", clone_col = "clone_id", 
                                     regionDefinition = IMGT_VDJ)
    seq_376 <- ExampleDb[376,]
    seq_376_reg_def <- makeRegion(juncLength = 60, 
                                  sequenceImgt = seq_376[,"germline_alignment"], 
                                  regionDefinition = IMGT_VDJ)
    
    clone_6940_reg <- getCloneRegion(clone_num = 6940, db = ExampleDb, 
                                     seq_col = "sequence_alignment", clone_col = "clone_id", 
                                     regionDefinition = IMGT_V_BY_REGIONS)
    seq_793 <- ExampleDb[793,]
    seq_793_reg_def <- makeRegion(juncLength = 84, 
                                  sequenceImgt = seq_793[,"germline_alignment"], 
                                  regionDefinition = IMGT_V_BY_REGIONS)
    
    expect_equal(as.numeric(seq_343[,"clone_id"]),2834)
    expect_equal(as.numeric(seq_343[,"junction_length"]),48)
    expect_that(clone_2834_reg, is_a("RegionDefinition"))
    expect_equal(clone_2834_reg,seq_343_reg_def)
    
    expect_equal(as.numeric(seq_376[,"clone_id"]),3227)
    expect_equal(as.numeric(seq_376[,"junction_length"]),60)
    expect_that(clone_3227_reg, is_a("RegionDefinition"))
    expect_equal(clone_3227_reg,seq_376_reg_def)
    
    expect_equal(as.numeric(seq_793[,"clone_id"]),6940)
    expect_equal(as.numeric(seq_793[,"junction_length"]),84)
    expect_that(clone_6940_reg, is_a("RegionDefinition"))
    expect_equal(clone_6940_reg,seq_793_reg_def)
})

