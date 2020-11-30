
test_that("makeClonesList", {
    load(file.path("..", "data-tests", "ExampleDb.rda")) 
    expect_equal(length(makeClonesList(ExampleDb, clone_col="CLONE"))[1],1199)
    expect_equal(makeClonesList(subset(ExampleDb,CLONE %in% c(30,125,115,467)), clone_col="CLONE"), 
                 as.list(c("30","115","125","467")))
})

test_that("makeChangeoCloneCurClone", {
    data("ExampleDb", package = "alakazam")
    clone_3128_db <- subset(ExampleDb, clone_id == 3128)
    clone_3128_changeO <- makeChangeoClone(clone_3128_db, seq="sequence_alignment", germ="germline_alignment")
    clone_3128_changeO_direct <- makeChangeoCloneCurClone(ExampleDb, cur_clone_num = 3128, 
                                                         germ="germline_alignment", 
                                                         junc_len = "junction_length",
                                                         seq="sequence_alignment")
    expect_equal(dim(clone_3128_db)[1], 100)
    expect_equal(clone_3128_changeO@junc_len, 60)
    expect_equal(clone_3128_changeO@data$collapse_count[37:40],c(7, 46, 2, 3))
    expect_equal(sum(clone_3128_changeO@data$collapse_count),dim(clone_3128_db)[1])
    expect_equal(clone_3128_changeO_direct@junc_len, 60)
    expect_equal(clone_3128_changeO_direct@v_gene,clone_3128_changeO@v_gene)
    expect_equal(clone_3128_changeO_direct@j_gene,clone_3128_changeO@j_gene)
    expect_equal(clone_3128_changeO_direct@germline,clone_3128_changeO@germline)
    expect_equal(clone_3128_changeO_direct@junc_len,clone_3128_changeO@junc_len)
    expect_equal(clone_3128_changeO_direct@clone,clone_3128_changeO@clone)
    # comparing first 3 columns, as rest of columns exist only in clone_3128_changeO_direct@data
    expect_equal(clone_3128_changeO_direct@data[,1:3],clone_3128_changeO@data[,1:3])
    
    expect_equal(clone_3128_changeO_direct@data$v_call,
                 replicate(length(clone_3128_changeO_direct@data$v_call),clone_3128_changeO@v_gene))
    expect_equal(clone_3128_changeO_direct@data$j_call,
                 replicate(length(clone_3128_changeO_direct@data$j_call),clone_3128_changeO@j_gene))
    expect_equal(clone_3128_changeO_direct@data$junction_length,
                 replicate(length(clone_3128_changeO_direct@data$junction_length),clone_3128_changeO@junc_len))
})

test_that("makeGraphCurClone", {
    data("ExampleDb", package="alakazam")
    clone_3141_obj <- makeChangeoCloneCurClone(db=ExampleDb, cur_clone_num=3141, 
                                              seq="sequence_alignment", 
                                              junc_len = "junction_length",
                                              germ="germline_alignment") 
    expect_equal(dim(clone_3141_obj@data)[1],17)
    expect_equal(sum(clone_3141_obj@data$collapse_count),44)
    dnapars_exec <- "~/dummy/phylip-3.69/dnapars"
    # dnapars_exec <- "c:\\Users\\milcat\\phylip-3.698\\exe\\dnapars.exe"
    if (file.access(dnapars_exec, mode=1) == -1) {
        expect_error(
            clone_3141_graph <- makeGraphCurClone(clone_3141_obj,dnapars_exec),
            "The file ~/dummy/phylip-3.69/dnapars cannot be executed"
        )
    } else {
        clone_3141_graph <- makeGraphCurClone(clone_3141_obj,dnapars_exec)
        expect_that(clone_3141_graph, is_a("igraph"))
        expect_equal(igraph::vcount(clone_3141_graph),18)
        expect_equal(igraph::ecount(clone_3141_graph),17)
    }
    
    clone_3184_obj <- makeChangeoCloneCurClone(db=ExampleDb, cur_clone_num=3184, 
                                              seq="sequence_alignment", 
                                              junc_len = "junction_length",
                                              germ="germline_alignment") 
    if (file.access(dnapars_exec, mode=1) == -1) {
        expect_error(
            clone_3184_graph <- makeGraphCurClone(clone_3184_obj,dnapars_exec),
            "The file ~/dummy/phylip-3.69/dnapars cannot be executed"
        )
    } else {
        clone_3184_graph <- makeGraphCurClone(clone_3184_obj,dnapars_exec)
        expect_true("Germline" %in% V(clone_3184_graph)$name)
        # Check that there are some Inferred sequences here:
        expect_true(any(grepl("Inferred", V(clone_3184_graph)$name)))
    }
})

test_that("makeGraphDf", {
    data("ExampleDb", package="alakazam")
    dnapars_exec <- "~/dummy/phylip-3.69/dnapars"
    #dnapars_exec <- "c:\\Users\\milcat\\phylip-3.698\\exe\\dnapars.exe"
    clone_3177_obj <- makeChangeoCloneCurClone(db=ExampleDb, cur_clone_num=3177, 
                                              seq="sequence_alignment", 
                                              junc_len = "junction_length",
                                              germ="germline_alignment") 
    
    if (file.access(dnapars_exec, mode=1) == -1) {
        expect_error(
            clone_3177_graph <- makeGraphCurClone(clone_3177_obj,dnapars_exec),
            "The file ~/dummy/phylip-3.69/dnapars cannot be executed"
        )
    } else {
        clone_3177_graph <- makeGraphCurClone(clone_3177_obj,dnapars_exec)
        clone_3177_graphDF <- makeGraphDf(curCloneGraph=clone_3177_graph,  
                                          curCloneObj=clone_3177_obj, objSeqId="sequence_id",
                                          objSeq="sequence")
        expect_equal(dim(clone_3177_graphDF), c(16,22))
        expect_equal(sum(clone_3177_graphDF$collapse_count, na.rm = T),30)
    } 
})

test_that("makeRegion", {
    data("ExampleDb", package="alakazam")
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

