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
