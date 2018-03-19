load(file.path("..", "data-tests", "ExampleTrees.rda"))

#### shmulateSeq ####

test_that("Test shmulateSeq", {
    
    cat("\nTest shmulateSeq\n")
    
    # Example input sequence
    input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"
     
    # Simulate using the default human S5F targeting model

    set.seed(56)
    output <- shmulateSeq(input_seq, numMutations = 6)

    expected <- "NAATTTGACGACACGGCCGTGGATGACTGTGCGAGAGATGCTTTA"   
    expect_equal(output, expected)
    
    set.seed(56)
    output <- shmulateSeq(input_seq, numMutations = 6, targetingModel = HH_S5F)
    expect_equal(output, expected)
    
    i <- strsplit(input_seq,"")[[1]]
    o <- strsplit(expected, "")[[1]]

    expect_equal(sum(i != o), 6)
})

#### shmulateTree ####

test_that("Test shmulateTree", { 
    
    cat("\nTest shmulateTree\n")
    
    input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"
    graph <- ExampleTrees[[17]]
    set.seed(7)
    tree <- shmulateTree(input_seq, graph, targetingModel = MK_RS5NF)
    expect_equal(tree$NAME, c("Inferred1", "GN5SHBT07JDYW5", "GN5SHBT03EP4KC",
                              "GN5SHBT01AKANC", "GN5SHBT01A3SFZ", "GN5SHBT08HUU7M",
                              "GN5SHBT04CEA6I", "GN5SHBT06IXJIH"))
    expect_equal(tree$DISTANCE, c(0, 2, 3, 4, 6, 3, 1, 1))
    expect_equal(tree$SEQUENCE, 
                 c( "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA", 
                    "NGATCTGACGACAGGGCCGTGTTTTACTGTGCGAGAGATAGTTTA", 
                    "NGATCTGACGACAGGGCCGTGTTATATTGTGCGAGAGATAATTTA", 
                    "NGATCTGACGACACGGCCGTATATTACTGTGCGAAATATAATTTA", 
                    "NGATCTGACGACAGGGCCTTATTCTACTGTACGAAAGATAATTTA",
                    "NGATCTGACGACAGGGCGGTATTATATTGTACGAGAGATAATTTA", 
                    "NGATCTGACGACAGGGCCGTATTTTACTGTGCGAGAGATAGTTTA", 
                    "NGATCTGACGACAGGGCCGTGTTTTACTGTGAGAGAGATAGTTTA"))
})
