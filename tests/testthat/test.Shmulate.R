load(file.path("..", "data-tests", "ExampleTrees.rda"))

test_that("Test shmulateSeq", {
    
    # Example input sequence
    input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"
     
    # Simulate using the default human S5F targeting model

    set.seed(56)
    output <- shmulateSeq(input_seq, num_muts = 6)

    expected <- "NAATTTGACGACACGGCCGTGGATGACTGTGCGAGAGATGCTTTA"   
    expect_equal(output, expected)
    
    set.seed(56)
    output <- shmulateSeq(input_seq, num_muts = 6, targeting_model = HS5FModel)
    expect_equal(output, expected)
    
    i <- strsplit(input_seq,"")[[1]]
    o <- strsplit(expected, "")[[1]]

    expect_equal(sum(i != o), 6)
})

test_that("Test shmulateTree", { 
    input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"
    graph <- ExampleTrees[[17]]
    set.seed(7)
    tree <- shmulateTree(input_seq, graph, targeting_model = MRS5NFModel)
    expect_equal(tree$distance, c(0, 2, 4, 3, 6, 1, 1, 3))
    expect_equal(tree$sequence, 
                 c("NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA", 
                   "NGATCTGACGACACGGCCGTGTATTACTGTGCGAAAGATAGTGTA", 
                   "NGATCTGACGACACGGCCGCATATTACTGTGCGAAAGATAATTTA", 
                   "NGATCTGACGACACGGCCGTGTATTACTGTGCGAAAAAGAATGTA",
                   "NGATCTGACGACACGGCCATATATTACCATACGAAAGATAGTATA", 
                   "NGATCTGACGACACGGCCGTGTATTACTGTGCAAAAGATAGTGTA", 
                   "NGATCTGACGACACGGCCGTATATTACTGTGCGAAAGATAGTGTA", 
                   "NGATCTGACGACACGACCGTGTATTACTGTGCGAAGAAGAATCTA"))
})