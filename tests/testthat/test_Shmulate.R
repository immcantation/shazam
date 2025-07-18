load(file.path("..", "data-tests", "ExampleTrees.rda"))

#### shmulateSeq ####
test_that("sampleMut throws error if vector of positions longer than sequence.", {
    sim_seq <- "ATA"
    sim_leng <- 3
    targetingModel<-HH_S5F
    targeting <- calculateTargeting(germlineSeq = sim_seq, targetingModel = targetingModel)
    targeting <- targeting[NUCLEOTIDES[1:4], ]
    targeting[is.na(targeting)] <- 0
    mutation_types <- computeMutationTypes(sim_seq)
    targeting[mutation_types=="stop"] <- 0
    positions <- numeric(4)
    expect_error(sampleMut(sim_leng, targeting, positions), 
                 regexp = "The vector of positions is longer than the length of the sequence")
})

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
    
    # numMutations must be a whole number
    expect_error(shmulateSeq("ATA", numMutations=1.5), regexp="whole number")
    
    # numMutations can't be larger than sequence length
    expect_error(shmulateSeq("ATA", numMutations=4), 
                 regexp="larger than the length of the sequence")
    
    # This will throw error, because G is trimmed
    expect_warning(expect_error(shmulateSeq("ATCG", numMutations=4), 
                 regexp="larger than the length of the sequence"),
                 regexp="Trimming sequence to last codon")
    
    # This will throw error, because only two positions used to add mutations
    expect_warning(expect_error(shmulateSeq("ATCG", numMutations=2, start=1, end=2), 
                 regexp="larger than the length of the sequence"),
                 regexp="Trimming sequence to last codon")
    
    # This will throw warning, because 4 positions requested to add mutations,
    # from 1 to 4, but nt 4 will be removed when trimming sequence to the last
    # codon
    expect_warning(sim_seq <- shmulateSeq("ATCG", numMutations=3, start=1, end=4), 
                 regexp="Trimming sequence to last codon")    
    expect_true(grepl("G$", sim_seq))
    expect_false(grepl("^ATC", sim_seq))
    expect_equal(nchar(sim_seq), nchar("ATCG"))
    
    # Throw error if unexpected codons found
    expect_error(shmulateSeq("ABC", numMutations=3), 
                 regexp="Unrecognized codons found")
    expect_warning(
        expect_error(shmulateSeq("atc", numMutations=3), 
                 regexp="Unrecognized codons found"),
        regexp = "shazam is case sensitive"
    )
})

test_that("Test shmulateSeq with frequency is accurate", {
    
    cat("\nTest shmulateSeq with frequency\n")
    
    # Example input sequence
    input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"
    freq <- 0.01
    n_reps <- 1000
    
    # Simulate using the default human S5F targeting model
    set.seed(56)
    
    num_muts <- 0
    initial_seq <- strsplit(input_seq,"")[[1]]
    for (i in 1:n_reps){
        output <- shmulateSeq(input_seq, numMutations = freq, frequency = T)
        o <- strsplit(output, "")[[1]]
        num_muts <- num_muts + sum(initial_seq != o)

    }
    expect_equal(num_muts, n_reps*nchar(input_seq)*freq, tolerance = 20)
})

#### shmulateTree ####

test_that("Test shmulateTree", { 
    
    cat("\nTest shmulateTree\n")
    
    input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"
    graph <- ExampleTrees[[17]]
    set.seed(7)
    tree <- shmulateTree(input_seq, graph, targetingModel = MK_RS5NF)
    expect_equal(tree$name, c("Inferred1", "GN5SHBT07JDYW5", "GN5SHBT03EP4KC",
                              "GN5SHBT01AKANC", "GN5SHBT01A3SFZ", "GN5SHBT08HUU7M",
                              "GN5SHBT04CEA6I", "GN5SHBT06IXJIH"))
    expect_equal(tree$distance, c(0, 2, 3, 4, 6, 3, 1, 1))
    expect_equal(tree$sequence, 
                 c( "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA", 
                    "NGATCTGACGACAGGGCCGTGTTTTACTGTGCGAGAGATAGTTTA", 
                    "NGATCTGACGACAGGGCCGTGTTATATTGTGCGAGAGATAATTTA", 
                    "NGATCTGACGACACGGCCGTATATTACTGTGCGAAATATAATTTA", 
                    "NGATCTGACGACAGGGCCTTATTCTACTGTACGAAAGATAATTTA",
                    "NGATCTGACGACAGGGCGGTATTATATTGTACGAGAGATAATTTA", 
                    "NGATCTGACGACAGGGCCGTATTTTACTGTGCGAGAGATAGTTTA", 
                    "NGATCTGACGACAGGGCCGTGTTTTACTGTGAGAGAGATAGTTTA"))
})
