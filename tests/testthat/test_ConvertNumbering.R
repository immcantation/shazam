#### ConvertNumbering ####
library(testthat)
test_that("Ensure numberings are consistent and stable", {
    ## Reproduce example
    imgt_seqs <- c("51", "23", "31", "18", "4", "90", "58", "59", "3")
    kabat_seqs <-c("*55", "82B", "91", "23")
    known_miss_imgt <- c("73", "110", "114", "115")
    known_arb_kabat <- c("G", "*")
    
    # Should be able to retain original numberings
    expect_equal(convertNumbering("IGH", "KABAT", "KABAT", kabat_seqs), kabat_seqs)
    expect_equal(convertNumbering("IGH", "IMGT", "IMGT", imgt_seqs), imgt_seqs)
    
    kabat_imgt_expected <- c("63", "93", "103", "24")
    imgt_kabat_expected <- c("46", "22", "*30", "17", "4", "81", "*52A", "*52B", "3")
    known_miss_expected <- c("G", "G", "G", "G")
    known_arb_expected <- c("NA", "NA")
    
    # Expected conversions
    expect_equal(convertNumbering("IGH", "KABAT", "IMGT", kabat_seqs),kabat_imgt_expected)
    expect_equal(convertNumbering("IGH", "IMGT", "KABAT", imgt_seqs),imgt_kabat_expected)
    expect_equal(convertNumbering("IGH", "IMGT", "KABAT", known_miss_imgt),known_miss_expected)
    expect_equal(convertNumbering("IGH", "KABAT", "IMGT", known_arb_kabat),known_arb_expected)
    
    # Missing values or incorrect formats
    issue_seqs <- c("55", "23", "52", "82D")
    expect_error(convertNumbering("IGH", "KABAT", "IMGT", issue_seqs), 
                 "Formatting of following characters does not match reference:  55, 82D", fixed=TRUE)
})

