computeBaselinePDF <- function(db){
  pdfs <- computeBayesianScore(db[,c("OBSERVED_R_CDR", "OBSERVED_S_CDR",
                                     "OBSERVED_R_FWR", "OBSERVED_S_FWR",
                                     "EXPECTED_R_CDR", "EXPECTED_S_CDR",
                                     "EXPECTED_R_FWR", "EXPECTED_S_FWR")],
                               test="Focused", max_sigma=20,length_sigma=4001)
  db$CDR_PDFs <- pdfs[[1]]
  db$FWR_PDFs <- pdfs[[2]]
  return(db)
}
