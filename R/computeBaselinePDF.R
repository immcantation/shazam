#' Compute the BASELINe (selection?) probability density function.
#'
#' \code{computeBaselinePDF} computes the BASELINe probability density functions for...
#'
#' @param    db  data.frame containing sequence data. \code{db} must contain the columns
#'               A, B, C, D.
#' @return   A modified \code{db} data.frame with added columns E, F, G.
#'
#' @references
#' Uduman M, Yaari G, Hershberg U, Stern JA, Shlomchik MJ, Kleinstein SH.
#'   Detecting selection in immunoglobulin sequences.
#'   Nucleic Acids Res. 2011 Jul;39(Web Server issue):W499-504.
#'
#' Yaari G, Uduman M, Kleinstein SH.
#'   Quantifying selection in high-throughput Immunoglobulin sequencing data sets.
#'   Nucleic Acids Res. 2012 Sep 1;40(17):e134.
#' @seealso  \code{\link{addObservedMutations}}, \code{\link{addExpectedFrequencies}}, \code{\link{calcGroupedBaseline}}.
#' @examples
#' # TODO
#' # Working example
#'
#' @export
computeBaselinePDF <- function(db,numbOfCores=1){

  availableCores <- getNumbOfCores()
  if(!(numbOfCores<=availableCores))numbOfCores=availableCores
  cluster <- makeCluster(numbOfCores, type = "SOCK")
  registerDoSNOW(cluster)
  numbOfSeqs <- nrow(db)
  pdfs <-
    foreach(i=icount(numbOfSeqs), .packages='shm', .combine=doparProgressBar(n=numbOfSeqs), .multicombine=TRUE) %dopar% {

      vecObsExpected <- as.numeric(db[i,c("OBSERVED_R_CDR", "OBSERVED_S_CDR", "OBSERVED_R_FWR", "OBSERVED_S_FWR",
                                          "EXPECTED_R_CDR", "EXPECTED_S_CDR", "EXPECTED_R_FWR", "EXPECTED_S_FWR")])

      #CDR
      binomialP = vecObsExpected[5]/sum(vecObsExpected[c(5,6,8)])
      binomialN = sum(vecObsExpected[c(1,2,4)])
      binomialX = vecObsExpected[1]
      cdrPDF <- calculate_bayes(x=binomialX,n=binomialN,p=binomialP,max_sigma=20,length_sigma=4001)

      #FWR
      binomialP = vecObsExpected[7]/sum(vecObsExpected[c(7,8,6)])
      binomialN = sum(vecObsExpected[c(3,4,2)])
      binomialX = vecObsExpected[3]
      fwrPDF <- calculate_bayes(x=binomialX,n=binomialN,p=binomialP,max_sigma=20,length_sigma=4001)

      list(CDR_PDF=cdrPDF,FWR_PDF=fwrPDF)
    }
  cat("\n")
  stopCluster(cluster)

  return(cbind(db,pdfs))
}

