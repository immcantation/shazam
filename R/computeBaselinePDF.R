#' Computes the BASELINe probability density function
#'
#' \code{computeBaselinePDF} computes the BASELINe probability density functions for a
#' set of sequences.
#'
#' @param    db     data.frame containing observed mutations and expected mutation 
#'                  frequencies for a set of sequence data. See Details for descriptions
#'                  of the required input columns.
#' @param    nproc  number of cores to distribute the operation over.
#' 
#' @return   (A modified \code{db} data.frame) Should convert to an object.
#'
#' @details
#' The input data.frame \code{db} must contain the following columns:
#' \itemize{
#'     \item  \code{OBSERVED_R_CDR}:  number of replacement mutations in the CDRs.
#'     \item  \code{OBSERVED_S_CDR}:  number of silent mutations in the CDRs.
#'     \item  \code{OBSERVED_R_FWR}:  number of replacement mutations in the FWRs.
#'     \item  \code{OBSERVED_S_FWR}:  number of silent mutations in the FWRs.
#'     \item  \code{EXPECTED_R_CDR}:  expected frequency of replacement mutations 
#'                                    in the CDRs.
#'     \item  \code{EXPECTED_S_CDR}:  expected frequency of silent mutations 
#'                                    in the CDRs.
#'     \item  \code{EXPECTED_R_FWR}:  expected frequency of replacement mutations 
#'                                    in the FWRs.
#'     \item  \code{EXPECTED_S_FWR}:  expected frequency of silent mutations 
#'                                    in the FWRs.
#' }
#' 
#' Brief explanation of how this works!!!
#' 
#' @references
#' \enumerate{
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499–504.
#'   \item  Yaari G, et al. Quantifying selection in high-throughput Immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134.
#'  }
#'  
#' @seealso  See \code{\link{getObservedMutations}} and \code{\link{getExpectedMutationFrequencies}}
#'           for generating the necessary mutation statistics. 
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add observed mutations and expected mutation frequencies to db
#' db <- addObservedMutations(db)
#' db <- addExpectedFrequencies(db)
#' 
#' # Calculate BASELINe PDFs
#' db_pdf <- computeBaselinePDF(db)
#' 
#' @export
computeBaselinePDF <- function(db, nproc=1){

  availableCores <- getnproc()
  if(!(nproc<=availableCores))nproc=availableCores
  cluster <- makeCluster(nproc, type = "SOCK")
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

