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
