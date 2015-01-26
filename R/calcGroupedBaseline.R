#' Calculated CDR and FWR BASELINe selection strengths
#'
#' \code{calcGroupedBaseline} groups the sequences by specified columns
#'
#' @param    db  data.frame containing sequence data. \code{db} must contain the columns CDR_PDFs and FWR_PDFs.
#' @return   A data.frame of the BASELINe sigma values and 95% Conf. Intervals for the CDR and FWR.
#'
#' @references    
#' Uduman M, Yaari G, Hershberg U, Stern JA, Shlomchik MJ, Kleinstein SH. 
#'   Detecting selection in immunoglobulin sequences. 
#'   Nucleic Acids Res. 2011 Jul;39(Web Server issue):W499-504.
#' 
#' Yaari G, Uduman M, Kleinstein SH. 
#'   Quantifying selection in high-throughput Immunoglobulin sequencing data sets. 
#'   Nucleic Acids Res. 2012 Sep 1;40(17):e134.
#' @seealso  \code{\link{addObservedMutations}}, \code{\link{addExpectedFrequencies}}, \code{\link{computeBaselinePDF}}.
#' @examples
#' # TODO
#' # Working example
#' 
#' @export
calcGroupedBaseline <- function(db, columnsToGroupBy){
  #df[,columnsToGroupBy] <- lapply(df[,columnsToGroupBy] , factor)
  db <- ddply(db, columnsToGroupBy, .fun=fun_calcGroupedBaseline)
  return(db)
}

# custom summary function
fun_calcGroupedBaseline <- function(i){
  columnName <- paste("CDR_PDFs",sep="")
  baselineInfo_CDR <- calcBayesOutputInfo(groupPosteriors(i[,"CDR_PDFs"]))
  baselineInfo_FWR <- calcBayesOutputInfo(groupPosteriors(i[,"FWR_PDFs"]))
  # make a dataframe
  d <- data.frame(
    CDR_Sigma = baselineInfo_CDR[1],
    CDR_CI95_Lower = baselineInfo_CDR[2],
    CDR_CI95_Upper = baselineInfo_CDR[3],
    FWR_Sigma = baselineInfo_FWR[1],
    FWR_CI95_Lower = baselineInfo_FWR[2],
    FWR_CI95_Upper = baselineInfo_FWR[3]
  )
  return(d)
}

