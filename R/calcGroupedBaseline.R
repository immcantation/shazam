#' Calculated CDR and FWR BASELINe selection strengths
#'
#' \code{calcGroupedBaseline} groups the sequences by specified columns
#'
#' @param    db                 data.frame containing sequence data. \code{db} must contain the columns CDR_PDFs and FWR_PDFs.
#' @param    columnsToGroupBy   The columns you want to group and obtain combined selection scores for.
#' @param    numbOfCores        The number of cores to distribute the function over.
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
calcGroupedBaseline <- function(db, columnsToGroupBy, numbOfCores=1){

  availableCores <- getNumbOfCores()
  runAsParallel <- FALSE
  if(!(numbOfCores<=availableCores))numbOfCores=availableCores
  if(numbOfCores>1){
    cluster <- makeCluster(numbOfCores, type = "SOCK")
    registerDoSNOW(cluster)
    clusterEvalQ(cluster, library(shm))
    runAsParallel <- TRUE
  }
  if( ('CDR_PDF' %in% colnames(db)) & ('FWR_PDF' %in% colnames(db))){
    db <- ddply(db, columnsToGroupBy, .fun=function_GroupBaseline, .parallel=runAsParallel)
  }else if( 'BASELINE_Sigma' %in% colnames(db) ){
    db <- ddply(db, columnsToGroupBy, .fun=function_ReGroupBaseline, .parallel=runAsParallel)
  }else{
    return(NULL)
  }

  stopCluster(cluster)
  return(db)
}


#' Calucaltes BASELINe selection strengths for previously grouped sequences
#'
#' \code{function_ReGroupBaseline} is a helper function that further groups the sequences that have
#'  been grouped using the \code{calcGroupedBaseline} (specifically, \code{function_GroupBaseline}).
#'
function_ReGroupBaseline <- function(i){
  cdr_groupN <- sum(!is.na(i[,"CDR_PDF"]))
  fwr_groupN <- sum(!is.na(i[,"FWR_PDF"]))
  cdr_groupPDF <- groupPosteriors(i[,"CDR_PDF"])
  fwr_groupPDF <- groupPosteriors(i[,"FWR_PDF"])

  baselineInfo_CDR <- calcBayesOutputInfo(cdr_groupPDF)
  baselineInfo_FWR <- calcBayesOutputInfo(fwr_groupPDF)

  # make a dataframe
  dCDR <- data.frame(
    Region="CDR",
    Sigma = baselineInfo_CDR[1],
    CI95_Lower = baselineInfo_CDR[2],
    CI95_Upper = baselineInfo_CDR[3],
    N=cdr_groupN,
    PDF=array(list(cdr_groupPDF))
  )
  dFWR <- data.frame(
    Region="FWR",
    Sigma = baselineInfo_FWR[1],
    CI95_Lower = baselineInfo_FWR[2],
    CI95_Upper = baselineInfo_FWR[3],
    N=fwr_groupN,
    PDF=array(list(fwr_groupPDF))
  )

  d <- rbind(dCDR, dFWR)
  return(d)
}



#' Calucaltes the BASELINe selection strengths for independent sequences grouped together
#'
#' \code{function_GroupBaseline} is a helper function that groups the sequences by specified columns,
#' and returns the Baseline Selection strengt, 95% confidence intervals and the grouped (convoluted)
#' probability density functions.
#'
function_GroupBaseline <- function(i){

  cdr_groupN <- sum(!is.na(i[,"CDR_PDF"]))
  fwr_groupN <- sum(!is.na(i[,"FWR_PDF"]))
  cdr_groupPDF <- groupPosteriors(i[,"CDR_PDF"])
  fwr_groupPDF <- groupPosteriors(i[,"FWR_PDF"])

  baselineInfo_CDR <- calcBayesOutputInfo(cdr_groupPDF)
  baselineInfo_FWR <- calcBayesOutputInfo(fwr_groupPDF)

  # make a dataframe
  dCDR <- data.frame(
    Region="CDR",
    Sigma = baselineInfo_CDR[1],
    CI95_Lower = baselineInfo_CDR[2],
    CI95_Upper = baselineInfo_CDR[3],
    N=cdr_groupN,
    PDF=array(list(cdr_groupPDF))
  )
  dFWR <- data.frame(
    Region="FWR",
    Sigma = baselineInfo_FWR[1],
    CI95_Lower = baselineInfo_FWR[2],
    CI95_Upper = baselineInfo_FWR[3],
    N=fwr_groupN,
    PDF=array(list(fwr_groupPDF))
  )

  d <- rbind(dCDR, dFWR)
  return(d)
}
# function_GroupBaseline <- function(i){
#   cdr_groupN <- sum(!is.na(i[,"CDR_PDF"]))
#   fwr_groupN <- sum(!is.na(i[,"FWR_PDF"]))
#   cdr_groupPDF <- groupPosteriors(i[,"CDR_PDF"])
#   fwr_groupPDF <- groupPosteriors(i[,"FWR_PDF"])
#
#   baselineInfo_CDR <- calcBayesOutputInfo(cdr_groupPDF)
#   baselineInfo_FWR <- calcBayesOutputInfo(fwr_groupPDF)
#
#   # make a dataframe
#   d <- data.frame(
#     CDR_Sigma = baselineInfo_CDR[1],
#     CDR_CI95_Lower = baselineInfo_CDR[2],
#     CDR_CI95_Upper = baselineInfo_CDR[3],
#     FWR_Sigma = baselineInfo_FWR[1],
#     FWR_CI95_Lower = baselineInfo_FWR[2],
#     FWR_CI95_Upper = baselineInfo_FWR[3],
#     CDR_N = cdr_groupN,
#     FWR_N = fwr_groupN
#   )
#
#   d <- cbind(d, CDR_PDF=array(list(cdr_groupPDF)), FWR_PDF=array(list(fwr_groupPDF) ))
#   return(d)
# }
