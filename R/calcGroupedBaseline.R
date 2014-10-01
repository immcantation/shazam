#' Calculated CDR and FWR BASELINe selection strengths
#'
#' Given a DB file contianing the columns "CDR_PDFs" and "FWR_PDFs".\cr
#' this function groups the sequences by specified columns

#' @param   db  a data.frame of the DB file.
#' @return  a data.frame of the BASELINe sigma values and 95% Conf. Intervals for the CDR and FWR
calcGroupedBaseline <- function(db,columnsToGroupBy){
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

