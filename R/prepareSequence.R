#' Prepare SEQUENCE_GAP column
#'
#' This function prepares the SEQUENCE_GAP column of the DB object. If the DB object uses a different 
#' name for the column or you need to prepare a different column please specify the column to prepare.
#'
#' @param   db          a data.frame of the DB file
#' @param   columnName  The name of the column to prepare.
#' @return  a data.frame of the DB file
#' @export
prepareSequence <- function(db, columnName = c("SEQUENCE_GAP", "GERMLINE_GAP_D_MASK") ) {
  db[,columnName] <- sapply(db[,columnName],toupper, simplify = "array", USE.NAMES=F)
  return(db)
}
