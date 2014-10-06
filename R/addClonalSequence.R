#' Adds columns to DB, containing the concensus clonal sequence
#'
#' This identifies the consensus clonal sequence to the DB file.\cr
#'

#' @param   db  a data.frame of the DB file.
#' @param   sequenceColumn  The name of the sequence column.
#' @param   germlineColumn  The name of the germline column.
#' @return  db  a data.frame of the DB file
#' @export
addClonalSequence <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK")  {
  db <- ddply(db, "CLONE", transform, SEQUENCE_GAP_CLONE=(collapseCloneTry(SEQUENCE_GAP,GERMLINE_GAP_D_MASK,readEnd))[[1]][1],.progress="text")
  return(db)
}
