#' Adds columns to DB, containing the numbers of observed mutations
#'
#' This adds the observed number of mutations to the DB file.\cr
#'

#' @param   db  a data.frame of the DB file.
#' @param   sequenceColumn  The name of the sequence column.
#' @param   germlineColumn  The name of the germline column.
#' @return  db  a data.frame of the DB file
#' @export
addObservedMutations <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK")  {
  db[,c("OBSERVED_R_CDR", "OBSERVED_S_CDR", "OBSERVED_R_FWR", "OBSERVED_S_FWR")] = t(sapply( 1:nrow(db), function(x){countMutations(db[x,sequenceColumn], db[x,germlineColumn])}, simplify="array"))
  return(db)
}
