#' List the numbers of observed mutations
#'
#' This lists the observed number of mutation.\cr
#'

#' @param   db  a data.frame of the DB file.
#' @param   sequenceColumn  The name of the sequence column.
#' @param   germlineColumn  The name of the germline column.
#' @return  list of mutations in each clone
#' @export
listObservedMutations <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK")  {
  listObsMutations <- lapply( 1:nrow(db), function(x){listMutations(db[x,sequenceColumn], db[x,germlineColumn])} )
  return(listObsMutations)
}
