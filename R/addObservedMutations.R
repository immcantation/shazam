#' Calculate observed mutations
#'
#' \code{addObservedMutations} calculates the observed number of mutations with the 
#' CD and FW regions of IMGT-gapped sequences. Is this count or frequency?
#' 
#' Methods.  
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  The name of the sequence column.
#' @param    germlineColumn  The name of the germline column.
#' @return   A modified \code{db} data.frame with observed mutations in the 
#'           OBSERVED_R_CDR (CDR replacement), OBSERVED_S_CDR (CDR silent), 
#'           OBSERVED_R_FWR, (FWR replacement, and OBSERVED_S_FWR (FWR silent) columns.
#' 
#' @seealso  \code{\link{addClonalSequence}}, \code{\link{addExpectedFrequencies}}.
#' @examples
#' # TODO
#' # Working example
#' 
#' @export
addObservedMutations <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK")  {
  pb <- txtProgressBar(min=1,max=nrow(db),width=20)
  cat("Progress: 0%      50%     100%\n")
  cat("          ")
  db[, c("OBSERVED_R_CDR", 
         "OBSERVED_S_CDR", 
         "OBSERVED_R_FWR", 
         "OBSERVED_S_FWR")] = t(sapply(1:nrow(db), function(x) { setTxtProgressBar(pb, x);
                                                                 countMutations(db[x,sequenceColumn], db[x,germlineColumn]) },
                                       simplify="array"))
  cat("\n")
  close(pb)
  return(db)
}
