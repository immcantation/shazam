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
#' @param    numbOfCores     The number of cores to distribute the function over.
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
addObservedMutations <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK", numbOfCores=1)  {
  numbOfSeqs <- nrow(db)
  if(numbOfCores==1){
    pb <- txtProgressBar(min=1,max=numbOfSeqs,width=20)
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
  }else{
    availableCores <- getNumbOfCores()
    if(!(numbOfCores<=availableCores))numbOfCores=availableCores
    cluster <- makeCluster(numbOfCores, type = "SOCK")
    registerDoSNOW(cluster)

    obsMutations <-
        foreach(i=icount(numbOfSeqs), .packages='shm', .combine=doparProgressBar(n=numbOfSeqs), .multicombine=TRUE) %dopar% {
          countMutations(db[i,sequenceColumn], db[i,germlineColumn])
        }

    stopCluster(cluster)


    db[, c("OBSERVED_R_CDR",
           "OBSERVED_S_CDR",
           "OBSERVED_R_FWR",
           "OBSERVED_S_FWR")] <- obsMutations
    cat("\n")

  }

  return(db)
}
