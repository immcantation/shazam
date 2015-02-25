#' Calculate mutation frequencies
#'
#' \code{addMutationFrequencies} calculates the mutation frequencies.
#'
#' Methods.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  The name of the sequence column.
#' @param    germlineColumn  The name of the germline column.
#' @param    nproc     The number of cores to distribute the function over.
#' @return   A modified \code{db} data.frame with observed mutations in the
#'           OBSERVED_R_CDR (CDR replacement), OBSERVED_S_CDR (CDR silent),
#'           OBSERVED_R_FWR, (FWR replacement, and OBSERVED_S_FWR (FWR silent) columns.
#'
#' @examples
#' # TODO
#' # Working example
#'
#' @export
addMutationFrequencies <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK", nproc=1)  {

  numbOfSeqs <- nrow(db)

  availableCores <- getnproc()
  if(!(nproc<=availableCores))nproc=availableCores
  cluster <- makeCluster(nproc, type = "SOCK")
  registerDoSNOW(cluster)

  muFreq <-
      foreach(i=icount(numbOfSeqs), .packages='shm', .combine=doparProgressBar(n=numbOfSeqs)) %dopar% {
        #Count numb of mutations
        inputSeq <-  db[i,sequenceColumn]
        totalMu <- sum(countMutations(inputSeq, db[i,germlineColumn]))
        #Get mu freq
        seq<-substring(inputSeq,1,312)
        seq<-gsub("[N.-]","", seq)
        seqLen <- nchar(seq)
        return(totalMu/seqLen)
      }
  cat("\n")
  stopCluster(cluster)

  db[,"MUTATION_FREQUENCY"] <- muFreq
  return(db)
}
