#' Calculate observed mutations
#'
#' \code{addObservedMutations} calculates the observed number of V-segment mutations 
#' within the framework (FW) and complementarity determining (CD) regions of IMGT-gapped 
#' nucleotide sequences. Mutation counts are appended to the input data.frame as 
#' additional columns.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn  name of the column containing IMGT-gapped germline sequences.
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with observed mutation counts for each 
#'           sequence listed in the following columns: 
#'           \itemize{
#'             \item  \code{OBSERVED_R_CDR}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{OBSERVED_S_CDR}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{OBSERVED_R_FWR}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{OBSERVED_S_FWR}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           
#' @details
#' How does this work?
#' 
#' @references
#' Which ones?
#' 
#' @seealso  See \code{\link{addExpectedFrequencies}} for calculating expected mutation
#'           frequencies.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add observed mutations to db
#' db_new <- addObservedMutations(db)
#' head(db_new[c(1, 15:18)])
#'
#' @export
addObservedMutations <- function(db, sequenceColumn="SEQUENCE_GAP", 
                                 germlineColumn="GERMLINE_GAP_D_MASK", nproc=1) {
  numbOfSeqs <- nrow(db)
  if(nproc == 1) {
    pb <- txtProgressBar(min=1, max=numbOfSeqs, width=20)
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
  } else {
    availableCores <- getnproc()
    if(!(nproc<=availableCores))nproc=availableCores
    cluster <- makeCluster(nproc, type = "SOCK")
    registerDoSNOW(cluster)

    obsMutations <-
        foreach(i=icount(numbOfSeqs), .packages='shm', .combine=doparProgressBar(n=numbOfSeqs), .multicombine=TRUE) %dopar% {
          countMutations(db[i,sequenceColumn], db[i,germlineColumn])
        }

    stopCluster(cluster)

    # Add observed mutations to db
    db[, c("OBSERVED_R_CDR",
           "OBSERVED_S_CDR",
           "OBSERVED_R_FWR",
           "OBSERVED_S_FWR")] <- obsMutations
    cat("\n")
  }

  return(db)
}
