#' Identifies clonal consensus sequences
#'
#' \code{addClonalSequence} identifies the consensus sequence of each clonal group in the input
#' data.frame.
#'
#' Details (methods) section. Paragraph 1.
#'
#' Paragraph 2.
#'
#' @param    db              data.frame containing sequence data.
#' @param    cloneColumn     column containing clonal cluster assignments.
#' @param    sequenceColumn  gapped sequence column name.
#' @param    germlineColumn  gapped germline column name.
#' @param    nproc           The number of cores to distribute the function over.
#' @return   A modified \code{db} data.frame with clonal consensus sequences in the
#'           SEQUENCE_GAP_CLONE column.
#'
#' @seealso  \code{\link{addObservedMutations}}, \code{\link{addExpectedFrequencies}}
#' @examples
#' # TODO
#' # Working example
#'
#' @export
addClonalSequence <- function(db, cloneColumn="CLONE", sequenceColumn="SEQUENCE_GAP",
                              germlineColumn="GERMLINE_GAP_D_MASK", nproc=1)  {
  #db <- ddply(db, "CLONE", transform,
  #            SEQUENCE_GAP_CLONE=(collapseCloneTry(SEQUENCE_GAP, GERMLINE_GAP_D_MASK, readEnd))[[1]][1],
  #            .progress="text")

  runAsParallel <- FALSE
  progressBar <- "text"
  if(nproc>1){
    cluster <- makeCluster(nproc, type = "SOCK")
    registerDoSNOW(cluster)
    clusterEvalQ(cluster, library(shm))
    runAsParallel <- TRUE
    progressBar <- "none"
  }


  db <- ddply(db, cloneColumn, here(mutate),
              SEQUENCE_GAP_CLONE=(collapseCloneTry(eval(parse(text=sequenceColumn)),
                                                   eval(parse(text=germlineColumn)),
                                                   readEnd))[[1]][1],
              .progress=progressBar, .parallel=runAsParallel)

  if(nproc>1)stopCluster(cluster)

  return(db)
}
