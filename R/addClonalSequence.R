#' Identifies clonal consensus sequences
#'
#' \code{addClonalSequence} identifies the consensus sequence of each clonal group and
#' appends a column to the input data.frame containing the clonal consensus for each 
#' sequence.
#'
#' @param    db              data.frame containing sequence data.
#' @param    cloneColumn     name of the column containing clonal cluster identifiers.
#' @param    sequenceColumn  name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn  name of the column containing IMGT-gapped germline sequences.
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with clonal consensus sequences in the
#'           SEQUENCE_GAP_CLONE column.
#'
#' @details
#' How does this work?
#' 
#' @references
#' Which ones?
#' 
#' @seealso  What uses this?
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add SEQUENCE_GAP_CLONE column to db
#' db_new <- addClonalSequence(db)
#' head(db_new[c(1, 15)])
#'
#' @export
addClonalSequence <- function(db, cloneColumn="CLONE", sequenceColumn="SEQUENCE_GAP",
                              germlineColumn="GERMLINE_GAP_D_MASK", nproc=1)  {
  # Start SNOW cluster and set progress bar
  if (nproc > 1) {
    runAsParallel <- TRUE
    progressBar <- "none"  
    cluster <- makeCluster(nproc, type="SOCK")
    registerDoSNOW(cluster)
    clusterEvalQ(cluster, library(shm))
  } else {
      runAsParallel <- FALSE
      progressBar <- "text"
  }

  # Generate clonal consensus sequences
  if (progressBar != "none") { cat("-> BUILDING CLONAL SEQUENCES\n") }
  db <- ddply(db, cloneColumn, here(mutate),
              SEQUENCE_GAP_CLONE=(collapseCloneTry(eval(parse(text=sequenceColumn)),
                                                   eval(parse(text=germlineColumn)),
                                                   VLENGTH))[[1]][1],
              .progress=progressBar, .parallel=runAsParallel)
    
  # Stop SNOW cluster
  if(nproc > 1) { stopCluster(cluster) }

  return(db)
}
