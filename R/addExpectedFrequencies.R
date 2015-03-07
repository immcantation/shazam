#' Determine expected mutation frequencies
#'
#' \code{addExpectedFrequencies} calculate the expected frequency of mutations for
#' each sequence. Expectations are calculated the framework (FW) and complementarity 
#' determining (CD) regions of IMGT-gapped nucleotide sequences. Expected mutation 
#' frequencies are appended to the input data.frame as additional columns.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn  name of the column containing IMGT-gapped germline sequences.
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with observed mutation counts for each 
#'           sequence listed in the following columns: 
#'           \itemize{
#'             \item  \code{EXPECTED_R_CDR}:  expected frequency of replacement mutations 
#'                                            in CDR1 and CDR2 of the V-segment.
#'             \item  \code{EXPECTED_S_CDR}:  expected frequency of silent mutations 
#'                                            in CDR1 and CDR2 of the V-segment.
#'             \item  \code{EXPECTED_R_FWR}:  expected frequency of replacement mutations 
#'                                            in FWR1, FWR2 and FWR3 of the V-segment.
#'             \item  \code{EXPECTED_S_FWR}:  expected frequency of silent mutations 
#'                                            in FWR1, FWR2 and FWR3 of the V-segment.
#'           }
#'
#' @details
#' How does this work?
#' 
#' @references
#' Which ones?
#'
#' @seealso  See \code{\link{addObservedMutations}} for counting observed mutations.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add observed mutations to db
#' db_new <- addExpectedFrequencies(db)
#' head(db_new[c(1, 15:18)])
#'
#' @export
addExpectedFrequencies <- function(db, sequenceColumn="SEQUENCE_GAP", 
                                   germlineColumn="GERMLINE_GAP_D_MASK", nproc=1) {
  if (nproc == 1) {
      facGL <- factor(db[,germlineColumn])
      facLevels = levels(facGL)
      cat("Computing mutabilities...\n")
      pb <- txtProgressBar(min=1,max=length(facLevels),width=20)
      cat("Progress: 0%      50%     100%\n")
      cat("          ")
      LisGLs_MutabilityU = lapply(1:length(facLevels),  function(x){
        setTxtProgressBar(pb, x);
        computeMutabilities(facLevels[x])
      })
      facIndex = match(facGL,facLevels)
      cat("\n")
      close(pb)
    
      LisGLs_Mutability = lapply(1:nrow(db),  function(x){
        cInput = rep(NA,nchar(db[x,sequenceColumn]))
        cInput[s2c(db[x,sequenceColumn])!="N"] = 1
        LisGLs_MutabilityU[[facIndex[x]]] * cInput
      })
    
      cat("Computing targeting...\n")
      pb <- txtProgressBar(min=1,max=nrow(db),width=20)
      cat("Progress: 0%      50%     100%\n")
      cat("          ")
      LisGLs_Targeting =  lapply(1:nrow(db),  function(x){
        setTxtProgressBar(pb, x);
        computeTargeting(db[x,germlineColumn],LisGLs_Mutability[[x]])
      })
      cat("\n")
      close(pb)
    
      cat("Computing mutation types...\n")
      pb <- txtProgressBar(min=1,max=nrow(db),width=20)
      cat("Progress: 0%      50%     100%\n")
      cat("          ")
      LisGLs_MutationTypes  = lapply(1:nrow(db),function(x){
        setTxtProgressBar(pb, x);
        computeMutationTypes(db[x,germlineColumn])
      })
      cat("\n")
      close(pb)
    
      cat("Computing expected frequencies of mutations...\n")
      pb <- txtProgressBar(min=1,max=nrow(db),width=20)
      cat("Progress: 0%      50%     100%\n")
      cat("          ")
      LisGLs_Exp = lapply(1:nrow(db),  function(x){
        setTxtProgressBar(pb, x);
        computeExpected(LisGLs_Targeting[[x]],LisGLs_MutationTypes[[x]])
      })
      cat("\n")
      close(pb)
    
      ul_LisGLs_Exp =  unlist(LisGLs_Exp)
      matExp <- matrix(ul_LisGLs_Exp,ncol=4,nrow=(length(ul_LisGLs_Exp)/4),byrow=T)
      matExp <- matExp/apply(matExp,1,sum,na.rm=T)
      db[,c("EXPECTED_R_CDR", "EXPECTED_S_CDR", "EXPECTED_R_FWR", "EXPECTED_S_FWR")] <- matExp
  } else {
    availableCores <- getnproc()
    if(!(nproc<=availableCores))nproc=availableCores
    cluster <- makeCluster(nproc, type = "SOCK")
    registerDoSNOW(cluster)

    nInterations <- nrow(db)
    matExp <-
    foreach(i=icount(nInterations), .packages='shm', .combine=doparProgressBar(n=nInterations)) %dopar% {
      # Calculate mutability
      seqMutabilities <- computeMutabilities(db[i,germlineColumn])
      # Only include non N positions
      cInput = rep(NA,nchar(db[i,sequenceColumn]))
      cInput[s2c(db[i,sequenceColumn])!="N"] = 1
      seqMutabilities <- seqMutabilities * cInput

      # Calculate targeting
      seqTargeting <- computeTargeting(db[i,germlineColumn],seqMutabilities)

      # Calculate mutation types
      seqMutationType <- computeMutationTypes(db[i,germlineColumn])

      #Expected Freq
      seqExpectedFreq <- computeExpected(seqTargeting, seqMutationType)
      seqExpectedFreq <- seqExpectedFreq/sum(seqExpectedFreq, na.rm=T)
    }
    cat("\n")
    stopCluster(cluster)

    db[,c("EXPECTED_R_CDR", "EXPECTED_S_CDR", "EXPECTED_R_FWR", "EXPECTED_S_FWR")] <- matExp
  }
  
  return(db)
}
