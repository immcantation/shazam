#' Determine expected mutation frequencies
#'
#' \code{addExpectedFrequencies} calculate the expected frequency of mutations for
#' each sequence.  Expectations are calculated for CD and FW regions and both
#' silent and replacement mutations.
#' 
#' Method details.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the sequence column.
#' @param    germlineColumn  name of the germline column.
#' @return   A modified \code{db} data.frame with expected frequencies in the 
#'           EXPECTED_R_CDR (CDR replacement), EXPECTED_S_CDR (CDR silent), 
#'           EXPECTED_R_FWR, (FWR replacement, and EXPECTED_S_FWR (FWR silent) columns.
#' 
#' @seealso  \code{\link{addObservedMutations}}
#' @examples
#' # TODO
#' # Working example
#' 
#' @export
addExpectedFrequencies <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK")  {

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
  return(db)

}
