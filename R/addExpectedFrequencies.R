#' Adds columns to DB, containing the expected frequencies of mutations
#'
#' This adds the expected frequencies number of mutations to the DB file.\cr
#'

#' @param   db  a data.frame of the DB file.
#' @param   sequenceColumn  The name of the sequence column.
#' @param   germlineColumn  The name of the germline column.
#' @return  db  a data.frame of the DB file
addExpectedFrequencies <- function(db, sequenceColumn="SEQUENCE_GAP", germlineColumn="GERMLINE_GAP_D_MASK")  {

  facGL <- factor(db[,germlineColumn])
  facLevels = levels(facGL)
  LisGLs_MutabilityU = lapply(1:length(facLevels),  function(x){
    computeMutabilities(facLevels[x])
  })
  facIndex = match(facGL,facLevels)

  LisGLs_Mutability = lapply(1:nrow(db),  function(x){
    cInput = rep(NA,nchar(db[x,sequenceColumn]))
    cInput[s2c(db[x,sequenceColumn])!="N"] = 1
    LisGLs_MutabilityU[[facIndex[x]]] * cInput
  })

  LisGLs_Targeting =  lapply(1:nrow(db),  function(x){
    computeTargeting(db[x,germlineColumn],LisGLs_Mutability[[x]])
  })

  LisGLs_MutationTypes  = lapply(1:nrow(db),function(x){
    computeMutationTypes(db[x,germlineColumn])
  })

  LisGLs_Exp = lapply(1:nrow(db),  function(x){
    computeExpected(LisGLs_Targeting[[x]],LisGLs_MutationTypes[[x]])
  })

  ul_LisGLs_Exp =  unlist(LisGLs_Exp)
  matExp <- matrix(ul_LisGLs_Exp,ncol=4,nrow=(length(ul_LisGLs_Exp)/4),byrow=T)
  matExp <- matExp/apply(matExp,1,sum,na.rm=T)
  db[,c("EXPECTED_R_CDR", "EXPECTED_S_CDR", "EXPECTED_R_FWR", "EXPECTED_S_FWR")] <- matExp
  return(db)

}
