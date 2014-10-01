#given a mat of targeting & it's corresponding mutationtypes returns
#a vector of Exp_RCDR,Exp_SCDR,Exp_RFWR,Exp_RFWR
computeExpected <- function(paramTargeting,paramMutationTypes){
  # Replacements
  RPos = which(paramMutationTypes=="R")
  #FWR
  Exp_R_FWR = sum(paramTargeting[ RPos[which(FWR_Nuc_Mat[RPos]==T)] ],na.rm=T)
  #CDR
  Exp_R_CDR = sum(paramTargeting[ RPos[which(CDR_Nuc_Mat[RPos]==T)] ],na.rm=T)
  # Silents
  SPos = which(paramMutationTypes=="S")
  #FWR
  Exp_S_FWR = sum(paramTargeting[ SPos[which(FWR_Nuc_Mat[SPos]==T)] ],na.rm=T)
  #CDR
  Exp_S_CDR = sum(paramTargeting[ SPos[which(CDR_Nuc_Mat[SPos]==T)] ],na.rm=T)

  return(c(Exp_R_CDR,Exp_S_CDR,Exp_R_FWR,Exp_S_FWR))
}
