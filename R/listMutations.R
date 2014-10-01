listMutations <- function(seqInput, seqGL) {
  if( is.na(c(seqInput, seqGL)) ) return(array(NA,4))
  seqI = s2c(seqInput)
  seqG = s2c(seqGL)
  matIGL = matrix(c(seqI,seqG),ncol=length(seqI),nrow=2,byrow=T)
  mutations <- analyzeMutations2NucUri(matIGL)
  return(mutations)
}
