listMutations <- function(seqInput, seqGL) {
  #if( is.na(c(seqInput, seqGL)) ) return(array(NA,4))
  if( is.na(c(seqInput, seqGL)) ) return(NA)
  seqI = s2c(seqInput)
  seqG = s2c(seqGL)
  matIGL = matrix(c(seqI,seqG),ncol=length(seqI),nrow=2,byrow=T)
  mutations <- analyzeMutations2NucUri(matIGL)
  mutations <- mutations[ !is.na(mutations) ]
  positions <- as.numeric(names(mutations))
  mutations <- mutations[positions<=VLENGTH]
  if(length(mutations)>0){
    return(mutations)
  }else{
    return(NA)
  }
}
