#' Adds columns to DB, containing the expected frequencies of mutations
#'
#' This adds the expected frequencies number of mutations to the DB file.\cr
#'

#' @param   param_strSeq  The name of the sequence column.
#' @param   param_vecMutabilities  The name of the germline column.
#' @return  matTargeting
computeMutationTypes <- function(param_strSeq){
  #cat(param_strSeq,"\n")
  lenSeq <- nchar(param_strSeq)
  trimmedSeq = trimToLastCodon(param_strSeq)
  lenTrimmedSeq = nchar(trimmedSeq)
  vecCodons = sapply({1:(lenSeq/3)}*3-2,function(x){substr(trimmedSeq,x,x+2)})
  matMutationTypes = matrix( NA ,ncol=lenSeq,nrow=4, byrow=F)
  matMutationTypes[,1:lenTrimmedSeq] = matrix( unlist(CODON_TABLE[,vecCodons]) ,ncol=lenTrimmedSeq,nrow=4, byrow=F)
  #matMutationTypes = CODON_TABLE[,vecCodons]
  #dimnames( matMutationTypes ) =  list(NUCLEOTIDES,1:(ncol(matMutationTypes)))
  matMutationTypes <- rbind(  matMutationTypes, matrix(NA,ncol=ncol(matMutationTypes),nrow=1))
  dimnames( matMutationTypes ) =  list(NUCLEOTIDES,s2c(param_strSeq))
  selfbases <- as.numeric(factor(colnames(matMutationTypes),levels=NUCLEOTIDES_FAC))
  tmp <- sapply(1:ncol(matMutationTypes),function(x){matMutationTypes[selfbases[x],x]<<-NA})
  if(nrow(matMutationTypes)==5)rownames(matMutationTypes)[5] <- "N"
  return(matMutationTypes)
}
