# Compute mutabilities of sequence based on the tri-nucleotide model
computeMutabilities <- function(paramSeq){
  seqLen = nchar(paramSeq)
  seqMutabilites = rep(NA,seqLen)
  
  gaplessSeq = gsub("-", "", paramSeq)
  gaplessSeqLen = nchar(gaplessSeq)
  gaplessSeqMutabilites = rep(NA,gaplessSeqLen)
  
  gaplessSeq <- paste("NN",gaplessSeq,"NN",sep="")
  gaplessSeqLen <- nchar(gaplessSeq)
  pos<- 3:(gaplessSeqLen-2)
  subSeq =  substr(rep(gaplessSeq,gaplessSeqLen-4),(pos-2),(pos+2))    
  gaplessSeqMutabilites[pos-2] = sapply(subSeq,function(x){ getMutability5(x) }, simplify=T)
  seqMutabilites[which(s2c(paramSeq)!="-")]<- gaplessSeqMutabilites
  return(seqMutabilites)
  
  
}