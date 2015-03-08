#' Adds columns to DB, containing the expected frequencies of mutations
#'
#' This adds the expected frequencies number of mutations to the DB file.
#'
#' @param   param_strSeq  The name of the sequence column.
#' @param   param_vecMutabilities  The name of the germline column.
#' @return  matTargeting
computeTargeting <- function(param_strSeq, param_vecMutabilities){
    NUCLEOTIDES <- c("A", "C", "G", "T", "N")
    seqLen <- nchar(param_strSeq)    
    
    seqsubstitution = matrix(NA,ncol=seqLen,nrow=5,dimnames=list(NUCLEOTIDES,s2c(param_strSeq)))
    gaplessSeq = gsub("\\.", "", param_strSeq)
    gaplessSeq <- paste("NN",gaplessSeq,"NN",sep="")
    gaplessSeqLen <- nchar(gaplessSeq)
    #gaplessSeqsubstitution = matrix(NA,ncol=gaplessSeqLen,nrow=4)

    pos<- 3:(gaplessSeqLen-2)
    subSeq =  substr(rep(gaplessSeq,gaplessSeqLen-4),(pos-2),(pos+2))
    gaplessSeqsubstitution = sapply(subSeq,function(x){ getTransistionProb5(x) }, simplify=T)
    seqsubstitution[,which(s2c(param_strSeq)!=".")] <- gaplessSeqsubstitution

    matTargeting <- sweep(seqsubstitution,2,param_vecMutabilities,`*`)
    return (matTargeting)

}
