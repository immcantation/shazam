#' Builds a mutability model
#'
#' This builds a 5-mer based mutability model.\cr
#'

#' @param   db  a data.frame of the DB file.
#' @param   model  The type of model to create "S" = silent, "RS" = Silent = Replacement mutations
#' @param   substitutionModel  matrix of fivemers and the substitution counts (from createSubstitutionModel)
#' @param   sequenceColumn  The name of the sequence column
#' @param   germlineColumn  The name of the germline column
#' @param   VCallColumn The name of the genotyped columnn
#' @param   multipleMutation 0 = Treat independently, 1 = Ignore codons with multiple mutations
#' @return  list of fivemers and the substitution counts
#' @export
createMutabilityModel <- function(db,
                                    model="RS",
                                    substitutionModel=substitutionModel,
                                    sequenceColumn="SEQUENCE_GAP",
                                    germlineColumn="GERMLINE_GAP_D_MASK",
                                    vCallColumn="V_CALL_GENOTYPED",
                                    multipleMutation=0)  {
  NUCLEOTIDES <- c("A", "C", "G", "T")

  mutations <- listObservedMutations(db)

  substitutionModel <- apply(substitutionModel,2,function(x)x/sum(x))
  template <- rep(0, 1024)
  names(template) <- words(5,NUCLEOTIDES)

  COUNT<-list()
  for(index in 1:length(mutations)){
    COUNT[[index]]<-template
    indexMutation <- mutations[[index]]
    if(!sum(is.na(indexMutation))){
      cSeq <-  s2c(db[index,sequenceColumn])
      cGL  <-  s2c(db[index,germlineColumn])
      positions <- as.numeric(names(indexMutation))
      positions <- positions[!is.na(positions)]
      for( position in  positions){
        wrd5<-substr(db[index,germlineColumn],position-2,position+2)
        codonNucs = getCodonPos(position)
        codonGL = cGL[codonNucs]
        codonSeq = cSeq[codonNucs]
        muCodonPos = {position-1}%%3+1
        seqAtMutation <- codonSeq[muCodonPos]
        glAtMutation <- codonGL[muCodonPos]
        if( !any(codonGL%in%c("N","-", ".")) & !any(codonSeq%in%c("N","-", ".")) ){
          if(multipleMutation==0){ # if independent mutations are included
            codonSeq <- codonGL
            codonSeq[muCodonPos] <- seqAtMutation
          }
          if(!length(grep("N",wrd5))){
            COUNT[[index]][wrd5]<- COUNT[[index]][wrd5]+1;
          }
        }
      }
    }
  }

  BG_COUNT<-list()
  for(index in 1:length(mutations)){
    cat(index,"\n")
    BG_COUNT[[index]]<-template
    sSeq <- gsub("\\.","",db[index,sequenceColumn])
    sGL <- gsub("\\.","",db[index,germlineColumn])
    cSeq <-  s2c(sSeq)
    cGL  <-  s2c(sGL)
    positions <- 3:(length(cGL)-2)
    for( position in  positions){
      wrd5<-substr(sGL,position-2,position+2)
      codonNucs = getCodonPos(position)
      codonGL = cGL[codonNucs]
      codonSeq = cSeq[codonNucs]
      muCodonPos = {position-1}%%3+1
      seqAtMutation <- codonSeq[muCodonPos]
      glAtMutation <- codonGL[muCodonPos]
      if( !any(codonGL%in%c("N","-")) & !any(codonSeq%in%c("N","-")) ){
        if(multipleMutation==0){ # if independent mutations are included
          codonSeq <- codonGL
          codonSeq[muCodonPos] <- seqAtMutation
        }
        codonPermutate <- matrix(rep(codonGL,3),ncol=3,byrow=T)
        codonPermutate[,muCodonPos] <- canMutateTo(glAtMutation)[-4]
        codonPermutate <- apply(codonPermutate,1,paste,collapse="")
        codonPermutate <- matrix( c( codonPermutate, rep(c2s(codonGL),3) ), ncol=2, byrow=F)
        muType <- mutationTypeOptimized(codonPermutate)
        if(!length(grep("N",wrd5)) & !length(grep("-",wrd5))){
          for(m in 1:3){
            if(muType[m]=="S"){
              BG_COUNT[[index]][wrd5]<- BG_COUNT[[index]][wrd5] + substitutionModel[substr(codonPermutate[m,1],muCodonPos,muCodonPos),wrd5];
            }
          }
        }
      }
    }
    BG_COUNT[[index]][BG_COUNT[[index]]==0]<-NA
  }

  Mutability<-list()
  for(i in 1:length(mutations)){
    Mutability[[i]]<-list()
    Mutability[[i]][[1]]<-COUNT[[i]]/BG_COUNT[[i]]
    Mutability[[i]][[1]]<-Mutability[[i]][[1]]/sum(Mutability[[i]][[1]],na.rm=TRUE)
    Mutability[[i]][[2]]<-length( mutations[[i]])
    Mutability[[i]][[3]]<-sum(COUNT[[i]],na.rm=TRUE)

  }
  return(Mutability)
  NUCLEOTIDES <- c("A", "C", "G", "T", "N")
}


