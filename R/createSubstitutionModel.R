#' Builds a substitution model
#'
#' This builds a 5-mer based substitution model.\cr
#'

#' @param   db  a data.frame of the DB file.
#' @param   model  The type of model to create "S" = silent, "RS" = Silent = Replacement mutations
#' @param   sequenceColumn  The name of the sequence column
#' @param   germlineColumn  The name of the germline column
#' @param   VCallColumn The name of the genotyped columnn
#' @param   multipleMutation 0 = Treat independently, 1 = Ignore codons with multiple mutations
#' @return  matrix of fivemers and the substitution counts
#' @export
createSubstitutionModel <- function(db,
                                    model="RS",
                                    sequenceColumn="SEQUENCE_GAP",
                                    germlineColumn="GERMLINE_GAP_D_MASK",
                                    vCallColumn="V_CALL_GENOTYPED",
                                    multipleMutation=0)  {

  NUCLEOTIDES <- c("A", "C", "G", "T")
  substitutionMatrix <- matrix(0,ncol=4,nrow=4,dimnames=list(NUCLEOTIDES,NUCLEOTIDES))
  substitutionList = list()
  WORDS <- words(4,NUCLEOTIDES)
  VH_FAMILIES <- paste(rep("VH",7),1:7,sep="")

  for(VH_FAMILY in VH_FAMILIES){
    substitutionList[[VH_FAMILY]] = list()
    for(WORD in WORDS){
      substitutionList[[VH_FAMILY]][[WORD]] = substitutionMatrix
    }
  }


  VH_FAMILES <- substr( db[,vCallColumn],
                        gregexpr("IGH",db[,vCallColumn])[[1]][1]+4,
                        gregexpr("IGH",db[,vCallColumn])[[1]][1]+4 )
  VH_FAMILES <- paste("VH",VH_FAMILES,sep="")

  mutations <- listObservedMutations(db)

  #Silent model
  if(model=="S"){

    for(index in 1:length(mutations)){
      cSeq <-  s2c(db[index,sequenceColumn])
      cGL  <-  s2c(db[index,germlineColumn])
      indexMutation <- mutations[[index]]
      VH_FAMILY <- VH_FAMILES[index]

      positions <- as.numeric(names(indexMutation))
      positions <- positions[positions<=readEnd]
      positions <- positions[!is.na(positions)]
      for( position in  positions){
        wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
        codonNucs = getCodonPos(position)
        codonGL = cGL[codonNucs]
        codonSeq = cSeq[codonNucs]
        muCodonPos = {position-1}%%3+1
        seqAtMutation <- codonSeq[muCodonPos]
        glAtMutation <- codonGL[muCodonPos]
        if( !any(codonGL=="N") & !any(codonSeq=="N") ){
          if(multipleMutation==0){ # if independent mutations are included
            codonSeq <- codonGL
            codonSeq[muCodonPos] <- seqAtMutation
          }
          codonPermutate <- matrix(rep(codonGL,3),ncol=3,byrow=T)
          codonPermutate[,muCodonPos] <- canMutateTo(glAtMutation)[-4]
          codonPermutate <- apply(codonPermutate,1,paste,collapse="")
          codonPermutate <- matrix( c( codonPermutate, rep(c2s(codonGL),3) ), ncol=2, byrow=F)
          muType <- mutationTypeOptimized(codonPermutate)
          if(!length(grep("N",wrd))){
            if( sum(muType=="S") == length(muType) ){
              substitutionList[[VH_FAMILY]][[wrd]][glAtMutation,seqAtMutation] <- (substitutionList[[VH_FAMILY]][[wrd]][glAtMutation,seqAtMutation] + 1)
            }
          }
        }
      }
    }


  }else{ #RS model (All mutations)

    for(index in 1:length(mutations)){
      #cat(index,"\n")
      cSeq <-  s2c(db[index,sequenceColumn])
      cGL  <-  s2c(db[index,germlineColumn])
      indexMutation <- mutations[[index]]
      VH_FAMILY <- VH_FAMILES[index]

      positions <- as.numeric(names(indexMutation))
      positions <- positions[positions<=readEnd]
      positions <- positions[!is.na(positions)]
      for( position in  positions){
        wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
        codonNucs = getCodonPos(position)
        codonGL = cGL[codonNucs]
        codonSeq = cSeq[codonNucs]
        muCodonPos = {position-1}%%3+1
        seqAtMutation <- codonSeq[muCodonPos]
        glAtMutation <- codonGL[muCodonPos]
        if( !any(codonGL=="N") & !any(codonSeq=="N") ){
          if(multipleMutation==0){ # if independent mutations are included
            codonSeq <- codonGL
            codonSeq[muCodonPos] <- seqAtMutation
          }
          if(!length(grep("N",wrd))){
            substitutionList[[VH_FAMILY]][[wrd]][glAtMutation,seqAtMutation] <- (substitutionList[[VH_FAMILY]][[wrd]][glAtMutation,seqAtMutation] + 1)
          }
        }
      }
    }

  }


  SAMPLE = "A"
  WORDS = words(4,NUCLEOTIDES)
  VH_FAMILIES <- paste(rep("VH",7),1:7,sep="")

  X1<-rep(WORDS,7)
  X2<-rep(VH_FAMILIES,each=256)
  arrNames <- paste(X2,X1,sep="_")
  rm(X2,X1)


  listSubstitution <- array(0,dim=c(256*7,4,4),dimnames=list(arrNames,NUCLEOTIDES,NUCLEOTIDES))


  #substitutionList
  for(VH_FAMILY in VH_FAMILIES){
    listSubstitution[paste(VH_FAMILY,WORDS,sep="_"),,]<-t(sapply(WORDS,function(WORD){substitutionList[[VH_FAMILY]][[WORD]]}))
  }


  M<-list()
  NAMES<-sapply(dimnames(listSubstitution)[[1]],function(x)strsplit(x,"_",fixed=TRUE)[[1]])
  CTS<-list()
  X1<-list()

  CT_FiveMers <-
  c( paste(substring(words(3,NUCLEOTIDES),1,1), "C",substring(words(3,NUCLEOTIDES),2,3),sep=""),
     paste(substring(words(3,NUCLEOTIDES),1,1),"G",substring(words(3,NUCLEOTIDES),2,3),sep=""),
     paste(substring(words(3,NUCLEOTIDES),1,1),"T",substring(words(3,NUCLEOTIDES),2,3),sep="")
     )

    for(CT in CT_FiveMers){
      CTS[[CT]]<-WORDS[sapply(WORDS,function(x)substring(x,1,4))== CT]
      INDEx<-NAMES[2,]%in%CTS[[CT]]
      M[[CT]]<- t(sapply(1:4,function(i)apply(listSubstitution[INDEx,i,],2,sum)))
      rownames(M[[CT]]) <- NUCLEOTIDES
    }
    M[["E"]]<-M[[paste(substring(words(3,NUCLEOTIDES),1,1),"C",substring(words(3,NUCLEOTIDES),2,3),sep="")[1]]]
    for(CT in c(paste(substring(words(3,NUCLEOTIDES),1,1),"C",substring(words(3,NUCLEOTIDES),2,3),sep=""),paste(substring(words(3,NUCLEOTIDES),1,1),"G",substring(words(3,NUCLEOTIDES),2,3),sep=""),paste(substring(words(3,NUCLEOTIDES),1,1),"T",substring(words(3,NUCLEOTIDES),2,3),sep=""))[-1]){
      M[["E"]]<- M[["E"]]+M[[CT]]
    }
    rownames(M[["E"]]) <- NUCLEOTIDES



  simplifivemer<-function(fivemer=M,FIVEMER="CCATT",Thresh=20){
    Nuc=substr(FIVEMER,3,3)
    Nei=paste(substr(FIVEMER,1,2),substr(FIVEMER,4,5),collapse="",sep="")

    if(sum(fivemer[[Nei]][Nuc,])>Thresh & sum(fivemer[[Nei]][Nuc,]==0)==1){
      return(fivemer[[Nei]][Nuc,]);
    }
    else{
      if(substring(Nei,2,2)!="A"){
        FIVE=fivemer[[Nei]][Nuc,]
        for(i in 1:4){
          for(j in 1:4){
            MutatedNeighbor=paste(NUCLEOTIDES[i],substring(Nei,2,3),NUCLEOTIDES[j],collapse="",sep="")
            FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
          }
        }
        return(FIVE)
      }
      if(substring(Nei,2,2)=="A"){
        FIVE=fivemer[[paste(substring(Nei,1,1),canMutateTo(substring(Nei,2,2))[1],substring(Nei,3,4),collapse="",sep="")]][Nuc,]
        for(i in 2){
          for(j in 2:3){
            MutatedNeighbor=paste(substring(Nei,1,(i-1)),canMutateTo(substring(Nei,i,i))[j],substring(Nei,(i+1),4),collapse="",sep="")
            FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
          }
        }
        if(sum(FIVE)>Thresh & sum(FIVE==0)==1){
          return(FIVE)
        }
        else{
          for(i in 1:3){
            for(j in 1:3){
              MutatedNeighbor=paste(canMutateTo(substring(Nei,1,1))[i],canMutateTo(substring(Nei,2,2))[j],substring(Nei,3,4),collapse="",sep="")
              FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
              MutatedNeighbor=paste(substring(Nei,1,1),canMutateTo(substring(Nei,2,2))[j],canMutateTo(substring(Nei,3,3))[i],substring(Nei,4,4),collapse="",sep="")
              FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
              MutatedNeighbor=paste(substring(Nei,1,1),canMutateTo(substring(Nei,2,2))[j],substring(Nei,3,3),canMutateTo(substring(Nei,4,4))[i],collapse="",sep="")
              FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
            }
          }
        }
        return(FIVE)
      }
    }
  }

  substitutionModel <-sapply(words(5,NUCLEOTIDES),function(x)simplifivemer(M,x))
  NUCLEOTIDES <- c("A", "C", "G", "T", "N")
  return(substitutionModel)
}

