CODON_TABLE <- as.data.frame(matrix(NA,ncol=64,nrow=12))

intCounter = 1
for(pOne in NUCLEOTIDES[1:4]){
  for(pTwo in NUCLEOTIDES[1:4]){
    for(pThree in NUCLEOTIDES[1:4]){
      codon = paste(pOne,pTwo,pThree,sep="")
      colnames(CODON_TABLE)[intCounter] =  codon
      intCounter = intCounter + 1
      CODON_TABLE[,codon] = mutationTypeOptimized(cbind(permutateAllCodon(codon),rep(codon,12)))
    }  
  }
}

chars = c("N","A","C","G","T", ".")
for(a in chars){
  for(b in chars){
    for(c in chars){
      if(a=="N" | b=="N" | c=="N"){ 
        #cat(paste(a,b,c),sep="","\n") 
        CODON_TABLE[,paste(a,b,c,sep="")] = rep(NA,12)
      }
    }  
  }
}

chars = c("N","A","C","G","T", ".")
for(a in chars){
  for(b in chars){
    for(c in chars){
      if(a=="." | b=="." | c=="."){ 
        #cat(paste(a,b,c),sep="","\n") 
        CODON_TABLE[,paste(a,b,c,sep="")] = rep(NA,12)
      }
    }  
  }
}

CODON_TABLE <- as.matrix(CODON_TABLE,nrow(CODON_TABLE),ncol(CODON_TABLE))
save(CODON_TABLE,file="data/CODON_TABLE.RData")
