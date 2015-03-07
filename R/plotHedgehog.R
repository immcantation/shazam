#' Plot Mutabilities in Hedgehog circle format
#'
#'
#' @param    mutabilityModel    matrix of mutabilities from createMutabilityModel()
#' @param    outputName         Filename for the output files
#' @param    outputPath         Path to save output
#' @export
plotHedgehog <- function(mutabilityModel, outputName=NULL, outputPath=NULL) {

  NUCLEOTIDES <- c("A","C","G","T")
  Mutability <- mutabilityModel


  x<-list()
  y<-list()
  Mu<-list()
  SAMPLE = "SAMPLE"
  SAMPLES = "SAMPLE"

  MutabilityMatrix<-sapply(Mutability,function(x)x[[1]][1:1024])
  MutabilityWeights<-sapply(Mutability,function(x)x[[2]][1])
  MutabilityMutations<-sapply(Mutability,function(x)x[[3]][1])
  Mu[[SAMPLE]]<-MutabilityMutations

  Template<-Mutability[[1]][[1]]
  Template[is.na(Template)]<-0
  Template[Template>0]<-0

  tmp<-Aggregate_Mutability(MutabilityMatrix,MutabilityWeights,Template)
  tmp[[1]][tmp[[1]]==0]<-NA
  y[[SAMPLE]]<-tmp
  x[[SAMPLE]]<-sort(tmp[[1]][tmp[[2]]>100])


  MutabilityMatrix<-sapply(y,function(x)x[[1]][1:1024])
  MutabilityWeights<-sapply(y,function(x)x[[2]][1:1024])

  L<-ncol(MutabilityMatrix)
  tmp9<-sapply(1:1024,function(i)weighted.mean(MutabilityMatrix[i,],MutabilityWeights[i,],na.rm=TRUE))
  MuNumbers<-apply(MutabilityWeights,1,sum)
  names(tmp9)<-names(Template)
  Mutability_ALL<-sort(tmp9[MuNumbers>500])
  tmp19<-tmp9/sum(tmp9[MuNumbers>500],na.rm=TRUE)*sum(!is.na(tmp9[MuNumbers>500]))

  MutabilityMatrixNorm<-apply(MutabilityMatrix,2,function(x)x/sum(x,na.rm=TRUE)*sum(!is.na(x)))


  tmp99<-sapply(1:1024,function(i)wt.sd(MutabilityMatrixNorm[i,],MutabilityWeights[i,]))
  Mutability_ALL_CI<-cbind(tmp19-tmp99,tmp19+tmp99)
  Mutability_ALL_CI[MuNumbers<=500,]<-NA
  rownames(Mutability_ALL_CI)<-names(Template)
  Mutability_ALL_CI[Mutability_ALL_CI[,1]<0,1]<-0


  NUCLEOTIDES <- c("G","T","A","C")
  tr<-NUCLEOTIDES
  names(tr)<-c("G","T","A","C")

  comp<-NUCLEOTIDES
  names(comp)<-c("T","G","C","A")

  tr_comp<-NUCLEOTIDES
  names(tr_comp)<-c("C","A","T","G")


  Mutability_ALL_AVE<-(sapply(words(5,NUCLEOTIDES),Fill_HOT))
  for(i in names(which(is.na(Mutability_ALL_AVE)))){
    Mutability_ALL_AVE[i]<-Fill_HOT(i,mutability=Mutability_ALL_AVE)
  }
  for(i in names(which((Mutability_ALL_AVE)<0.000001))){
    Mutability_ALL_AVE[i]<-Fill_HOT(i,mutability=Mutability_ALL_AVE)
  }

  Mutability_ALL_AVE[is.nan(Mutability_ALL_AVE)]=0 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AC added!!
  Mutability_ALL_AVE<-sort(Mutability_ALL_AVE)


  Mutability_CI1_AVE<-(sapply(words(5,NUCLEOTIDES),function(x)Fill_HOT(x,mutability=Mutability_ALL_CI[!is.na(Mutability_ALL_CI[,1]),1])))

  for(i in names(which(is.na(Mutability_CI1_AVE)))){
    Mutability_CI1_AVE[i]<-Fill_HOT(i,mutability=Mutability_CI1_AVE)
  }
  for(i in names(which((Mutability_CI1_AVE)<0.000001))){
    Mutability_CI1_AVE[i]<-Fill_HOT(i,mutability=Mutability_CI1_AVE)
  }
  Mutability_CI1_AVE<-sort(Mutability_CI1_AVE)


  Mutability_CI2_AVE<-(sapply(words(5,NUCLEOTIDES),function(x)Fill_HOT(x,mutability=Mutability_ALL_CI[!is.na(Mutability_ALL_CI[,2]),2])))

  for(i in names(which(is.na(Mutability_CI2_AVE)))){
    Mutability_CI2_AVE[i]<-Fill_HOT(i,mutability=Mutability_CI2_AVE)
  }
  for(i in names(which((Mutability_CI2_AVE)<0.000001))){
    Mutability_CI2_AVE[i]<-Fill_HOT(i,mutability=Mutability_CI2_AVE)
  }
  Mutability_CI2_AVE<-sort(Mutability_CI2_AVE)


  MUTABILITY<-Mutability_ALL_AVE
  MUTABILITY<-MUTABILITY/sum(MUTABILITY)*1024

  OUT_CIRCLE = file.path(outputPath, paste0(outputName, ".pdf"))
  OUT_MUT_CSV = file.path(outputPath, paste0(outputName, ".csv"))


  WORDS = words(4,NUCLEOTIDES)
  VH_FAMILIES <- paste(rep("VH",7),1:7,sep="")
  ISOTYPES = paste(rep("IGH",4),c("A","D","G","M"),sep="")

  X1<-rep(WORDS,7*1*length(SAMPLES))
  X2<-rep(VH_FAMILIES,each=256,times=1*length(SAMPLES))
  X4<-rep(SAMPLES,each=256*7*1)

  arrNames <- paste(X4,X2,X1,sep="_")
  rm( X4,X2,X1)


  WRC = c(grep("[AT][GA]C..", words(5, s2c("AGCT")), perl=TRUE, value=TRUE),
          grep("..G[CT][AT]", words(5, s2c("AGCT")), perl=TRUE, value=TRUE))

  SYC = c(grep("[CG][CT]C..", words(5, s2c("AGCT")), perl=TRUE, value=TRUE),
          grep("..G[GA][CG]", words(5, s2c("AGCT")), perl=TRUE, value=TRUE))

  WT = c(grep(".[AT]A..", words(5, s2c("AGCT")), perl=TRUE, value=TRUE),
         grep("..T[AT].", words(5, s2c("AGCT")), perl=TRUE, value=TRUE))

  writeLines("Fivemer Mutability Source lower25 upper25", OUT_MUT_CSV)
  MUT4FILE<-cbind(MUTABILITY,rep("Inferred",1024))
  MUT4FILE[rownames(MUT4FILE)%in%names(Mutability_ALL),2]<-"Measured"
  MUT4FILE<-cbind(MUT4FILE,(Mutability_CI1_AVE[rownames(MUT4FILE)]),(Mutability_CI2_AVE[rownames(MUT4FILE)]))
  write.table(MUT4FILE, OUT_MUT_CSV, append=TRUE,sep=" ",col.names=FALSE)

  Mutability_ALL_NA<-Mutability_ALL_AVE
  Mutability_ALL_NA[!names(Mutability_ALL_AVE)%in%names(Mutability_ALL)]<-NA

  Mutability_ALL_NA_NORM<-Mutability_ALL_NA/sum(Mutability_ALL_NA,na.rm=TRUE)*sum(!is.na(Mutability_ALL_NA))


  # SCALES=c(135,135,85,135)
  SCALES=c(85,85,85,85)/85*90
  N<-1024
  Yshift = 1.5

  pdf(OUT_CIRCLE,30,10)
  par(mar=c(0,0,1.3,0))
  cols=c(rgb(1,1,0.5),rgb(0.5,0.5,1),rgb(0.5,1,0.5),rgb(1,0.5,0.5))
  plot(-5:5,xlim=c(-30,30),ylim=c(-9,11),pch=NA,xlab="",ylab="",axes=FALSE,main="",cex.main=2)
  for(NUC in NUCLEOTIDES[c(1,4,2,3)]){
    ##Colored rings
    X=c(-25.5+5,-12.+5,-1+5,11.5+5)[match(NUC,NUCLEOTIDES[c(1,4,2,3)])]*1.2 ## position of the center of the circle
    ##Outer ring with mutability values
    par(xpd=TRUE)
    if(match(NUC,NUCLEOTIDES[c(1,4,2,3)])%in%c(1,3)){##If Nucleotide is A or C

      #!!!!!!!!!!!!!   ark(xshift= X, yshift= Yshift,level=1,R=0.4,start=pi/2*(i-1),theta=2*pi,col=NA,lwd=1,border=1) ##inner circle for 5'/3'
      for(i in 1:4)ark(xshift= X, yshift= Yshift,level=1,R=0.8,start=pi/2*(i-1),theta=pi/2,col=cols[i],lwd=0.001,ADD=0.4) ##Ring #1
      for(i in 1:4)for(j in 1:4) ark(xshift= X, yshift= Yshift,level=2,R=0.8,start=pi/2*(i-1)+pi/2*(j-1)/4,theta=pi/2/4,col=cols[j],lwd=0.001,ADD=0.4)  ##Ring #2
      for(i in 1:4)for(j in 1:4) ark(xshift= X, yshift= Yshift,level=3,R=0.8,start=0,theta=pi*2,col=cols[match(NUC,NUCLEOTIDES[])],lwd=0.001,ADD=0.4) ##Ring #3
      for(i in 1:4)for(j in 1:4)for(k in 1:4) ark(xshift= X, yshift= Yshift,level=4,R=0.8,start=pi/2*(i-1)+pi/2*(j-1)/4+pi/2*(k-1)/4/4,theta=pi/2/4/4,col=cols[k],lwd=0.001,ADD=0.4) ##Ring #4
      for(i in 1:4)for(j in 1:4)for(k in 1:4) for(l in 1:4) ark(xshift= X, yshift= Yshift,level=5,R=0.8,start=pi/2*(i-1)+pi/2*(j-1)/4+pi/2*(k-1)/4/4+pi/2*(l-1)/4/4/4,theta=pi/2/4/4/4,col=cols[l],lwd=0.001,ADD=0.4) ##Ring #5
      ##Text inside rings
      for(i in 1:4)text(X+rt2xy(0.4+.4,pi/4+pi/2*(i-1))[1,1],rt2xy(0.4+.4,pi/4+pi/2*(i-1))[1,2]+Yshift,NUCLEOTIDES[i]) ##Ring #1
      for(i in 1:4)for(j in 1:4)text(X+rt2xy(0.4+1.2,pi/4/4+pi/2/4*((i-1)*(4)+j-1))[1,1],Yshift+rt2xy(0.4+1.2,pi/4/4+pi/2/4*((i-1)*4+j-1))[1,2],NUCLEOTIDES[j]) ##Ring #2
      for(i in 1:4)for(j in 1:4)text(X+rt2xy(0.4+2.,pi/4+pi/2*(i-1))[1,1],Yshift+rt2xy(0.4+2.,pi/4+pi/2*(i-1))[1,2],NUC)  ##Ring #3
      for(i in 1:4)for(j in 1:4)for(k in 1:4)text(X+rt2xy(0.4+2.8,pi/4/4/4+pi/2/4/4*((i-1)*(16)+(j-1)*4+k-1))[1,1],Yshift+rt2xy(0.4+2.8,pi/4/4/4+pi/2/4/4*((i-1)*16+(j-1)*4+k-1))[1,2],NUCLEOTIDES[k],cex=0.5)  ##Ring #4
      for(i in 1:4)for(j in 1:4)for(k in 1:4)for(l in 1:4)text(X+rt2xy(0.4+3.6,pi/4/4/4/4+pi/2/4/4/4*((i-1)*(16)*4+(j-1)*4*4+(k-1)*4+l-1))[1,1],Yshift+rt2xy(0.4+3.6,pi/4/4/4/4+pi/2/4/4/4*((i-1)*16*4+(j-1)*4*4+(k-1)*4+l-1))[1,2],NUCLEOTIDES[l],cex=0.25)  ##Ring #5
      ##define colors
      COLS<-rep(8,length(words(5,NUCLEOTIDES)[substring(words(5,NUCLEOTIDES),3,3)==NUC]))
      COLS[words(5,NUCLEOTIDES)[substring(words(5,NUCLEOTIDES),3,3)==NUC]%in%WRC]<-2
      COLS[words(5,NUCLEOTIDES)[substring(words(5,NUCLEOTIDES),3,3)==NUC]%in%SYC]<-4
      COLS[words(5,NUCLEOTIDES)[substring(words(5,NUCLEOTIDES),3,3)==NUC]%in%WT]<-3
      text(X,11.5,paste("NN",NUC,"NN",sep=""),cex=2) ## title for the panel


      for(i in 1:4)for(j in 1:4)for(k in 1:4)for(l in 1:4) ark(xshift= X, yshift= Yshift,level=6,R=0.8,scale=(Mutability_ALL_NA_NORM[words(5,NUCLEOTIDES)[substring(words(5,NUCLEOTIDES),3,3)==NUC]]*1)[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+l],start=pi/2*(i-1)+pi/2*(j-1)/4+pi/2*(k-1)/4/4+pi/2*(l-1)/4/4/4,theta=pi/2/4/4/4,col=COLS[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+(l)],lwd=0.001,ADD=0.4)##Plot mutabilities

      for(i in 1:4)for(j in 1:4)for(k in 1:4)for(l in 1:4) arkErr(xshift= X, yshift= Yshift,level=6,R=0.8,scale=c((Mutability_ALL_CI[words(5,NUCLEOTIDES)[substring(words(5,NUCLEOTIDES),3,3)==NUC],1]*1)[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+l],(Mutability_ALL_CI[words(5,NUCLEOTIDES)[substring(words(5,NUCLEOTIDES),3,3)==NUC],2]*1)[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+l]),start=1.5*pi/2/4/4/4/4+pi/2*(i-1)+pi/2*(j-1)/4+pi/2*(k-1)/4/4+pi/2*(l-1)/4/4/4,theta=pi/2/4/4/4/4,col=1,lwd=0.001,ADD=0.4)##Plot Error bars
      text(X+(0.0),Yshift+0.0,"5\'",cex=2)
      text(X+(-3.1),6.8,"3\'",cex=2)
    }

    if(match(NUC,NUCLEOTIDES[c(1,4,2,3)])%in%c(2,4)){##If Nucleotide is T or G


      #!!!!!!!!!!! ark(xshift= X, yshift= Yshift,level=1,R=0.4,start=pi/2*(i-1),theta=2*pi,col=NA,lwd=1,border=1) ##inner circle for 5'/3'
      for(i in 1:4)ark(xshift= X, yshift= Yshift,level=1,R=0.8,start=pi/2*(i-1),theta=pi/2,col=cols[c(4,3,2,1)][i],lwd=0.001,ADD=0.4) ##Ring #1
      for(i in 1:4)for(j in 1:4) ark(xshift= X, yshift= Yshift,level=2,R=0.8,start=pi/2*(i-1)+pi/2*(j-1)/4,theta=pi/2/4,col=cols[c(4,3,2,1)][j],lwd=0.001,ADD=0.4)  ##Ring #2
      for(i in 1:4)for(j in 1:4) ark(xshift= X, yshift= Yshift,level=3,R=0.8,start=0,theta=pi*2,col=cols[c(4,3,2,1)][match(NUC,NUCLEOTIDES[c(4,3,2,1)])],lwd=0.001,ADD=0.4) ##Ring #3
      for(i in 1:4)for(j in 1:4)for(k in 1:4) ark(xshift= X, yshift= Yshift,level=4,R=0.8,start=pi/2*(i-1)+pi/2*(j-1)/4+pi/2*(k-1)/4/4,theta=pi/2/4/4,col=cols[c(4,3,2,1)][k],lwd=0.001,ADD=0.4) ##Ring #4
      for(i in 1:4)for(j in 1:4)for(k in 1:4) for(l in 1:4) ark(xshift= X, yshift= Yshift,level=5,R=0.8,start=pi/2*(i-1)+pi/2*(j-1)/4+pi/2*(k-1)/4/4+pi/2*(l-1)/4/4/4,theta=pi/2/4/4/4,col=cols[c(4,3,2,1)][l],lwd=0.001,ADD=0.4) ##Ring #5
      ##Text inside rings
      for(i in 1:4)text(X+rt2xy(0.4+.4,pi/4+pi/2*(i-1))[1,1],Yshift+rt2xy(0.4+.4,pi/4+pi/2*(i-1))[1,2],NUCLEOTIDES[c(4,3,2,1)][i]) ##Ring #1
      for(i in 1:4)for(j in 1:4)text(X+rt2xy(0.4+1.2,pi/4/4+pi/2/4*((i-1)*(4)+j-1))[1,1],Yshift+rt2xy(0.4+1.2,pi/4/4+pi/2/4*((i-1)*4+j-1))[1,2],NUCLEOTIDES[c(4,3,2,1)][j]) ##Ring #2
      for(i in 1:4)for(j in 1:4)text(X+rt2xy(0.4+2.,pi/4+pi/2*(i-1))[1,1],Yshift+rt2xy(0.4+2.,pi/4+pi/2*(i-1))[1,2],NUC)  ##Ring #3
      for(i in 1:4)for(j in 1:4)for(k in 1:4)text(X+rt2xy(0.4+2.8,pi/4/4/4+pi/2/4/4*((i-1)*(16)+(j-1)*4+k-1))[1,1],Yshift+rt2xy(0.4+2.8,pi/4/4/4+pi/2/4/4*((i-1)*16+(j-1)*4+k-1))[1,2],NUCLEOTIDES[c(4,3,2,1)][k],cex=0.5)  ##Ring #4
      for(i in 1:4)for(j in 1:4)for(k in 1:4)for(l in 1:4)text(X+rt2xy(0.4+3.6,pi/4/4/4/4+pi/2/4/4/4*((i-1)*(16)*4+(j-1)*4*4+(k-1)*4+l-1))[1,1],Yshift+rt2xy(0.4+3.6,pi/4/4/4/4+pi/2/4/4/4*((i-1)*16*4+(j-1)*4*4+(k-1)*4+l-1))[1,2],NUCLEOTIDES[c(4,3,2,1)][l],cex=0.25)  ##Ring #5
      ##define colors
      COLS<-rep(8,length(words(5,NUCLEOTIDES[c(4,3,2,1)])[substring(words(5,NUCLEOTIDES[c(4,3,2,1)]),3,3)==NUC]))
      COLS[words(5,NUCLEOTIDES[c(4,3,2,1)])[substring(words(5,NUCLEOTIDES[c(4,3,2,1)]),3,3)==NUC]%in%WRC]<-2
      COLS[words(5,NUCLEOTIDES[c(4,3,2,1)])[substring(words(5,NUCLEOTIDES[c(4,3,2,1)]),3,3)==NUC]%in%SYC]<-4
      COLS[words(5,NUCLEOTIDES[c(4,3,2,1)])[substring(words(5,NUCLEOTIDES[c(4,3,2,1)]),3,3)==NUC]%in%WT]<-3
      text(X,11.5,paste("NN",NUC,"NN",sep=""),cex=2) ## title for the panel



      for(i in 1:4)for(j in 1:4)for(k in 1:4)for(l in 1:4) ark(xshift= X, yshift= Yshift,level=6,R=0.8,scale=(Mutability_ALL_NA_NORM[words(5,NUCLEOTIDES[c(4,3,2,1)])[substring(words(5,NUCLEOTIDES[c(4,3,2,1)]),3,3)==NUC]]*1)[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+l],start=pi/2*(l-1)+pi/2*(k-1)/4+pi/2*(j-1)/4/4+pi/2*(i-1)/4/4/4,theta=pi/2/4/4/4,col=COLS[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+(l)],lwd=0.001,ADD=0.4)##Plot mutabilities

      for(i in 1:4)for(j in 1:4)for(k in 1:4)for(l in 1:4) arkErr(xshift= X, yshift= Yshift,level=6,R=0.8,scale=c((Mutability_ALL_CI[words(5,NUCLEOTIDES[c(4,3,2,1)])[substring(words(5,NUCLEOTIDES[c(4,3,2,1)]),3,3)==NUC],1]*1)[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+l],(Mutability_ALL_CI[words(5,NUCLEOTIDES[c(4,3,2,1)])[substring(words(5,NUCLEOTIDES[c(4,3,2,1)]),3,3)==NUC],2]*1)[(i-1)*4*4*4+(j-1)*4*4+(k-1)*4+l]),start=1.5*pi/2/4/4/4/4+pi/2*(l-1)+pi/2*(k-1)/4+pi/2*(j-1)/4/4+pi/2*(i-1)/4/4/4,theta=pi/2/4/4/4/4,col=1,lwd=0.001,ADD=0.4)##Plot Error bars
      text(X+(0.0),Yshift+0.0,"3\'",cex=2)
      text(X+(-3.1),6.8,"5\'",cex=2)
    }
  }
  legend('topleft',c('WRC/GYW','WA/TW','SYC/GRS','Neutral'),col=c(2,3,4,8),pch=15,text.col=c(2,3,4,8),box.col=NA,cex=2)

  dev.off()
  return(NULL)
}



Aggregate_Mutability<-function(MutabilityMatrix,MutabilityWeights,Template){
  L<-ncol(MutabilityMatrix)
  tmp<-sapply(1:1024,function(i)weighted.mean(MutabilityMatrix[i,],MutabilityWeights,na.rm=TRUE))
  MuNumbers<-sapply(1:1024,function(i)sum(MutabilityWeights[!is.na(MutabilityMatrix[i,])]))
  names(tmp)<-names(Template)
  return(list(tmp,MuNumbers))
}



Fill_HOT<-function(FIVEMER,mutability=Mutability_ALL){
  if(FIVEMER%in%names(mutability)){
    if(!is.na(mutability[[FIVEMER]])){
      if(mutability[[FIVEMER]]>=0.0){
        return(mutability[[FIVEMER]])
      }
    }
  }
  return(0)
}

rt2xy<-function(r,theta){
  x=r*cos(theta)
  y=r*sin(theta)
  return(cbind(x,y))
}


arkErr<-function(level=4,R=1,N=2048,theta=pi,start=pi/4,col=2,lwd=1,scale=c(R,R),xshift=0,yshift=0,ADD=0.0){
  tmp<-rt2xy(ADD+(level-1)*R+scale[1],seq(start,start+theta,length.out=N*(theta)/2/pi))
  x1<-tmp[,1]
  y1<-tmp[,2]
  tmp<-rt2xy(ADD+(level-1)*R+scale[2],seq(start,start+theta,length.out=N*(theta)/2/pi))
  x2<-tmp[,1]
  y2<-tmp[,2]
  polygon(xshift+c(x1,rev(x2)),yshift+c(y1,rev(y2)),col=col,border=col,lwd=lwd)
}



ark<-function(level=4,R=1,N=2048,theta=pi,start=pi/4,col=2,lwd=1,scale=R,xshift=0,yshift=0,ADD=0,border=col){
  tmp<-rt2xy(ADD+(level-1)*R,seq(start,start+theta,length.out=N*(theta)/2/pi))
  x1<-tmp[,1]
  y1<-tmp[,2]
  tmp<-rt2xy(ADD+(level-1)*R+scale,seq(start,start+theta,length.out=N*(theta)/2/pi))
  x2<-tmp[,1]
  y2<-tmp[,2]
  polygon(xshift+c(x1,rev(x2)),yshift+c(y1,rev(y2)),col=col,border=col,lwd=lwd)
}

ark3v<-function(level=4,R=1,N=2048,theta=pi,start=pi/4,col=c(2,1,3),col2=1,lwd=1,scale=c(R1,R2,R3),xshift=0,yshift=0,ADD=0,border=col){
  index=1
  for(S in c(scale[1]+scale[2]+scale[3],scale[1]+scale[2],scale[1])){
    tmp<-rt2xy(ADD+(level-1)*R,seq(start,start+theta,length.out=N*(theta)/2/pi))
    x1<-tmp[,1]
    y1<-tmp[,2]
    tmp<-rt2xy(ADD+(level-1)*R+S,seq(start,start+theta,length.out=N*(theta)/2/pi))
    x2<-tmp[,1]
    y2<-tmp[,2]
    polygon(xshift+c(x1,rev(x2)),yshift+c(y1,rev(y2)),col=col[4-index],border=col2,lwd=lwd)
    index=index+1
  }
}
