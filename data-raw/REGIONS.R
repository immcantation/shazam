# IMGT_V_NO_CDR3 <- factor( c( rep("FWR", 78), rep("CDR", 36),  rep("FWR", 51), rep("CDR", 30), rep("FWR", 117) ),
#                             levels = c("FWR", "CDR")
#                         )
#                     
# save(IMGT_V_NO_CDR3, file="data/IMGT_V_NO_CDR3.RData")


region <- "1:26:38:55:65:104:-"
region <- as.numeric(strsplit(region,":")[[1]])
# FWR/CDR boundaries
flagTrim <- F
if( is.na(region[7])){
  flagTrim <- T
  region[7]<-region[6]
}
readStart = min(region,na.rm=T)
readEnd = max(region,na.rm=T)
if(readStart>1){
  region = region - (readStart - 1)
}
region_Nuc = c( (region[1]*3-2) , (region[2:7]*3) )
region_Cod = region

readStart = (readStart*3)-2
readEnd = (readEnd*3)

FWR_Nuc <- c( rep(TRUE,(region_Nuc[2])),
              rep(FALSE,(region_Nuc[3]-region_Nuc[2])),
              rep(TRUE,(region_Nuc[4]-region_Nuc[3])),
              rep(FALSE,(region_Nuc[5]-region_Nuc[4])),
              rep(TRUE,(region_Nuc[6]-region_Nuc[5])),
              rep(FALSE,(region_Nuc[7]-region_Nuc[6]))
)
CDR_Nuc <- (1-FWR_Nuc)
CDR_Nuc <- as.logical(CDR_Nuc)

  FWR_Nuc_Mat <- matrix( rep(FWR_Nuc,6), ncol=length(FWR_Nuc), nrow=5, byrow=T)
  CDR_Nuc_Mat <- matrix( rep(CDR_Nuc,6), ncol=length(CDR_Nuc), nrow=5, byrow=T)

FWR_Codon <- c( rep(TRUE,(region[2])),
                rep(FALSE,(region[3]-region[2])),
                rep(TRUE,(region[4]-region[3])),
                rep(FALSE,(region[5]-region[4])),
                rep(TRUE,(region[6]-region[5])),
                rep(FALSE,(region[7]-region[6]))
)
CDR_Codon <- (1-FWR_Codon)
CDR_Codon <- as.logical(CDR_Codon)


save(FWR_Nuc_Mat,file="data/FWR_Nuc_Mat.RData")
save(CDR_Nuc_Mat,file="data/CDR_Nuc_Mat.RData")

