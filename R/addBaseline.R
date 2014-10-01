addBASELINe <- function(db){
  baselineInfo <- t(sapply(1:nrow(db), function(i){fun_addBASELINe(db,i)}))
  colnames(baselineInfo) <- paste(c("CDR","FWR"),rep(c("Sigma", "CI95_Lower", "CI95_Upper"),each=2),sep="_")
  db <- cbind(db, baselineInfo)
  return(db)
}

# custom summary function
fun_addBASELINe <- function(db,i){
  baselineInfo_CDR <- calcBayesOutputInfo(groupPosteriors(db[i, "CDR_PDFs"]))
  baselineInfo_FWR <- calcBayesOutputInfo(groupPosteriors(db[i, "FWR_PDFs"]))
  # make a dataframe
  d <- array(c(
    CDR_Sigma = baselineInfo_CDR[1],
    CDR_CI95_Lower = baselineInfo_CDR[2],
    CDR_CI95_Upper = baselineInfo_CDR[3],
    FWR_Sigma = baselineInfo_FWR[1],
    FWR_CI95_Lower = baselineInfo_FWR[2],
    FWR_CI95_Upper = baselineInfo_FWR[3]
  ))
  return(d)
}


