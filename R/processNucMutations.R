#' Count mutations by region
#'
#' This function counts the number of mutations in a sequence and groups them by region.
#'
#' @param   mutations  an array of mutations. The name indicates the position and the values are the kinds of mutations.
#' @return  retArrary  an array of R_CDR, S_CDR, R_FWR and S_FWR
processNucMutations <- function(mutations){
  facmutations <- factor(mutations, levels=c("R","S") )
  tableMutations <- table(REGIONS[as.numeric(names(mutations))],facmutations)
  return( c(tableMutations["CDR","R"], tableMutations["CDR","S"], tableMutations["FWR","R"], tableMutations["FWR","S"]) )
}