#' Count mutations by region
#'
#' This function counts the number of mutations in a sequence and groups them by region.
#'
#' @param    mutations  an array of mutations. The name indicates the position and the values are the kinds of mutations.
#' @return   an array of R_CDR, S_CDR, R_FWR and S_FWR
processNucMutations <- function(mutations) {
    regions <- factor(c(rep("FWR", 78), 
                        rep("CDR", 36), 
                        rep("FWR", 51), 
                        rep("CDR", 30), 
                        rep("FWR", 117)), 
                      levels=c("FWR", "CDR"))
    facmutations <- factor(mutations, levels=c("R", "S"))
    tableMutations <- table(regions[as.numeric(names(mutations))], facmutations)
    return(c(tableMutations["CDR", "R"], 
             tableMutations["CDR", "S"], 
             tableMutations["FWR", "R"], 
             tableMutations["FWR", "S"]))
}