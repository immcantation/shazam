#' Count number of mutations in sequence
#'
#' This function counts the number of mutations in a sequence.
#'
#' @param   seqInput  Input sequence
#' @param   seqGL  Germline sequence
#' @return  array of observed mutations
#' 
#' @export
countMutations <- function(seqInput, seqGL) {
    if( is.na(c(seqInput, seqGL)) ) return(array(NA,4))
    seqI = s2c(seqInput)
    seqG = s2c(seqGL)
    matIGL = matrix(c(seqI,seqG),ncol=length(seqI),nrow=2,byrow=T)
    mutations <- analyzeMutations2NucUri(matIGL)

    if(is.na(mutations)){
      return(array(0,4))
    }else{
      return(processNucMutations(mutations))
    }
}
