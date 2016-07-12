# SHMulate

#' @include shazam.R
NULL


#' Calculate targeting probability along sequence
#'
#' @param germlineSeq     input sequence
#' @param targetingModel  model underlying SHM to be used
#'
#' @return matrix of probabilities of each nucleotide being mutated to any other nucleotide.
#'
#' @details
#' Applies the targeting model to the input sequence to determine the probability
#' that a given nucleotide at a given position will mutate to any of the other
#' nucleotides. This is calculated for every position and every possible mutation.
calcTargeting <- function(germlineSeq, targetingModel=shazam::HS5FModel) {
  #* counts on constant variable NUCLEOTIDES (ACGTN-.)
  
  s_germlineSeq <- germlineSeq
  
  # Removing IMGT gaps (they should come in threes)
  # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
  gaplessSeq <- gsub("\\.\\.\\.", "XXX", s_germlineSeq)
  #If there is a single gap left convert it to an N
  gaplessSeq <- gsub("\\.", "N", gaplessSeq)
  
  # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
  s_germlineSeq <- gsub("XXX", "...", gaplessSeq)
  c_germlineSeq <- s2c(s_germlineSeq)
  
  # Matrix to hold targeting values for each position in c_germlineSeq
  gsTargeting <- matrix(NA, ncol=nchar(s_germlineSeq),
                        nrow=length(NUCLEOTIDES[1:5]), # ATGCN
                        dimnames=list(NUCLEOTIDES[1:5], c_germlineSeq))
  
  # Now remove the IMGT gaps so that the correct 5mers can be made to calculate
  # targeting. e.g.
  # GAGAAA......TAG yields: "GAGAA" "AGAAA" "GAAAT" "AAATA" "AATAG"
  # (because the IMGT gaps are NOT real gaps in sequence!!!)
  gaplessSeq <- gsub("\\.\\.\\.", "", s_germlineSeq)
  gaplessSeqLen <- nchar(gaplessSeq)
  
  # Slide through 5-mers and look up targeting
  gaplessSeq <- paste("NN", gaplessSeq, "NN", sep="")
  gaplessSeqLen <- nchar(gaplessSeq)
  pos<- 3:(gaplessSeqLen-2)
  subSeq =  substr(rep(gaplessSeq, gaplessSeqLen-4), (pos-2), (pos+2))
  
  gsTargeting_gapless <- sapply(subSeq, 
                                function(x){ targetingModel@targeting[NUCLEOTIDES[1:5], x] })
  
  gsTargeting[, c_germlineSeq!="."] <- gsTargeting_gapless
  
  # Set self-mutating targeting values to be NA
  mutatingToSelf <- colnames(gsTargeting)
  # print(mutatingToSelf[!(mutatingToSelf %in% NUCLEOTIDES[1:4])])
  mutatingToSelf[!(mutatingToSelf %in% NUCLEOTIDES[1:4])] <- "N" # ATGC
  
  tmp <- sapply(1:ncol(gsTargeting), function(pos){
    gsTargeting[mutatingToSelf[pos], pos] <<- NA }) # must use <<- (<- won't work)
  
  gsTargeting[!is.finite(gsTargeting)] <- NA
  
  return(gsTargeting[NUCLEOTIDES[1:4],]) # ATGC
}

#' Create all codons one mutation away from input codon.
#'
#' @param codon   starting codon to which mutations are added.
#'
#' @return a vector of codons.
#'
#' @details
#' All codons one mutation away from the input codon are generated.
allCodonMuts <- function(codon) {
  #* counts on constant variable NUCLEOTIDES (ACGTN-.)
  
  codon_char <- s2c(codon)
  matCodons <- t(array(codon_char, dim=c(3,12)))
  matCodons[1:4, 1] <- NUCLEOTIDES[1:4] # ATGC
  matCodons[5:8, 2] <- NUCLEOTIDES[1:4] # ATGC
  matCodons[9:12,3] <- NUCLEOTIDES[1:4] # ATGC
  return(apply(matCodons,1,c2s))
}


#' Asses if a mutation is R or S based on codon information.
#'
#' @param codonFrom     starting codon
#' @param codonTo       codon with mutation
#'
#' @return type of mutation, one of "S", "R", "Stop", or NA.
mutationType <- function(codonFrom, codonTo) {
  #* counts on constant variable AMINO_ACIDS
  
  if (is.na(codonFrom) | is.na(codonTo) | is.na(AMINO_ACIDS[codonFrom]) | is.na(AMINO_ACIDS[codonTo])) {
    mutationType <- NA
  } else {
    mutationType <- "S"
    if(AMINO_ACIDS[codonFrom] != AMINO_ACIDS[codonTo]) {
      mutationType <- "R"
    }
    if(AMINO_ACIDS[codonFrom]=="*" | AMINO_ACIDS[codonTo]=="*") {
      mutationType <- "Stop"
    }
  }
  return(mutationType)
}

