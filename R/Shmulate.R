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


#' Compute the mutations types
#'
#' @param seq   sequence for which to compute mutation types
#'
#' @return matrix of mutation types for each position in the sequence.
#'
#' @details
#' For each position in the input sequence, use the codon table to
#' determine what types of mutations are possible. Returns matrix
#' of all possible mutations and corresponding types.
computeMutationTypes <- function(seq){
  #* counts on constant variable CODON_TABLE, NUCLEOTIDES (ACTGN-.)
  #* caution: this breaks down if length of seq is not a multiple of 3
  
  leng_seq <- nchar(seq)
  try(if( (leng_seq %%3 !=0) ) stop("length of input sequence must be a multiple of 3"))
  
  codons <- sapply(seq(1, leng_seq, by=3), function(x) {substr(seq,x,x+2)})
  mut_types <- matrix(unlist(CODON_TABLE[, codons]), ncol=leng_seq, nrow=4, byrow=F)
  dimnames(mut_types) <-  list(NUCLEOTIDES[1:4], 1:leng_seq)
  return(mut_types)
}

#' Find encompassing codon
#'
#' @param nuc_pos    position for which codon is to be found
#' @param frame      reading frame in which to determine codon
#'
#' @return vector of positions of codon encompassing input position.
#'
#' @details
#' Given a nuclotide position, find the positions of the three nucleotides
#' that encompass the codon in the given reading frame of the sequence.
# e.g. nuc 86 is part of nucs 85,86,87
getCodonPos <- function(nuc_pos, frame=0) {
  codon_num <- ( ceiling((nuc_pos + frame) / 3) ) * 3
  codon <- (codon_num-2):codon_num
  return(codon)
}

#' Pick a position to mutate
#'
#' @param sim_leng      length of sequence in which mutation is being simulated
#' @param targeting     probabilities of each position in the sequence being mutated
#' @param positions     vector of positions which have already been mutated
#'
#' @return list of position being mutated and updated vector of mutated positions.
#'
#' @details
#' Sample positions in the sequence to mutate given targeting probability
#' until a new position is selected. This new position is then added to the
#' vector of mutated positions and returned.
sampleMut <- function(sim_leng, targeting, positions) {
  pos <- 0
  # Sample mutations until new position is selected
  while (pos %in% positions) {
    # Randomly select a mutation
    mut <- sample(1:(4*sim_leng), 1, replace=F, prob=as.vector(targeting))
    pos <- ceiling(mut/4)
  }
  return(list(mut=mut, pos=pos))
}