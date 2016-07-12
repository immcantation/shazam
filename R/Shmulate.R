# SHMulate

#' @include MutationProfiling.R
#' @include shazam.R
NULL


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

