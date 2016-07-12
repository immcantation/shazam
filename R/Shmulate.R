# SHMulate

#' @include MutationProfiling.R
#' @include shazam.R
NULL


#### shmulation ####


#' Simulate mutations in a single sequence
#'
#' @param input_seq    sequence in which mutations are to be introduced
#' @param num_muts     number of mutations to be introduced into input sequence
#'
#' @return mutated sequence.
#'
#' @details
#' Generates mutations in sequence one by one while updating targeting
#' probability of each position after each mutation.
#' @export
shmulateSeq <- function(input_seq, num_muts) {
  #* counts on constant variables CODON_TABLE, NUCLEOTIDES (ACTGN-.)
  
  # Trim sequence to last codon (getCodonPos from MutationProfiling.R)
  if(getCodonPos(nchar(input_seq))[3] > nchar(input_seq)) {
    sim_seq <- substr(input_seq, 1, getCodonPos(nchar(input_seq))[1]-1)
  } else {
    sim_seq <- input_seq
  }
  sim_seq <- gsub("\\.", "-", sim_seq)
  sim_leng <- nchar(sim_seq)
  stopifnot((sim_leng %% 3)==0)
  
  # Calculate possible mutations (given codon table)
  mutation_types <- computeMutationTypes(sim_seq)
  
  # Calculate probabilities of mutations at each position given targeting
  # from MutationProfiling.R; includes a N row
  targeting <- calculateTargeting(germlineSeq = sim_seq) 
  # keep only ACGT rows
  targeting <- targeting[NUCLEOTIDES[1:4], ] 
  # set NA to 0
  targeting[is.na(targeting)] <- 0 
  # Make probability of stop codon 0
  targeting[mutation_types=="Stop"] <- 0
  
  # Initialize counters
  total_muts <- 0
  positions <- numeric(num_muts)
  
  while(total_muts < num_muts) {
    # Get position to mutate and update counters
    mutpos <- sampleMut(sim_leng, targeting, positions)
    total_muts <- total_muts + 1
    positions[total_muts] <- mutpos$pos
    
    # Implement mutation in simulation sequence
    mut_nuc <- 4 - (4*mutpos$pos - mutpos$mut)
    sim_char <- s2c(sim_seq)
    sim_char[mutpos$pos] <- NUCLEOTIDES[mut_nuc]
    sim_seq <- c2s(sim_char)
    
    # Update targeting
    lower <- max(mutpos$pos-4, 1)
    upper <- min(mutpos$pos+4, sim_leng)
    targeting[, lower:upper] <- calculateTargeting(germlineSeq = 
                                                     substr(sim_seq, lower, upper))[NUCLEOTIDES[1:4], ]
    targeting[is.na(targeting)] <- 0
    
    # Update possible mutations
    lower <- getCodonPos(lower)[1]
    upper <- getCodonPos(upper)[3]
    mutation_types[, lower:upper] <- computeMutationTypes(substr(sim_seq, lower, upper))
    # Make probability of stop codon 0
    if(any(mutation_types[, lower:upper]=="Stop", na.rm=T)) {
      targeting[, lower:upper][mutation_types[, lower:upper]=="Stop"] <- 0
    }
  }
  return(sim_seq)
}


#### helper functions ####

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
