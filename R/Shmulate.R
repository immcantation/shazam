# SHMulate

#' @include MutationProfiling.R
#' @include Shazam.R
NULL


#### shmulation ####

#' Simulate mutations in a single sequence
#'
#' @param input_seq    sequence in which mutations are to be introduced
#' @param num_muts     number of mutations to be introduced into \code{input_seq}
#'
#' @return A mutated sequence.
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


#' Simulate sequences to populate a tree
#'
#' \code{shmulateTree} returns a set of simulated sequences generated from an input sequence and an
#' igraph object. The input sequence is used to replace the founder node of the \code{igraph} lineage
#' tree and sequences are simulated with mutations corresponding to edge weights in the tree.
#' Sequences will not be generated for groups of nodes that are specified to be excluded.
#'
#' @param input_seq   sequence in which mutations are to be introduced.
#' @param graph       \code{igraph} object with vertex annotations whose edges are to be recreated.
#' @param field       annotation field to use for both unweighted path length exclusion and
#'                    consideration as a founder node. If \code{NULL} do not exclude any nodes.
#' @param exclude     vector of annotation values in the given field to exclude from potential
#'                    founder set. If \code{NULL} do not exclude any nodes. Has no effect if \code{field=NULL}.
#' @param jun_frac    fraction of characters in the junction region to add proportional number
#'                    of trunk mutations to the sequence.
#'
#' @return A \code{data.frame} of simulated sequences.
#' @export
shmulateTree <- function(input_seq, graph, field=NULL, exclude=NULL, jun_frac=NULL) {
  # Determine founder (mrca) of lineage tree
  # specify alakazam as opposed to ape::getMRCA
  founder_df <- alakazam::getMRCA(graph, path="distance", root="Germline", 
                                  field=field, exclude=exclude)
  
  # Get adjacency matrix
  adj <- get.adjacency(graph, attr="weight", sparse=F)
  # Get names of nodes for which sequences are not to be returned
  skip_names <- c()
  if (!is.null(field)) {
    g <- get.vertex.attribute(graph, name = field)
    g_names <- get.vertex.attribute(graph, name = 'name')
    skip_names <- g_names[g %in% exclude]
  }
  
  # Create data.frame to hold simulated sequences
  sim_tree <- data.frame('name'=founder_df$NAME[1],
                         'sequence'=input_seq, 'distance'=0,
                         stringsAsFactors=F)
  parent_nodes <- founder_df$NAME[1]
  nchil <- sum(adj[parent_nodes,]>0)
  
  # Add trunk mutations proportional to fraction of sequence in junction
  if(!is.null(jun_frac)) {
    adj[parent_nodes,] <- round(adj[parent_nodes,]*(1+jun_frac))
  }
  while(nchil>0) {
    new_parents <- c()
    # Loop through parent-children combos
    for(p in parent_nodes) {
      children <- colnames(adj)[adj[p,]>0]
      for(ch in children) {
        # Add child to new parents
        new_parents <- union(new_parents, ch)
        # Simulate sequence for that edge
        seq <- shmulateSeq(sim_tree$sequence[sim_tree$name==p], adj[p, ch])
        new_node <- data.frame('name'=ch, 'sequence'=seq, 'distance'=adj[p, ch])
        # Update output data.frame (bind_rows from dplyr)
        sim_tree <- bind_rows(sim_tree, new_node)
      }
    }
    # Re-calculate number of children
    parent_nodes <- new_parents
    nchil <- sum(adj[parent_nodes,]>0)
  }
  # Remove sequences that are to be excluded
  sim_tree <- subset(sim_tree, !(name %in% skip_names))
  return(sim_tree)
}


#### helper functions ####

#' Compute the mutations types
#'
#' @param seq   sequence for which to compute mutation types
#'
#' @return A \code{matrix} of mutation types for each position in the sequence.
#'
#' @details
#' For each position in the input sequence, use \code{CODON_TABLE} to
#' determine what types of mutations are possible. Returns \code{matrix}
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
#' @return A \code{list} of position being mutated and updated vector of mutated positions.
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
