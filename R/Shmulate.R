# SHMulate

#' @include MutationProfiling.R
#' @include Shazam.R
NULL


#### SHMulation ####

#' Simulate mutations in a single sequence
#'
#' Generates random mutations in a sequence iteratively using a targeting model.
#' Targeting probabilities at each position are updated after each iteration.
#' 
#' @param    sequence        sequence string in which mutations are to be introduced.
#'                           Accepted alphabet: \code{\{A, T, G, C, N, .\}}. Note
#'                           that \code{-} is not accepted.
#' @param    numMutations    number of mutations to be introduced into \code{sequence}.
#' @param    targetingModel  5-mer \link{TargetingModel} object to be used for computing 
#'                           probabilities of mutations at each position. Defaults to
#'                           \link{HH_S5F}.
#'                           
#' @return   A string defining the mutated sequence.
#' 
#' @details
#' If the input \code{sequence} has a non-triplet overhang at the end, it will be trimmed
#' to the last codon. For example, \code{ATGCATGC} will be trimmed to \code{ATGCAT}.
#' 
#' Mutations are not introduced to positions in the input \code{sequence} that contain 
#' \code{.} or \code{N}.
#' 
#' @seealso  See \link{shmulateTree} for imposing mutations on a lineage tree. 
#'           See \link{HH_S5F} and \link{MK_RS5NF} for predefined 
#'           \link{TargetingModel} objects.
#' 
#' @examples
#' # Define example input sequence
#' sequence <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATA.TTTA"
#' 
#' # Simulate using the default human 5-mer targeting model
#' shmulateSeq(sequence, numMutations=6)
#' 
#' @export
shmulateSeq <- function(sequence, numMutations, targetingModel=HH_S5F) {
    #* counts on constant variables CODON_TABLE, NUCLEOTIDES (ACTGN-.)
    
    # Check targeting model
    if (!is(targetingModel, "TargetingModel")) {
        stop(deparse(substitute(targetingModel)), " is not a valid TargetingModel object")
    }

    # Trim sequence to last codon (getCodonPos from MutationProfiling.R)
    if(getCodonPos(stri_length(sequence))[3] > stri_length(sequence)) {
        sim_seq <- stri_sub(str=sequence, from=1, 
                            to=getCodonPos(stri_length(sequence))[1]-1)
    } else {
        sim_seq <- sequence
    }
    
    sim_leng <- stri_length(sim_seq)
    stopifnot((sim_leng %% 3)==0)
    
    # Calculate possible mutations (given codon table)
    mutation_types <- computeMutationTypes(sim_seq)
    
    # Calculate probabilities of mutations at each position given targeting
    # from MutationProfiling.R; includes a N row
    
    # Columns corresponding to "N" and "." positions will have NA across all rows
    # These get converted to a probability of 0, ensuring that sampleMut() will
    # never choose these positions
    targeting <- calculateTargeting(germlineSeq = sim_seq, targetingModel = targetingModel) 
    # keep only ACGT rows
    targeting <- targeting[NUCLEOTIDES[1:4], ] 
    # set NA to 0
    targeting[is.na(targeting)] <- 0 
    # Make probability of stop codon 0
    targeting[mutation_types=="Stop"] <- 0
    
    # Initialize counters
    total_muts <- 0
    positions <- numeric(numMutations)
    
    while(total_muts < numMutations) {
        # Get position to mutate and update counters
        mutpos <- sampleMut(sim_leng, targeting, positions)
        total_muts <- total_muts + 1
        positions[total_muts] <- mutpos$pos
        
        # Implement mutation in simulation sequence
        mut_nuc <- 4 - (4*mutpos$pos - mutpos$mut)
        stri_sub(str=sim_seq, from=mutpos$pos, to=mutpos$pos) <- NUCLEOTIDES[mut_nuc]
        
        # Update targeting
        lower <- max(mutpos$pos-4, 1)
        upper <- min(mutpos$pos+4, sim_leng)
        targeting[, lower:upper] <- calculateTargeting(germlineSeq = stri_sub(str = sim_seq, 
                                                                              from = lower, 
                                                                              to = upper),
                                                       targetingModel = targetingModel)[NUCLEOTIDES[1:4], ]
        targeting[is.na(targeting)] <- 0
        
        # Update possible mutations
        lower <- getCodonPos(lower)[1]
        upper <- getCodonPos(upper)[3]
        mutation_types[, lower:upper] <- computeMutationTypes(stri_sub(str = sim_seq, 
                                                                       from = lower, 
                                                                       to = upper))
        # Make probability of stop codon 0
        if(any(mutation_types[, lower:upper]=="Stop", na.rm=T)) {
            targeting[, lower:upper][mutation_types[, lower:upper]=="Stop"] <- 0
        }
    }
    # sanity check: length of sim_seq should remain unchanged after simulation
    stopifnot(sim_leng==stri_length(sim_seq))
    return(sim_seq)
}


#' Simulate mutations in a lineage tree
#'
#' \code{shmulateTree} returns a set of simulated sequences generated from an input 
#' sequence and a lineage tree. The input sequence is used to replace the most recent 
#' common ancestor (MRCA) node of the \code{igraph} object defining the lineage tree. 
#' Sequences are then simulated with mutations corresponding to edge weights in the tree. 
#' Sequences will not be generated for groups of nodes that are specified to be excluded.
#'
#' @param    sequence        string defining the MRCA sequence to seed mutations from.
#' @param    graph           \code{igraph} object defining the seed lineage tree, with 
#'                           vertex annotations, whose edges are to be recreated.
#' @param    targetingModel  5-mer \link{TargetingModel} object to be used for computing 
#'                           probabilities of mutations at each position. Defaults to
#'                           \link{HH_S5F}.
#' @param    field           annotation to use for both unweighted path length exclusion 
#'                           and consideration as the MRCA node. If \code{NULL} do not 
#'                           exclude any nodes.
#' @param    exclude         vector of annotation values in \code{field} to exclude from 
#'                           potential MRCA set. If \code{NULL} do not exclude any nodes.
#'                           Has no effect if \code{field=NULL}.
#' @param    junctionWeight  fraction of the nucleotide sequence that is within the 
#'                           junction region. When specified this adds a proportional 
#'                           number of mutations to the immediate offspring nodes of the 
#'                           MRCA. Requires a value between 0 and 1. If \code{NULL} then 
#'                           edge weights are unmodified from the input \code{graph}.
#'
#' @return   A \code{data.frame} of simulated sequences with columns:
#'           \itemize{
#'             \item \code{NAME}:      name of the corresponding node in the input 
#'                                     \code{graph}.  
#'             \item \code{SEQUENCE}:  mutated sequence.
#'             \item \code{DISTANCE}:  Hamming distance of the mutated sequence from 
#'                                     the seed \code{sequence}.
#'           }
#' 
#' @seealso  See \link{shmulateSeq} for imposing mutations on a single sequence. 
#'           See \link{HH_S5F} and \link{MK_RS5NF} for predefined 
#'           \link{TargetingModel} objects.
#' 
#' @examples
#' # Load example lineage and define example MRCA
#' data(ExampleTrees, package="alakazam")
#' graph <- ExampleTrees[[17]]
#' sequence <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"
#' 
#' # Simulate using the default human 5-mer targeting model
#' shmulateTree(sequence, graph)
#' 
#' # Simulate using the mouse 5-mer targeting model
#' # Exclude nodes without a sample identifier
#' # Add 20% mutation rate to the immediate offsprings of the MRCA
#' shmulateTree(sequence, graph, targetingModel=MK_RS5NF,
#'              field="SAMPLE", exclude=NA, junctionWeight=0.2)
#'  
#' @export
shmulateTree <- function(sequence, graph, targetingModel=HH_S5F,
                         field=NULL, exclude=NULL, junctionWeight=NULL) {
    ## DEBUG
    # targetingModel=HH_S5F; field=NULL; exclude=NULL; junctionWeight=NULL
    
    # Check targeting model
    if (!is(targetingModel, "TargetingModel")) {
        stop(deparse(substitute(targetingModel)), " is not a valid TargetingModel object")
    }
    
    # Determine MRCA of lineage tree
    mrca_df <- alakazam::getMRCA(graph, path="distance", root="Germline", 
                                 field=field, exclude=exclude)
    
    # Get adjacency matrix
    adj <- as_adjacency_matrix(graph, attr="weight", sparse=FALSE)
    # Get names of nodes for which sequences are not to be returned
    skip_names <- c()
    if (!is.null(field)) {
        g <- vertex_attr(graph, name=field)
        g_names <- vertex_attr(graph, name="name")
        skip_names <- g_names[g %in% exclude]
    }
    
    # Create data.frame to hold simulated sequences
    # this will include a row for Germline
    sim_tree <- data.frame(matrix(NA, ncol=3, nrow=length(V(graph)),
                           dimnames=list(NULL, c("NAME", "SEQUENCE", "DISTANCE"))))
    sim_tree$NAME <- vertex_attr(graph, name="name")
    
    # remove row for Germline
    sim_tree <- sim_tree[-which(sim_tree$NAME=="Germline"), ]
        
    parent_nodes <- mrca_df$NAME[1]
    nchild <- sum(adj[parent_nodes, ] > 0)
    
    sim_tree$SEQUENCE[which(sim_tree$NAME==parent_nodes)] <- sequence
    sim_tree$DISTANCE[which(sim_tree$NAME==parent_nodes)] <- 0
    
    # Add mutations to the immediate offsprings of the MRCA
    # Number of mutations added is proportional to fraction of sequence in junction
    if(!is.null(junctionWeight)) {
        adj[parent_nodes, ] <- round(adj[parent_nodes, ] * (1 + junctionWeight))
    }
    
    while(nchild > 0) {
        new_parents <- c()
        # Loop through parent-children combos
        for(p in parent_nodes) {
            children <- colnames(adj)[adj[p, ] > 0]
            for(ch in children) {
                # Add child to new parents
                new_parents <- union(new_parents, ch)
                # Simulate sequence for that edge
                seq <- shmulateSeq(sequence=sim_tree$SEQUENCE[sim_tree$NAME == p], 
                                   numMutations=adj[p, ch],
                                   targetingModel=targetingModel)
                # Update output data.frame
                chRowIdx = which(sim_tree$NAME==ch)
                sim_tree$SEQUENCE[chRowIdx] <- seq
                sim_tree$DISTANCE[chRowIdx] <- adj[p, ch]
            }
        }
        
        # Re-calculate number of children
        parent_nodes <- new_parents
        nchild <- sum(adj[parent_nodes, ] > 0)
    }
    
    # Remove sequences that are to be excluded
    sim_tree <- sim_tree[!(sim_tree$NAME %in% skip_names), ]
    # Remove NAs
    # e.g. if node B is an offspring of node A, and node A has been excluded
    # then node B will have $SEQUENCE and $DISTANCE of NAs
    sim_tree <- sim_tree[!is.na(sim_tree$SEQUENCE), ]
    
    return(sim_tree)
}


#### Helper functions ####

# Compute the mutations types
#
# For each position in the input sequence, use \code{CODON_TABLE} to
# determine what types of mutations are possible. Returns \code{matrix}
# of all possible mutations and corresponding types.
#
# @param   inputSeq   sequence for which to compute mutation types
# @return  A \code{matrix} of mutation types for each position in the sequence.
computeMutationTypes <- function(inputSeq){
    #* counts on constant variable CODON_TABLE, NUCLEOTIDES (ACTGN-.)
    #* caution: this breaks down if length of seq is not a multiple of 3
    
    leng_seq <- stri_length(inputSeq)
    try(if( (leng_seq %%3 !=0) ) stop("length of input sequence must be a multiple of 3"))
    
    codons <- sapply(seq(1, leng_seq, by=3), function(x) {substr(inputSeq,x,x+2)})
    mut_types <- matrix(unlist(CODON_TABLE[, codons]), ncol=leng_seq, nrow=4, byrow=F)
    dimnames(mut_types) <-  list(NUCLEOTIDES[1:4], 1:leng_seq)
    return(mut_types)
}


# Pick a position to mutate
#
# Sample positions in the sequence to mutate given targeting probability
# until a new position is selected. This new position is then added to the
# vector of mutated positions and returned.
#
# @param   sim_leng   length of sequence in which mutation is being simulated
# @param   targeting  probabilities of each position in the sequence being mutated
# @param   positions  vector of positions which have already been mutated
#
# @return   A \code{list} of position being mutated and updated vector of mutated positions.
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
