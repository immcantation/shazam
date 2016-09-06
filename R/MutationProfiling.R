# Mutation profiling

#' @include Shazam.R
NULL

#### Clonal consensus building functions ####

#' Constructs effective clonal sequences
#'
#' \code{collapseClones} creates an effective sequence for each clonal 
#' group and appends columns to the input \code{data.frame} containing the effective 
#' sequence and germline for each clone.
#'
#' @param   db                  \code{data.frame} containing sequence data.
#' @param   cloneColumn         \code{character} name of the column containing clonal 
#'                              identifiers.
#' @param   sequenceColumn      \code{character} name of the column containing input 
#'                              sequences.
#' @param   germlineColumn      \code{character} name of the column containing germline 
#'                              sequences.
#' @param   expandedDb          \code{logical} indicating whether or not to return the 
#'                              expanded \code{db}, containing all the sequences (as opposed
#'                              to returning just one sequence per clone collapsed by )
#' @param   regionDefinition    \link{RegionDefinition} object defining the regions
#'                              and boundaries of the Ig sequences.
#' @param   nonTerminalOnly     \code{logical} indicating whether to include mutations
#'                              at the leaves.
#' @param   nproc               Number of cores to distribute the operation over. If the 
#'                              \code{cluster} has already been set earlier, then pass the 
#'                              \code{cluster}. This will ensure that it is not reset.
#'                              
#' 
#' @return   A modified \code{db} with clonal consensus sequences added 
#'           in the following columns:
#'           \itemize{
#'             \item \code{CLONAL_SEQUENCE}:  effective sequence for the clone.
#'             \item \code{CLONAL_GERMLINE}:  germline sequence for the clone.
#'                                            Generally, this will be unchanged from
#'                                            the data in \code{germlineColumn}, but
#'                                            may be truncated when the input sequence
#'                                            is truncated due to inconsistencies 
#'                                            in the lengths of the input sequences or
#'                                            \code{regionDefinition} limits.
#'           }
#'
#' @details
#' 
#' For sequences identified to be part of the same clone, an effective clonal sequence, 
#' representative of all mutations in a clone, is constructed using a stochastic approach. 
#' Each position in th effective sequence is determined by a weighted sampling 
#' of each mutated non-ambiguous base (excluding "N", "." or "-" characters) from all 
#' the sequences in the clone. For example, in a clone with 5 sequences that have "C" 
#' at position 1, and 5 sequences with "T" at this same position, the effective sequence 
#' will have a "C" 50\% and "T" 50\% of the time it is called.
#' 
#' Non-terminal branch mutations are defined as the set of mutations that occur on 
#' branches of the lineage tree that are not connected to a leaf. For computational 
#' efficiency, the set of non-terminal branch mutations is approximated as those that are
#' shared between more than one sequence in a clone. In this case the terminal branch 
#' mutations are filtered out.
#' 
#' @seealso
#' See \link{IMGT_SCHEMES} for a set of predefined \link{RegionDefinition} objects.
#' 
#' @examples
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
#' 
#' # Build clonal consensus for the full sequence
#' clones <- collapseClones(db, nproc=1)
#'                          
#' # Build clonal consensus for V-region only 
#' # Return the same number of rows as the input
#' clones <- collapseClones(db, regionDefinition=IMGT_V_NO_CDR3, 
#'                          expandedDb=TRUE, nproc=1)
#' 
#' @export
collapseClones <- function(db, 
                           cloneColumn="CLONE", 
                           sequenceColumn="SEQUENCE_IMGT",
                           germlineColumn="GERMLINE_IMGT_D_MASK",
                           expandedDb=FALSE,
                           regionDefinition=NULL,
                           nonTerminalOnly=FALSE,
                           nproc=1) {
    ## DEBUG
    # cloneColumn="CLONE"; sequenceColumn="SEQUENCE_IMGT"; germlineColumn="GERMLINE_IMGT_D_MASK"
    # expandedDb=FALSE; regionDefinition=NULL; nonTerminalOnly=FALSE; nproc=1
    
    # Hack for visibility of data.table and foreach index variables
    idx <- yidx <- .I <- NULL

    # Check for valid columns
    check <- checkColumns(db, c(cloneColumn, sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    
    if (class(expandedDb) != "logical") {
        stop ("expandedDb must be TRUE or FALSE.")
    }
    
    if (class(nonTerminalOnly) != "logical") {
        stop ("nonTerminalOnly must be TRUE or FALSE.")
    }
    
    # Convert clone identiers to strings
    db[[cloneColumn]] <- as.character(db[[cloneColumn]])
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of CLONES, instead of doing dplyr
    dt <- data.table(db)
    setkeyv(dt, cloneColumn)
    # Get the group indexes
    #dt <- dt[, list(yidx=list(.I)), by=list(CLONE)]
    dt <- dt[, list(yidx=list(.I)), by=cloneColumn]
    groups <- dt[, yidx]
    lenGroups <- length(groups)
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else {
        if (nproc != 0) { 
            #cluster <- makeCluster(nproc, type="SOCK") 
            cluster <- parallel::makeCluster(nproc, type= "PSOCK")
        }
        parallel::clusterExport(cluster, 
                                list('db', 'sequenceColumn', 'germlineColumn', 'cloneColumn',
                                     'regionDefinition', 'calcClonalConsensus',
                                     'groups', 'c2s', 's2c', 'words', 'translate'), 
                                envir=environment() )
        registerDoParallel(cluster)
    }
    
    # Printing status to console
    cat("Collapsing clonal sequences...\n")
    
    #cons_list <- foreach(idx=iterators::icount(lenGroups), .combine=c, .verbose=FALSE, 
    cons_mat <- foreach(idx=iterators::icount(lenGroups), .combine="cbind", 
                        .verbose=FALSE, .errorhandling='pass') %dopar% {
        calcClonalConsensus(inputSeq=db[[sequenceColumn]][groups[[idx]]],
                            germlineSeq=db[[germlineColumn]][groups[[idx]]],
                            regionDefinition=regionDefinition, 
                            nonTerminalOnly=nonTerminalOnly)
    }

    # Stop cluster
    if(nproc > 1) { parallel::stopCluster(cluster) }
    
    # Build return data.frame
    if (expandedDb) { 
        # Fill all rows with the consensus sequence
        clone_id <-  unique(db[[cloneColumn]])
        clone_index <- match(db[[cloneColumn]], clone_id)
        cons_db <- db
        cons_db$CLONAL_SEQUENCE <- cons_mat[1, clone_index]
        cons_db$CLONAL_GERMLINE <- cons_mat[2, clone_index]
    } else {
        # Return only the first low of each clone
        clone_id <-  unique(db[[cloneColumn]])
        clone_index <- match(clone_id, db[[cloneColumn]])
        cons_db <- db[clone_index, ]
        cons_db$CLONAL_SEQUENCE <- cons_mat[1, ]
        cons_db$CLONAL_GERMLINE <- cons_mat[2, ]
    }
    
    return(cons_db)
}



# Helper function for calcDBClonalConsensus
# 
# @param   inputSeq         character vector of observed sequences
# @param   germlineSeq         character vector of germline sequences
# @param   regionDefinition    \link{RegionDefinition} object defining the regions
#                              and boundaries of the Ig sequences.
# @param   nonTerminalOnly     \code{logical} indicating whether to include mutations
#                              at the leaves.
# @return  A named vector length two with "input" and "germline" consensus sequences
calcClonalConsensus <- function(inputSeq, germlineSeq, regionDefinition=NULL, 
                                nonTerminalOnly=FALSE) {
    ## DEBUG
    # inputSeq=db$SEQUENCE_IMGT[4:6]; germlineSeq=db$GERMLINE_IMGT_D_MASK[4:6]; regionDefinition=NULL; nonTerminalOnly=FALSE
    
    # If only one sequence in clone, return it
    if (length(inputSeq) == 1) {
        returnSeq <- c("input"=inputSeq, "germline"=germlineSeq)
        return(returnSeq)
    }
    
    # Find length of shortest input sequence
    # This is used to trim all the sequences to that length
    # or, if a regionDefinition is passed, then only analyze till the end of the defined length
    len_inputSeq <- stri_length(inputSeq)
    len_shortest <- min(len_inputSeq, na.rm=TRUE)
    if(!is.null(regionDefinition)) {
        len_shortest <- min(len_shortest, regionDefinition@seqLength, na.rm=TRUE)
    }
    
    if (class(nonTerminalOnly) != "logical") {
        stop ("nonTerminalOnly must be TRUE or FALSE.")
    }
    
    #Find the length of the longest germline sequence
    len_germlineSeq <- stri_length(germlineSeq)
    germlineSeq <- germlineSeq[(which.max(len_germlineSeq))]
    
    # Identify the consensus sequence
    inputSeq <- unique(inputSeq)
    charInputSeqs <- sapply(inputSeq, function(x) { s2c(x)[1:len_shortest] })
    charGLSeq <- s2c(germlineSeq)
    matClone <- sapply(1:len_shortest, function(i) {
        # Identify the nucleotides (in seqs and germline) at the current position
        posNucs = unique(charInputSeqs[i,])
        posGL = charGLSeq[i]
        error = FALSE
        
        # If the current position is a gap in both germline and the sequence,
        # return a gap
        if(posGL %in% c("-", ".") & sum(!(posNucs%in%c("-", ".", "N", "n")))==0 ){
            return(c(".", error))
        }
        
        # If all the sequences in the clone have the same nucleotide at the current
        # position, return the value at the current positions
        if(length(posNucs)==1) {
            return(c(posNucs[1], error))
        } else {         
            if("N"%in%posNucs) { error=TRUE }
            
            # If the current nucleotide matches germline, return germline 
            if(sum(!posNucs[!posNucs%in%c("N", "n")] %in% posGL) == 0) {
                return(c(posGL, error))
            }else{
                #return( c(sample(posNucs[posNucs!="N"],1),error) )
                
                # If we look at all nodes (including terminal nodes), sample a nucleotide from the possible
                # nucleotides in the clonal sequences at this position
                if(!nonTerminalOnly){
                    return( c(sample(charInputSeqs[i, !charInputSeqs[i, ] %in% c("N", "n") & charInputSeqs[i, ] != posGL], 1), error))
                }else{
                    
                    # If we look at only non-terminal nodes, we only sample the nucleotides that appear more 
                    # than once (this is a quick approximation)
                    posNucs = charInputSeqs[i,!charInputSeqs[i,]%in% c("N", "n") & charInputSeqs[i,]!=posGL]
                    posNucsTable = table(posNucs)
                    if(sum(posNucsTable > 1) == 0) {
                        return(c(posGL,error))
                    } else {
                        return(c(sample(posNucs[posNucs %in% names(posNucsTable)[posNucsTable > 1]], 1), error))
                    }
                }
                
            }
        }
        if (error==TRUE) { warning("Error while attempting to collapse by clone!") }
    })
    
    returnSeq <- c("input"=stri_join(matClone[1,], collapse=""), 
                   "germline"=stri_join(charGLSeq[1:len_shortest], collapse=""))
    
    return (returnSeq)
}



#### Mutation counting functions ####

#' Calculate observed numbers of mutations
#'
#' \code{observedMutations} calculates the observed number of mutations for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      \code{character} name of the column containing input 
#'                               sequences.
#' @param    germlineColumn      \code{character} name of the column containing 
#'                               the germline or reference sequence.
#' @param    frequency           \code{logical} indicating whether or not to calculate
#'                               mutation frequencies. Default is \code{FALSE}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences. If NULL, mutations 
#'                               are counted for entire sequence.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    nproc               number of cores to distribute the operation over. If the 
#'                               cluster has already been set the call function with 
#'                               \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                               \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with observed mutation counts for each 
#'           sequence listed. The columns names are dynamically created based on the
#'           regions in the \code{regionDefinition}. For example, when using the
#'           \link{IMGT_V_NO_CDR3} definition, which defines positions for CDR and
#'           FWR, the following columns are added:
#'           \itemize{
#'             \item  \code{OBSERVED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{OBSERVED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{OBSERVED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{OBSERVED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           
#' @details
#' Mutation count are determined by comparing the input sequences (in the column specified 
#' by \code{sequenceColumn}) to the germline sequence (in the column specified by 
#' \code{germlineColumn}). 
#' 
#' The mutations are binned as either replacement (R) or silent (S) across the different 
#' regions of the sequences as defined by \code{regionDefinition}. Typically, this would 
#' be the framework (FWR) and complementarity determining (CDR) regions of IMGT-gapped 
#' nucleotide sequences. Mutation counts are appended to the input \code{db} as 
#' additional columns.
#' 
#' @seealso  
#' \link{calcObservedMutations} is called by this function to get the number of mutations 
#' in each sequence grouped by the \link{RegionDefinition}. 
#' See \link{IMGT_SCHEMES} for a set of predefined \link{RegionDefinition} objects.
#' See \link{expectedMutations} for calculating expected mutation frequencies.
#'           
#' 
#' @examples
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
#'
#' # Calculate mutation frequency over the entire sequence
#' db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
#'                             germlineColumn="GERMLINE_IMGT_D_MASK",
#'                             frequency=TRUE,
#'                             nproc=1)
#'
#' # Count of V-region mutations split by FWR and CDR
#' # With mutations only considered replacement if charge changes
#' db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
#'                             germlineColumn="GERMLINE_IMGT_D_MASK",
#'                             regionDefinition=IMGT_V_NO_CDR3,
#'                             mutationDefinition=CHARGE_MUTATIONS,
#'                             nproc=1)
#'                      
#' @export
observedMutations <- function(db, 
                              sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              frequency=FALSE,
                              regionDefinition=NULL,
                              mutationDefinition=NULL,
                              nproc=1) {
    # Hack for visibility of data.table and foreach index variables
    idx <- NULL
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if (nproc > 1) {        
        cluster <- parallel::makeCluster(nproc, type = "PSOCK")
        parallel::clusterExport(cluster, list('db', 'sequenceColumn', 'germlineColumn', 
                                              'regionDefinition', 'frequency',
                                              'calcObservedMutations','s2c','c2s','NUCLEOTIDES',
                                              'getCodonPos','getContextInCodon','mutationType',
                                              'translateCodonToAminoAcid','AMINO_ACIDS','binMutationsByRegion',
                                              'collapseMatrixToVector'), 
                                envir=environment())
        registerDoParallel(cluster)
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Printing status to console
    cat("Calculating observed number of mutations...\n")
    
    # Identify all the mutations in the sequences
    numbOfSeqs <- nrow(db)
    observedMutations_list <-
        foreach(idx=iterators::icount(numbOfSeqs)) %dopar% {
            calcObservedMutations(db[idx, sequenceColumn], 
                                  db[idx, germlineColumn],
                                  frequency=frequency,
                                  regionDefinition=regionDefinition,
                                  mutationDefinition=mutationDefinition)
        }
    
    # Convert list of mutations to data.frame
    if (!is.null(regionDefinition)) {
        labels_length <- length(regionDefinition@labels)
    } else{
        #labels_length=1
        labels_length <- length(makeNullRegionDefinition()@labels)
    }
    observed_mutations <- do.call( rbind, lapply(observedMutations_list, function(x) { 
        length(x) <- labels_length 
        return(x)
    }))
    
    
    sep <- "_"
    if (ncol(observed_mutations) > 1) sep <- "_"
    observed_mutations[is.na(observed_mutations)] <- 0
    if (frequency == TRUE) {
        colnames(observed_mutations) <- paste("MU_FREQ", colnames(observed_mutations), sep=sep)
    } else {
        colnames(observed_mutations) <- paste("OBSERVED", colnames(observed_mutations), sep=sep)
    }
    
    # Properly shutting down the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, observed_mutations)
    return(db_new)    
}


#' Count the number of observed mutations in a sequence.
#'
#' \code{calcObservedMutations} determines all the mutations in a given input seqeunce compared
#' to its germline sequence.
#'
#' @param    inputSeq            input sequence.
#' @param    germlineSeq         germline sequence.
#' @param    frequency           \code{logical} indicating whether or not to calculate
#'                               mutation frequencies. The denominator used is the number of bases
#'                               that are non-N in both the input and the germline sequences.
#'                               Default is \code{FALSE}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences. Note, only the part of
#'                               sequences defined in \code{regionDefinition} are analyzed.
#'                               If NULL, mutations are counted for entire sequence.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    returnRaw           return the positions of point mutations and their corresponding
#'                               mutation types (as opposed to counts of mutations; hence "raw"); 
#'                               as well as the number of non-N bases used as the denominator when
#'                               calculating frequency. Default if \code{FALSE}.                               
#'                               
#' @return   For \code{returnRaw=FALSE}, an \code{array} with the number of replacement (R) 
#'           and silent(S) mutations. For \code{returnRaw=TRUE}, a list containing a data 
#'           frame (\code{$pos}) whose columns (\code{position}, \code{type}, and \code{region}) 
#'           indicate the position, mutation type (R or S), and region of each mutation; and a 
#'           vector (\code{$nonN}) indicating the number of non-N bases in regions defined by
#'           \code{regionDefinition}.
#'           
#' @details
#' Each mutation is considered independently in the germline context. Note, only the part of 
#' \code{inputSeq} defined in \code{regionDefinition} is analyzed. For example, when using 
#' the default \link{IMGT_V_NO_CDR3} definition, then mutations in positions beyond 
#' 312 will be ignored. 
#' 
#' Note that only replacement (R) and silent (S) mutations are included in the 
#' results. Stop mutations and mutations such as the case in which NNN in the germline
#' sequence is observed as NNC in the input sequence are excluded. In other words,
#' a result that is \code{NA} or zero indicates absence of R and S mutations, not 
#' necessarily all types of mutations, such as the excluded ones mentioned above.
#' 
#' @seealso  See \link{observedMutations} for counting the number of observed mutations 
#' in a \code{data.frame}.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' in_seq <- ExampleDb[100, "SEQUENCE_IMGT"]
#' germ_seq <-  ExampleDb[100, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Identify all mutations in the sequence
#' ex1_raw = calcObservedMutations(in_seq, germ_seq, returnRaw=T)
#' # Count all mutations in the sequence
#' ex1_count = calcObservedMutations(in_seq, germ_seq, returnRaw=F)
#' ex1_freq = calcObservedMutations(in_seq, germ_seq, returnRaw=F, frequency=T)
#' # Compare this with ex1_count
#' table(ex1_raw$pos$region, ex1_raw$pos$type)
#' # Compare this with ex1_freq
#' table(ex1_raw$pos$region, ex1_raw$pos$type) / ex1_raw$nonN
#' 
#' # Identify only mutations the V segment minus CDR3
#' ex2_raw = calcObservedMutations(in_seq, germ_seq, 
#'                                 regionDefinition=IMGT_V_NO_CDR3, returnRaw=T)
#' # Count only mutations the V segment minus CDR3
#' ex2_count = calcObservedMutations(in_seq, germ_seq, 
#'                                   regionDefinition=IMGT_V_NO_CDR3, returnRaw=F)
#' ex2_freq = calcObservedMutations(in_seq, germ_seq, 
#'                                  regionDefinition=IMGT_V_NO_CDR3, returnRaw=F,
#'                                  frequency=T)
#' # Compare this with ex2_count
#' table(ex2_raw$pos$region, ex2_raw$pos$type)                                 
#' # Compare this with ex2_freq
#' table(ex2_raw$pos$region, ex2_raw$pos$type) / ex2_raw$nonN                                        
#' 
#' # Identify mutations by change in hydropathy class
#' ex3_raw = calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V_NO_CDR3,
#'                                 mutationDefinition=HYDROPATHY_MUTATIONS, returnRaw=T)
#' # Count mutations by change in hydropathy class
#' ex3_count = calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V_NO_CDR3,
#'                                   mutationDefinition=HYDROPATHY_MUTATIONS, returnRaw=F)
#' ex3_freq = calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V_NO_CDR3,
#'                                  mutationDefinition=HYDROPATHY_MUTATIONS, returnRaw=F, 
#'                                  frequency=T)
#' # Compre this with ex3_count
#' table(ex3_raw$pos$region, ex3_raw$pos$type)                                        
#' # Compare this with ex3_freq
#' table(ex3_raw$pos$region, ex3_raw$pos$type) / ex3_raw$nonN                                        
#'                                 
#' @export
calcObservedMutations <- function(inputSeq, germlineSeq, frequency=FALSE,
                                  regionDefinition=NULL, mutationDefinition=NULL,
                                  returnRaw=FALSE) {
    # Assign mutation definition
    aminoAcidClasses <- if (is.null(mutationDefinition)) { NULL } else { mutationDefinition@classes }
        
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
    germlineSeq <- gsub("\\.\\.\\.", "XXX", germlineSeq)
    #If there is a single gap left convert it to an N
    germlineSeq <- gsub("\\.", "N", germlineSeq)
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    germlineSeq <- gsub("XXX", "...", germlineSeq)
    
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
    inputSeq <- gsub("\\.\\.\\.", "XXX", inputSeq)
    #If there is a single gap left convert it to an N
    inputSeq <- gsub("\\.", "N", inputSeq)
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    inputSeq <- gsub("XXX", "...", inputSeq)    
    
    # Trim the input and germline sequence to the shortest
    len_inputSeq <- nchar(inputSeq)
    len_germlineSeq <- nchar(germlineSeq)
    
    # If a regionDefinition is passed,
    # then only analyze till the end of the defined length
    if(!is.null(regionDefinition)) {
        rdLength  <- regionDefinition@seqLength
    } else {
        rdLength <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        # Create full sequence RegionDefinition object
        regionDefinition <- makeNullRegionDefinition(rdLength)
    }
    len_shortest <- min(c(len_inputSeq, len_germlineSeq, rdLength), na.rm=TRUE)
    
    c_inputSeq <- s2c(inputSeq)[1:len_shortest]
    c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
    
    # If the sequence and germline (which now should be the same length) is shorter
    # than the rdLength, pad it with Ns
    if(len_shortest<rdLength){
        fillWithNs <- array("N",rdLength-len_shortest)
        c_inputSeq <- c( c_inputSeq, fillWithNs)
        c_germlineSeq <- c( c_germlineSeq, fillWithNs)
    }
    
    mutations_array_raw <- NA
    mutations_array <- NA
    mutations = (c_germlineSeq != c_inputSeq) & (c_germlineSeq%in%NUCLEOTIDES[1:5]) & (c_inputSeq%in%NUCLEOTIDES[1:5])
    if (sum(mutations) > 0) {
        # The nucleotide positions of the mutations
        mutations_pos <- which(mutations==TRUE)
        # For every mutations_pos, extract the entire codon from germline
        mutations_pos_codons <- array(sapply(mutations_pos,getCodonPos))
        c_germlineSeq_codons <- c_germlineSeq[mutations_pos_codons]
        # For every mutations_pos, extract the codon from germline (without other mutations 
        # at the same codon, if any).
        c_inputSeq_codons <- array(sapply(mutations_pos, function(x) {
            seqP = c_germlineSeq[getCodonPos(x)]
            seqP[getContextInCodon(x)] = c_inputSeq[x]
            return(seqP) }))
        # split the string of codons into vector of codons
        c_germlineSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_germlineSeq_codons)), " ")[[1]]
        c_inputSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_inputSeq_codons)), " ")[[1]]
        
        # Determine whether the mutations are R or S
        mutations_array_raw <- apply(rbind(c_germlineSeq_codons, c_inputSeq_codons), 2, 
                                     function(x) { mutationType(c2s(x[1]), c2s(x[2]), aminoAcidClasses=aminoAcidClasses) })
        names(mutations_array_raw) = mutations_pos
        # remove NAs (arisen from cases such as NNN [germline] and NNC [input])
        mutations_array_raw <- mutations_array_raw[!is.na(mutations_array_raw)]
        
        if(length(mutations_array_raw)==0){
            # NA if mutations_array_raw contains all NAs and they have all been removed
            mutations_array_raw <- NA
            mutations_array <- NA    
        } else {
            # mutation types other than "R" and "S" (e.g. "Stop") will be removed by binMutationsByRegion,
            # and stored in mutations_array
            mutations_array <- binMutationsByRegion(mutations_array_raw, regionDefinition)
            if (frequency==TRUE) {
                tempNames <- sapply(regionDefinition@labels, function(x) { substr(x, 1, nchar(x)-2) })
                # Subset boundaries to only non N bases (in both seq and gl)
                boundaries <- regionDefinition@boundaries[c_inputSeq%in%NUCLEOTIDES[1:4] &  
                                                          c_germlineSeq%in%NUCLEOTIDES[1:4]]
                # Freq = numb of mutations / numb of non N bases (in both seq and gl)
                denoms <- sapply(tempNames, function(x) { sum(boundaries==x) })
                mutations_array <- mutations_array/denoms
            }
            # for consistency, manually remove non-"R"/"S" from mutations_array_raw
            # i.e. not counting mutations such as "Stop"
            mutations_array_raw = mutations_array_raw[mutations_array_raw %in% c("R", "S")]
        }        
    }
    
    # return positions of point mutations and their mutation types ("raw")
    if (returnRaw){
      # number of non-N bases (in both seq and gl)
      nonN.regions <- unique(sapply(regionDefinition@labels, function(x) { substr(x, 1, nchar(x)-2) }))
      nonN.boundaries <- regionDefinition@boundaries[c_inputSeq%in%NUCLEOTIDES[1:4] &  
                                                       c_germlineSeq%in%NUCLEOTIDES[1:4]]
      nonN.denoms <- sapply(nonN.regions, function(x) { sum(nonN.boundaries==x) })
      
      if (length(mutations_array_raw) == sum(is.na(mutations_array_raw))) {
        # if mutations_array_raw is NA, or 
        # if mutations_array_raw is empty due to all mutations being "Stop" and hence removed
        # avoid is.na(mutations_array_raw) to avoid warning in case mutations_array_raw is a vector
        return(list(pos=mutations_array_raw, nonN=nonN.denoms))
      } else {
        # df indicating position, mutation type (R or S), and region of each mutation
        rawDf = data.frame(as.numeric(names(mutations_array_raw)))
        rawDf = cbind(rawDf,
                      as.character(mutations_array_raw), # as.character to remove names of the vector
                      as.character(regionDefinition@boundaries[as.numeric(names(mutations_array_raw))]),
                      stringsAsFactors=F)
        colnames(rawDf) = c("position", "type", "region")
      return(list(pos=rawDf, nonN=nonN.denoms))
      }
    } else {
    # return counts of each mutation type  
      return(mutations_array)
    }
}


# Aggregate mutations by region
#
# \code{binMutationsByRegion} takes an array of observed mutations (e.g. from 
# \code{calcObservedMutations}) and bins them by the different regions defined in the 
# \code{regionDefinition}.
#
# @param   mutationsArray     \code{array} containing the mutations (R/S) with the names
#                             indicating the nucleotide positions of the mutations.                             
# @param   regionDefinition   \link{RegionDefinition} object defining the regions
#                             and boundaries of the Ig sequences.
# 
# @return An \code{array} of R/S mutations binned across all the unique regions defined
# by \code{regionDefinition}.
# 
# @details
# Note, only the part of sequences defined in \code{regionDefinition} are analyzed.
# For example, if the default \link{IMGT_V_NO_CDR3} definition is used, then mutations
# in positions beyond 312 will be ignored.
# 
# @seealso  
# See \link{observedMutations} for identifying and counting the 
# number of observed mutations.
# This function is also used in \link{calcObservedMutations}.
# 
# @examples
# # Generate a random mutation array
# numbOfMutations <- sample(3:10, 1) 
# posOfMutations <- sort(sample(330, numbOfMutations))
# mutation_types <- sample(c("R","S"), length(posOfMutations), replace=TRUE)
# mutations_array <- array(mutation_types, dimnames=list(posOfMutations))
# 
# # Random mutations
# binMutationsByRegion(mutations_array, regionDefinition=NULL)
# binMutationsByRegion(mutations_array, regionDefinition=IMGT_V_NO_CDR3)
binMutationsByRegion <- function(mutationsArray, 
                                 regionDefinition=NULL) {
    # Create full sequence RegionDefinition object 
    # The seqLength will be the largest index of a mutation
    if (is.null(regionDefinition)) {
        regionDefinition <- makeNullRegionDefinition(max(as.numeric(names(mutationsArray))))
    }
    
    # Make a factor of R/S
    mutatedPositions <- as.numeric(names(mutationsArray))
    mutations <- array(NA,  dim=regionDefinition@seqLength)
    mutations[mutatedPositions] <- mutationsArray
    mutations <- mutations[1:regionDefinition@seqLength]
    mutations <- factor(mutations,levels=c("R", "S"))
    
    mutations_region_counts <- collapseMatrixToVector(table(regionDefinition@boundaries, mutations))
    
    sortingOrder <- match(regionDefinition@labels, names(mutations_region_counts))
    mutations_region_counts <- mutations_region_counts[sortingOrder]
    return(mutations_region_counts)
}


#### Sliding window approach ####
#' Sliding window approach towards filtering a single sequence
#'
#' \code{slideWindowSeq} determines whether an input sequence contains equal to or more than 
#' a given number of mutations in a given length of consecutive nucleotides (a "window") 
#' when compared to a germline sequence.
#' 
#' @param    inputSeq            input sequence.
#' @param    germlineSeq         germline sequence.
#' @param    mutThresh           threshold on the number of mutations in \code{windowSize} 
#'                               consecutive nucleotides. Must be between 1 and \code{windowSize} 
#'                               inclusive. 
#' @param    windowSize          length of consecutive nucleotides. Must be at least 2.
#'                               
#' @return  \code{TRUE} if there are equal to or more than \code{mutThresh} number of mutations
#'          in any window of \code{windowSize} consecutive nucleotides (i.e. the sequence should
#'          be filtered); \code{FALSE} if otherwise.
#' 
#' @seealso  \link{calcObservedMutations} is called by \code{slideWindowSeq} to identify observed 
#'           mutations. See \link{slideWindowDb} for applying the sliding window approach on a 
#'           \code{data.frame}. See \link{slideWindowTune} for parameter tuning for \code{mutThresh}
#'           and \code{windowSize}.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' in_seq <- ExampleDb[100, "SEQUENCE_IMGT"]
#' germ_seq <-  ExampleDb[100, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Determine if in_seq has 6 or more mutations in 10 consecutive nucleotides
#' slideWindowSeq(inputSeq=in_seq, germlineSeq=germ_seq, mutThresh=6, windowSize=10)
#'                                 
#' @export
slideWindowSeq = function(inputSeq, germlineSeq, mutThresh, windowSize){
  # identify all R and S mutations in input sequence
  inputMut = calcObservedMutations(inputSeq=inputSeq, germlineSeq=germlineSeq, returnRaw=T)$pos
  
  # extract positions of mutations
  # inputMut must either be NA (no observed mutation) or a df
  # avoid is.na (in case inputMut is a data frame then will get multiple T/F values and hence warning)
  if (!is.data.frame(inputMut)) {
    inputMutPos = NA
  } else {
    inputMutPos = inputMut$position
  }

  # call helper
  return(slideWindowSeqHelper(mutPos=inputMutPos, mutThresh=mutThresh, windowSize=windowSize))
}


# Helper for sliding window approach towards filtering sequences
#
# @param    mutPos              a numeric vector containing positions of point mutations as 
#                               returned in the \code{position} column by 
#                               \code{calcObserverdMutations()} with \code{returnRaw=TRUE}. 
#                               Can be \code{NA}, in which case the returned value will be 
#                               \code{FALSE}.
# @param    mutThresh           threshold on the number of mutations in \code{windowSize} 
#                               consecutive nucleotides. Must be between 1 and \code{windowSize} 
#                               inclusive.
# @param    windowSize          length of consecutive nucleotides. Must be at least 2.
#
# @return   \code{TRUE} if there are equal to or more than \code{mutThresh} number of mutations
#           in any window of \code{windowSize} consecutive nucleotides; \code{FALSE} if otherwise.
#
slideWindowSeqHelper = function(mutPos, mutThresh, windowSize){
  # check preconditions
  stopifnot(mutThresh>=1 & mutThresh<=windowSize & windowSize>=2)
  
  if (length(mutPos)==1 && is.na(mutPos)) {
    # use && instead of & to short-circuit in case length(mutPos)!=1 (otherwise warning)
    return(F)
  } else {
    # general idea:
    # only need to check windows containing mutations (as opposed to every possible window)
    for (i in 1:length(mutPos)){
      # get window limits
      lower = mutPos[i]
      upper = lower + windowSize - 1
      # how many mutations fall within current window
      windowCount = sum(mutPos>=lower & mutPos <= upper)
      # return as soon as a window has >= mutThresh mutations
      if (windowCount >= mutThresh) {return(T)}
    }
    return(F)
  }
}


#' Sliding window approach towards filtering sequences in a \code{data.frame}
#'
#' \code{slideWindowDb} determines whether each input sequence in a \code{data.frame} 
#' contains equal to or more than a given number of mutations in a given length of 
#' consecutive nucleotides (a "window") when compared to their respective germline 
#' sequence.
#' 
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    mutThresh           threshold on the number of mutations in \code{windowSize} 
#'                               consecutive nucleotides. Must be between 1 and \code{windowSize} 
#'                               inclusive. 
#' @param    windowSize          length of consecutive nucleotides. Must be at least 2.
#'                               
#' @return   a logical vector. The length of the vector matches the number of input sequences in 
#'           \code{db}. Each entry in the vector indicates whether the corresponding input sequence
#'           should be filtered based on the given parameters.
#' 
#' @seealso  See \link{slideWindowSeq} for applying the sliding window approach on a single sequence. 
#'           See \link{slideWindowTune} for parameter tuning for \code{mutThresh} and \code{windowSize}.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' 
#' # Apply the sliding window approach on a subset of ExampleDb
#' slideWindowDb(db = ExampleDb[1:10, ], mutThresh=6, windowSize=10)
#' 
#' @export
slideWindowDb = function(db, sequenceColumn = "SEQUENCE_IMGT", 
                         germlineColumn = "GERMLINE_IMGT_D_MASK",
                         mutThresh, windowSize){
  db.filter = sapply(1:nrow(db), function(i){slideWindowSeq(inputSeq = db[i, sequenceColumn],
                                                            germlineSeq = db[i, germlineColumn],
                                                            mutThresh = mutThresh,
                                                            windowSize = windowSize)})
  return(db.filter)
}


#' Parameter tuning for sliding window approach
#'
#' Apply \link{slideWindowDb} over a search grid made of combinations of \code{mutThresh} and 
#' \code{windowSize} to help with picking a pair of values for these parameters. Parameter 
#' tuning can be performed by choosing a combination that gives a reasonable number of 
#' filtered/remaining sequences. 
#' 
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    dbMutList           if supplied, this should be a list of \code{data.frame}s returned 
#'                               as \code{$pos} of the nested list produced by 
#'                               \link{calcObservedMutations} with \code{returnRaw=TRUE}; otherwise, 
#'                               \link{calcObservedMutations} is called on columns \code{sequenceColumn}
#'                               and \code{germlineColumn} of \code{db}. Default is \code{NULL}. 
#' @param    mutThreshRange      range of threshold on the number of mutations in \code{windowSize} 
#'                               consecutive nucleotides to try. Must be between 1 and 
#'                               maximum \code{windowSizeRange} inclusive. 
#' @param    windowSizeRange     range of length of consecutive nucleotides to try. The lower end
#'                               must be at least 2.
#' @param    verbose             whether to print out messages indicating current progress. Default
#'                               is \code{TRUE}.              
#'                               
#' @return   a list of logical matrices. Each matrix corresponds to a \code{windowSize} in 
#'           \code{windowSizeRange}. Each column in a matrix corresponds to a \code{mutThresh} in
#'           \code{mutThreshRange}.
#' 
#' @details  If, in a given combination of \code{mutThresh} and \code{windowSize}, \code{mutThresh} 
#'           is greater than \code{windowSize}, \code{NA}s will be returned for that particular
#'           combination. A message indicating that the combination has been "skipped" will be 
#'           printed if \code{verbose=TRUE}.
#'           
#'           If \link{calcObservedMutations} was previously run on \code{db} and saved, supplying
#'           \code{$pos} from the saved result as \code{dbMutList} could save time by skipping a
#'           second call of \link{calcObservedMutations}. This could be helpful especially when 
#'           \code{db} is large.
#' 
#' @seealso  \link{slideWindowDb} is called on \code{db} for tuning. See \link{slideWindowTunePlot} 
#'           for visualization. See \link{calcObservedMutations} for generating \code{dbMutList}.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' 
#' # Try out thresholds of 2-4 mutations in window sizes of 7-9 nucleotides 
#' # on a subset of ExampleDb
#' # In this case, all combinations are legal and thus computed
#' slideWindowTune(db = ExampleDb[1:5, ], mutThreshRange = 2:4, windowSizeRange = 7:9)
#' 
#' # In the following case, illegal combinations are skipped, returning NAs                                  
#' slideWindowTune(db = ExampleDb[1:5, ], mutThreshRange = 2:4, windowSizeRange = 2:4)
#'                                                             
#' # Run calcObservedMutations separately and skip calling it again in slideWindowTune
#' exDbMutList = sapply(1:5, 
#'                    function(i){
#'                      calcObservedMutations(inputSeq = ExampleDb[i, "SEQUENCE_IMGT"],
#'                                            germlineSeq = ExampleDb[i, "GERMLINE_IMGT_D_MASK"],
#'                                            returnRaw = TRUE)$pos})
#' slideWindowTune(db = ExampleDb[1:5, ], dbMutList = exDbMutList, 
#'                 mutThreshRange = 2:4, windowSizeRange = 2:4)                                                 
#'                                                            
#'                                                                                                                                                                                                                                                                      
#' @export
slideWindowTune = function(db, sequenceColumn = "SEQUENCE_IMGT", 
                           germlineColumn = "GERMLINE_IMGT_D_MASK",
                           dbMutList = NULL,
                           mutThreshRange, windowSizeRange, verbose=TRUE){
  # check preconditions
  stopifnot(!is.null(db))
  stopifnot( min(mutThreshRange)>=1 & 
             max(mutThreshRange)<=max(windowSizeRange) &
             min(windowSizeRange)>=2 )
  
  
  # get positions of R/S mutations for sequences in db
  # do this here and then call slideWindowSeqHelper (so it's done only once)
  # instead of calling slideWindowDb which does this every time it is called
  if (is.null(dbMutList)) {
    inputMutList = sapply(1:nrow(db), 
                          function(i){
                            calcObservedMutations(inputSeq=db[i, sequenceColumn],
                                                  germlineSeq=db[i, germlineColumn],
                                                  returnRaw=T)$pos})    
  } else {
    if (verbose) {cat("dbMutList supplied; skipped calling calcObservedMutations()\n")}
    inputMutList = dbMutList
  }
  
  inputMutPosList = lapply(inputMutList, 
                           function(x){
                             if (!is.data.frame(x)) {
                               return(NA) 
                             } else {
                               return(x$position)
                             }
                           })
    
  # apply slideWindow on combinations of windowSize and mutThresh
  for (size in windowSizeRange) {
    if (verbose) {cat(paste0("now computing for windowSize = ", size, "\n"))}
    
    for (thresh in mutThreshRange) {
      if (thresh <= size){
        if (verbose) {cat(paste0(">>> mutThresh = ", thresh, "\n"))}
        # apply slideWindow using current pair of parameters
        cur.logical = unlist(lapply(inputMutPosList,
                                    slideWindowSeqHelper,
                                    mutThresh = thresh, windowSize = size))
      } else {
        if (verbose) {cat(paste0(">>> mutThresh = ", thresh, " > windowSize = ", 
                                 size, " (skipped)\n"))}
        # NA if skipped
        cur.logical = rep(NA, nrow(db))
      }
      # store results for each thresh as a column in a logical matrix
      if (thresh == mutThreshRange[1]) {
        cur.mtx = matrix(data=cur.logical, nrow=length(cur.logical))
      } else {
        cur.mtx = cbind(cur.mtx, cur.logical)
      }
    }
    colnames(cur.mtx) = as.character(mutThreshRange)
    
    # store results for each size (and threshes under that size) as a logical matrix in a list
    if (size == windowSizeRange[1]) {
      cur.list = list(cur.mtx)
    } else {
      cur.list = c(cur.list, list(cur.mtx))
    }
  }
  names(cur.list) = as.character(windowSizeRange)
  
  return(cur.list)
}


#' Visualize parameter tuning for sliding window approach
#'
#' Visualize results from \link{slideWindowTune}
#' 
#' @param    tuneList            a list of logical matrices returned by \link{slideWindowTune}.
#' @param    plotFiltered        whether to plot the number of filtered sequences (as opposed to
#'                               the number of remaining sequences). Default is \code{TRUE}.
#' @param    percentage          whether to plot on the y-axis the percentage of filtered sequences
#'                               (as opposed to the absolute number). Default is \code{FALSE}.                             
#' @param    jitter.x            whether to jitter x-axis values. Default is \code{FALSE}.                               
#' @param    jitter.x.amt        amount of jittering to be applied on x-axis values if 
#'                               \code{jitter.x=TRUE}. Default is 0.1.
#' @param    jitter.y            whether to jitter y-axis values. Default is \code{FALSE}.
#' @param    jitter.y.amt        amount of jittering to be applied on y-axis values if 
#'                               \code{jitter.y=TRUE}. Default is 0.1.                               
#' @param    pchs                point types to pass on to \link{plot}.
#' @param    ltys                line types to pass on to \link{plot}.
#' @param    cols                colors to pass on to \link{plot}.                             
#' @param    plotLegend          whether to plot legend. Default is \code{TRUE}.
#' @param    legendPos           position of legend to pass on to \link{legend}. Default is
#'                               \code{"topright"}.
#' @param    legendHoriz         whether to make legend horizontal. Default is \code{FALSE}.
#' @param    legendCex           numeric values by which legend should be magnified relative to 1.
#' 
#' @details  For each \code{windowSize}, the numbers of sequences filtered or remaining after applying
#'           the sliding window approach are plotted on the y-axis against thresholds on the number of
#'           mutations in a window on the x-axis.
#'           
#'           When plotting, a user-defined \code{amount} of jittering can be applied on values plotted
#'           on either axis or both axes via adjusting \code{jitter.x}, \code{jitter.y}, 
#'           \code{jitter.x.amt} and \code{jitter.y.amt}. This may be help with visually distinguishing
#'           lines for different window sizes in case they are very close or identical to each other. 
#'           If plotting percentages (\code{percentage=TRUE}) and using jittering on the y-axis values 
#'           (\code{jitter.y=TRUE}), it is strongly recommended that \code{jitter.y.amt} be set very
#'           small (e.g. 0.01). 
#'           
#'           \code{NA} for a combination of \code{mutThresh} and \code{windowSize} where 
#'           \code{mutThresh} is greater than \code{windowSize} will not be plotted. 
#' 
#' @seealso  See \link{slideWindowTune} for how to get \code{tuneList}. See \link{jitter} for 
#'           use of \code{amount} of jittering.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' 
#' # Try out thresholds of 2-4 mutations in window sizes of 3-5 nucleotides 
#' # on a subset of ExampleDb
#' tuneList = slideWindowTune(db = ExampleDb[1:10, ], 
#'                            mutThreshRange = 2:4, windowSizeRange = 3:5,
#'                            verbose = FALSE)
#'
#' # Visualize
#' # Plot numbers of sequences filtered without jittering y-axis values
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
#'                     plotFiltered=T, jitter.y=F)
#'                     
#' # Notice that some of the lines overlap
#' # Jittering could help
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
#'                     plotFiltered=T, jitter.y=T)
#'                     
#' # Plot numbers of sequences remaining instead of filtered
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
#'                     plotFiltered=F, jitter.y=T, legendPos="bottomright")
#'                     
#' # Plot percentages of sequences filtered with a tiny amount of jittering
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
#'                     plotFiltered=T, percentage=T, 
#'                     jitter.y=T, jitter.y.amt=0.01)
#'                                                             
#' @export
slideWindowTunePlot = function(tuneList, plotFiltered = TRUE, percentage = FALSE,
                               jitter.x = FALSE, jitter.x.amt = 0.1,
                               jitter.y = FALSE, jitter.y.amt = 0.1,
                               pchs = 1, ltys = 2, cols = 1,
                               plotLegend = TRUE, legendPos = "topright", 
                               legendHoriz = FALSE, legendCex = 1){
  
  # invert (!) tuneList if plotting retained sequences
  ylab.part.2 = "filtered"
  if (!plotFiltered) {
    tuneList = lapply(tuneList, function(x){!x})
    ylab.part.2 = "remaining"}
  
  # if number of pchs/ltys/cols provided does not match number of lines expected
  # expand into vector with repeating values (otherwise legend would break)
  if (length(pchs)!=length(tuneList)) {pchs = rep(pchs, length(tuneList))}
  if (length(ltys)!=length(tuneList)) {ltys = rep(ltys, length(tuneList))}
  if (length(cols)!=length(tuneList)) {cols = rep(cols, length(tuneList))}
  
  # tabulate tuneList (and if applicable convert to percentage)
  plotList = lapply(tuneList, colSums)
  if (percentage) {plotList = lapply(plotList, function(x){x/nrow(tuneList[[1]])})}
  
  # get x-axis values (i.e. mutThreshRange; colnames of matrix in tuneList with most columns)
  #threshes = as.numeric(colnames(tuneList[[which.max(lapply(lapply(tuneList, colnames), length))]]))
  threshes = as.numeric(colnames(tuneList[[1]]))
  
  # plot for first window size
  x1 = threshes
  if (jitter.x) {x1 = jitter(x1, amount=jitter.x.amt)}
  y1 = plotList[[1]]
  if (jitter.y) {y1 = jitter(y1, amount=jitter.y.amt)}
  
  if (percentage) {
    ylab.part.1 = "Percentage of sequences"
    # ylim
    ylim.padding = abs(diff(range(plotList, na.rm=T)))*0.01
    ylims = c(max(0, min(range(plotList, na.rm=T)) - ylim.padding), 
              min(1, max(range(plotList, na.rm=T)) + ylim.padding) )

  } else {
    ylab.part.1 = "Number of sequences"
    # ylim: non-negative lower limit; upper limit slight above max tabulated sum
    ylims = c( max(0, min(range(plotList, na.rm=T)) - max(1, jitter.y.amt) ), 
               max(range(plotList, na.rm=T)) + max(1, jitter.y.amt) )
  }
  
  plot(x1, # mutThreshRange on x-axis
       y1, # tabulated sums in plotList on y-axis
       ylim = ylims,
       # xlim: +/- jitter.x.amt*2 to accommodate for amount of jittering on x-axis
       xlim = c(min(threshes)-jitter.x.amt*2, max(threshes+jitter.x.amt*2)),
       xaxt="n", xlab="Threshold on number of mutations",
       ylab=paste(ylab.part.1, ylab.part.2),
       cex.lab=1.5, cex.axis=1.5, type="b", lwd=1.5,
       pch=pchs[1], lty=ltys[1], col=cols[1])
  axis(side=1, at=threshes, cex.axis=1.5)
  
  # plot for the rest of the window sizes
  for (i in 1:length(plotList)){
    if (i>=2) {
      
      xi = threshes
      if (jitter.x) {xi = jitter(xi, amount=jitter.x.amt)}
      yi = plotList[[i]]
      if (jitter.y) {yi = jitter(yi, amount=jitter.y.amt)}
      
      points(xi, yi, type='b', lwd=1.5,
             pch=pchs[i], lty=ltys[i], col=cols[i])
    }
  }
  
  # add legend
  if (plotLegend) {
    legend(legendPos, legend = c("Window Size", names(tuneList)),
           horiz = legendHoriz, cex = legendCex,
           pch=c(NA, pchs), lty=c(NA, ltys), col=c(NA, cols))
  }
}


#### Expected frequencies calculating functions ####

#' Calculate expected mutation frequencies
#'
#' \code{expectedMutations} calculates the expected mutation frequencies for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      \code{character} name of the column containing input 
#'                               sequences.
#' @param    germlineColumn      \code{character} name of the column containing 
#'                               the germline or reference sequence.
#' @param    targetingModel      \link{TargetingModel} object. Default is \link{HS5FModel}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    nproc               \code{numeric} number of cores to distribute the operation
#'                               over. If the cluster has already been set the call function with 
#'                               \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                               \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with expected mutation frequencies 
#'           for each region defined in \code{regionDefinition}.
#'          
#'           The columns names are dynamically created based on the regions in  
#'           \code{regionDefinition}. For example, when using the \link{IMGT_V_NO_CDR3}
#'           definition, which defines positions for CDR and FWR, the following columns are
#'           added:  
#'           \itemize{
#'             \item  \code{EXPECTED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{EXPECTED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{EXPECTED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{EXPECTED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           
#' @details
#' Only the part of the sequences defined in \code{regionDefinition} are analyzed. 
#' For example, when using the \link{IMGT_V_NO_CDR3} definition, mutations in
#' positions beyond 312 will be ignored.
#' 
#' @seealso  
#' \link{calcExpectedMutations} is called by this function to calculate the expected 
#' mutation frequencies. See \link{observedMutations} for getting observed 
#' mutation counts. See \link{IMGT_SCHEMES} for a set of predefined 
#' \link{RegionDefinition} objects.
#' 
#' @examples
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
#'
#' # Calculate expected mutations over V region
#' db_exp <- expectedMutations(db,
#'                             sequenceColumn="SEQUENCE_IMGT",
#'                             germlineColumn="GERMLINE_IMGT_D_MASK",
#'                             regionDefinition=IMGT_V_NO_CDR3,
#'                             nproc=1)
#' 
#' # Calculate hydropathy expected mutations over V region
#' db_exp <- expectedMutations(db,
#'                            sequenceColumn="SEQUENCE_IMGT",
#'                            germlineColumn="GERMLINE_IMGT_D_MASK",
#'                            regionDefinition=IMGT_V_NO_CDR3,
#'                            mutationDefinition=HYDROPATHY_MUTATIONS,
#'                            nproc=1)
#'
#' @export
expectedMutations <- function(db, 
                              sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              targetingModel=HS5FModel,
                              regionDefinition=NULL,
                              mutationDefinition=NULL,
                              nproc=1) {
    # Hack for visibility of data.table and foreach index variables
    idx <- NULL
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc(), na.rm=T)
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if (nproc > 1) {        
        cluster <- parallel::makeCluster(nproc, type = "PSOCK")
        parallel::clusterExport(cluster, list('db', 'sequenceColumn', 'germlineColumn', 
                                              'regionDefinition','targetingModel',
                                              'calcExpectedMutations','calculateTargeting',
                                              's2c','c2s','NUCLEOTIDES','HS5FModel',
                                              'calculateMutationalPaths','CODON_TABLE'),
                                envir=environment() )
        registerDoParallel(cluster)
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    
    # Printing status to console
    cat("Calculating the expected frequencies of mutations...\n")
    
    # Calculate targeting for each sequence (based on the germline)
    # Should be a 5 x N matrix where N in the number of nucleotides defined by
    # the regionDefinition
    numbOfSeqs <- nrow(db)
    
    targeting_list <-
        foreach (idx=iterators::icount(numbOfSeqs)) %dopar% {
            calcExpectedMutations(germlineSeq=db[idx, germlineColumn],
                                  inputSeq=db[idx, sequenceColumn],
                                  targetingModel=targetingModel,
                                  regionDefinition=regionDefinition,
                                  mutationDefinition=mutationDefinition)
        }
    
    # Convert list of expected mutation freq to data.frame
    if (is.null(regionDefinition)) {
        labels_length <- length(makeNullRegionDefinition()@labels)
    } else {
        labels_length <- length(regionDefinition@labels)
    }
    expectedMutationFrequencies <- do.call(rbind, lapply(targeting_list, function(x) { 
        length(x) <- labels_length 
        return(x) })) 
    
    expectedMutationFrequencies[is.na(expectedMutationFrequencies)] <- 0
    colnames(expectedMutationFrequencies) <- paste0("EXPECTED_", colnames(expectedMutationFrequencies))
    
    # Properly shutting down the cluster
    if(nproc>1){ parallel::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, expectedMutationFrequencies)
    return(db_new)    
    
}


#' Calculate expected mutation frequencies of a sequence
#'
#' \code{calcExpectedMutations} calculates the expected mutation
#' frequencies of a given sequence. This is primarily a helper function for
#' \link{expectedMutations}. 
#'
#' @param    germlineSeq         germline (reference) sequence.
#' @param    inputSeq            input (observed) sequence. If this is not \code{NULL}, 
#'                               then \code{germlineSeq} will be processed to be the same
#'                               same length as \code{inputSeq} and positions in 
#'                               \code{germlineSeq} corresponding to positions with Ns in 
#'                               \code{inputSeq} will also be assigned an N. 
#' @param    targetingModel      \link{TargetingModel} object. Default is \link{HS5FModel}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#'                               
#' @return   A \code{numeric} vector of the expected frequencies of mutations in the 
#'           regions in the \code{regionDefinition}. For example, when using the default 
#'           \link{IMGT_V_NO_CDR3} definition, which defines positions for CDR and 
#'           FWR, the following columns are calculated:
#'           \itemize{
#'              \item  \code{EXPECTED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                             CDR2 of the V-segment.
#'              \item  \code{EXPECTED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                             of the V-segment.
#'              \item  \code{EXPECTED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                             FWR2 and FWR3 of the V-segment.
#'              \item  \code{EXPECTED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                             FWR3 of the V-segment.
#'            }
#'           
#' @details
#' \code{calcExpectedMutations} calculates the expected mutation frequencies of a 
#' given sequence and its germline. 
#' 
#' Note, only the part of the sequences defined in \code{regionDefinition} are analyzed. 
#' For example, when using the default \link{IMGT_V_NO_CDR3} definition, mutations in
#' positions beyond 312 will be ignored.
#' 
#' @seealso  \link{expectedMutations} calls this function.
#' To create a custom \code{targetingModel} see \link{createTargetingModel}.
#' See \link{calcObservedMutations} for getting observed mutation counts.
#' 
#' @examples
#' # Load example data
#' data(ExampleDb, package="alakazam")
#' 
#' # Use first entry in the exampled data for input and germline sequence
#' in_seq <- ExampleDb[1, "SEQUENCE_IMGT"]
#' germ_seq <-  ExampleDb[1, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Identify all mutations in the sequence
#' calcExpectedMutations(in_seq, germ_seq)
#' 
#' # Identify only mutations the V segment minus CDR3
#' calcExpectedMutations(in_seq, germ_seq, regionDefinition=IMGT_V_NO_CDR3)
#' 
#' # Define mutations based on hydropathy
#' calcExpectedMutations(in_seq, germ_seq, regionDefinition=IMGT_V_NO_CDR3,
#'                       mutationDefinition=HYDROPATHY_MUTATIONS)
#' 
#' @export
calcExpectedMutations <- function(germlineSeq,
                                  inputSeq=NULL,
                                  targetingModel=HS5FModel,
                                  regionDefinition=NULL,
                                  mutationDefinition=NULL) {
    # Assign codon table
    codonTable <- if (is.null(mutationDefinition)) { CODON_TABLE } else { mutationDefinition@codonTable }
    
    # Get targeting
    targeting <- calculateTargeting(germlineSeq=germlineSeq, 
                                    inputSeq=inputSeq,
                                    targetingModel=targetingModel,
                                    regionDefinition=regionDefinition)
    
    # Determine the mutations paths (i.e. determine R and S mutation frequencies)
    mutationalPaths <- calculateMutationalPaths(germlineSeq=c2s(colnames(targeting)), 
                                                regionDefinition=regionDefinition,
                                                codonTable=codonTable)
    
    typesOfMutations <- c("R", "S")
    mutationalPaths[!(mutationalPaths %in% typesOfMutations)] <- NA
    
    if (is.null(regionDefinition)) {
        rdLength <- max(nchar(inputSeq), nchar(germlineSeq), na.rm=TRUE)
        regionDefinition <- makeNullRegionDefinition(rdLength)
    }
    listExpectedMutationFrequencies <- list()
    for(region in regionDefinition@regions){
        for(typeOfMutation in typesOfMutations){
            region_mutation <- paste(region, typeOfMutation, sep="_")    
            
            targeting_region <- targeting[1:4, regionDefinition@boundaries %in% region]
            mutationalPaths_region <- mutationalPaths[, regionDefinition@boundaries[1:ncol(mutationalPaths)] %in% region]
            targeting_typeOfMutation_region <- sum(targeting_region[mutationalPaths_region == typeOfMutation], 
                                                   na.rm=TRUE)
            
            listExpectedMutationFrequencies[[region_mutation]] <- targeting_typeOfMutation_region
            
        }
    }
    expectedMutationFrequencies <- unlist(listExpectedMutationFrequencies)
    expectedMutationFrequencies[!is.finite(expectedMutationFrequencies)] <- NA
    expectedMutationFrequencies <- expectedMutationFrequencies/sum(expectedMutationFrequencies, na.rm=TRUE)
    return(expectedMutationFrequencies)    
}


calculateTargeting <- function(germlineSeq,
                               inputSeq=NULL,
                               targetingModel=HS5FModel,
                               regionDefinition=NULL) {
    
    # If an inputSequence is passed then process the germlineSequence
    # to be the same legth, mask germlineSequence with Ns where inputSequence is also N
    # If not needed then  you may skip this step by passing in inputSequence=NULL 
    # (which is default). 
    if(!is.null(inputSeq)){    
        # Trim the input and germline sequence to the shortest
        len_inputSeq <- nchar(inputSeq)
        len_germlineSeq <- nchar(germlineSeq)
        # If a regionDefinition is passed,
        # then only analyze till the end of the defined length
        if(!is.null(regionDefinition)){
            length_regionDefinition  <- regionDefinition@seqLength
        } else{
            length_regionDefinition <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        }
        len_shortest <- min( c(len_inputSeq,len_germlineSeq,length_regionDefinition),  na.rm=TRUE)
        
        c_inputSeq <- s2c(inputSeq)[1:len_shortest]
        c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
        
        # If the sequence and germline (which now should be the same length) is shorter
        # than the length_regionDefinition, pad it with Ns
        if(len_shortest<length_regionDefinition){
            fillWithNs <- array("N",length_regionDefinition-len_shortest)
            c_inputSeq <- c( c_inputSeq, fillWithNs)
            c_germlineSeq <- c( c_germlineSeq, fillWithNs)
        }
        
        # Mask germline with Ns where input sequence has Ns
        c_germlineSeq[ c_inputSeq=="N" |  !c_inputSeq%in%c(NUCLEOTIDES[1:5],".") ] = "N"    
        s_germlineSeq <- c2s(c_germlineSeq)
    }else{
        s_germlineSeq <- germlineSeq
        c_germlineSeq <- s2c(s_germlineSeq)
    }
    
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
    gaplessSeq <- gsub("\\.\\.\\.", "XXX", s_germlineSeq)
    #If there is a single gap left convert it to an N
    gaplessSeq <- gsub("\\.", "N", gaplessSeq)
    
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    s_germlineSeq <- gsub("XXX", "...", gaplessSeq)
    c_germlineSeq <- s2c(s_germlineSeq)
    # Matrix to hold targeting values for each position in c_germlineSeq
    germlineSeqTargeting <- matrix(NA, 
                                   ncol=nchar(s_germlineSeq), 
                                   nrow=length(NUCLEOTIDES[1:5]),
                                   dimnames=list(NUCLEOTIDES[1:5], c_germlineSeq))
    
    # Now remove the IMGT gaps so that the correct 5mers can be made to calculate
    # targeting. e.g.
    # GAGAAA......TAG yields: "GAGAA" "AGAAA" "GAAAT" "AAATA" "AATAG"
    # (because the IMGT gaps are NOT real gaps in sequence!!!)
    gaplessSeq <- gsub("\\.\\.\\.", "", s_germlineSeq)
    gaplessSeqLen <- nchar(gaplessSeq)
    
    #Slide through 5-mers and look up targeting
    gaplessSeq <- paste("NN", gaplessSeq, "NN", sep="")
    gaplessSeqLen <- nchar(gaplessSeq)
    pos<- 3:(gaplessSeqLen - 2)
    subSeq =  substr(rep(gaplessSeq, gaplessSeqLen - 4), (pos - 2), (pos + 2))
    germlineSeqTargeting_gapless <- targetingModel@targeting[,subSeq]
#     germlineSeqTargeting_gapless <- sapply(subSeq, function(x) { 
#         targetingModel@targeting[, x] })
    
    germlineSeqTargeting[, c_germlineSeq != "."] <- germlineSeqTargeting_gapless  
    
    # Set self-mutating targeting values to be NA
    mutatingToSelf <- colnames(germlineSeqTargeting)
    mutatingToSelf[!(mutatingToSelf %in% NUCLEOTIDES[1:5])] <- "N"
#     # TODO: What's with this <<- business?
#     # TODO: I think this is assigning NA to all self-mutations, which are already NA
#     sapply(1:ncol(germlineSeqTargeting), function(pos) { germlineSeqTargeting[mutatingToSelf[pos], pos] <<- NA })
    
    germlineSeqTargeting[!is.finite(germlineSeqTargeting)] <- NA
    return(germlineSeqTargeting)
}

calculateMutationalPaths <- function(germlineSeq,
                                     inputSeq=NULL,
                                     regionDefinition=NULL,
                                     codonTable=NULL) {    
    # Set codon table if required
    if (is.null(codonTable)) { codonTable <- CODON_TABLE }
    
    # If an inputSequence is passed then process the germlineSequence
    # to be the same length, mask germlineSequence with Ns where inputSequence is also N
    # If this function is being called after running calculateTargeting you may skip
    # this step by passing in inputSequence=NULL (which is default). This way you save
    # some processing time.
    if(!is.null(inputSeq)){    
        # Trim the input and germline sequence to the shortest
        len_inputSeq <- nchar(inputSeq)
        len_germlineSeq <- nchar(germlineSeq)
        # If a regionDefinition is passed,
        # then only analyze till the end of the defined length
        if(!is.null(regionDefinition)){
            length_regionDefinition  <- regionDefinition@seqLength
        } else{
            length_regionDefinition <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        }
        len_shortest <- min( c(len_inputSeq,len_germlineSeq,length_regionDefinition),  na.rm=TRUE)
        
        c_inputSeq <- s2c(inputSeq)[1:len_shortest]
        c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
        
        # If the sequence and germline (which now should be the same length) is shorter
        # than the length_regionDefinition, pad it with Ns
        if(len_shortest<length_regionDefinition){
            fillWithNs <- array("N",length_regionDefinition-len_shortest)
            c_inputSeq <- c( c_inputSeq, fillWithNs)
            c_germlineSeq <- c( c_germlineSeq, fillWithNs)
        }
        
        # Mask germline with Ns where input sequence has Ns
        c_germlineSeq[c_inputSeq=="N" |  !c_inputSeq %in% c(NUCLEOTIDES[1:5], ".") ] = "N"    
        s_germlineSeq <- c2s(c_germlineSeq)
    } else {
        s_germlineSeq <- germlineSeq
        c_germlineSeq <- s2c(s_germlineSeq)
    }
    

    s_germlineSeq_len <- nchar(s_germlineSeq)    
    vecCodons <- sapply({1:(s_germlineSeq_len/3)}*3 - 2, function(x) { substr(s_germlineSeq, x, x + 2) })
    vecCodons[!vecCodons %in% colnames(codonTable)] <- "NNN"
    matMutationTypes = matrix(codonTable[, vecCodons], nrow=4, byrow=F,
                              dimnames=list(NUCLEOTIDES[1:4], c_germlineSeq[ 1:(3 * length(vecCodons))]))
    
    return(matMutationTypes)
}

#### Additional helper functions ####
# Given a nuclotide position, returns the codon number
# e.g. nuc 86  = codon 29
getCodonNumb <- function(nucPos){
  return( ceiling(nucPos/3) )
}

# Given a codon, returns all the nuc positions that make the codon
getCodonNucs <- function(codonNumb){
  getCodonPos(codonNumb*3)
}

# Given a nucleotide postions return the position in the codon
getContextInCodon <- function(nucPos){
  return( {nucPos-1}%%3+1 )
}

# Given a nuclotide position, returns the pos of the 3 nucs that made the codon
# e.g. nuc 86 is part of nucs 85,86,87
getCodonPos <- function(nucPos) {
  codonNum =  (ceiling(nucPos / 3)) * 3
  return ((codonNum - 2):codonNum)
}

# Translate codon to amino acid
translateCodonToAminoAcid <- function(Codon) {
  return (AMINO_ACIDS[Codon])
}

# Given two codons, tells you if the mutation is R or S (based on your definition)
#
# @param   codonFrom         starting codon
# @param   codonTo           ending codon
# @param   aminoAcidClasses  vector of amino acid trait classes
#                            if NULL then R or S is determined by amino acid identity
# @return  Mutation type as "R" (replacement), "S" (silent), "Stop" (stop) or NA (input is NA).
#
# @examples
# # Without classes
# shazam:::mutationType("TTT", "TTC")
# shazam:::mutationType("TTT", "TTA")
# shazam:::mutationType("TTT", "TGA")
#
# # With classes
# classes <- HYDROPATHY_MUTATIONS@classes
# shazam:::mutationType("TTT", "TTC", aminoAcidClasses=classes)
# shazam:::mutationType("TTT", "TTA", aminoAcidClasses=classes)
# shazam:::mutationType("TTT", "TCT", aminoAcidClasses=classes)
# shazam:::mutationType("TTT", "TGA", aminoAcidClasses=classes)
mutationType <- function(codonFrom, codonTo, aminoAcidClasses=NULL) {
    # codonFrom="TTT"; codonTo="TTA"
    # codonFrom="TTT"; codonTo="TGA"
    
    # Translate codons
    aaFrom <- translateCodonToAminoAcid(codonFrom)
    aaTo <- translateCodonToAminoAcid(codonTo)
    
    # If any codon is NA then return NA
    if (any(is.na(c(codonFrom, codonTo, aaFrom, aaTo)))) { 
        return(NA) 
    }
    
    # If any amino acid is Stop then return "Stop"
    if (any(c(aaFrom, aaTo) == "*")) { 
        return("Stop") 
    }
    
    if (is.null(aminoAcidClasses)) {
        # Check for exact identity if no amino acid classes are specified
        mutation <- if (aaFrom == aaTo) { "S" } else { "R" }
    } else {
        # Check for amino acid class identity if classes are specified
        mutation <- if (aminoAcidClasses[aaFrom] == aminoAcidClasses[aaTo]) { "S" } else { "R" }
    }
    return(mutation)
}