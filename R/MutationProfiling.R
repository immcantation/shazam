# Mutation profiling
# 
# @author     Mohamed Uduman, Gur Yaari
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.03.08

#' @include shm.R
NULL


#### Clonal Consensus building functions ####

#' Identifies clonal consensus sequences
#'
#' \code{calcDBClonalConsensus} identifies the consensus sequence of each clonal 
#' group and appends a column to the input \code{data.frame} containing the clonal 
#' consensus for each sequence.
#'
#' @param   db                  \code{data.frame} containing sequence data.
#' @param   cloneColumn         \code{character} name of the column containing clonal 
#'                              identifiers.
#' @param   sequenceColumn      \code{character} name of the column containing input 
#'                              sequences.
#' @param   germlineColumn      \code{character} name of the column containing germline 
#'                              sequences.
#' @param   collapseByClone     \code{logical} indicating whether or not to collapse the 
#'                              \code{db} by the \code{cloneColumn}.
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences.
#' @param   nproc               Number of cores to distribute the operation over. If the 
#'                              \code{cluster} has already been set earlier, then pass the 
#'                              \code{cluster}. This will ensure that it is not reset.
#'                              
#' 
#' @return   A modified \code{db} data.frame with clonal consensus sequences in the
#'           CLONAL_CONSENSUS_SEQUENCE column.
#'
#' @details
#' For sequences identified to be part of the same clone, this function defines an 
#' effective sequence that will be representative for all mutations in the clone. Each 
#' position in this consensus (or effective) sequence is created by a weighted sampling 
#' of each mutated base (and non "N", "." or "-" characters) from all the sequences in 
#' the clone. 
#' 
#' For example, in a clone with 5 sequences that have a C at position 1, and 5 sequences
#' with a T at this same position, the consensus sequence will have a C 50\%  and T 50\% 
#' of the time it is called.
#' 
#' The function returns an updated \code{db} that collpases all the sequences by clones 
#' defined in the \code{cloneColumn} column argument.
#' 
#' Non-terminal branch mutations are defined as the set of mutations that occur on 
#' branches of the lineage tree that are not connected to a leaf. For computational 
#' efficiency, the set of non-terminal branch mutations is approximated as those that are
#' shared between more than one sequence in a clone. In this case the terminal branch 
#' mutations are filtered out.
#' 
#' This function can be parallelized if \code{db} contains thousands of sequences. 
#' Specify the number of cores available using the \code{nproc} parameter.
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#' 
#' # Run calcDBClonalConsensus
#' db_new <- calcDBClonalConsensus(db, 
#'                              cloneColumn="CLONE", 
#'                              sequenceColumn="SEQUENCE_IMGT",
#'                              germlineColumn="GERMLINE_IMGT_D_MASK",
#'                              collapseByClone=TRUE)
#' 
#' @export
calcDBClonalConsensus <- function(db, 
                               cloneColumn="CLONE", 
                               sequenceColumn="SEQUENCE_IMGT",
                               germlineColumn="GERMLINE_IMGT_D_MASK",
                               collapseByClone=TRUE,
                               regionDefinition=IMGT_V_NO_CDR3,
                               nproc=1) {
    # Check for valid columns
    check <- checkColumns(db, c(cloneColumn, sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    
    db[,cloneColumn] <- as.numeric(db[,cloneColumn])
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of CLONES, instead of doing dplyr
    dt <- data.table(db)
    setkeyv(dt, cloneColumn )
    # Get the group indexes
    dt <- dt[ , list( yidx = list(.I) ) , by = list(CLONE) ]
    groups <- dt[,yidx]
    lenGroups <- length(groups)
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else {
        if(nproc != 0) { cluster <- makeCluster(nproc, type="SOCK") }
        snow::clusterExport( cluster, list('db', 
                                     'sequenceColumn', 'germlineColumn', "cloneColumn",
                                     'regionDefinition', 
                                     'groups', 'c2s', 's2c', 'words', 'translate'), 
                       envir=environment() )
        snow::clusterEvalQ(cluster, library(shm))
        registerDoSNOW(cluster)
    }
    
    # Printing status to console
    cat("Building clonal consensus sequences...\n")
    
    list_ClonalConsensus <-
      foreach(i=iterators::icount(lenGroups), .combine=c, .verbose=FALSE, .errorhandling='pass') %dopar% {
        calcClonalConsensus(inputSeq = db[groups[[i]],sequenceColumn],
                           germlineSeq = db[groups[[i]],germlineColumn],
                           regionDefinition = regionDefinition)
      }
    
    # Stop SNOW cluster
    if(nproc > 1) { snow::stopCluster(cluster) }
    
    # If collapseByClone is TRUE then collapse the db by clones
    if(collapseByClone==TRUE){ 
        uniqueCloneIDs <-  unique(db[,cloneColumn])
        indexOfFirstOccurenceOfClone <- match(uniqueCloneIDs, db[,cloneColumn])
        db_ClonalConsensus <- db[indexOfFirstOccurenceOfClone, ]
        db_ClonalConsensus$CLONAL_CONSENSUS_SEQUENCE <- unlist(list_ClonalConsensus)
    }else{
      # Match the ClonalConsensus to all the sequences in the clone
      vec_ClonalConsensus <- unlist(list_ClonalConsensus)
      expanded_ClonalConsensus <- tapply(db[,cloneColumn],1:nrow(db),function(x){return(vec_ClonalConsensus[x])})
      db_ClonalConsensus <- db
      db_ClonalConsensus$CLONAL_CONSENSUS_SEQUENCE <- unlist(expanded_ClonalConsensus)
    }
    
    return(db_ClonalConsensus)
}



# Helper function for calcDBClonalConsensus
calcClonalConsensus <- function(inputSeq, germlineSeq, 
                               regionDefinition=NULL, 
                               nonTerminalOnly=0){    
    
    # Find length of shortest input sequence
    # This is used to trim all the sequencesto that length
    # or, if a regionDefinition is passed, then only analyze till the end of the defined length
    len_inputSeq <- sapply(inputSeq, function(x){nchar(x)})
    len_shortest <- min(len_inputSeq, na.rm=TRUE)
    if(!is.null(regionDefinition)){len_shortest <- min(len_shortest, regionDefinition@seqLength, na.rm=TRUE)}        
    
    #Find the length of the longest germline sequence
    len_germlineSeq <- sapply(germlineSeq, function(x){nchar(x)})
    len_longest <- max(len_germlineSeq, na.rm=TRUE)
    germlineSeq <- germlineSeq[(which(len_longest==len_germlineSeq))[1]]    
    
    # Identify the consensus sequence
    # TODO: Figure out the T/F
    charInputSeqs <- sapply(inputSeq, function(x){ s2c(x)[1:len_shortest]})
    charGLSeq <- s2c(germlineSeq)
    matClone <- sapply(1:len_shortest, function(i){
        
        # Identify the nucleotides (in seqs and germline) at the current position
        posNucs = unique(charInputSeqs[i,])
        posGL = charGLSeq[i]
        error = FALSE
        
        # If the current position is a gap in both germline and the sequence,
        # return a gap
        if(posGL=="-" & sum(!(posNucs%in%c("-","N")))==0 ){
            return(c("-",error))
        }
        
        # If all the sequences in the clone have the same nucleotide at the current
        # position, return the value at the current positions
        if(length(posNucs)==1)
            return(c(posNucs[1],error))
        else{         
            # if the nucleotides at the current position are not all the same
            
            # TODO: The error message is not doing anything currently... 
            if("N"%in%posNucs){
                error=TRUE
            }
            
            # If the current nucleotide matches germline, return germline 
            if(sum(!posNucs[posNucs!="N"]%in%posGL)==0){
                return( c(posGL,error) )
            }else{
                #return( c(sample(posNucs[posNucs!="N"],1),error) )
                
                # If we look at all nodes (including terminal nodes), sample a nucleotide from the possible
                # nucleotides in the clonal sequences at this position
                if(nonTerminalOnly==0){
                    return( c(sample(charInputSeqs[i,charInputSeqs[i,]!="N" & charInputSeqs[i,]!=posGL],1),error) )
                }else{
                    
                    # If we look at only non-terminal nodes, we only sample the nucleotides that appear more 
                    # than once (this is a quick approximation)
                    posNucs = charInputSeqs[i,charInputSeqs[i,]!="N" & charInputSeqs[i,]!=posGL]
                    posNucsTable = table(posNucs)
                    if(sum(posNucsTable>1)==0){
                        return( c(posGL,error) )
                    }else{
                        return( c(sample( posNucs[posNucs%in%names(posNucsTable)[posNucsTable>1]],1),error) )
                    }
                }
                
            }
        }
        if(error==TRUE){warning("Error while attempting to collapse by clone!")}
    })
    
    return( c2s(matClone[1,]) )
}



#### Mutation counting functions ####

#' Calculate observed numbers of mutations
#'
#' \code{calcDBObservedMutations} calculates the observed number of mutations for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                 \code{data.frame} containing sequence data.
#' @param    sequenceColumn     \code{character} name of the column containing input 
#'                              sequences.
#' @param    germlineColumn     \code{character} name of the column containing 
#'                              the germline or reference sequence.
#' @param    frequency          \code{logical} indicating whether or not to calculate
#'                              mutation frequencies. Default is \code{FALSE}.
#' @param    regionDefinition   \link{RegionDefinition} object defining the regions
#'                              and boundaries of the Ig sequences. 
#' @param    nproc              number of cores to distribute the operation over. If the 
#'                              cluster has already been set the call function with 
#'                              \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                              \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with observed mutation counts for each 
#'           sequence listed. The columns names are dynamically created based on the
#'           regions in the \code{regionDefinition}. For example, when using the default
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
#' @seealso  \link{calcObservedMutations} is called by this function to get the list of mutations 
#'           in each sequence. \link{binMutationsByRegion} is called by this function to 
#'           aggregate the mutations by the \code{regionDefinition}. 
#'           See \link{calcDBExpectedMutations} for calculating expected mutation frequencies.
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#'
#' # Run calcDBObservedMutations()
#' db_new <- calcDBObservedMutations(db,
#'                                sequenceColumn="SEQUENCE_IMGT",
#'                                germlineColumn="GERMLINE_IMGT_D_MASK",
#'                                frequency=TRUE,
#'                                regionDefinition=IMGT_V_NO_CDR3,
#'                                nproc=1)
#'                      
#' @export
calcDBObservedMutations <- function(db, 
                                   sequenceColumn="SEQUENCE_IMGT",
                                   germlineColumn="GERMLINE_IMGT_D_MASK",
                                   frequency=FALSE,
                                   regionDefinition=IMGT_V_NO_CDR3,
                                   nproc=1) {
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    if(frequency==TRUE){  
      regionDefinition=NULL  
    }
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- snow::makeCluster(nproc, type = "SOCK")
        snow::clusterExport( cluster, list('db', 'sequenceColumn', 'germlineColumn', 'regionDefinition', 'frequency'), envir=environment() )
        snow::clusterEvalQ( cluster, library("shm") )
        registerDoSNOW(cluster)
    } else if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Printing status to console
    cat("Calculating observed number of mutations...\n")
    
    # Identify all the mutations in the sequences
    # observedMutations helper function returns a list (1 element per sequence)
    # containing an array of mutations (s or R) and the array labels indicate
    # the nucleotide position of the mutations.
    numbOfSeqs <- nrow(db)
    observedMutations_list <-
        foreach( i=iterators::icount(numbOfSeqs) ) %dopar% {
            calcObservedMutations( db[i,sequenceColumn], 
                            db[i,germlineColumn],
                            frequency,
                            regionDefinition,
                            binByRegions=TRUE)
        }
    
    # Convert list of mutations to data.frame
    if(!is.null(regionDefinition)) {
      labels_length <- length(regionDefinition@labels)
    } else{
      labels_length=1
    }
    observed_mutations <- do.call(rbind, lapply(observedMutations_list, function(x){ 
        length(x) <- labels_length 
        return(x)
    })) 
    
    sep <- ""
    if(ncol(observed_mutations)>1) sep <- "_"
    observed_mutations[is.na(observed_mutations)] <- 0
    if(frequency==TRUE){
      colnames(observed_mutations) <- paste0( paste0("MU_FREQ",sep), colnames(observed_mutations))
    }else{
      colnames(observed_mutations) <- paste0(paste0("OBSERVED",sep), colnames(observed_mutations))
    }
    
    # Properly shutting down the cluster
    if(nproc>1){ snow::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, observed_mutations)
    return(db_new)    
}



#' Count the number of observed mutations in a sequence.
#'
#' \code{calcObservedMutations} determines all the mutations in a given input seqeunce compared
#' to its germline sequence.
#'
#' @param    inputSeq          input sequence.
#' @param    germlineSeq       germline sequence.
#' @param    frequency          \code{logical} indicating whether or not to calculate
#'                              mutation frequencies. Default is \code{FALSE}.
#' @param    regionDefinition  \link{RegionDefinition} object defining the regions
#'                             and boundaries of the Ig sequences. Note, only the part of
#'                             sequences defined in \code{regionDefinition} are analyzed.
#'                             This argument is required if \code{binByRegions = TRUE}.                      
#' @param    binByRegions      if \code{TRUE} then aggregate the mutations by 
#'                             the regions defined in \code{regionDefinition}.
#' @return   An \code{array} of the mutations, replacement (R) or silent(S), with the 
#'           names indicatng the nucleotide postion of the mutations in the sequence.
#'           
#' @details
#' Each mutation is considered independently in its codon context. Note, only the part of 
#' \code{inputSeq} defined in \code{regionDefinition} is analyzed. For example, when using 
#' the default \link{IMGT_V_NO_CDR3} definition, then mutations in positions beyond 
#' 312 will be ignored.
#' 
#' @seealso  See \link{calcDBObservedMutations} for counting the number of observed mutations. 
#'           See \link{binMutationsByRegion} for aggregation of mutations by region. 
#' 
#' @examples
#' dbPath <- system.file("extdata", "Influenza.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#' 
#' # Extracting the first entry in the sample db to use for input and germline sequences.
#' inputSeq <- db[1, "SEQUENCE_IMGT"]
#' germlineSeq <-  db[1, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Identify all mutations in the sequence
#' mutations <- calcObservedMutations(inputSeq, germlineSeq)
#' 
#' #Identify only mutations the V segment minus CDR3
#' mutations <- calcObservedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3)
#'  
#' @export
calcObservedMutations <- function(inputSeq, 
                           germlineSeq,
                           frequency=FALSE,
                           regionDefinition=NULL, 
                           binByRegions=FALSE) {
    
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
    
    mutations_array <- NA
    mutations = (c_germlineSeq != c_inputSeq) & (c_germlineSeq%in%NUCLEOTIDES[1:5]) & (c_inputSeq%in%NUCLEOTIDES[1:5])
    if(sum(mutations)>0){
        # The nucleotide positions of the mutations
        mutations_pos <- which(mutations==TRUE)
        # For every mutations_pos, extract the entire codon from germline
        mutations_pos_codons <- array(sapply(mutations_pos,getCodonPos))
        c_germlineSeq_codons <- c_germlineSeq[mutations_pos_codons]
        # For every mutations_pos, extract the codon from germline (without other mutations 
        # at the same codon, if any).
        c_inputSeq_codons <- array(sapply(mutations_pos, function(x){
            seqP = c_germlineSeq[getCodonPos(x)]
            seqP[getContextInCodon(x)] = c_inputSeq[x]
            return(seqP)}))
        # split the string of codons into vector of codons
        c_germlineSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_germlineSeq_codons)), " ")[[1]]
        c_inputSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_inputSeq_codons)), " ")[[1]]
        
        # Determine whether the mutations are R or S
        mutations_array <- apply(rbind(c_germlineSeq_codons , c_inputSeq_codons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
        names(mutations_array) = mutations_pos
        mutations_array<- mutations_array[!is.na(mutations_array)]
        if(length(mutations_array)==sum(is.na(mutations_array))){
            mutations_array<-NA    
        }else{ #If there are mutations present proceed to aggregate (if requested)
          if(frequency==TRUE){
            nonNLength <-length( c_inputSeq[ c_inputSeq%in%NUCLEOTIDES[1:4]] )
            nMutations <- length(mutations_array)
            mutations_array <- nMutations/nonNLength
            if(nonNLength==0) mutations_array <- 0
          }else{
            if(binByRegions & !is.null(regionDefinition)){ 
              mutations_array <- binMutationsByRegion(mutations_array,regionDefinition)
            }
          }
        }        
    }    
    return(mutations_array)
}




#' Aggregate mutations by region
#'
#' \code{binMutationsByRegion} takes an array of observed mutations (e.g. from 
#' \code{\link{calcObservedMutations}}) and bins them by the different regions defined in the 
#' \code{regionDefinition}.
#'
#' @param   mutations_array    \code{array} containing the mutations (R/S) with the names
#'                             indicating the nucleotide positions of the mutations.                             
#' @param   regionDefinition   \code{\link{RegionDefinition}} object defining the regions
#'                             and boundaries of the Ig sequences.
#' 
#' @return An \code{array} of R/S mutations binned across all the unique regions defined
#' by \code{regionDefinition}.
#' 
#' @details
#' Note, only the part of sequences defined in \code{regionDefinition} are analyzed.
#' For example, if the default \link{IMGT_V_NO_CDR3} definition is used, then mutations
#' in positions beyond 312 will be ignored.
#' 
#' @seealso  
#' See \code{\link{calcDBObservedMutations}} for identifying and counting the 
#' number of observed mutations.
#' This function is also used in \code{\link{calcObservedMutations}}.
#' 
#' @examples
#' # Generate a random mutation array
#' numbOfMutations <- sample(3:10, 1) 
#' posOfMutations <- sort(sample(330, numbOfMutations))
#' mutationTypes <- sample(c("R","S"), length(posOfMutations), replace=TRUE)
#' mutations_array <- array(mutationTypes, dimnames=list(posOfMutations))
#' 
#' # Random mutations
#' binMutationsByRegion(mutations_array, regionDefinition=IMGT_V_NO_CDR3)
#' 
#' @export
binMutationsByRegion <- function(mutations_array, 
                                 regionDefinition=IMGT_V_NO_CDR3) {
    # Make a factor of R/S
    mutatedPositions <- as.numeric(names(mutations_array))
    mutations <- array(NA,  dim=regionDefinition@seqLength)
    mutations[mutatedPositions] <- mutations_array
    mutations <- mutations[1:regionDefinition@seqLength]
    mutations <- factor(mutations,levels=c("R","S"))
    
    mutations_region_counts <-  
        collapseMatrixToVector( table(regionDefinition@boundaries, mutations) )
    
    sortingOrder <- match(regionDefinition@labels, names(mutations_region_counts))
    mutations_region_counts <- mutations_region_counts[sortingOrder]
    return(mutations_region_counts)
}



#### Expected frequencies calculating functions ####

#' Calculate expected mutation frequencies
#'
#' \code{calcDBExpectedMutations} calculates the expected mutation frequencies for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                \code{data.frame} containing sequence data.
#' @param    sequenceColumn    \code{character} name of the column containing input 
#'                             sequences.
#' @param    germlineColumn    \code{character} name of the column containing 
#'                             the germline or reference sequence.
#' @param    targetingModel    \link{TargetingModel} object. Default is \link{HS5FModel}.
#' @param    regionDefinition  \link{RegionDefinition} object defining the regions
#'                             and boundaries of the Ig sequences.
#' @param    nproc             \code{numeric} number of cores to distribute the operation
#'                             over. If the cluster has already been set the call function with 
#'                             \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                             \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with expected mutation frequencies 
#'           for each region defined in \code{regionDefinition}.
#'          
#'           The columns names are dynamically created based on the regions in  
#'           \code{regionDefinition}. For example, when using the default \link{IMGT_V_NO_CDR3}
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
#' For example, when using the default \link{IMGT_V_NO_CDR3} definition, mutations in
#' positions beyond 312 will be ignored.
#' 
#' @seealso  \link{calcExpectedMutations} is called by this function to calculate the
#' expected mutation frequencies. See \link{calcDBObservedMutations} for getting observed 
#' mutation counts.
#' 
#' @examples
#' # Load example data
#' library("alakazam")
#' dbPath <- system.file("extdata", "Influenza.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#' # Subset data for demo purposes
#' db <- db[1:10, ]
#'
#' # Run calcDBExpectedMutations()
#' db <- calcDBExpectedMutations(db,
#'                               sequenceColumn="SEQUENCE_IMGT",
#'                               germlineColumn="GERMLINE_IMGT_D_MASK",
#'                               regionDefinition=IMGT_V_NO_CDR3,
#'                               nproc=1)
#'
#' @export
calcDBExpectedMutations <- function(db, 
                                    sequenceColumn="SEQUENCE_IMGT",
                                    germlineColumn="GERMLINE_IMGT_D_MASK",
                                    targetingModel=HS5FModel,
                                    regionDefinition=IMGT_V_NO_CDR3,
                                    nproc=1) {
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
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
    if(nproc>1){        
        cluster <- snow::makeCluster(nproc, type = "SOCK")
        snow::clusterExport( cluster, list('db', 'sequenceColumn', 'germlineColumn', 
                                     'regionDefinition','targetingModel'), envir=environment() )
        snow::clusterEvalQ( cluster, library("shm") )
        registerDoSNOW(cluster)
    } else if( nproc==1 ) {
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
        foreach( i=iterators::icount(numbOfSeqs) ) %dopar% {
          calcExpectedMutations( germlineSeq = db[i,germlineColumn],
                                inputSeq = db[i,sequenceColumn],
                                targetingModel = HS5FModel,
                                regionDefinition = regionDefinition)
        }
    
    # Convert list of expected mutation freq to data.frame
    labels_length <- length(regionDefinition@labels)
    expectedMutationFrequencies <- do.call(rbind, lapply(targeting_list, function(x){ 
        length(x) <- labels_length 
        return(x)
    })) 
    
    expectedMutationFrequencies[is.na(expectedMutationFrequencies)] <- 0
    colnames(expectedMutationFrequencies) <- paste0("EXPECTED_", colnames(expectedMutationFrequencies))
    
    # Properly shutting down the cluster
    if(nproc>1){ snow::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, expectedMutationFrequencies)
    return(db_new)    
    
}


#' Calculate expected mutation frequencies of a sequence
#'
#' \code{calcExpectedMutations} calculates the expected mutation
#' frequencies of a given sequence. This is primarily a helper function for
#' \link{calcDBExpectedMutations}. 
#'
#' @param    germlineSeq       germline (reference) sequence.
#' @param    inputSeq          input (observed) sequence. If this is not \code{NULL}, 
#'                             then \code{germlineSeq} will be processed to be the same
#'                             same length as \code{inputSeq} and positions in 
#'                             \code{germlineSeq} corresponding to positions with Ns in 
#'                             \code{inputSeq} will also be assigned an N. 
#' @param    targetingModel    \link{TargetingModel} object. Default is \link{HS5FModel}.
#' @param    regionDefinition  \link{RegionDefinition} object defining the regions
#'                             and boundaries of the Ig sequences.
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
#' @seealso  \link{calcDBExpectedMutations} calls this function.
#' To create a custom \code{targetingModel} see \link{createTargetingModel}.
#' See \link{calcObservedMutations} for getting observed mutation counts.
#' 
#' @examples
#' dbPath <- system.file("extdata", "Influenza.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#' 
#' # Extracting the first entry in the sample db to use for input and germline sequences.
#' inputSeq <- db[1, "SEQUENCE_IMGT"]
#' germlineSeq <-  db[1, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Identify all mutations in the sequence
#' expectedFreq <- calcExpectedMutations(inputSeq, germlineSeq)
#' 
#' # Identify only mutations the V segment minus CDR3
#' expectedFreq <- calcObservedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3)
#'
#' @export
calcExpectedMutations <- function(germlineSeq,
                                  inputSeq=NULL,
                                  targetingModel=HS5FModel,
                                  regionDefinition=IMGT_V_NO_CDR3){
    
    
    targeting <- 
        calculateTargeting(germlineSeq = germlineSeq, 
                           inputSeq = inputSeq,
                           targetingModel = targetingModel,
                           regionDefinition = regionDefinition)
    
    # Determine the mutations paths (i.e. determine R and S mutation frequencies)
    mutationalPaths <- 
        calculateMutationalPaths( germlineSeq = c2s(colnames(targeting)), 
                                  regionDefinition = regionDefinition)
    
    typesOfMutations <- c("R","S")
    mutationalPaths[!(mutationalPaths%in%typesOfMutations)] <- NA
    
    listExpectedMutationFrequencies <- list()
    for(region in regionDefinition@regions){
        for(typeOfMutation in typesOfMutations){
            region_mutation <- paste(region,typeOfMutation,sep="_")    
            
            targeting_typeOfMutation_region <- 
                sum(targeting[ regionDefinition@boundaries%in%region & 
                                   mutationalPaths%in%typeOfMutation ], na.rm=TRUE )
            
            listExpectedMutationFrequencies[[region_mutation]] <- targeting_typeOfMutation_region
            
        }
    }
    expectedMutationFrequencies <- unlist(listExpectedMutationFrequencies)
    expectedMutationFrequencies[!is.finite(expectedMutationFrequencies)] <- NA
    expectedMutationFrequencies <- expectedMutationFrequencies/sum(expectedMutationFrequencies,na.rm=TRUE)
    return(expectedMutationFrequencies)    
}


calculateTargeting <- function(germlineSeq,
                               inputSeq=NULL,
                               targetingModel=HS5FModel,
                               regionDefinition=IMGT_V_NO_CDR3) {
    
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
                                   dimnames=list( NUCLEOTIDES[1:5], c_germlineSeq))
    
    # Now remove the IMGT gaps so that the correct 5mers can be made to calculate
    # targeting. e.g.
    # GAGAAA......TAG yields: "GAGAA" "AGAAA" "GAAAT" "AAATA" "AATAG"
    # (because the IMGT gaps are NOT real gaps in sequence!!!)
    gaplessSeq <- gsub("\\.\\.\\.", "", s_germlineSeq)
    gaplessSeqLen <- nchar(gaplessSeq)
    
    #Slide through 5-mers and look up targeting
    gaplessSeq <- paste("NN",gaplessSeq,"NN",sep="")
    gaplessSeqLen <- nchar(gaplessSeq)
    pos<- 3:(gaplessSeqLen-2)
    subSeq =  substr(rep(gaplessSeq,gaplessSeqLen-4),(pos-2),(pos+2))    
    germlineSeqTargeting_gapless <- sapply(subSeq,function(x){ 
        targetingModel@targeting[,x] 
    })
    
    germlineSeqTargeting[,c_germlineSeq!="."] <- germlineSeqTargeting_gapless  
    
    # Set self-mutating targeting values to be NA
    mutatingToSelf <- colnames(germlineSeqTargeting)
    mutatingToSelf[ !(mutatingToSelf%in%NUCLEOTIDES[1:5])  ] <- "N"
    tmp <- sapply( 1:ncol(germlineSeqTargeting), function(pos){ germlineSeqTargeting[ mutatingToSelf[pos],pos ] <<- NA })
    
    germlineSeqTargeting[!is.finite(germlineSeqTargeting)] <- NA
    return(germlineSeqTargeting)
}

calculateMutationalPaths <- function(germlineSeq,
                                     inputSeq=NULL,
                                     regionDefinition=IMGT_V_NO_CDR3) {    
    
    # If an inputSequence is passed then process the germlienSequence
    # to be the same legth, mask germlienSequence with Ns where inputSequence is also N
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
        c_germlineSeq[ c_inputSeq=="N" |  !c_inputSeq%in%c(NUCLEOTIDES[1:5],".") ] = "N"    
        s_germlineSeq <- c2s(c_germlineSeq)
    }else{
        s_germlineSeq <- germlineSeq
        c_germlineSeq <- s2c(s_germlineSeq)
    }
    
    s_germlineSeq_len <- nchar(s_germlineSeq)    
    vecCodons = sapply({1:(s_germlineSeq_len/3)}*3-2,function(x){substr(s_germlineSeq,x,x+2)})
    vecCodons[!vecCodons %in% colnames(CODON_TABLE)] = "NNN"
    matMutationTypes = matrix( CODON_TABLE[,vecCodons], nrow=4, byrow=F,
                               dimnames=list(NUCLEOTIDES[1:4], c_germlineSeq))
    
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
  codonNum =  (ceiling(nucPos/3))*3
  return ((codonNum-2):codonNum)
}

# Translate codon to amino acid
translateCodonToAminoAcid <- function(Codon) {
  return (AMINO_ACIDS[Codon])
}

# Given two codons, tells you if the mutation is R or S (based on your definition)
mutationType <- function(codonFrom, codonTo, testID=1) {
  if (testID == 4) {
    if (is.na(codonFrom) | is.na(codonTo) | is.na(translateCodonToAminoAcid(codonFrom)) | is.na(translateCodonToAminoAcid(codonTo)) ){
      return(NA)
    } else {
      mutationType = "S"
      if( translateAminoAcidToTraitChange(translateCodonToAminoAcid(codonFrom)) != translateAminoAcidToTraitChange(translateCodonToAminoAcid(codonTo)) ){
        mutationType = "R"
      }
      if(translateCodonToAminoAcid(codonTo)=="*" | translateCodonToAminoAcid(codonFrom)=="*"){
        mutationType = "Stop"
      }
      return(mutationType)
    }
  } else if (testID == 5) {
    if (is.na(codonFrom) | is.na(codonTo) | is.na(translateCodonToAminoAcid(codonFrom)) | is.na(translateCodonToAminoAcid(codonTo)) ){
      return(NA)
    } else {
      if (codonFrom==codonTo) {
        mutationType = "S"
      } else {
        codonFrom = s2c(codonFrom)
        codonTo = s2c(codonTo)
        mutationType = "Stop"
        nucOfI = codonFrom[which(codonTo!=codonFrom)]
        if(nucOfI=="C"){
          mutationType = "R"
        }else if(nucOfI=="G"){
          mutationType = "S"
        }
      }
      return(mutationType)
    }
  } else {
    if (is.na(codonFrom) | is.na(codonTo) | is.na(translateCodonToAminoAcid(codonFrom)) | is.na(translateCodonToAminoAcid(codonTo)) ){
      return(NA)
    } else {
      mutationType = "S"
      if( translateCodonToAminoAcid(codonFrom) != translateCodonToAminoAcid(codonTo) ){
        mutationType = "R"
      }
      if(translateCodonToAminoAcid(codonTo)=="*" | translateCodonToAminoAcid(codonFrom)=="*"){
        mutationType = "Stop"
      }
      return(mutationType)
    }
  }
}