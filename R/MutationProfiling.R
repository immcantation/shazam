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
#' \code{getClonalConsensus} identifies the consensus sequence of each clonal 
#' group and appends a column to the input data.frame containing the clonal consensus
#' for each sequence.
#'
#' @param    db              data.frame containing sequence data.
#' @param    cloneColumn     name of the column containing clonal cluster identifiers.
#' @param    sequenceColumn  name of the column containing input sequences.
#' @param    germlineColumn  name of the column containing germline sequences.
#' @param    collapseByClone if TRUE, collapse the \code{db} by the \code{cloneColumn}.
#' @param    trimSequence    optional argument to trim the consensus sequence to a 
#'                           specified length (speicifed in number of nucleotides). 
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with clonal consensus sequences in the
#'           CLONAL_CONSENSUS_SEQUENCE column.
#'
#' @details
#' For seqeunces identified to be part of the same clone. this function defines an 
#' effective sequence that will be representative for all mutations in the clone. Each 
#' position in this consensus (or effective) sequence is created by a weighted sampling 
#' of each mutated base (and non "N", '.' or '-') from all the sequences in the clone. 
#' 
#' For example, in a clone with 5 sequences that have a C at position 1, and 5 sequences
#' with a T at this same position, the consensus sequence will have a C 50% of the time
#' and a T 50% of the time. Thus, the results of this function can change somewhat every 
#' time it is called.
#' 
#' NOTE: The function returns an updated ChangeODB data.frame that collpases all the 
#' sequences by clones defined in the \code{cloneColumn} column passed as a parameter.
#' 
#' Non-terminal branch mutations are defined as the set of mutations that occur on 
#' branches of the lineage tree that are not connected to a leaf. For computational 
#' efficiency, the set of non-terminal branch mutations is approximated as those that are
#' shared between more than one sequence in a clone. In this case the terminal branch 
#' mutations are filtered out.
#' 
#' This function can be paralellized if \code{db} contains thousands of sequences. 
#' Specify the number of cores/CPUS available using the \code{nproc} parameter.
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#' 
#' # Run getClonalConsensus
#' db_new <- getClonalConsensus(db, 
#'                              cloneColumn="CLONE", 
#'                              sequenceColumn="SEQUENCE_IMGT",
#'                              germlineColumn="GERMLINE_IMGT_D_MASK",
#'                              collapseByClone=TRUE)
#'                              
#' @export
getClonalConsensus <- function(db, 
                               cloneColumn="CLONE", 
                               sequenceColumn="SEQUENCE_IMGT",
                               germlineColumn="GERMLINE_IMGT_D_MASK",
                               collapseByClone=TRUE,
                               trimSequence=NULL,
                               nproc=1) {
    
    db[,cloneColumn] <- as.numeric(db[,cloneColumn])
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if (nproc > 1) {        
        runAsParallel <- TRUE # Flag used in ddply to indicate whether to run as paralell
        progressBar <- "none" # Flag used in ddply to indicate whether to display a progress bar
        # Register clusters   
        cluster <- makeCluster(nproc, type="SOCK")
        clusterExport(cluster, list('trimSequence', 'sequenceColumn', 'germlineColumn', 'cloneColumn'), envir=environment())
        clusterEvalQ(cluster, library(shm))
        clusterEvalQ(cluster, library(seqinr))
        registerDoSNOW(cluster)
    } else {
        # If running on single CPU, then set flags 
        runAsParallel <- FALSE # Flag used in ddply to indicate whether to run as paralell
        progressBar <- "text"  # Flag used in ddply to indicate whether to display a progress bar
    }
    
    # Printing status to console
    cat("Building clonal consensus sequences...\n")
    
    # Subset the db by clones and calling the clonalConsensus helper function
    # which identifies the consensus seqeunce for the clone.
    db_ClonalConsensus <- 
        ddply(db, cloneColumn, here(mutate),
              CLONAL_CONSENSUS_SEQUENCE=(clonalConsensus(inputSeq = eval(parse(text=sequenceColumn)),
                                                         glSeq = eval(parse(text=germlineColumn)),
                                                         trimSequence = trimSequence)),
              .progress=progressBar, 
              .parallel=runAsParallel)
    
    # Stop SNOW cluster
    if(nproc > 1) { stopCluster(cluster) }
    
    # If collapseByClone is TRUE then collapse the db by clones
    if(collapseByClone){ 
        uniqueCloneIDs <-  unique(db_ClonalConsensus[,cloneColumn])
        indexOfFirstOccurenceOfClone <- match(uniqueCloneIDs, db_ClonalConsensus[,cloneColumn])
        db_ClonalConsensus <- db_ClonalConsensus[indexOfFirstOccurenceOfClone, ]
    }
    
    return(db_ClonalConsensus)
}



# Helper function for getClonalConsensus
clonalConsensus <- function(inputSeq, glSeq, trimSequence=NULL, nonTerminalOnly=0){
    
    # Since this function is called using ddply, each argument will be a vector
    # the length of the number of sequences that comprises the clone.
    # The trimSequence and nonTerminalOnly arguments will just be repeated, so 
    # only the first element is considered
    #trimSequence <- trimSequence[1]
    #nonTerminalOnly <- nonTerminalOnly[1]
    
    # Find length of shortest input sequence
    # This is used to trim all the sequencesto that length
    # or the length of , if specified, trimSequence which ever is shorter
    len_inputSeq <- sapply(inputSeq, function(x){nchar(x)})
    len_shortest <- min(len_inputSeq, na.rm=TRUE)
    if(!is.null(trimSequence)){len_shortest <- min(len_shortest, trimSequence, na.rm=TRUE)}
    
    #Find the length of the longest germline sequence
    len_glSeq <- sapply(glSeq, function(x){nchar(x)})
    len_longest <- max(len_glSeq, na.rm=TRUE)
    glSeq <- glSeq[(which(len_longest==len_glSeq))[1]]
    
    
    # Identify the consensus sequence
    # TODO: Figure out the T/F
    charInputSeqs <- sapply(inputSeq, function(x){ s2c(x)[1:len_shortest]})
    charGLSeq <- s2c(glSeq)
    matClone <- sapply(1:len_shortest, function(i){
        posNucs = unique(charInputSeqs[i,])
        posGL = charGLSeq[i]
        error = FALSE
        if(posGL=="-" & sum(!(posNucs%in%c("-","N")))==0 ){
            return(c("-",error))
        }
        if(length(posNucs)==1)
            return(c(posNucs[1],error))
        else{
            if("N"%in%posNucs){
                error=TRUE
            }
            if(sum(!posNucs[posNucs!="N"]%in%posGL)==0){
                return( c(posGL,error) )
            }else{
                #return( c(sample(posNucs[posNucs!="N"],1),error) )
                if(nonTerminalOnly==0){
                    return( c(sample(charInputSeqs[i,charInputSeqs[i,]!="N" & charInputSeqs[i,]!=posGL],1),error) )
                }else{
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
    })
    return( c2s(matClone[1,]) )
}



#### Mutation counting functions ####

#' Calculate observed mutations
#'
#' \code{getObservedMutations} calculates the observed number of mutations for each sequence
#' in a given input data.frame (\code{db}). The mutations are determined by comparing the 
#' input sequences (in the column specified by \code{sequenceColumn}) to the germline 
#' seqeunce (in the column specified by \code{germlineColumn}). The mutations are 
#' binned or aggregated by replacement(R) or silent(S) across the different regions of the
#' seqeunce as defined in the \code{regionDefinition}. Typically, this would be the 
#' framework (FWR) and complementarity determining (CDR) regions of IMGT-gapped 
#' nucleotide sequences. Mutation counts are appended to the input data.frame as 
#' additional columns.
#'
#' @param   db                  data.frame containing sequence data.
#' @param   sequenceColumn      name of the column containing sample/input sequences.
#' @param   germlineColumn      name of the column containing germline sequences.
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences. Note, only the part of
#'                              sequences defined in \code{regionDefinition} are analyzed.
#'                              Any mutations outside the definition will be ignored. E.g.
#'                              If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                              used, then mutations in positions greater than 312 will not
#'                              be counted.
#' @param   nproc               number of cores to distribute the operation over.
#' 
#' @return  A modified \code{db} data.frame with observed mutation counts for each 
#'           sequence listed. The columns names are dynamically created based on the
#'           regions in the \code{regionDefinition}. E.g. For the default
#'           \code{\link{IMGT_V_NO_CDR3}} definition, which defines positions for CDR and
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
#' \code{getObservedMutations} calculates the observed number of mutations for each sequence
#' in a given input data.frame (\code{db}). The mutations are determined by comparing the 
#' input sequences (in the column specified by \code{sequenceColumn}) to the germline 
#' sequence (in the column specified by \code{germlineColumn}). The mutations are 
#' binned or aggregated by replacement(R) or silent(S) across the different regions of the
#' seqeunce as defined in the \code{regionDefinition}. Typically, this would be the 
#' framework (FWR) and complementarity determining (CDR) regions of IMGT-gapped 
#' nucleotide sequences. Mutation counts are appended to the input data.frame as 
#' additional columns.
#' 
#' @seealso  
#' \code{\link{countMutations}} is called by this function to get the list of
#' mutations in each sequence.
#' \code{\link{binMutationsByRegion}} is called by this function to aggregate the mutations
#' by the \code{regionDefinition}.
#' 
#' Also see \code{\link{addExpectedFrequencies}} for calculating expected mutation
#'           frequencies.
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#' # Subset data for demo purposes
#' db <- db[1:10,]
#'
#' #Run getObservedMutations()
#' db_new <- getObservedMutations(db,
#'                      sequenceColumn="SEQUENCE_IMGT",
#'                      germlineColumn="GERMLINE_IMGT_D_MASK",
#'                      regionDefinition=IMGT_V_NO_CDR3,
#'                      nproc=1)
#'                      
#' @export
getObservedMutations <- function(db, 
                                 sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 regionDefinition=IMGT_V_NO_CDR3,
                                 nproc=1) {
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- makeCluster(nproc, type = "SOCK")
        registerDoSNOW(cluster)
        clusterExport( cluster, list('db', 'sequenceColumn', 'germlineColumn', 'regionDefinition')  )
        clusterEvalQ( cluster, library("shm") )
    }else{
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Identify all the mutations in the sequences
    # observedMutations helper function returns a list (1 element per sequence)
    # containing an array of mutations (s or R) and the array labels indicate
    # the nucleotide position of the mutations.
    numbOfSeqs <- nrow(db)
    observedMutations_list <-
        foreach( i=icount(numbOfSeqs) ) %dopar% {
            countMutations( db[i,sequenceColumn], 
                            db[i,germlineColumn], 
                            regionDefinition,
                            binByRegions=TRUE)
        }
    
    # Convert list of mutations to data.frame
    labels_length <- length(regionDefinition@labels)
    observed_mutations <- do.call(rbind, lapply(observedMutations_list, function(x){ 
                                        length(x) <- labels_length 
                                        return(x)
                                    })) 
    
    observed_mutations[is.na(observed_mutations)] <- 0
    colnames(observed_mutations) <- paste0("OBSERVED_", colnames(observed_mutations))
    # Properly shutting down the cluster
    stopCluster(cluster)
    
    # Bind the observed mutations to db
    db_new <- cbind(db, observed_mutations)
    return(db_new)    
}



#' Count the number of observed mutations in a given sequence and its germline.
#'
#' \code{countMutations} determines all the mutations in a given input seqeunce and its
#'  germline seqeunce
#'
#' @param   inputSeq            the input sequence
#' @param   germlineSeq         the germline sequence
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences. Note, only the part of
#'                              sequences defined in \code{regionDefinition} are analyzed.
#'                              Any mutations outside the definition will be ignored. E.g.
#'                              If the \code{\link{IMGT_V_NO_CDR3}} definition is used, 
#'                              then mutations in positions greater than 312 will not
#'                              be counted.
#' @param   binByRegions        if TRUE, then aggregate the mutations by the regions
#'                              defined in \code{regionDefinition}, which also needs to be
#'                              passed. If not, binByRegions is ignored.
#' @return   an \code{array} of the mutations (replacement (R) or silent(S)) with the 
#'              names indicatng the nucleotide postion of the mutations in the sequence.
#'              Note. Each mutation is considered independently in its codon context.
#'              If \code{binByRegions}=\code{TRUE}, then the mutations are binned by the regions
#'              defined by \code{regionDefinition} and an \code{array} of R/S mutations
#'              across all the unique regions is returned (See \code{\link{binMutationsByRegion}}). 
#'           
#' @details
#' Counts all the mutations observed in the input sequence. 
#' Note. Each mutation is considered independently in its codon context.
#' 
#' @seealso  
#' See \code{\link{getObservedMutations}} for indentifying and counting the 
#' numer of observed mutations.
#' 
#' See \code{\link{binMutationsByRegion}} for 
#' 
#' @examples
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#' 
#' # Extracting the first entry in the sample db to use for input and germline seqeucnes.
#' inputSeq <- db[1,"SEQUENCE_IMGT"]
#' glSeq <-  db[1,"GERMLINE_IMGT_D_MASK"]
#' 
#' #Identify all mutations in the sequence
#' mutations <- countMutations(inputSeq, glSeq)
#' 
#' #Identify only mutations the V segment minus CDR3
#' mutations <- countMutations(inputSeq, glSeq, regionDefinition=IMGT_V_NO_CDR3)
#'  
#' @export
countMutations <- function(inputSeq, 
                           glSeq, 
                           regionDefinition=NULL, 
                           binByRegions=FALSE) {
    
    # Trim the input and germline sequence to the shortest
    len_inputSeq <- nchar(inputSeq)
    len_glSeq <- nchar(glSeq)
    # If a regionDefinition is passed,
    # then only analyze till the end of the defined length
    if(!is.null(regionDefinition)){
        length_regionDefinition  <- regionDefinition@seqLength
    } else{
        length_regionDefinition <- max(len_inputSeq, len_glSeq, na.rm=TRUE)
    }
    len_shortest <- min( c(len_inputSeq,len_glSeq,length_regionDefinition),  na.rm=TRUE)
    
    c_inputSeq = s2c(inputSeq)[1:len_shortest]
    c_glSeq = s2c(glSeq)[1:len_shortest]
    
    # If the sequence and germline (which now should be the same length) is shorter
    # than the length_regionDefinition, pad it with Ns
    if(len_shortest<length_regionDefinition){
        fillWithNs <- array("N",region_length-len_shortest)
        c_inputSeq <- c( c_inputSeq, fillWithNs)
        c_glSeq <- c( c_glSeq, fillWithNs)
    }
    
    mutations_array <- NA
    mutations = c_glSeq != c_inputSeq
    if(sum(mutations)>0){
        # The nucleotide positions of the mutations
        mutations_pos <- which(mutations==TRUE)
        # For every mutations_pos, extract the entire codon from germline
        mutations_pos_codons <- array(sapply(mutations_pos,getCodonPos))
        c_glSeq_codons <- c_glSeq[mutations_pos_codons]
        # For every mutations_pos, extract the codon from germline (without other mutations 
        # at the same codon, if any).
        c_inputSeq_codons <- array(sapply(mutations_pos, function(x){
            seqP = c_glSeq[getCodonPos(x)]
            seqP[getContextInCodon(x)] = c_inputSeq[x]
            return(seqP)}))
        # split the string of codons into vector of codons
        c_glSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_glSeq_codons)), " ")[[1]]
        c_inputSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_inputSeq_codons)), " ")[[1]]
        
        # Determine whether the mutations are R or S
        mutations_array <- apply(rbind(c_glSeq_codons , c_inputSeq_codons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
        names(mutations_array) = mutations_pos
        mutations_array<- mutations_array[!is.na(mutations_array)]
        if(length(mutations_array)==sum(is.na(mutations_array))){
            mutations_array<-NA    
        }else{
            if(binByRegions & !is.null(regionDefinition)){ 
                mutations_array <- binMutationsByRegion(mutations_array,regionDefinition)
            }
        }        
    }    
    return(mutations_array)
}




#' Bin (i.e. aggregate) mutations (e.g. R or S) by the defined regions (e.g. CDR or FWR)
#'
#' \code{binMutationsByRegion} takes an array of observed mutations (e.g. from 
#' \code{\link{countMutations}}) and bins them by the different regions defined in the 
#' \code{regionDefinition}.
#'
#' @param   mutations_array    array containing the mutations (R/S) with the names
#'                              indicatign the nucleotide positions of the mutations.
#'                              See \code{\link{observedMutations}} for more information.
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences. Note, only the part of
#'                              sequences defined in \code{regionDefinition} are analyzed.
#'                              Any mutations outside the definition will be ignored. E.g.
#'                              If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                              used, then mutations in positions greater than 312 will not
#'                              be counted.
#' 
#' @return an \code{array} of R/S mutations binned across all the unique regions, defined
#' by \code{regionDefinition} is returned.
#' 
#' @seealso  
#' See \code{\link{getObservedMutations}} for indentifying and counting the 
#' numer of observed mutations in a \code{db}.
#' This function is also used in \code{\link{countMutations}}.
#' 
#' @examples
#' # Generate a sampel mutations_array 
#' numbOfMutations <- sample(3:10,1) #Random (between 3-10) number of mutations
#' posOfMutations <- sort(sample(330,numbOfMutations)) #Random positions of mutations
#' mutationTypes <- sample( c("R","S"), length(posOfMutations), replace=TRUE)
#' mutations_array <- array( mutationTypes, dimnames=list(posOfMutations) )
#' binMutationsByRegion(mutations_array, regionDefinition=IMGT_V_NO_CDR3)
#' 
#' @export
binMutationsByRegion <- function( mutations_array, 
                                  regionDefinition=IMGT_V_NO_CDR3){
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




#' Calculate mutation frequencies
#'
#' \code{addMutationFrequencies} calculates mutation frequencies of Ig sequences.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn  name of the column containing IMGT-gapped germline sequences.
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with the mutation frequency of each sequence
#'           in the appended MUTATION_FREQUENCY column.
#'           
#' @details
#' How does this work?
#' 
#' @seealso  See \code{\link{addObservedMutations}} and \code{\link{addExpectedFrequencies}} 
#'           for calculating CDR and FWR observed mutations and expected mutation 
#'           frequencies, respectively.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add mutation frequencies to db
#' db_new <- addMutationFrequencies(db)
#' head(db_new[c(1, 15)])
#'
#' @export
addMutationFrequencies <- function(db, sequenceColumn="SEQUENCE_GAP", 
                                   germlineColumn="GERMLINE_GAP_D_MASK", nproc=1)  {
    
    numbOfSeqs <- nrow(db)
    
    availableCores <- getnproc()
    if(!(nproc<=availableCores))nproc=availableCores
    cluster <- makeCluster(nproc, type = "SOCK")
    registerDoSNOW(cluster)
    
    muFreq <-
        foreach(i=icount(numbOfSeqs), .packages='shm', .combine=doparProgressBar(n=numbOfSeqs)) %dopar% {
            #Count numb of mutations
            inputSeq <-  db[i,sequenceColumn]
            totalMu <- sum(countMutations(inputSeq, db[i,germlineColumn]))
            #Get mu freq
            seq<-substring(inputSeq,1,312)
            seq<-gsub("[N.-]","", seq)
            seqLen <- nchar(seq)
            return(totalMu/seqLen)
        }
    cat("\n")
    stopCluster(cluster)
    
    db[,"MUTATION_FREQUENCY"] <- muFreq
    return(db)
}


#' Determine expected mutation frequencies
#'
#' \code{addExpectedFrequencies} calculate the expected frequency of mutations for
#' each sequence. Expectations are calculated the framework (FW) and complementarity 
#' determining (CD) regions of IMGT-gapped nucleotide sequences. Expected mutation 
#' frequencies are appended to the input data.frame as additional columns.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn  name of the column containing IMGT-gapped germline sequences.
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with observed mutation counts for each 
#'           sequence listed in the following columns: 
#'           \itemize{
#'             \item  \code{EXPECTED_R_CDR}:  expected frequency of replacement mutations 
#'                                            in CDR1 and CDR2 of the V-segment.
#'             \item  \code{EXPECTED_S_CDR}:  expected frequency of silent mutations 
#'                                            in CDR1 and CDR2 of the V-segment.
#'             \item  \code{EXPECTED_R_FWR}:  expected frequency of replacement mutations 
#'                                            in FWR1, FWR2 and FWR3 of the V-segment.
#'             \item  \code{EXPECTED_S_FWR}:  expected frequency of silent mutations 
#'                                            in FWR1, FWR2 and FWR3 of the V-segment.
#'           }
#'
#' @details
#' How does this work?
#' 
#' @references
#' Which ones?
#'
#' @seealso  See \code{\link{addObservedMutations}} for counting observed mutations.
#'
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add observed mutations to db
#' db_new <- addExpectedFrequencies(db)
#' head(db_new[c(1, 15:18)])
#'
#' @export
addExpectedFrequencies <- function(db, sequenceColumn="SEQUENCE_GAP", 
                                   germlineColumn="GERMLINE_GAP_D_MASK", nproc=1) {
    if (nproc == 1) {
        facGL <- factor(db[,germlineColumn])
        facLevels = levels(facGL)
        cat("Computing mutabilities...\n")
        pb <- txtProgressBar(min=1,max=length(facLevels),width=20)
        cat("Progress: 0%      50%     100%\n")
        cat("          ")
        LisGLs_MutabilityU = lapply(1:length(facLevels),  function(x){
            setTxtProgressBar(pb, x);
            computeMutabilities(facLevels[x])
        })
        facIndex = match(facGL,facLevels)
        cat("\n")
        close(pb)
        
        LisGLs_Mutability = lapply(1:nrow(db),  function(x){
            cInput = rep(NA,nchar(db[x,sequenceColumn]))
            cInput[s2c(db[x,sequenceColumn])!="N"] = 1
            LisGLs_MutabilityU[[facIndex[x]]] * cInput
        })
        
        cat("Computing targeting...\n")
        pb <- txtProgressBar(min=1,max=nrow(db),width=20)
        cat("Progress: 0%      50%     100%\n")
        cat("          ")
        LisGLs_Targeting =  lapply(1:nrow(db),  function(x){
            setTxtProgressBar(pb, x);
            computeTargeting(db[x,germlineColumn],LisGLs_Mutability[[x]])
        })
        cat("\n")
        close(pb)
        
        cat("Computing mutation types...\n")
        pb <- txtProgressBar(min=1,max=nrow(db),width=20)
        cat("Progress: 0%      50%     100%\n")
        cat("          ")
        LisGLs_MutationTypes  = lapply(1:nrow(db),function(x){
            setTxtProgressBar(pb, x);
            computeMutationTypes(db[x,germlineColumn])
        })
        cat("\n")
        close(pb)
        
        cat("Computing expected frequencies of mutations...\n")
        pb <- txtProgressBar(min=1,max=nrow(db),width=20)
        cat("Progress: 0%      50%     100%\n")
        cat("          ")
        LisGLs_Exp = lapply(1:nrow(db),  function(x){
            setTxtProgressBar(pb, x);
            computeExpected(LisGLs_Targeting[[x]],LisGLs_MutationTypes[[x]])
        })
        cat("\n")
        close(pb)
        
        ul_LisGLs_Exp =  unlist(LisGLs_Exp)
        matExp <- matrix(ul_LisGLs_Exp,ncol=4,nrow=(length(ul_LisGLs_Exp)/4),byrow=T)
        matExp <- matExp/apply(matExp,1,sum,na.rm=T)
        db[,c("EXPECTED_R_CDR", "EXPECTED_S_CDR", "EXPECTED_R_FWR", "EXPECTED_S_FWR")] <- matExp
    } else {
        availableCores <- getnproc()
        if(!(nproc<=availableCores))nproc=availableCores
        cluster <- makeCluster(nproc, type = "SOCK")
        registerDoSNOW(cluster)
        
        nInterations <- nrow(db)
        matExp <-
            foreach(i=icount(nInterations), .packages='shm', .combine=doparProgressBar(n=nInterations)) %dopar% {
                # Calculate mutability
                seqMutabilities <- computeMutabilities(db[i,germlineColumn])
                # Only include non N positions
                cInput = rep(NA,nchar(db[i,sequenceColumn]))
                cInput[s2c(db[i,sequenceColumn])!="N"] = 1
                seqMutabilities <- seqMutabilities * cInput
                
                # Calculate targeting
                seqTargeting <- computeTargeting(db[i,germlineColumn],seqMutabilities)
                
                # Calculate mutation types
                seqMutationType <- computeMutationTypes(db[i,germlineColumn])
                
                #Expected Freq
                seqExpectedFreq <- computeExpected(seqTargeting, seqMutationType)
                seqExpectedFreq <- seqExpectedFreq/sum(seqExpectedFreq, na.rm=T)
            }
        cat("\n")
        stopCluster(cluster)
        
        db[,c("EXPECTED_R_CDR", "EXPECTED_S_CDR", "EXPECTED_R_FWR", "EXPECTED_S_FWR")] <- matExp
    }
    
    return(db)
}
# List mutations
listMutations <- function(seqInput, seqGL) {
    #if( is.na(c(seqInput, seqGL)) ) return(array(NA,4))
    if (is.na(seqInput) | is.na(seqGL)) { return(NA) }
    seqI = s2c(seqInput)
    seqG = s2c(seqGL)
    matIGL = matrix(c(seqI, seqG), ncol=length(seqI), nrow=2, byrow=T)
    mutations <- analyzeMutations2NucUri(matIGL)
    mutations <- mutations[!is.na(mutations)]
    positions <- as.numeric(names(mutations))
    mutations <- mutations[positions <= VLENGTH]
    if (length(mutations) > 0) {
        return(mutations)
    } else {
        return(NA)
    }
}


#' List the numbers of observed mutations
#'
#' This lists the observed number of mutation.
#'
#' @param   db  a data.frame of the DB file.
#' @param   sequenceColumn  The name of the sequence column.
#' @param   germlineColumn  The name of the germline column.
#' 
#' @return  list of mutations in each clone
#' 
#' @export
listObservedMutations <- function(db, sequenceColumn="SEQUENCE_IMGT", 
                                  germlineColumn="GERMLINE_IMGT_D_MASK")  {
    mutations <- mapply(listMutations, db[, sequenceColumn], db[, germlineColumn], 
                        USE.NAMES=FALSE)
    return(mutations)
}
