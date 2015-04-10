# Mutation profiling
# 
# @author     Mohamed Uduman, Gur Yaari
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.03.08

#' @include shm.R
NULL

#### Consensus functions ####

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
#' @param    trimSequence    optional argument that trims the consensus sequence to 
#'                           the specified length (speicifed in number of nucleotides). 
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
#' library(alakazam)
#' 
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # For every clone identify  SEQUENCE_GAP_CLONE column to db
#' db_new <- getClonalConsensus(db)
#' head(db_new[c(1, 15)])
#'
#' @export
getClonalConsensus <- function(db, 
                               cloneColumn="CLONE", 
                               sequenceColumn="SEQUENCE_IMGT",
                               germlineColumn="GERMLINE_IMGT_D_MASK",
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
    df_ClonalConsensus <- 
        ddply(db, cloneColumn, here(summarize),
              CLONAL_CONSENSUS_SEQUENCE=(clonalConsensus(inputSeq = eval(parse(text=sequenceColumn)),
                                                         glSeq = eval(parse(text=germlineColumn)),
                                                         trimSequence)),
              .progress=progressBar, 
              .parallel=runAsParallel)
    
    # Stop SNOW cluster
    if(nproc > 1) { stopCluster(cluster) }
    
    #db_collapsed_by_clone <- db[ match(df_ClonalConsensus$CLONE, db[,cloneColumn]), ] 
    #db_collapsed_by_clone$CLONAL_CONSENSUS_SEQUENCE <- df_ClonalConsensus$CLONAL_CONSENSUS_SEQUENCE
    
    return(df_ClonalConsensus)
}

# Helper function for getClonalConsensus
# Given 
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
#' \code{addObservedMutations} calculates the observed number of V-segment mutations 
#' within the framework (FW) and complementarity determining (CD) regions of IMGT-gapped 
#' nucleotide sequences. Mutation counts are appended to the input data.frame as 
#' additional columns.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn  name of the column containing IMGT-gapped germline sequences.
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with observed mutation counts for each 
#'           sequence listed in the following columns: 
#'           \itemize{
#'             \item  \code{OBSERVED_R_CDR}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{OBSERVED_S_CDR}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{OBSERVED_R_FWR}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{OBSERVED_S_FWR}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           
#' @details
#' How does this work?
#' 
#' @references
#' Which ones?
#' 
#' @seealso  See \code{\link{addExpectedFrequencies}} for calculating expected mutation
#'           frequencies.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add observed mutations to db
#' db_new <- addObservedMutations(db)
#' head(db_new[c(1, 15:18)])
#'
#' @export
addObservedMutations <- function(db, sequenceColumn="SEQUENCE_GAP", 
                                 germlineColumn="GERMLINE_GAP_D_MASK", nproc=1) {
    numbOfSeqs <- nrow(db)
    if(nproc == 1) {
        pb <- txtProgressBar(min=1, max=numbOfSeqs, width=20)
        cat("Progress: 0%      50%     100%\n")
        cat("          ")
        db[, c("OBSERVED_R_CDR",
               "OBSERVED_S_CDR",
               "OBSERVED_R_FWR",
               "OBSERVED_S_FWR")] = t(sapply(1:nrow(db), function(x) { setTxtProgressBar(pb, x);
                                                                       countMutations(db[x,sequenceColumn], db[x,germlineColumn]) },
                                             simplify="array"))
        cat("\n")
        close(pb)
    } else {
        availableCores <- getnproc()
        if(!(nproc<=availableCores))nproc=availableCores
        cluster <- makeCluster(nproc, type = "SOCK")
        registerDoSNOW(cluster)
        
        obsMutations <-
            foreach(i=icount(numbOfSeqs), .packages='shm', .combine=doparProgressBar(n=numbOfSeqs), .multicombine=TRUE) %dopar% {
                countMutations(db[i,sequenceColumn], db[i,germlineColumn])
            }
        
        stopCluster(cluster)
        
        # Add observed mutations to db
        db[, c("OBSERVED_R_CDR",
               "OBSERVED_S_CDR",
               "OBSERVED_R_FWR",
               "OBSERVED_S_FWR")] <- obsMutations
        cat("\n")
    }
    
    return(db)
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

