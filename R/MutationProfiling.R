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
#' \code{addClonalSequence} identifies the consensus sequence of each clonal group and
#' appends a column to the input data.frame containing the clonal consensus for each 
#' sequence.
#'
#' @param    db              data.frame containing sequence data.
#' @param    cloneColumn     name of the column containing clonal cluster identifiers.
#' @param    sequenceColumn  name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn  name of the column containing IMGT-gapped germline sequences.
#' @param    nproc           number of cores to distribute the operation over.
#' 
#' @return   A modified \code{db} data.frame with clonal consensus sequences in the
#'           SEQUENCE_GAP_CLONE column.
#'
#' @details
#' How does this work?
#' 
#' @references
#' Which ones?
#' 
#' @seealso  What uses this?
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Add SEQUENCE_GAP_CLONE column to db
#' db_new <- addClonalSequence(db)
#' head(db_new[c(1, 15)])
#'
#' @export
addClonalSequence <- function(db, cloneColumn="CLONE", sequenceColumn="SEQUENCE_GAP",
                              germlineColumn="GERMLINE_GAP_D_MASK", nproc=1)  {
    # Start SNOW cluster and set progress bar
    if (nproc > 1) {
        runAsParallel <- TRUE
        progressBar <- "none"  
        cluster <- makeCluster(nproc, type="SOCK")
        registerDoSNOW(cluster)
        clusterEvalQ(cluster, library(shm))
    } else {
        runAsParallel <- FALSE
        progressBar <- "text"
    }
    
    # Generate clonal consensus sequences
    if (progressBar != "none") { cat("-> BUILDING CLONAL SEQUENCES\n") }
    db <- ddply(db, cloneColumn, here(mutate),
                SEQUENCE_GAP_CLONE=(collapseCloneTry(eval(parse(text=sequenceColumn)),
                                                     eval(parse(text=germlineColumn)),
                                                     VLENGTH))[[1]][1],
                .progress=progressBar, .parallel=runAsParallel)
    
    # Stop SNOW cluster
    if(nproc > 1) { stopCluster(cluster) }
    
    return(db)
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
