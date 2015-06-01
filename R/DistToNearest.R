# Generates distance to nearest neighbor by using S5F mutability model
#
# @author     Namita Gupta, Gur Yaari, Mohamed Uduman
# @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @date       2015.03.09

#' @include shm.R
NULL


# Returns a 5-mer sliding window of given sequence
#
# @param   strSequence   The sequence string
# @return  An array of 5-mer sliding windows
#
# @examples
# slidingArrayOf5mers("ACGTNACGTNACGTN")
slidingArrayOf5mers <- function(strSequence){
    seqLength <- nchar(strSequence)
    return( substr( rep(strSequence,seqLength-4), 1:(seqLength-4), 5:seqLength ) )
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1                first nucleotide sequence, broken into 5mers.
# @param    seq2                second nucleotide sequence, broken into 5mers.
# @param    targeting_model     targeting model.
# @param    normalize           The method of normalization. Default is "none".
#                               "length" = normalize distance by length of junction.
# @return   distance between two sequences.
#
# @examples
# seq1 = c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
# seq2 = c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
#
# distSeq5mers(seq1, seq2, HS5FModel)
distSeq5mers <- function(seq1, seq2, targeting_model, 
                        normalize=c("none" ,"length", "mutations")) {
  # Evaluate choices
  normalize <- match.arg(normalize)
  
  # Get distance from targeting model
  targeting_dist <- calcTargetingDistance(targeting_model)
  
  # Compute length of sequence (for normalization, if specified)
  juncLength <- length(seq1)

  # Compute distance only on fivemers that have mutations
  fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
  fivemersWithNonNuc <- ( !is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & !is.na(match(substr(seq2,3,3),c("A","C","G","T"))) )
  fivemersWithMu <- fivemersWithMu & fivemersWithNonNuc
  seq1 <- seq1[fivemersWithMu]
  seq2 <- seq2[fivemersWithMu]
  
  # Number of mutations (for normalization, if specified)
  numbOfMutation <- sum(fivemersWithMu)

  avgDist <- NA
  a <- tryCatch({
    if (length(seq1)==1){
      seq1_to_seq2 <- targeting_dist[substr(seq2,3,3),seq1]
      seq2_to_seq1 <- targeting_dist[substr(seq1,3,3),seq2]
    } else {
      seq1_to_seq2 <- sum( diag(targeting_dist[substr(seq2,3,3),seq1]) )
      seq2_to_seq1 <- sum( diag(targeting_dist[substr(seq1,3,3),seq2]) )
    }
    avgDist <- mean(c(seq1_to_seq2, seq2_to_seq1))
  },error = function(e){
    return(NA)
  })

  # Normalize distances
  if (normalize == "length") { 
      avgDist <- avgDist/juncLength
  } else if (normalize == "mutations") { 
      avgDist <- avgDist/numbOfMutation 
  }
  
  return(avgDist)
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1          first nucleotide sequence.
# @param    seq2          second nucleotide sequence.
# @param    normalize     The method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   distance between two sequences.
#
# @examples
# seq1 = c("A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T")
# seq2 = c("A", "C", "G", "A", "A", "C", "G", "T", "A", "C", "G", "T")
# 
# distSeqM1N(seq1, seq2)
distSeqM1N <- function(seq1, seq2, normalize=c("none" ,"length", "mutations")) {
  # Evaluate choices
  normalize <- match.arg(normalize)
  
  # Compute length of sequence (for normalization, if specified)
  juncLength <- length(seq1)
  
  # Compute distance only on positions that have mutations
  withMu <- seq1 != seq2
  withNonNuc <- ( !is.na(match(seq1, c("A","C","G","T"))) & !is.na(match(seq2, c("A","C","G","T"))) )
  withMu <- withMu & withNonNuc
  seq1 <- seq1[withMu]
  seq2 <- seq2[withMu]
  
  # Number of mutations (for normalization, if specified)
  numbOfMutation <- sum(withMu)
  
  dist <- NA
  a <- tryCatch({
    if (length(seq1)==1){
      dist <- M1NDistance[seq2,seq1]
    } else {
      dist <- sum( diag(M1NDistance[seq2,seq1]) )
    }
  },error = function(e){
    return(NA)
  })
  
  # Normalize distances
  if (normalize == "length") { 
    dist <- dist/juncLength
  } else if (normalize == "mutations") { 
    dist <- dist/numbOfMutation 
  }
  
  return(dist)
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1          first nucleotide sequence.
# @param    seq2          second nucleotide sequence.
# @param    model         DNA (ham) or amino acid (aa) hamming distance model
# @param    normalize     The method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   distance between two sequences.
#
# @examples
# seq1 = c("ATG-C")
# seq2 = c("AT--C")
# 
# distSeqHam(seq1, seq2)
distSeqHam <- function(seq1, seq2, model=c("ham","aa"),
                       normalize=c("none" ,"length", "mutations")) {
  # Evaluate choices
  model <- match.arg(model)
  normalize <- match.arg(normalize)
  
  # Get character distance matrix
  if (model == "ham") {
    # Calculate distance
    dist_mat <- getDNADistMatrix(gap=0)
    dist <- getSeqDistance(seq1, seq2, dist_mat=dist_mat)
    
    # Compute length of sequence (for normalization, if specified)
    juncLength <- nchar(seq1)
  } else if (model == "aa") {
    
    # Translate sequences
    seq1 <- strsplit(tolower(gsub("[-.]","N",seq1)), "")[[1]]
    seq2 <- strsplit(tolower(gsub("[-.]","N",seq2)), "")[[1]]
    aa1 <- translate(seq1, ambiguous=T)
    aa2 <- translate(seq2, ambiguous=T)
    
    # Calculate distance
    dist_mat <- getAADistMatrix()
    # Calculate distance
    dist <- getSeqDistance(aa1, aa2, dist_mat=dist_mat)
    
    # Compute length of sequence (for normalization, if specified)
    juncLength <- length(aa1)
  }
  
  # Normalize distances
  if (normalize == "length") { 
    dist <- dist/juncLength
  } else if (normalize == "mutations") { 
    numbOfMutation <- dist
    dist <- dist/numbOfMutation 
  }
  
  return(dist)
}


# Given an array of nucleotide sequences, find the pairwise distances
# 
# @param   arrJunctions   character vector of nucleotide sequences.
# @param   model          SHM targeting model.
# @param   normalize      The method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
# @return  A matrix of pairwise distances between junction sequences.
# 
# @details
# needs method details
# 
# @seealso needs links
# 
# @examples
# # TODO
# # working example
getPairwiseDistances <- function(arrJunctions, targeting_model, 
                                 normalize=c("none" ,"length", "mutations")) {
  # Initial checks
  normalize <- match.arg(normalize)
  
  # Convert junctions to uppercase
  arrJunctions <- toupper(arrJunctions)
  # Convert gaps to Ns
  arrJunctions <- gsub('[-.]', 'N', arrJunctions, fixed=T)
  # Add 'NN' to front and end of each sequence for fivemers
  arrJunctions <- as.vector(sapply(arrJunctions, function(x){ paste("NN", x, "NN", sep="") }))

  numbOfJunctions<-length(arrJunctions)

  #Junctions are broken in to 5-mers based on a sliding window (of one) and placed in matrix
  #Each column is a junction
  #E.g. junctions 1234567, ABCDEFG, JKLMNOP becomes:
  # 12345   ABCDE   JKLMN
  # 23456   BCDEF   KLMNO
  # 34567   CDEFG   LMNOP
  .matSeqSlidingFiveMer <- sapply(arrJunctions, function(x) { slidingArrayOf5mers(x) }, simplify="matrix")

  # Compute pairwise distance between all sequences' fivemers (by column)
  matDistance <-
    sapply(1:numbOfJunctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfJunctions, function(j) {
      distSeq5mers(.matSeqSlidingFiveMer[,i],
                   .matSeqSlidingFiveMer[,j],
                   targeting_model,
                   normalize=normalize)
    })))
  # Make distance matrix symmetric
  matDistance <- matDistance + t(matDistance)
  return(matDistance)
}


# Given an array of junctions, generate distance array for pairwise distances with
# 0 if junction is non-unique and return this array along with unique junctions that 
# are only observed once
findUniqueJunctions <- function(arrJunctions) {
  # Initialize array of distances
  arrJunctionsDist <- rep(NA,length(arrJunctions))
  
  # Filter unique junctions
  arrJunctionsUnique <- unique(arrJunctions)
  
  # Map indices of unique to its non-unique in the original arrJunctions
  indexJunctions <- match(arrJunctions, arrJunctionsUnique)
  
  # Identify junctions with multiple non-unique sequences and set its distances to 0
  indexJunctionsCounts <- table(indexJunctions)
  indexRepeated <- as.numeric(names(indexJunctionsCounts)[indexJunctionsCounts>1])
  indexRepeated <- indexJunctions%in%indexRepeated
  arrJunctionsDist[ indexRepeated ] <- 0
  names(arrJunctionsDist) <- arrJunctions
  
  # Subset unique junctions to those that are only observed once
  arrJunctionsUnique <- arrJunctionsUnique[indexJunctionsCounts==1]
  
  return(list('arrJunctionsDist'=arrJunctionsDist, 
              'arrJunctionsUnique'=arrJunctionsUnique))
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param    arrJunctions  character vector of junction sequences.
# @param    targeting_model     targeting model
# @param    normalize     method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   A vector of distances to the closest sequence.
# @examples
# arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#                    "ACGAACGTATGT", "ACGAACGTATGC",
#                    "ACGAACGTATCC","AAAAAAAAAAAA")
# getClosestBy5mers(arrJunctions, HS5FModel, normalize="none" )
getClosestBy5mers <- function(arrJunctions, targeting_model, 
                                 normalize=c("none" ,"length", "mutations")) {
  # Initial checks
  normalize <- match.arg(normalize)
  
  # Find unique sequences and return distance array with mapping
  l <- findUniqueJunctions(arrJunctions)
  arrJunctionsDist <- l$arrJunctionsDist
  arrJunctionsUnique <- l$arrJunctionsUnique

  # Compute distances between junctions
  numbOfUniqueJunctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJunctions)
  if (numbOfUniqueJunctions>1){
    # Calculate symmetric distance matrix
    matDistance <- getPairwiseDistances(arrJunctionsUnique, targeting_model, normalize)
    # Find minimum distance for each sequence
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJunctions, function(i){ min(matDistance[-i,i]) })
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }

  # Fill the distances for unique sequences
  arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist
  return(round(arrJunctionsDist,4))
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param    arrJunctions  character vector of junction sequences.
# @param    normalize     method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   A vector of distances to the closest sequence.
# @examples
# arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#                    "ACGAACGTATGT", "ACGAACGTATGC",
#                    "ACGAACGTATCC","AAAAAAAAAAAA")
# getClosestM1N(arrJunctions, normalize="none" )
getClosestM1N <- function(arrJunctions, normalize=c("none" ,"length", "mutations")) {
  # Initial checks
  normalize <- match.arg(normalize)
  
  # Find unique sequences and return distance array with mapping
  l <- findUniqueJunctions(arrJunctions)
  arrJunctionsDist <- l$arrJunctionsDist
  arrJunctionsUnique <- l$arrJunctionsUnique
  
  # Compute distances between junctions
  numbOfUniqueJunctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJunctions)
  if (numbOfUniqueJunctions>1){
    # Calculate symmetric distance matrix
    charDf <- ldply(strsplit(arrJunctionsUnique, ''))
    matDistance <-
      sapply(1:numbOfUniqueJunctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfUniqueJunctions, function(j) {
        distSeqM1N(charDf[i,], charDf[j,], normalize=normalize)
      })))
    matDistance <- matDistance + t(matDistance)
    # Find minimum distance for each sequence
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJunctions, function(i){ min(matDistance[-i,i]) })
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }
  
  # Fill the distances for unique sequences
  arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist
  return(round(arrJunctionsDist,4))
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param    arrJunctions  character vector of junction sequences.
# @param    model         DNA (ham) or amino acid (aa) hamming distance model
# @param    normalize     method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   A vector of distances to the closest sequence.
# @examples
# arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#                    "ACGAACGTATGT", "ACGAACGTATGC",
#                    "ACGAACGTATCC","AAAAAAAAAAAA")
# getClosestM1N(arrJunctions, normalize="none" )
getClosestHam <- function(arrJunctions, model=c("ham","aa"),
                          normalize=c("none" ,"length", "mutations")) {
  # Initial checks
  model <- match.arg(model)
  normalize <- match.arg(normalize)
  
  # Find unique sequences and return distance array with mapping
  l <- findUniqueJunctions(arrJunctions)
  arrJunctionsDist <- l$arrJunctionsDist
  arrJunctionsUnique <- l$arrJunctionsUnique
  
  # Compute distances between junctions
  numbOfUniqueJunctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJunctions)
  if (numbOfUniqueJunctions>1){
    # Calculate symmetric distance matrix
    matDistance <-
      sapply(1:numbOfUniqueJunctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfUniqueJunctions, function(j) {
        distSeqHam(arrJunctionsUnique[i], arrJunctionsUnique[j], model=model, normalize=normalize)
      })))
    matDistance <- matDistance + t(matDistance)
    # Find minimum distance for each sequence
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJunctions, function(i){ min(matDistance[-i,i]) })
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }
  
  # Fill the distances for unique sequences
  arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist
  return(round(arrJunctionsDist,4))
}

#' Distance to nearest neighbor
#'
#' Get distance of every sequence to its nearest sequence sharing same V gene, J gene, and
#' sequence length.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing nucleotide sequences to compare. 
#'                           Also used to determine sequence length for grouping.
#' @param    vCallColumn     name of the column containing the V-segment allele calls.
#' @param    jCallColumn     name of the column containing the J-segment allele calls.
#' @param    model           underlying SHM model, which must be one of 
#'                           \code{c("m1n", "ham", "aa", "m3n", "hs5f")}.
#'                           See Details for further information.
#' @param    normalize       method of normalization. The default is "none". If the "length" 
#'                           method is chosen, then distance is divided by the length of the 
#'                           junction sequence.
#' @param    first           if \code{TRUE} only the first call of the gene assignments is used.
#'                           If \code{FALSE} the union of ambiguous gene assignments is used to 
#'                           group all sequences with any overlapping gene calls.
#' @param    nproc           number of cores to distribute the function over.
#'
#' @return   Returns a modified \code{db} data.frame with nearest neighbor distances in the 
#'           \code{DIST_NEAREST} column.
#'
#' @details
#' The distance to nearest neighbor can be used to estimate a threshold for assigning Ig
#' sequences to clonal groups. A histogram of the resulting vector is often bimodal, 
#' with the ideal threshold being a value that separates the two modes.
#' 
#' "hs5f" and "m3n" use distance derived from the \link{HS5FModel} and \link{M3NModel} respectively 
#' using \link{calcTargetingDistance}. "m1n" uses \link{M1NDistance} to calculate distances.
#' "ham" uses a nucleotide hamming distance matrix from \link{getDNADistMatrix}, with 
#' gaps being zero. "aa" uses an amino acid hamming distance matrix from \link{getAADistMatrix}.
#' 
#' @references
#' \enumerate{
#'   \item  Smith DS, et al. Di- and trinucleotide target preferences of somatic 
#'            mutagenesis in normal and autoreactive B cells. 
#'            J Immunol. 1996 156:2642-52. 
#'   \item  Glanville J, Kuo TC, von Budingen H-C, et al. 
#'            Naive antibody gene-segment frequencies are heritable and unaltered by 
#'            chronic lymphocyte ablation. 
#'            Proc Natl Acad Sci USA. 2011 108(50):20066-71.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4:358.
#'  }
#'  
#' @seealso  See \link{calcTargetingDistance} for generating nucleotide distance matrices 
#'           from a \link{TargetingModel} object. See \link{M1NDistance}, \link{M3NModel}, 
#'           \link{HS5FModel}, \link{getDNADistMatrix}, and \link{getAADistMatrix}
#'           for individual model details.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#'
#' # Use genotyped V assignments and HS5F model
#' dist_hs5f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", model="hs5f", first=FALSE)
#' hist(dist_hs5f$DIST_NEAREST, breaks=100, xlim=c(0, 60))
#' 
#' # Use M1N model and normalize by junction length
#' dist_m1n <- distToNearest(db, model="m1n", first=FALSE, normalize="length")
#' hist(dist_m1n$DIST_NEAREST, breaks=25, xlim=c(0, 1))
#'
#' @export
distToNearest <- function(db, sequenceColumn="JUNCTION", vCallColumn="V_CALL", 
                          jCallColumn="J_CALL", model=c("m1n", "ham", "aa", "m3n", "hs5f"), 
                          normalize=c("none", "length"), 
                          first=TRUE, nproc=1) {
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, vCallColumn, jCallColumn))
    if (check != TRUE) { stop(check) }
    
    # Get targeting model
    if (model == "hs5f") {
        targeting_model <- HS5FModel
    } else if (model == "m3n") {
        targeting_model <- M3NModel            
    }
    
    # Parse V and J columns to get gene
    # cat("V+J Column parsing\n")
    if (first) {
        db$V <- getGene(db[, vCallColumn])
        db$J <- getGene(db[, jCallColumn])
    } else {
        db$V1 <- getGene(db[, vCallColumn], first=FALSE)
        db$J1 <- getGene(db[, jCallColumn], first=FALSE)
        db$V <- db$V1
        db$J <- db$J1
        # Reassign V genes to most general group of genes
        for(ambig in unique(db$V1[grepl(',', db$V1)])) {
            for(g in strsplit(ambig, split=',')[[1]]) {
                db$V[grepl(g, db$V1)] = ambig
            }
        }
        # Reassign J genes to most general group of genes
        for(ambig in unique(db$J1[grepl(',',db$J1)])) {
            for(g in strsplit(ambig, split=',')) {
                db$J[grepl(g, db$J1)] = ambig
            }
        }
    }
    
    # Create new column for distance to nearest neighbor
    db$DIST_NEAREST <- rep(NA, nrow(db))
    db$ROW_ID <- 1:nrow(db)
    db$L <- nchar(db[, sequenceColumn])
    
    # Create cluster of nproc size and export namespaces
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else {
        if(nproc != 0) { cluster <- makeCluster(nproc, type="SOCK") }
        cluster <- makeCluster(nproc, type = "SOCK")
        registerDoSNOW(cluster)
        clusterEvalQ(cluster, library(shm))
    }
    
    # Calculate distance to nearest neighbor
    # cat("Calculating distance to nearest neighbor\n")
    
    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of V J L, instead of doing dplyr
    dt <- data.table(db)
    # Get the group indexes
    dt <- dt[, list( yidx = list(.I) ) , by = list(V,J,L) ]
    groups <- dt[,yidx]
    lenGroups <- length(groups)
    
    # Export groups to the clusters
    if (nproc>1) { clusterExport(cluster, list("db", 
                                               "groups", 
                                               "sequenceColumn"), envir=environment()) }
    
    if (model %in% c("hs5f", "m3n")) {
        # Export targeting model to processes
        if (nproc>1) { clusterExport(cluster, list("targeting_model"), envir=environment()) }    
        list_db <-
            foreach(i=icount(lenGroups), .errorhandling='pass') %dopar% {
                db_group <- db[groups[[i]],]
                db_group$DIST_NEAREST <-
                    getClosestBy5mers( db[groups[[i]],sequenceColumn],
                                       targeting_model=targeting_model,
                                       normalize=normalize )
                return(db_group)
            }    
    } else if (model == "m1n") {    
        list_db <-
            foreach(i=icount(lenGroups), .errorhandling='pass') %dopar% {
                db_group <- db[groups[[i]],]
                db_group$DIST_NEAREST <-
                    getClosestM1N( db[groups[[i]],sequenceColumn],
                                   normalize=normalize )
                return(db_group)
            } 
    } else if (model %in% c("ham", "aa")) {    
        list_db <-
            foreach(i=icount(lenGroups), .errorhandling='pass') %dopar% {
                db_group <- db[groups[[i]],]
                db_group$DIST_NEAREST <-
                    getClosestHam( db[groups[[i]],sequenceColumn],
                                   normalize=normalize )
                return(db_group)
            }        
    }
    
    # Convert list from foreach into a db data.frame
    db <- do.call(plyr::rbind.fill, list_db)
    db <- db[order(db[,"ROW_ID"]),]
    
    # Stop the cluster
    if( nproc>1) { stopCluster(cluster) }
    
    return(db[, !(names(db) %in% c("V", "J", "L", "ROW_ID", "V1", "J1"))])
}
