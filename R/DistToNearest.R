# Generates distance to nearest neighbor by using S5F mutability model
#
# @author     Namita Gupta, Gur Yaari, Mohamed Uduman
# @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @date       2015.03.08


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param   seq1          first nucleotide sequence.
# @param   seq2          second nucleotide sequence.
# @param   sub_model     substitution model.
# @param   mut_model     mutability model.
# @param   normalize     The method of normalization. Default is "none".
#                        "length" = normalize distance by length of junction.
#                        "mutations" = normalize distance by number of mutations in junction.
#                        If a numeric value is passed, then the computed distance is divided by the given value.
# @return  distance between two sequences.
#
# @examples
# seq1 = c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
# seq2 = c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
#
# distSeqFast(seq1, seq2, HS5FModel@substitution, HS5FModel@mutability)
distSeqFast <- function(seq1, seq2, sub_model, mut_model, 
                        normalize=c("none" ,"length", "mutations")) {
  # Evaluate choices
  normalize <- match.arg(normalize)
  
  #Compute length of sequence (for normalization, if specified)
  juncLength <- length(seq1)

  #Compute distance only on fivemers that have mutations
  fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
  fivemersWithNonNuc <- ( !is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & !is.na(match(substr(seq2,3,3),c("A","C","G","T"))) )
  fivemersWithMu <- fivemersWithMu & fivemersWithNonNuc
  #Number of mutations (for normalization, if specified)
  numbOfMutation <- sum(fivemersWithMu)

  seq1 <- seq1[fivemersWithMu & fivemersWithNonNuc]
  seq2 <- seq2[fivemersWithMu & fivemersWithNonNuc]

  avgDist <- NA
  a <- tryCatch({
    if(length(seq1)==1){
      seq1_to_seq2 <- sub_model[substr(seq2,3,3),seq1] * mut_model[seq1]
      seq2_to_seq1 <- sub_model[substr(seq1,3,3),seq2] * mut_model[seq2]
    }else{
      seq1_to_seq2 <- sum( diag(sub_model[substr(seq2,3,3),seq1]) *  mut_model[seq1] )
      seq2_to_seq1 <- sum( diag(sub_model[substr(seq1,3,3),seq2]) *  mut_model[seq2] )
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


#' Given an array of junction sequences, find the pairwise distances
#'
#' @param   arrJunctions   character vector of junction sequences.
#' @param   model          name of SHM targeting model.
#' @param   normalize     The method of normalization. Default is "none".
#'                        "length" = normalize distance by length of junction.
#'                        "mutations" = normalize distance by number of mutations in junction.
#'                        If a numeric value is passed, then the computed distance is divided by the given value.
#' @return  A matrix of pairwise distances between junction sequences.
#' 
#' @details
#' needs method details
#' 
#' @seealso needs links
#' 
#' @examples
#' # TODO
#' # working example
#' 
#' @export
getPairwiseDistances <- function(arrJunctions, model=c("hs5f", "m3n"), 
                                 normalize=c("none" ,"length", "mutations")) {
  # Initial checks
  model <- match.arg(model)
  normalize <- match.arg(normalize)
    
  # Define targeting model
  if (model == "hs5f") {
      mut_model <- HS5FModel@mutability
      sub_model <- HS5FModel@substitution
  } else if (model == "m3n") {
      mut_model <- M3NModel@mutability
      sub_model <- M3NModel@substitution              
  }
  
  # Convert junctions to uppercase
  arrJunctions <- toupper(arrJunctions)
  # Convert gaps to Ns
  arrJunctions <- gsub('.', 'N', arrJunctions, fixed=T)
  # Add 'NN' to front and end of each sequence for fivemers
  arrJunctions <- as.vector(sapply(arrJunctions, function(x){ paste("NN", x, "NN", sep="") }))

  numbOfJuctions<-length(arrJunctions)

  #Junctions are broken in to 5-mers based on a sliding window (of one) and placed in matrix
  #Each column is a junction
  #E.g. junctions 1234567, ABCDEFG, JKLMNOP becomes:
  # 12345   ABCDE   JKLMN
  # 23456   BCDEF   KLMNO
  # 34567   CDEFG   LMNOP
  matSeqSlidingFiveMer <- sapply(arrJunctions, function(x) { slidingArrayOf5mers(x) }, simplify="matrix")

  # Compute pairwise distance between all sequences' fivemers (by column)
  matDistance <-
    sapply(1:numbOfJuctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfJuctions, function(j) {
      distSeqFast(matSeqSlidingFiveMer[,i],
                  matSeqSlidingFiveMer[,j],
                  sub_model,
                  mut_model,
                  normalize=normalize)
    })))
  # Make distance matrix symmetric
  matDistance <- matDistance + t(matDistance)
  return(matDistance)
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param   arrJunctions  character vector of junction sequences.
# @param   sub_model          substitution model.
# @param   mut_model           mutability model.
# @param   normalize     The method of normalization. Default is "none".
#                        "length" = normalize distance by length of junction.
#                        "mutations" = normalize distance by number of mutations in junction.
#                        If a numeric value is passed, then the computed distance is divided by the given value.
# @return  A vector of distances to the closest sequence.
# @examples
# arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#                    "ACGAACGTATGT", "ACGAACGTATGC",
#                    "ACGAACGTATCC","AAAAAAAAAAAA")
# getDistanceToClosest(arrJunctions, HS5FModel@substitution, HS5FModel@mutability, normalize="none" )
getDistanceToClosest <- function(arrJunctions, sub_model, mut_model, 
                                 normalize=c("none" ,"length", "mutations")) {
  # Initial checks
  normalize <- match.arg(normalize)
  
  #Initialize array of distances
  arrJunctionsDist <- rep(NA,length(arrJunctions))

  #Filter unique junctions
  arrJunctionsUnique <- unique(arrJunctions)

  #Map indexes of unique to its non-unique in the original arrJunctions
  indexJunctions <- match(arrJunctions, arrJunctionsUnique)

  #Identify junctions with multiple non-unique sequences and set its distances to 0
  indexJunctionsCounts <- table(indexJunctions)
  indexRepeated <- as.numeric(names(indexJunctionsCounts)[indexJunctionsCounts>1])
  indexRepeated <- indexJunctions%in%indexRepeated
  arrJunctionsDist[ indexRepeated ] <- rep(0,sum(indexRepeated))
  names(arrJunctionsDist) <- arrJunctions

  #Compute distances between junctions
  numbOfUniqueJuctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJuctions)
  if(numbOfUniqueJuctions>1){
    arrJunctionsUnique <- toupper(arrJunctionsUnique)
    arrJunctionsUnique <- gsub('.', '-', arrJunctionsUnique, fixed=TRUE)
    arrJunctionsUnique <- as.vector(sapply(arrJunctionsUnique,function(x){paste("NN",x,"NN",sep="")}))

    #Junctions are broken in to 5-mers based on a sliding window (of one) and placed in matrix
    #Each column is a junction
    #E.g. junctions 1234567, ABCDEFG, JKLMNOP becomes:
    # 12345   ABCDE   JKLMN
    # 23456   BCDEF   KLMNO
    # 34567   CDEFG   LMNOP
    matSeqSlidingFiveMer <- sapply(arrJunctionsUnique, function(x) { slidingArrayOf5mers(x) },simplify="matrix")

    # Compute pairwise distance between all sequences' fivemers (by column)
    matDistance <-
      sapply(1:numbOfUniqueJuctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfUniqueJuctions, function(j) {
        distSeqFast( matSeqSlidingFiveMer[,i],
                       matSeqSlidingFiveMer[,j],
                       sub_model,
                       mut_model,
                       normalize=normalize)
      })))

    matDistance <- matDistance + t(matDistance)
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJuctions, function(i){ min(matDistance[-i,i]) })
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }

  #Fill the distances for the sequences that are unique
  arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist[indexJunctionsCounts==1]
  return(round(arrJunctionsDist,4))
}


#' Distance to nearest neighbor
#'
#' Get distance of every sequence to its nearest sequence sharing same V gene, J gene, and
#' sequence length.
#'
#' @param   db          \code{data.frame} which must have the following columns: V_CALL and J_CALL.
#' @param   seq         the column containing nucleotide sequences to compare. Also used to determine
#'                      sequence length for grouping.
#' @param   genotyped   logical indicating whether \code{db} is genotyped; if genotyped is \code{TRUE},
#'                      \code{db} must have the column V_CALL_GENOTYPED.
#' @param   model       SHM targeting model; must be one of c("hs5f", "m3n"). See Details for further information.
#' @param   normalize   The method of normalization. Default is "none".
#'                      "length" = normalize distance by length of junction.
#'                      "mutations" = normalize distance by number of mutations in junction.
#'                      If a numeric value is passed, then the computed distance is divided by the given value.
#' @param   first       if \code{TRUE} only the first call the gene assignment is used;
#'                      if \code{FALSE} the union of ambiguous gene assignments is used to group all sequences with
#'                      any of those gene calls.
#' @param   nproc       The number of cores to distribute the function over.
#'
#' @return  Returns a modified \code{db} data.frame with nearest neighbor distances in the 
#'          DIST_NEAREST column.
#'
#' @details
#' Needs method details.
#' 
#' @references
#' \enumerate{
#'   \item  Smith DS, et al. Di- and trinucleotide target preferences of somatic 
#'            mutagenesis in normal and autoreactive B cells. 
#'            J Immunol. 1996 156:2642–52. 
#'   \item  Glanville J, Kuo TC, von Büdingen H-C, et al. Naive antibody gene-segment 
#'            frequencies are heritable and unaltered by chronic lymphocyte ablation. 
#'            Proc Natl Acad Sci USA. 2011 108(50):20066–71.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'  
#' @seealso needs links
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#'
#' # Calculate distance to nearest using genotyped V assignments and 5-mer model
#' dist <- distToNearest(db, vcall="V_CALL_GENOTYPED", model="hs5f", first=FALSE)
#' hist(dist$DIST_NEAREST, breaks=50, xlim=c(0, 15))
#'
#' @export
distToNearest <- function(db, seq="JUNCTION", vcall="V_CALL", jcall="J_CALL", 
                          model=c("hs5f", "m3n"), 
                          normalize=c("none" ,"length", "mutations"), 
                          first=TRUE, nproc=1) {
  # Initial checks
  model <- match.arg(model)
  normalize <- match.arg(normalize)
  
  if(!is.data.frame(db)) { stop('Must submit a data frame') }

  if (model == "hs5f") {
      mut_model <- HS5FModel@mutability
      sub_model <- HS5FModel@substitution
  } else if (model == "m3n") {
      mut_model <- M3NModel@mutability
      sub_model <- M3NModel@substitution              
  }

  # Parse V and J columns to get gene
  # cat("V+J Column parsing\n")
  if(first) {
    db$V <- getGene(db[, vcall])
    db$J <- getGene(db[, jcall])
  } else {
    db$V <- getGene(db[, vcall], first=FALSE)
    db$J <- getGene(db[, jcall], first=FALSE)
    # Reassign V genes to most general group of genes
    for(ambig in unique(db$V[grepl(',', db$V)])) {
      for(g in strsplit(ambig, split=',')) {
        db$V2[db$V==g] = ambig
      }
    }
    # Reassign J genes to most general group of genes
    for(ambig in unique(db$J[grepl(',',db$J)])) {
      for(g in strsplit(ambig, split=',')) {
        db$J[db$J==g] = ambig
      }
    }
  }

  # Create new column for distance to nearest neighbor
  db$DIST_NEAREST <- rep(NA, nrow(db))
  db$ROW_ID <- 1:nrow(db)
  db$L <- nchar(db[, seq])
  
  # Create cluster of nproc size and export namespaces
  cluster <- makeCluster(nproc, type = "SOCK")
  registerDoSNOW(cluster)
  clusterExport(cluster, list("sub_model", "mut_model"), envir=environment())
  clusterEvalQ(cluster, library(shm))
  
  # Calculate distance to nearest neighbor
  # cat("Calculating distance to nearest neighbor\n")
  db <- arrange(ddply(db, .(V, J, L), function(piece) 
    mutate(piece, DIST_NEAREST=getDistanceToClosest(eval(parse(text=seq)),
                                                    sub_model=sub_model,
                                                    mut_model=mut_model,
                                                    normalize=normalize)),
    .parallel=TRUE),
    ROW_ID)
  
  # Stop the cluster
  stopCluster(cluster)
  
  return(db[, !(names(db) %in% c("V", "J", "L", "ROW_ID"))])
}
