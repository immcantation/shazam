# Generates distance to nearest neighbor by using S5F mutability model
#
# @author     Namita Gupta, Gur Yaari, Mohamed Uduman
# @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @date       2014.11.24


# Load targeting model data from the shm package
#
# @param   model   string defining name of the model to load.
#                  One of "hs5f" or "m3n".
# @return  a list containing the substitution (subs) and mutability (mut) models
#' @export
loadModel <- function(model) {
  if(model == "hs5f") { data_file = "HS5F_Targeting.RData" }
  else if(model == "m3n") { data_file = "MTri_Targeting.RData" }
  else { stop("Are you sure you know what you're doing?\n") }

  tmp_env <- new.env()
  load(system.file("extdata", data_file, package="shm"), envir=tmp_env)
  Targeting <- get("Targeting", envir=tmp_env)
  rm(tmp_env)

  return(list(subs=Targeting[["Substitution"]], mut=Targeting[["Mutability"]]))
}


#' Get distance between two sequences of same length, broken by a sliding window of 5mers
#'
#' @param   seq1          first nucleotide sequence.
#' @param   seq2          second nucleotide sequence.
#' @param   subs          substitution model.
#' @param   mut           mutability model.
#' @param   normalize     The method of normalization. Default is none.
#'                        'juncLength' = normalize distance by length of junction.
#'                        'juncMutation' = normalize distance by number of mutations in junction.
#'                        If a numeric value is passed, then the computed distance is divided by the given value.
#' @return  distance between two sequences.
#'
#' @examples
#' seq1 = c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
#' seq2 = c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
#'
#' model="hs5f"
#' model_data <- loadModel(model)
#' dist_seq_fast(seq1, seq2, model_data[["subs"]], model_data[["mut"]], normalize="none")
#' @export
dist_seq_fast <- function(seq1, seq2, subs, mut, normalize="none") {
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
      seq1_to_seq2 <- subs[substr(seq2,3,3),seq1] * mut[seq1]
      seq2_to_seq1 <- subs[substr(seq1,3,3),seq2] * mut[seq2]
    }else{
      seq1_to_seq2 <- sum( diag(subs[substr(seq2,3,3),seq1]) *  mut[seq1] )
      seq2_to_seq1 <- sum( diag(subs[substr(seq1,3,3),seq2]) *  mut[seq2] )
    }
    avgDist <- mean(c(seq1_to_seq2, seq2_to_seq1))
  },error = function(e){
    return(NA)
  })

  #Normalizing distance
  if(normalize!="none"){
    if(normalize=="juncLength")avgDist<-avgDist/juncLength
    if(normalize=="juncMutation")avgDist<-avgDist/numbOfMutation
    if(is.numeric(normalize))avgDist<-avgDist/normalize
  }
  return(avgDist)
}


#' Given an array of junction sequences, find the pairwise distances
#'
#' @param   arrJunctions   character vector of junction sequences.
#' @param   model          name of SHM targeting model.
#' @param   normalize     The method of normalization. Default is none.
#'                        'juncLength' = normalize distance by length of junction.
#'                        'juncMutations' = normalize distance by number of mutations in junction.
#'                        If a numeric value is passed, then the computed distance is divided by the given value.
#' @return  A matrix of pairwise distances between junction sequences.
#' @export
getPairwiseDistances <- function(arrJunctions, model="hs5f", normalize="none") {

  # Load targeting model
  model_data <- loadModel(model)
  # Convert junctions to uppercase
  arrJunctions <- toupper(arrJunctions)
  # Convert gaps to Ns
  arrJunctions <- gsub('.', 'N', arrJunctions, fixed=T)
  # Add 'NN' to front and end of each sequence for fivemers
  arrJunctions <- as.vector(sapply(arrJunctions, function(x){paste("NN",x,"NN",sep="")}))

  numbOfJuctions<-length(arrJunctions)

  #Junctions are broken in to 5-mers based on a sliding window (of one) and placed in matrix
  #Each column is a junction
  #E.g. junctions 1234567, ABCDEFG, JKLMNOP becomes:
  # 12345   ABCDE   JKLMN
  # 23456   BCDEF   KLMNO
  # 34567   CDEFG   LMNOP
  matSeqSlidingFiveMer <- sapply(arrJunctions, function(x) { slidingArrayOf5mers(x) },simplify="matrix")

  # Compute pairwise distance between all sequences' fivemers (by column)
  matDistance <-
    sapply(1:numbOfJuctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfJuctions, function(j) {
      dist_seq_fast( matSeqSlidingFiveMer[,i],
                     matSeqSlidingFiveMer[,j],
                     model_data[["subs"]],
                     model_data[["mut"]],
                     normalize=normalize)
    })))
  # Make distance matrix symmetric
  matDistance <- matDistance + t(matDistance)
  return(matDistance)
}


#' Given an array of junction sequences, find the distance to the closest sequence
#'
#' @param   arrJunctions  character vector of junction sequences.
#' @param   subs          substitution model.
#' @param   mut           mutability model.
#' @param   normalize     The method of normalization. Default is none.
#'                        'juncLength' = normalize distance by length of junction.
#'                        'juncMutations' = normalize distance by number of mutations in junction.
#'                        If a numeric value is passed, then the computed distance is divided by the given value.
#' @return  A vector of distances to the closest sequence.
#' @examples
#' arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#'                    "ACGAACGTATGT", "ACGAACGTATGC",
#'                    "ACGAACGTATCC","AAAAAAAAAAAA")
#' model_data <- loadModel(model="hs5f")
#' subs <- model_data[['subs']]
#' mut <- model_data[['mut']]
#' getDistanceToClosest(arrJunctions, subs, muts, normalize="none" )
#' @export
getDistanceToClosest <- function(arrJunctions, subs, mut, normalize="none") {

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
        dist_seq_fast( matSeqSlidingFiveMer[,i],
                       matSeqSlidingFiveMer[,j],
                       model_data[["subs"]],
                       model_data[["mut"]],
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
#' hs5f model is the SHM targeting model from Yaari, G., et al. Frontiers in Immunology, 2013.
#' m3n model uses the SHM substitution matrix found in Smith, D., et al. J. Immunol., 1996.
#'
#' @param   db          \code{data.frame} which must have the following columns: V_CALL and J_CALL.
#' @param   seq         the column containing nucleotide sequences to compare. Also used to determine
#'                      sequence length for grouping.
#' @param   genotyped   logical indicating whether \code{db} is genotyped; if genotyped is \code{TRUE},
#'                      \code{db} must have the column V_CALL_GENOTYPED.
#' @param   first       if \code{TRUE} only the first call the gene assignment is used;
#'                      if \code{FALSE} the union of ambiguous gene assignments is used to group all sequences with
#'                      any of those gene calls.
#' @param   model       SHM targeting model; must be one of c("hs5f", "m3n"). See Details for further information.
#' @param   vector      if \code{TRUE} return a numeric vector of only the distances; if \code{FALSE} return the
#'                      entire input data.frame with a DIST_NEAREST column added.
#' @param   nproc       The number of cores to distribute the function over.
#' @param   normalize   The method of normalization. Default is none.
#'                      'juncLength' = normalize distance by length of junction.
#'                      'juncMutations' = normalize distance by number of mutations in junction.
#'                      If a numeric value is passed, then the computed distance is divided by the given value.
#'
#' @return  If \code{vector=TRUE} returns a numeric vector of distances of each sequence to its nearest neighbor.
#'          If \code{vector=FALSE} returns a modified \code{db} data.frame with nearest neighbor distances
#'          in the DIST_NEAREST column.
#'
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#'
#' # Calculate distance to nearest
#' dist <- distToNearest(df, genotyped=TRUE, first=FALSE, vector=TRUE)
#' hist(dist, breaks=50, xlim=c(0, 15))
#'
#' @export
distToNearest <- function(db, seq="JUNCTION", genotyped=FALSE, first=TRUE, model="hs5f",
                          vector=FALSE, nproc=1, normalize="none") {
  # Initial check
  if(!is.data.frame(db)) { stop('Must submit a data frame') }

  # Define V and J columns
  if(genotyped) {
    v_col <- "V_CALL_GENOTYPED"
  } else {
    v_col <- "V_CALL"
  }
  j_col <- "J_CALL"

  # Parse V and J columns to get gene
  # cat("V+J Column parsing\n")
  if(first) {
    db$V <- getGene(db[,v_col])
    db$J <- getGene(db[,j_col])
  } else {
    db$V <- getGene(db[,v_col], first=FALSE)
    db$J <- getGene(db[,j_col], first=FALSE)
    # Reassign V genes to most general group of genes
    for(ambig in unique(db$V[grepl(',',db$V)])) {
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

  # Load targeting model
  # cat("Loading Targeting Model\n")
  model_data <- loadModel(model)

  # Create new column for distance to nearest neighbor
  db$DIST_NEAREST <- rep(NA, nrow(db))
  db$ROW_ID <- 1:nrow(db)
  db$L <- nchar(db[, seq])
  
  # Create cluster of nproc size and export namespaces
  cluster <- makeCluster(nproc, type = "SOCK")
  registerDoSNOW(cluster)
  clusterExport(cluster, list('model_data'), envir=environment())
  clusterEvalQ(cluster, library(shm))
  
  # Calculate distance to nearest neighbor
  # cat("Calculating distance to nearest neighbor\n")
  db <- arrange(ddply(db, .(V, J, L), function(piece) 
    mutate(piece, DIST_NEAREST=getDistanceToClosest(eval(parse(text=seq)),
                                                    subs=model_data[['subs']],
                                                    mut=model_data[['mut']],
                                                    normalize=normalize)),
    .parallel=TRUE),
    ROW_ID)
  
  # Stop the cluster
  stopCluster(cluster)
  
  # Determine output
  if (vector) {
    return(db$DIST_NEAREST)
  } else {
    return(db[, !(names(db) %in% c("V", "J", "L", "ROW_ID"))])
  }
}
