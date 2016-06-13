# Generates distance to nearest neighbor

#' @include shazam.R
NULL


# Returns a 5-mer sliding window of given sequence
#
# @param   strSequence   The sequence string
# @return  An array of 5-mer sliding windows
#
# @examples
# slidingArrayOf5mers("ACGTNACGTNACGTN")
slidingArrayOf5mers <- function(strSequence){
    seqLength <- stri_length(strSequence)
    return(substr(rep(strSequence, seqLength - 4), 1:(seqLength - 4), 5:seqLength))
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1                first nucleotide sequence, broken into 5mers.
# @param    seq2                second nucleotide sequence, broken into 5mers.
# @param    targeting_model     targeting model.
# @param    normalize           The method of normalization. Default is "none".
#                               "length" = normalize distance by length of junction.
# @param    symmetry            if model is hs5f, distance between seq1 and seq2 is either 
#                               the average (avg) of seq1->seq2 and seq2->seq1 or the 
#                               minimum (min).
# @return   distance between two sequences.
#
# @examples
# seq1 <- c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", 
#          "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
# seq2 <- c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", 
#          "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
#
# shazam:::dist5Mers(seq1, seq2, HS5FModel)
dist5Mers <- function(seq1, seq2, targeting_model, 
                         normalize=c("none" ,"length", "mutations"),
                         symmetry=c("avg","min")) {
    # Evaluate choices
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    
    # Get distance from targeting model
    targeting_dist <- calcTargetingDistance(targeting_model)
    
    # Compute length of sequence (for normalization, if specified)
    juncLength <- length(seq1)
    
    # Compute distance only on fivemers that have mutations
    fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
    #fivemersWithNonNuc <- (!is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & 
    #                       !is.na(match(substr(seq2,3,3),c("A","C","G","T"))))
    #fivemersWithMu <- fivemersWithMu & fivemersWithNonNuc
    seq1 <- seq1[fivemersWithMu]
    seq2 <- seq2[fivemersWithMu]
    
    # Number of mutations (for normalization, if specified)
    numbOfMutation <- sum(fivemersWithMu)
    
    dist <- NA
    tryCatch({
        if (length(seq1)==1){
            seq1_to_seq2 <- targeting_dist[substr(seq2, 3, 3), seq1]
            seq2_to_seq1 <- targeting_dist[substr(seq1, 3, 3), seq2]
        } else {
            seq1_to_seq2 <- sum(diag(targeting_dist[substr(seq2, 3, 3), seq1]))
            seq2_to_seq1 <- sum(diag(targeting_dist[substr(seq1, 3, 3), seq2]))
        }
        if (symmetry == "avg") {
            dist <- mean(c(seq1_to_seq2, seq2_to_seq1))
        } else if (symmetry == "min") {
            dist <- min(c(seq1_to_seq2, seq2_to_seq1))
        }
    },error = function(e){
        warning(e)
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
# @param    model         DNA (ham) or amino acid (aa) hamming distance model or
#                         mouse (m1n) or human (hs1f) single nucleotide distance model
# @param    normalize     The method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   distance between two sequences.
#
# @examples
# seq <- c("ATG-C", "AT--C", "ATC-C", "ATQ-C")
# shazam:::dist1Mers(seq[1], seq[2])
# shazam:::dist1Mers(seq[1], seq[3])
# shazam:::dist1Mers(seq[1], seq[4])
dist1Mers <- function(seq1, seq2, model=c("ham","aa","m1n","hs1f"),
                       normalize=c("none" ,"length", "mutations")) {
    # TODO: No longer used.
    # Evaluate choices
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    
    # Get character distance matrix
    if (model == "ham") {
        dist_mat <- getDNAMatrix(gap=0)
    } else if (model == "m1n") {
        dist_mat <- M1NDistance
    } else if (model == "hs1f") {
        dist_mat <- HS1FDistance
    } else if (model == "aa") {
        
        # Translate sequences
        seq1 <- strsplit(tolower(gsub("[-.]", "N", seq1)), "")[[1]]
        seq2 <- strsplit(tolower(gsub("[-.]", "N", seq2)), "")[[1]]
        seq1 <- translate(seq1, ambiguous=T)
        seq2 <- translate(seq2, ambiguous=T)
        
        dist_mat <- getAAMatrix()
    }
    
    # Calculate distance
    dist <- tryCatch(seqDist(seq1, seq2, dist_mat=dist_mat),
                     error=function(e) { warning(e); return(NA) })
    
    # Normalize distances
    if (normalize == "length") { 
        dist <- dist / sum(stri_length(seq1))
    } else if (normalize == "mutations") {
        dist <- dist / sum(strsplit(seq1, "")[[1]] != strsplit(seq2, "")[[1]])
    }
    
    return(dist)
}


# Given an array of nucleotide sequences, find the pairwise distances
# 
# @param   arrJunctions   character vector of nucleotide sequences.
# @param   model          SHM targeting model.
# @param   normalize      The method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
# @param    symmetry      if model is hs5f, distance between seq1 and seq2 is either the
#                         average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
# @return  A matrix of pairwise distances between junction sequences.
# 
# @details
# needs method details
# 
# @seealso needs links
# 
# @examples
# # working example
getPairwiseDistances <- function(arrJunctions, targeting_model, 
                                 normalize=c("none" ,"length", "mutations"),
                                 symmetry=c("avg","min")) {
    # Initial checks
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    
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
    .matSeqSlidingFiveMer <- sapply(arrJunctions, function(x) { slidingArrayOf5mers(x) }, 
                                    simplify="matrix")
    
    # Compute pairwise distance between all sequences' fivemers (by column)
    matDistance <-
        sapply(1:numbOfJunctions, function(i) c(rep.int(0,i-1), 
                                                sapply(i:numbOfJunctions, function(j) {
            dist5Mers(.matSeqSlidingFiveMer[,i],
                         .matSeqSlidingFiveMer[,j],
                         targeting_model,
                         normalize=normalize,
                         symmetry=symmetry)
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
    # TODO: arrJunctionsDist is always NA? Not needed.
    arrJunctionsDist <- rep(NA, length(arrJunctions))
    
    # Filter unique junctions
    arrJunctionsUnique <- unique(arrJunctions)
    
    # Map indices of unique to its non-unique in the original arrJunctions
    indexJunctions <- match(arrJunctions, arrJunctionsUnique)
    
    # Identify junctions with multiple non-unique sequences and set its distances to 0
    indexJunctionsCounts <- table(indexJunctions)
    indexRepeated <- as.numeric(names(indexJunctionsCounts)[indexJunctionsCounts>1])
    indexRepeated <- indexJunctions%in%indexRepeated
    # arrJunctionsDist[ indexRepeated ] <- 0
    names(arrJunctionsDist) <- arrJunctions
    
    # Subset unique junctions to those that are only observed once
    #arrJunctionsUnique <- arrJunctionsUnique[indexJunctionsCounts==1]
    
    return(list('arrJunctionsDist'=arrJunctionsDist, 
                'arrJunctionsUnique'=arrJunctionsUnique,
                'indexJunctionsCounts'=indexJunctionsCounts))
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param    arrJunctions  character vector of junction sequences.
# @param    targeting_model     targeting model
# @param    normalize     method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in 
#                         junction.
# @param    symmetry      if model is hs5f, distance between seq1 and seq2 is either the
#                         average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
# @param    crossGroups   Columns for grouping to calculate distances across groups 
#                         (self vs others)
# @param    mst           if true, return comma-separated branch lengths from minimum 
#                         spanning tree
#
# @return   A vector of distances to the closest sequence.
# 
# @examples
# arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#                    "ACGAACGTATGT", "ACGAACGTATGC",
#                    "ACGAACGTATCC","AAAAAAAAAAAA")
# shazam:::getClosestBy5Mers(arrJunctions, HS5FModel, normalize="none")
getClosestBy5Mers <- function(arrJunctions, targeting_model, 
                              normalize=c("none" ,"length", "mutations"),
                              symmetry=c("avg","min"),
                              crossGroups=NULL, mst=FALSE) {
    # Initial checks
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    
    ## If crossGroup requested, but only one group found, return NA
    if (!is.null(crossGroups) & length(unique(crossGroups))<2) {
        arrJunctionsDist <- rep(NA, length(arrJunctions))
        return (arrJunctionsDist)
    }
    
    # Find unique sequences and return distance array with mapping
    l <- findUniqueJunctions(arrJunctions)
    arrJunctionsDist <- l$arrJunctionsDist
    arrJunctionsUnique <- l$arrJunctionsUnique
    # indexJunctionsCounts <- l$indexJunctionsCounts
    
    # Compute distances between junctions
    numbOfUniqueJunctions <- length(arrJunctionsUnique)
    arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJunctions)
    if (numbOfUniqueJunctions>1){
        # Calculate symmetric distance matrix
        matDistance <- getPairwiseDistances(arrJunctionsUnique, targeting_model, 
                                            normalize, symmetry)
        # Find minimum distance for each sequence
        if (is.null(crossGroups)) {
            if(!mst) {
                arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJunctions, function(i) { 
                    ## Return smaller value greater than 0
                    ## If all 0, return NA
                    gt0 <- which(matDistance[,i]>0)
                    if (length(gt0)==0) return (NA)
                    min(matDistance[,i][gt0]) 
                    min(matDistance[,i][matDistance[,i]>0]) 
                })
                names(arrUniqueJunctionsDist) <- arrJunctionsUnique
                arrUniqueJunctionsDist <- round(arrUniqueJunctionsDist, 4)
            } else {
                # Get adjacency matrix of minimum spanning tree
                adj <- ape::mst(matDistance)
                arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJunctions, 
                                                 function(i) { 
                                                     # Get value(s) from mst branches
                                                     gt0 <- which(adj[,i]==1)
                                                     # If none (broken mst!), return NA
                                                     if (length(gt0)==0) return (NA)
                                                     # If multiple values, comma-join
                                                     paste0(round(matDistance[,i][gt0], 4), 
                                                            collapse=",")
                                                 })
                names(arrUniqueJunctionsDist) <- arrJunctionsUnique
            }
        } else {
            arrJunctionsDist <- sapply(1:length(arrJunctions),function(i) {
                # Identify sequences to be considered when finding minimum
                # cross distance
                thisGroup <- crossGroups[i]
                otherGroups <-  which(crossGroups!=thisGroup)
                otherJunctions <- unique(arrJunctions[otherGroups])
                otherJunctionsUniqueIdx <- match(otherJunctions,arrJunctionsUnique)
                thisJunctionUniqueIdx <- match(arrJunctions[i],arrJunctionsUnique)
                matDistanceRow <- matDistance[thisJunctionUniqueIdx,otherJunctionsUniqueIdx]
                #min(range(matDistanceRow[matDistanceRow>0],finite=T))
                gt0 <- which(matDistanceRow>0)
                if (length(gt0)==0) return (NA)
                min(matDistanceRow[gt0])
            })  
            names(arrJunctionsDist) <- arrJunctions
            return (round(arrJunctionsDist,4))
        }
    }
    
    # Fill the distances for unique sequences
    # arrJunctionsDist[is.na(arrJunctionsDist)] <- 
    # arrUniqueJunctionsDist[indexJunctionsCounts==1]
    arrJunctionsDist <- arrUniqueJunctionsDist[match(names(arrJunctionsDist), 
                                                     names(arrUniqueJunctionsDist))]
    return(arrJunctionsDist)
}


# Given an array of sequences, find the distance to the closest sequence
#
# @param    sequences     character vector of sequences.
# @param    model         DNA (ham) or amino acid (aa) hamming distance model or
#                         mouse (m1n) or human (hs1f) single nucleotide distance model
# @param    normalize     method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in 
#                         junction.
# @param    crossGroups   column for grouping to calculate distances across groups 
#                         (self vs others)
# @param    mst           if true, return comma-separated branch lengths from minimum 
#                         spanning tree
#
# @return   A vector of distances to the closest sequence.
#
# @examples
# sequences <- c("ACGTACGTACGT", "ACGAACGTACGT", "ACGAACGTATGT", "ACGAACGTATGC",
#                "ACGAACGTATCC", "AAAAAAAAAAAA", "A-GAACGTATCC", "AAAAAA---AAA")
# shazam:::getClosestBy1Mers(sequences, model="ham", normalize="none")
# shazam:::getClosestBy1Mers(sequences, model="aa", normalize="none")
# shazam:::getClosestBy1Mers(sequences, model="ham", normalize="length")
# shazam:::getClosestBy1Mers(sequences, model="aa", normalize="length")
getClosestBy1Mers <- function(sequences, model=c("ham","aa","m1n","hs1f"),
                              normalize=c("none" ,"length", "mutations"),
                              crossGroups=NULL, mst=FALSE) {
    ## DEBUG
    # sequences <- c("ACGTACGTACGT", "ACGAACGTACGT", "AAAAAAAAAAAA", "A-AAAA---AAA")
    # model="aa"; normalize="length"; crossGroups=NULL; mst=FALSE
    
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    
    ## If crossGroup requested, but only one group found, return NA
    if (!is.null(crossGroups) & length(unique(crossGroups)) < 2) {
        seq_dist <- rep(NA, length(sequences))
        return (seq_dist)
    }
    
    # Find unique sequences and return distance array with mapping
    x <- findUniqueJunctions(sequences)
    ## TODO: x$arrJunctionsDist is always NAs. Seems pointless.
    seq_dist <- x$arrJunctionsDist
    seq_uniq <- x$arrJunctionsUnique
    names(seq_uniq) <- seq_uniq
    
    # Compute distances between sequences
    n_uniq <- length(seq_uniq)
    seq_uniq_dist <- rep(NA, n_uniq)
    if (n_uniq > 1) {
        if (model == "ham") {
            dist_mat <- getDNAMatrix(gap=0)
        } else if (model == "m1n") {
            dist_mat <- M1NDistance
        } else if (model == "hs1f") {
            dist_mat <- HS1FDistance
        } else if (model == "aa") {
            # Translate sequences
            seq_uniq <- setNames(alakazam::translateDNA(seq_uniq), seq_uniq)
            dist_mat <- getAAMatrix()
        }
        
        # Check for length mismatches
        seq_length <-  unique(stri_length(seq_uniq))
        if (length(seq_length) > 1) {
            stop("Unexpected. Different sequence lengths found.")
        }
        
        uniq_mat <- pairwiseDist(seq_uniq, dist_mat=dist_mat)
        
        ## DEBUG
        # cat("\n-> seq_uniq:\n")
        # print(seq_uniq)
        # cat("\n-> uniq_mat (raw):\n")
        # print(uniq_mat)

        # Normalize distances
        if (normalize == "length") { 
            uniq_mat <- uniq_mat / seq_length
        } else if (normalize == "mutations") {
            #dist <- dist/sum(strsplit(seq1,"")[[1]] != strsplit(seq2,"")[[1]])
            stop("Sorry! nomalize=mutations is not available.")
        }
        
        ## DEBUG
        # cat("\n-> seq_length:\n")
        # print(seq_length)
        # cat("\n-> uniq_mat (normalized):\n")
        # print(uniq_mat)
        
    } else {
        return(seq_dist)
    }
    
    # Find minimum distance for each sequence
    if (is.null(crossGroups)) {
        if(!mst) {
            # Return smaller value greater than 0
            # If all 0, return NA
            .dmin <- function(i) { 
                x <- uniq_mat[, i]
                gt0 <- which(x > 0)
                if (length(gt0) != 0) { min(x[gt0]) } else { NA }
            }
            
            ## TODO: Could be an apply over columns
            seq_uniq_dist <- setNames(sapply(1:n_uniq, .dmin), names(seq_uniq))
        } else {
            # Get adjacency matrix of minimum spanning tree
            adj <- ape::mst(uniq_mat)
            
            # TODO: This could be cleaner
            # Get value(s) from mst branches
            # If none (broken mst!), return NA
            # If multiple values, comma-join
            .dmst <- function(i) { 
                gt0 <- which(adj[, i] == 1)
                if (length(gt0) != 0) { 
                    stri_join(round(uniq_mat[, i][gt0], 4), collapse=",") 
                } else {
                    NA
                }
            }
            
            ## TODO: Could be an apply over columns
            seq_uniq_dist <- setNames(sapply(1:n_uniq, .dmst), names(seq_uniq))
        }
        
        # Define return distance vector
        seq_dist <- seq_uniq_dist[match(names(seq_dist), names(seq_uniq_dist))]
        
        ## DEBUG
        # cat("\n-> seq_uniq_dist:\n")
        # print(seq_uniq_dist)
        # cat("\n-> seq_dist:\n")
        # print(seq_dist)
    } else {
        # Identify sequences to be considered when finding minimum
        # cross distance
        .dcross <- function(i) {
            this_group <- crossGroups[i]
            other_groups <-  which(crossGroups != this_group)
            other_seq <- unique(sequences[other_groups])
            other_idx <- match(other_seq, seq_uniq)
            this_idx <- match(sequences[i], seq_uniq)
            r <- uniq_mat[this_idx, other_idx]
            gt0 <- which(r > 0)
            
            if (length(gt0) != 0) { min(r[gt0]) } else { NA }
        }
        
        # Define return distance vector
        seq_dist <- setNames(sapply(1:length(sequences), .dcross), sequences)
    }
    
    return(round(seq_dist, 4))
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
#'                           \code{c("m1n", "ham", "aa", "hs5f")}.
#'                           See Details for further information.
#' @param    normalize       method of normalization. The default is \code{"length"}, which 
#'                           divides the distance by the length of the sequence group. If 
#'                           \code{"none"} then no normalization if performed
#' @param    symmetry        if model is hs5f, distance between seq1 and seq2 is either the
#'                           average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
#' @param    first           if \code{TRUE} only the first call of the gene assignments 
#'                           is used. if \code{FALSE} the union of ambiguous gene 
#'                           assignments is used to group all sequences with any 
#'                           overlapping gene calls.
#' @param    nproc           number of cores to distribute the function over.
#' @param    fields          additional fields to use for grouping
#' @param    cross           columns for grouping to calculate distances across groups 
#'                           (self vs others)
#' @param    mst             if true, return comma-separated branch lengths from minimum 
#'                           spanning tree
#'
#' @return   Returns a modified \code{db} data.frame with nearest neighbor distances in the 
#'           \code{DIST_NEAREST} column if \code{crossGrups=NULL} or in the 
#'           \code{CROSS_DIST_NEAREST} column if \code{crossGroups} was specified.
#'
#' @details
#' The distance to nearest neighbor can be used to estimate a threshold for assigning Ig
#' sequences to clonal groups. A histogram of the resulting vector is often bimodal, 
#' with the ideal threshold being a value that separates the two modes.
#' 
#' "hs5f" use distance derived from the \link{HS5FModel}
#' using \link{calcTargetingDistance}. "hs1f" and "m1n" use \link{HS1FDistance} and 
#' \link{M1NDistance} to calculate distances respectively. "ham" uses a nucleotide 
#' hamming distance matrix from \link[alakazam]{getDNAMatrix}, with gaps being zero. 
#' "aa" uses an amino acid hamming distance matrix from \link[alakazam]{getAAMatrix}.
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
#'           from a \link{TargetingModel} object. See \link{M1NDistance}, 
#'           \link{HS5FModel}, \link[alakazam]{getDNAMatrix}, and \link[alakazam]{getAAMatrix}
#'           for individual model details.
#' 
#' @examples
#' # Subset data for demo purposes
#' db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
#'              BARCODE %in% c("RL016","RL018","RL019","RL021"))
#' 
#' # Use genotyped V assignments, HS1F model, and normalize by junction length
#' dist_hs1f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
#'                            model="hs1f", first=FALSE, normalize="length")
#'                            
#' # Plot histogram of non-NA distances
#' p1 <- ggplot(data=subset(dist_hs1f, !is.na(DIST_NEAREST))) + theme_bw() + 
#'     ggtitle("Distance to nearest: hs1f") + xlab("distance") +
#'     geom_histogram(aes(x=DIST_NEAREST), binwidth=0.025, 
#'                    fill="steelblue", color="white")
#' plot(p1)
#'
#' @export
distToNearest <- function(db, sequenceColumn="JUNCTION", vCallColumn="V_CALL", 
                          jCallColumn="J_CALL", model=c("hs1f", "m1n", "ham", "aa", "hs5f"), 
                          normalize=c("length", "none"), symmetry=c("avg","min"),
                          first=TRUE, nproc=1, fields=NULL, cross=NULL, mst=FALSE) {
    # Hack for visibility of data.table and foreach index variables
    idx <- yidx <- .I <- NULL
    
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    
    # Check for valid columns
    columns <- c(sequenceColumn, vCallColumn, jCallColumn,fields,cross)
    columns <- columns[!is.null(columns)]
    
    check <- checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    # Convert case check for invalid characters
    db[, sequenceColumn] <- toupper(db[[sequenceColumn]])
    #check <- grepl("[^ACGTN]", db[[sequenceColumn]], perl=TRUE)
    #if (any(check)) {
    #  stop("Invalid sequence characters in the ", sequenceColumn, " column.")
    #}
    
    # Get targeting model
    if (model == "hs5f") {
        targeting_model <- HS5FModel
    }
    
    # Parse V and J columns to get gene
    # cat("V+J Column parsing\n")
    if (first) {
        db$V <- getGene(db[[vCallColumn]])
        db$J <- getGene(db[[jCallColumn]])
    } else {
        db$V1 <- getGene(db[[vCallColumn]], first=FALSE)
        db$J1 <- getGene(db[[jCallColumn]], first=FALSE)
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
            for(g in strsplit(ambig, split=',')[[1]]) {
                db$J[grepl(g, db$J1)] = ambig
            }
        }
    }
    
    # Create new column for distance to nearest neighbor
    db$TMP_DIST_NEAREST <- rep(NA, nrow(db))
    db$ROW_ID <- 1:nrow(db)
    db$L <- stri_length(db[[sequenceColumn]])
    
    # Create cluster of nproc size and export namespaces
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if( nproc==1 ) {
        # If needed to run on a single core/cpu then, register DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    # Calculate distance to nearest neighbor
    # cat("Calculating distance to nearest neighbor\n")
    
    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of V J L, instead of doing dplyr
    dt <- data.table(db)
    # Get the group indexes
    group_cols <- c("V","J","L")
    if (!is.null(fields)) {
        group_cols <- append(group_cols,fields)
    }
    dt <- dt[, list( yidx = list(.I) ) , by = group_cols ]
    groups <- dt[,yidx]
    lenGroups <- length(groups)
    #     
    #     
    #     # Get the group indices
    #     db$group_indices <- db %>% 
    #         dplyr::group_indices_(.dots=group_cols)
    #     
    
    # Export groups to the clusters
    export_functions <- list("db",
                             "groups", 
                             "cross",
                             "sequenceColumn", 
                             "model",
                             "normalize",
                             "symmetry",
                             "getClosestBy1Mers", 
                             "HS1FDistance",
                             "calcTargetingDistance",
                             "findUniqueJunctions",
                             "getPairwiseDistances")
    if (nproc > 1) { parallel::clusterExport(cluster, export_functions, 
                                             envir=environment()) }
    
    if (model %in% c("hs5f")) {
        # Export targeting model to processes
        if (nproc > 1) { parallel::clusterExport(cluster, list("targeting_model",
                                                               "getClosestBy5Mers"), 
                                                 envir=environment()) }    
        list_db <-
            foreach(idx=iterators::icount(lenGroups), .errorhandling='pass') %dopar% {
                db_group <- db[groups[[idx]],]
                crossGroups <- NULL
                if (!is.null(cross)) {
                    crossGroups <- db_group %>% dplyr::group_indices_(.dots=cross)
                }
                arrSeqs <- as.vector(unlist(db[groups[[idx]], sequenceColumn]))
                db_group$TMP_DIST_NEAREST <-
                    getClosestBy5Mers(arrSeqs,
                                      targeting_model=targeting_model,
                                      normalize=normalize,
                                      symmetry=symmetry,
                                      crossGroups=crossGroups,
                                      mst=mst)
                return(db_group)
            }    
    } else if (model %in% c("ham", "aa", "m1n", "hs1f")) {    
        list_db <-
            foreach(idx=iterators::icount(lenGroups), .errorhandling='pass') %dopar% {
                db_group <- db[groups[[idx]], ]
                crossGroups <- NULL
                if (!is.null(cross)) {
                    crossGroups <- db_group %>% dplyr::group_indices_(.dots=cross)
                }
                arrSeqs <-  as.vector(unlist(db[groups[[idx]], sequenceColumn]))
                db_group$TMP_DIST_NEAREST <-
                    getClosestBy1Mers(arrSeqs,
                                      model=model,
                                      normalize=normalize,
                                      crossGroups=crossGroups,
                                      mst=mst)
                return(db_group)
            }        
    }
    
    ## DEBUG
    # print(list_db)
    # for (i in 1:length(list_db)) {
    #     cat(i, ": ", nrow(list_db[[i]]), "\n", sep="")
    # }

    # Convert list from foreach into a db data.frame
    db <- dplyr::bind_rows(list_db)
    db <- db[order(db$ROW_ID), ]
    
    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    if (!is.null(cross)) {
        db$CROSS_DIST_NEAREST <- db$TMP_DIST_NEAREST
    } else {
        db$DIST_NEAREST <- db$TMP_DIST_NEAREST
    }
    return(db[, !(names(db) %in% c("V", "J", "L", "ROW_ID", "V1", "J1","TMP_DIST_NEAREST"))])
}
