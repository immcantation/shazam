# Generates distance to nearest neighbor

#' @include Shazam.R
NULL

#### Classes ####

#' Output of the \code{gmm} method of findThreshold
#' 
#' \code{GmmThreshold} contains output from the \code{gmm} method \link{findThreshold}. 
#' It includes parameters of two Gaussian fits and threshold cut.
#'
#' @slot   x            input distance vector with NA or infinite values removed.
#' @slot   model        first-second fit functions. 
#' @slot   cutoff       type of threshold cut.
#' @slot   a1           mixing weight of the first curve.
#' @slot   b1           second parameter of the first curve. Either the mean of a Normal 
#'                      distribution or shape of a Gamma distribution.
#' @slot   c1           third parameter of the first curve. Either the standard deviation of a 
#'                      Normal distribution or scale of a Gamma distribution.
#' @slot   a2           mixing weight of the second curve.
#' @slot   b2           second parameter of the second curve. Either the mean of a Normal 
#'                      distribution or shape of a Gamma distribution.
#' @slot   c2           third parameter of the second curve. Either the standard deviation 
#'                      of a Normal distribution or scale of a Gamma distribution.
#' @slot   loglk        log-likelihood of the fit.
#' @slot   threshold    threshold.
#' @slot   sensitivity  sensitivity.
#' @slot   specificity  specificity.
#' @slot   pvalue       p-value from Hartigans' dip statistic (HDS) test. 
#'                      Values less than 0.05 indicate significant bimodality.
#'
#' @seealso      \link{findThreshold}
#'
#' @name         GmmThreshold-class
#' @rdname       GmmThreshold-class
#' @aliases      GmmThreshold
#' @exportClass  GmmThreshold
setClass("GmmThreshold",
         slots=c(x="numeric",
                 model = "character",
                 cutoff = "character",
                 a1="numeric", 
                 b1="numeric", 
                 c1="numeric", 
                 a2="numeric", 
                 b2="numeric", 
                 c2="numeric",
                 loglk="numeric",
                 threshold="numeric",
                 sensitivity="numeric", 
                 specificity="numeric",
                 pvalue="numeric"))

#' Output of the \code{dens} method of findThreshold
#' 
#' \code{DensityThreshold} contains output from the \code{dens} method \link{findThreshold}. 
#'
#' @slot   x          input distance vector with NA or infinite values removed.
#' @slot   bandwidth  bandwidth value fit during density estimation.
#' @slot   xdens      x-axis (distance value) vector for smoothed density estimate.
#' @slot   ydens      y-axis (density) vector for smoothed density estimate.
#' @slot   threshold  distance threshold that separates two modes of the input distribution.
#'
#' @seealso      \link{findThreshold}
#'
#' @name         DensityThreshold-class
#' @rdname       DensityThreshold-class
#' @aliases      DensityThreshold
#' @exportClass  DensityThreshold
setClass("DensityThreshold",
         slots=c(x="numeric",
                 bandwidth="numeric",
                 xdens="numeric",
                 ydens="numeric",
                 threshold="numeric"))


#### Methods ####

#' @param    x    GmmThreshold object
#' 
#' @rdname   GmmThreshold-class
#' @aliases  GmmThreshold-method
#' @export
setMethod("print", c(x="GmmThreshold"), function(x) { print(x@threshold) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotGmmThreshold}.
#' 
#' @rdname   GmmThreshold-class
#' @aliases  GmmThreshold-method
#' @export
setMethod("plot", c(x="GmmThreshold", y="missing"),
          function(x, y, ...) { plotGmmThreshold(x, ...) })

#' @param    x    DensityThreshold object
#' 
#' @rdname   DensityThreshold-class
#' @aliases  DensityThreshold-method
#' @export
setMethod("print", c(x="DensityThreshold"), function(x) { print(x@threshold) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotDensityThreshold}.
#' 
#' @rdname   DensityThreshold-class
#' @aliases  DensityThreshold-method
#' @export
setMethod("plot", c(x="DensityThreshold", y="missing"),
          function(x, y, ...) { plotDensityThreshold(x, ...) })


#### Distance to Nearest ####

# Returns a 5-mer sliding window of given sequence
#
# @param   sequence   sequence string
# @return  An array of 5-mer sliding windows
#
# @examples
# window5Mers("ACGTNACGTNACGTN")
window5Mers <- function(sequence) {
    n <- stri_length(sequence)
    w <- substr(rep(sequence, n - 4), 1:(n - 4), 5:n)
    
    return(w)
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1               first nucleotide sequence, broken into 5mers.
# @param    seq2               second nucleotide sequence, broken into 5mers.
# @param    targetingDistance  targeting distance obtained from a targeting model
#                              with the function \code{calcTargetingDistance}.
# @param    symmetry           if model is hs5f, distance between seq1 and seq2 is either 
#                              the average (avg) of seq1->seq2 and seq2->seq1 or the 
#                              minimum (min).
# @return   distance between two sequences.
#
# @examples
# seq1 <- c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", 
#          "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
# seq2 <- c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", 
#          "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
# targeting_distance <- calcTargetingDistance(HH_S5F)
# shazam:::dist5Mers(seq1, seq2, targeting_distance)
dist5Mers <- function(seq1, seq2, targetingDistance, 
                      symmetry=c("avg", "min", "raw")) {
    # Evaluate choices
    symmetry <- match.arg(symmetry)
    
    # Get distance from targeting model
    #targeting_dist <- calcTargetingDistance(targetingModel)
    
    
    # Check all characters in seq1 and seq2 are valid,
    # found in the targetingModel distance matrix
    validChars <- rownames(targetingDistance)
    allChars <- unique(strsplit(paste(c(seq1, seq2), collapse=""), "")[[1]])
    invalidChars <- allChars[allChars %in% validChars == F]
    if (length(invalidChars) > 0 ) {
        stop(paste0("Character not found in targeting_dist: ", paste(invalidChars, collapse=", ")))
    }
    
    # Compute distance only on fivemers that have mutations
    fivemersWithMu <- substr(seq1, 3, 3) != substr(seq2, 3, 3)
    #fivemersWithNonNuc <- (!is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & 
    #                       !is.na(match(substr(seq2,3,3),c("A","C","G","T"))))
    #fivemersWithMu <- fivemersWithMu & fivemersWithNonNuc
    seq1 <- seq1[fivemersWithMu]
    seq2 <- seq2[fivemersWithMu]
    
    # Number of mutations (for normalization, if specified)
    #numbOfMutation <- sum(fivemersWithMu)
    
    dist <- NA
    tryCatch({
        if (length(seq1)==1){
            seq1_to_seq2 <- targetingDistance[substr(seq2, 3, 3), seq1]
            seq2_to_seq1 <- targetingDistance[substr(seq1, 3, 3), seq2]
        } else {
            seq1_to_seq2 <- diag(targetingDistance[substr(seq2, 3, 3), seq1])
            seq2_to_seq1 <- diag(targetingDistance[substr(seq1, 3, 3), seq2])
        }
        if (symmetry == "avg") {
            dist <- sum(apply(cbind(seq1_to_seq2, seq2_to_seq1), 1, mean))
        } else if (symmetry == "min") {
            dist <- sum(apply(cbind(seq1_to_seq2, seq2_to_seq1), 1, min))
        } else if (symmetry == "raw")  {
            dist <- c(seq1_to_seq2, seq2_to_seq1)
        }
    }, error = function(e) {
        warning(e)
        return(NA)
    })
    
    return(dist)
}


# Given an array of nucleotide sequences, find the pairwise distances
# 
# @param   sequences          character vector of nucleotide sequences.
# @param   targetingDistance  targeting distance obtained from a targeting model
#                             with the function `calcTargetingDistance`
# @param   symmetry           if model is hs5f, distance between seq1 and seq2 is either the
#                             average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
#
# @return  A matrix of pairwise distances between junction sequences.
pairwise5MerDist <- function(sequences, 
                             targetingDistance, 
                             symmetry=c("avg", "min")) {
    # get names
    seq_names <- names(sequences)

    # Initial checks
    symmetry <- match.arg(symmetry)
    
    # Convert junctions to uppercase
    sequences <- toupper(sequences)
    # Convert gaps to Ns
    sequences <- gsub('[-.]', 'N', sequences, fixed=T)
    # Add 'NN' to front and end of each sequence for fivemers
    sequences <- as.vector(sapply(sequences, function(x){ paste("NN", x, "NN", sep="") }))
    
    n_seq <- length(sequences)
    
    #Junctions are broken in to 5-mers based on a sliding window (of one) and placed in matrix
    #Each column is a junction
    #E.g. junctions 1234567, ABCDEFG, JKLMNOP becomes:
    # 12345   ABCDE   JKLMN
    # 23456   BCDEF   KLMNO
    # 34567   CDEFG   LMNOP
    .matSeqSlidingFiveMer <- sapply(sequences, function(x) { window5Mers(x) }, 
                                    simplify="matrix")
    
    # Compute pairwise distance between all sequences' fivemers (by column)
    .dist <- function(i) {
        c(rep.int(0, i - 1), 
          sapply(i:n_seq, function(j) { dist5Mers(.matSeqSlidingFiveMer[,i],
                                                  .matSeqSlidingFiveMer[,j],
                                                  targetingDistance,
                                                  symmetry=symmetry) }))
    }
    dist_mat <- sapply(1:n_seq, .dist)

    # Make distance matrix symmetric
    dist_mat <- dist_mat + t(dist_mat)
    
    # assign names
    if (!is.null(seq_names)) {
        rownames(dist_mat) <- seq_names
        colnames(dist_mat) <- seq_names
    }
    
    return(dist_mat)
}

# Given an array of nucleotide sequences and a vector indices (a subset of array of nucleotide sequences), 
# find the pairwise distances
# 
# @param   sequences          character vector of nucleotide sequences.
# @paramindx                  numeric vector of subsamples indices
# @param   targetingDistance  targeting distance obtained from a targeting model
#                             with the function `calcTargetingDistance`
# @param   symmetry           if model is hs5f, distance between seq1 and seq2 is either the
#                             average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
#
# @return  A non-square matrix of pairwise distances between junction sequences.
nonsquare5MerDist <- function(sequences, indx,  targetingDistance, symmetry=c("avg", "min")) {
    
    # get names
    seq_names <- names(sequences)
    
    # Initial checks
    symmetry <- match.arg(symmetry)
    
    # Convert junctions to uppercase
    sequences <- toupper(sequences)
    # Convert gaps to Ns
    sequences <- gsub('[-.]', 'N', sequences, fixed=T)
    # Add 'NN' to front and end of each sequence for fivemers
    sequences <- as.vector(sapply(sequences, function(x){ paste("NN", x, "NN", sep="") }))
    
    n_seq <- length(sequences)
    
    #Junctions are broken in to 5-mers based on a sliding window (of one) and placed in matrix
    #Each column is a junction
    #E.g. junctions 1234567, ABCDEFG, JKLMNOP becomes:
    # 12345   ABCDE   JKLMN
    # 23456   BCDEF   KLMNO
    # 34567   CDEFG   LMNOP
    .matSeqSlidingFiveMer <- sapply(sequences, function(x) { window5Mers(x) }, 
                                    simplify="matrix")
    
    # # Compute pairwise distance between all sequences' fivemers (by column)
    # .dist <- function(i) {
    #     d <- c(rep.int(0, i - 1), 
    #       sapply(i:n_seq, function(j) { dist5Mers(.matSeqSlidingFiveMer[,i],
    #                                               .matSeqSlidingFiveMer[,j],
    #                                               targetingDistance,
    #                                               symmetry=symmetry) }))
    # }
    dist_mat <- matrix(NA, nrow=n_seq, ncol=n_seq)
    diag(dist_mat) <- 0
    indx <- sort(indx)
    for (i in 1:n_seq) {
        if (!(i %in% indx)) next
        for (j in 1:n_seq) {
            if (!is.na(dist_mat[i,j])) next
            dist_mat[i,j] = dist5Mers(.matSeqSlidingFiveMer[,i], .matSeqSlidingFiveMer[,j], 
                                      targetingDistance, symmetry=symmetry)
            dist_mat[j,i] = dist_mat[i,j]
        }
    }
    sub_dist_mat <- dist_mat[indx,]
    
    # assign names
    if (!is.null(seq_names)) {
        rownames(sub_dist_mat) <- seq_names[indx]
        colnames(sub_dist_mat) <- seq_names
    }
    return(sub_dist_mat)
}


# Subset to unique sequences
# 
# @param    sequences  character vector of sequences
#
# @return   Named vector of unique sequences, with names as the sequence itself.
findUniqSeq <- function(sequences) {
    seq_uniq <- unique(sequences)
    names(seq_uniq) <- seq_uniq
    
    return(seq_uniq)
}

# Get chars in the distance model
# 
# @param    model  
#
# @return    vector of unique chars in the distance model
# @examples 
# getCharsInModel("hh_s1f")
getCharsInModel <- function(model) {
    if (model == "ham") {
        chars <- colnames(getDNAMatrix(gap=0))
    } else if (model == "aa") {
        chars <- colnames(getAAMatrix())
    } else if (model == "hh_s1f") {
        chars <- colnames(HH_S1F_Distance)
    } else if (model == "hh_s5f") {
        chars <-rownames(HH_S5F@targeting)
    } else if (model == "mk_rs1nf") {
        chars <- colnames(MK_RS1NF_Distance)
    } else if (model == "mk_rs5nf") {
        chars <-rownames(MK_RS1NF@targeting)
    } else if (model == "hs1f_compat") {
        chars <- colnames(HS1F_Compat)
    } else if (model == "m1n_compat") {
        chars <- colnames(M1N_Compat)
    }
  
    return(chars)
}

# Validate the sequence
# 
# @param    seq
# @param    validChars
#
# @return    TRUE is all the character in the sequence are found in validChars; 
#            FALSE otherwise
# @examples 
# allValidChars("ATCG", getCharsInModel("hh_s1f"))
# allValidChars("ATCG.", getCharsInModel("hh_s1f"))
# allValidChars("ATCGJ", getCharsInModel("hh_s1f"))
allValidChars <- function(seq, validChars) {
    all(unique(strsplit(seq, "")[[1]]) %in% validChars)
}
    
# Given an array of sequences, find the distance to the closest sequence
#
# @param    sequences      character vector of sequences.
# @param    model          5-mer or 1-mer distance model
# @param    normalize      method of normalization. Default is "none".
#                          "len" = normalize distance by length of junction.
#                          "mut" = normalize distance by number of mutations in 
#                          junction.
# @param    symmetr        if model is hs5f or mrs5nf, distance between seq1 and seq2 is either the
#                          average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
# @param    crossGroups    column for grouping to calculate distances across groups 
#                          (self vs others).
# @param    mst            if true, return comma-separated branch lengths from minimum 
#                          spanning tree.
#
# @return   A vector of distances to the closest sequence.
#
# @examples
# sequences <- c("ACGTACGTACGT", "ACGAACGTACGT", "ACGAACGTATGT", "ACGAACGTATGC",
#                "ACGAACGTATCC", "AAAAAAAAAAAA", "A-GAACGTATCC", "AAAAAA---AAA")
# shazam:::nearestDist(sequences, model="ham", normalize="none")
# shazam:::nearestDist(sequences, model="aa", normalize="none")
# shazam:::nearestDist(sequences, model="ham", normalize="len")
# shazam:::nearestDist(sequences, model="aa", normalize="len")
nearestDist <- function(sequences, model=c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", 
                                          "hs1f_compat", "m1n_compat"),
                       normalize=c("none", "len", "mut"),
                       symmetry=c("avg", "min"),
                       crossGroups=NULL, mst=FALSE,
                       subsample=NULL) {
    ## DEBUG
    # sequences <- c("ACGTACGTACGT", "ACGAACGTACGT", "AAAAAAAAAAAA", "A-AAAA---AAA")
    # model="aa"; normalize="len"; crossGroups=NULL; mst=FALSE
    
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    
    ## If crossGroup requested, but only one group found, return NA
    if (!is.null(crossGroups) & length(unique(crossGroups)) < 2) {
        seq_dist <- rep(NA, length(sequences))
        return (seq_dist)
    }
    
    # Find unique sequences
    seq_uniq <- findUniqSeq(sequences)
    n_uniq <- length(seq_uniq)
    
    # corresponding crossGroups values for seq_uniq
    if (!is.null(crossGroups)) {
        stopifnot( all.equal(sequences[match(seq_uniq, sequences)], 
                             seq_uniq, check.attributes=FALSE) )
        crossGroups_uniq <- crossGroups[match(seq_uniq, sequences)]
    }
    
    # Initialize return vector and computation vector
    seq_dist <- setNames(rep(NA, length(sequences)), sequences)
    seq_uniq_dist <- rep(NA, n_uniq)
    
    # Compute distances between sequences
    if (n_uniq > 1) {
        # Check for length mismatches
        seq_length <-  unique(stri_length(seq_uniq))
        if (length(seq_length) > 1) {
            stop("Unexpected. Different sequence lengths found.")
        }
        # check subSampling
        subSampling <- all(!is.null(subsample), subsample < n_uniq)
        if (subSampling) indx <- sample(x=1:n_uniq, size=subsample, replace=FALSE, prob=NULL)
        # corresponding subsampling of crossGroups_uniq
        if (subSampling & !is.null(crossGroups)) {
            crossGroups_uniq_sub <- crossGroups_uniq[indx]
        }
        # Get distance matrix
        if (model == "ham") {
            if (subSampling) {
                dist_mat <- nonsquareDist(seq_uniq, indx, dist_mat=getDNAMatrix(gap=0))
            } else {
                dist_mat <- pairwiseDist(seq_uniq, dist_mat=getDNAMatrix(gap=0))
            }
        } else if (model == "aa") {
            seq_uniq <- setNames(alakazam::translateDNA(seq_uniq), seq_uniq)
            if (subSampling) {
                dist_mat <- nonsquareDist(seq_uniq, indx, dist_mat=getAAMatrix())
            } else {
                dist_mat <- pairwiseDist(seq_uniq, dist_mat=getAAMatrix())
            }
        } else if (model == "hh_s1f") {
            if (subSampling) {
                dist_mat <- nonsquareDist(seq_uniq, indx, dist_mat=HH_S1F_Distance)
            } else {
                dist_mat <- pairwiseDist(seq_uniq, dist_mat=HH_S1F_Distance)
            }
        } else if (model == "mk_rs1nf") {
            if (subSampling) {
                dist_mat <- nonsquareDist(seq_uniq, indx, dist_mat=MK_RS1NF_Distance)
            } else {
                dist_mat <- pairwiseDist(seq_uniq, dist_mat=MK_RS1NF_Distance)
            }
        } else if (model == "hh_s5f") {
            if (subSampling) {
                dist_mat <- nonsquare5MerDist(seq_uniq, indx, HH_S5F_Distance, symmetry=symmetry)
            } else {
                dist_mat <- pairwise5MerDist(seq_uniq, HH_S5F_Distance, symmetry=symmetry)
            }
        } else if (model == "mk_rs5nf") {
            if (subSampling) {
                dist_mat <- nonsquare5MerDist(seq_uniq, indx, MK_RS5NF_Distance, symmetry=symmetry)
            } else {
                dist_mat <- pairwise5MerDist(seq_uniq, MK_RS5NF_Distance, symmetry=symmetry)
            }
        } else if (model == "hs1f_compat") {
            if (subSampling) {
                dist_mat <- nonsquareDist(seq_uniq, indx, dist_mat=HS1F_Compat)
            } else {
                dist_mat <- pairwiseDist(seq_uniq, dist_mat=HS1F_Compat)
            }
        } else if (model == "m1n_compat") {
            if (subSampling) {
                dist_mat <- nonsquareDist(seq_uniq, indx, dist_mat=M1N_Compat)
            } else {
                dist_mat <- pairwiseDist(seq_uniq, dist_mat=M1N_Compat)
            }
        }                
        ## DEBUG
        # cat("\n-> seq_uniq:\n")
        # print(seq_uniq)
        # cat("\n-> dist_mat (raw):\n")
        # print(dist_mat)

        # Normalize distances
        if (normalize == "len") { 
            dist_mat <- dist_mat / seq_length
        } else if (normalize == "mut") {
            #dist <- dist/sum(strsplit(seq1,"")[[1]] != strsplit(seq2,"")[[1]])
            stop('Sorry! nomalize="mut" is not available.')
        }
        
        ## DEBUG
        # cat("\n-> seq_length:\n")
        # print(seq_length)
        # cat("\n-> dist_mat (normalized):\n")
        # print(dist_mat)
        
    } else {
        return(seq_dist)
    }
    
    # Find minimum distance for each sequence
    if (is.null(crossGroups)) {
        if(!mst) {
            # Return smaller value greater than 0
            # If all 0, return NA
            .dmin <- function(i) { 
                x <- dist_mat[, i]
                gt0 <- which(x > 0)
                if (length(gt0) != 0) { min(x[gt0]) } else { NA }
            }
            
            ## TODO: Could be an apply over columns
            seq_uniq_dist <- setNames(sapply(1:n_uniq, .dmin), names(seq_uniq))
        } else {
            # Get adjacency matrix of minimum spanning tree
            adj <- ape::mst(dist_mat)
            
            # TODO: This could be cleaner
            # Get value(s) from mst branches
            # If none (broken mst!), return NA
            # If multiple values, comma-join
            .dmst <- function(i) { 
                gt0 <- which(adj[, i] == 1)
                if (length(gt0) != 0) { 
                    stri_join(round(dist_mat[, i][gt0], 4), collapse=",") 
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
            #cat(i,"\n")
            this_group <- crossGroups[i]
            other_groups <-  which(crossGroups != this_group)
            other_seq <- unique(sequences[other_groups])
            other_idx <- match(other_seq, seq_uniq)
            this_idx <- match(sequences[i], seq_uniq)
            
            stopifnot( all.equal( other_seq, seq_uniq[other_idx] , check.attributes=FALSE ) )
            stopifnot( all.equal( sequences[i], seq_uniq[this_idx] , check.attributes=FALSE ) )
            
            # the next two checks may not always be true
            # this happens when all the out-group sequences are identical to the in-group sequences
            #stopifnot( all( crossGroups_uniq[other_idx] != this_group ) )
            #stopifnot( crossGroups_uniq[this_idx] == this_group )
            
            if (subSampling) {
                # When there is subsampling, nonsquareDist returns a non-n-by-n matrix 
                # This matrix has fewers than n rows, and exactly n cols
                # For each unique sequence, look for its cross-group distances in its column, 
                #     NOT in its row (because there will be fewer than n rows)
                
                # dist_mat rows correspond to seq_uniq[indx]
                # (indx itself is wrt seq_uniq)
                # (other_idx is also wrt seq_uni)
                
                # which other_seq are included in the subsampled seqs represented by
                #       the available rows in dist_mat?
                # wrt dist_mat
                other_avail_wrt_dist_mat <- which(indx %in% other_idx)
                
                if (length(other_avail_wrt_dist_mat)>0) {
                    # the next two checks may not always be true
                    # this happens when all the out-group sequences are identical to the in-group sequences
                    #stopifnot(all( crossGroups_uniq_sub[other_avail_wrt_dist_mat] != this_group ))
                    #stopifnot(all( crossGroups_uniq_sub[-other_avail_wrt_dist_mat] == this_group ))
                    
                    r <- dist_mat[other_avail_wrt_dist_mat, this_idx]
                } else {
                    stopifnot(all( crossGroups_uniq_sub == this_group ))
                    return(NA)
                }
            } else {
                # without subsampling
                # dist_mat is a n-by-n matrix
                stopifnot( all(other_idx <= nrow(dist_mat) ) ) 
                r <- dist_mat[other_idx, this_idx]
            }
            
            gt0 <- which(r > 0)
            
            if (length(gt0) != 0) { return(min(r[gt0])) } else { return(NA) }
        }
        
        # Define return distance vector
        seq_dist <- setNames(sapply(1:length(sequences), .dcross), sequences)
    }
    
    return(round(seq_dist, 4))
}


#' Distance to nearest neighbor
#'
#' Get non-zero distance of every heavy chain (\code{IGH}) sequence (as defined by 
#' \code{sequenceColumn}) to its nearest sequence in a partition of heavy chains sharing the same 
#' V gene, J gene, and junction length (VJL), or in a partition of single cells with heavy chains
#' sharing the same heavy chain VJL combination, or of single cells with heavy and light chains 
#' sharing the same heavy chain VJL and light chain VJL combinations.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing the junction for grouping and for calculating
#'                           nearest neighbot distances. Note that while both heavy and light chain junctions
#'                           may be used for VJL grouping, only the heavy chain junction is used to calculate 
#'                           distances.
#' @param    vCallColumn     name of the column containing the V-segment allele calls.
#' @param    jCallColumn     name of the column containing the J-segment allele calls.
#' @param    model           underlying SHM model, which must be one of 
#'                           \code{c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "hs1f_compat", "m1n_compat")}.
#'                           See Details for further information.
#' @param    normalize       method of normalization. The default is \code{"len"}, which 
#'                           divides the distance by the length of the sequence group. If 
#'                           \code{"none"} then no normalization if performed.
#' @param    symmetry        if model is hs5f, distance between seq1 and seq2 is either the
#'                           average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
#' @param    first           if \code{TRUE} only the first call of the gene assignments 
#'                           is used. if \code{FALSE} the union of ambiguous gene 
#'                           assignments is used to group all sequences with any 
#'                           overlapping gene calls.
#' @param    VJthenLen       a Boolean value specifying whether to perform partitioning as a 2-stage
#'                           process. If \code{TRUE}, partitions are made first based on V and J
#'                           annotations, and then further split based on junction lengths corresponding 
#'                           to \code{sequenceColumn}. If \code{FALSE}, perform partition as a 1-stage 
#'                           process during which V annotation, J annotation, and junction length are used 
#'                           to create partitions simultaneously. Defaults to \code{TRUE}.
#' @param    nproc           number of cores to distribute the function over.
#' @param    fields          additional fields to use for grouping.
#' @param    cross           character vector of column names to use for grouping to calculate 
#'                           distances across groups. Meaning the columns that define self versus others.
#' @param    mst             if \code{TRUE}, return comma-separated branch lengths from minimum 
#'                           spanning tree.
#' @param    subsample       number of sequences to subsample for speeding up pairwise-distance-matrix calculation. 
#'                           Subsampling is performed without replacement in each VJL group of heavy chain sequences. 
#'                           If \code{subsample} is larger than the unique number of heavy chain sequences in each 
#'                           VJL group, then the subsampling process is ignored for that group. For each heavy chain
#'                           sequence in \code{db}, the reported \code{dist_nearest} is the distance to the closest
#'                           heavy chain sequence in the subsampled set for the VJL group. If \code{NULL} no 
#'                           subsampling is performed.
#' @param    progress        if \code{TRUE} print a progress bar.
#' @param    cellIdColumn    name of the character column containing cell IDs. Only applicable and required for 
#'                           single-cell mode.
#' @param    locusColumn     name of the character column containing locus information. Only applicable and 
#'                           required for single-cell mode.
#' @param    groupUsingOnlyIGH    use only heavy chain (\code{IGH}) sequences for VJL grouping, disregarding 
#'                                light chains. Only applicable and required for single-cell mode. 
#'                                Default is \code{TRUE}. Also see \link[alakazam]{groupGenes}.                   
#' @param    keepVJLgroup         a Boolean value specifying whether to keep in the output the the column 
#'                                column indicating grouping based on VJL combinations. Only applicable for
#'                                1-stage partitioning (i.e. \code{VJthenLen=FALSE}). Also see 
#'                                \link[alakazam]{groupGenes}.
#' 
#' @return   Returns a modified \code{db} data.frame with nearest neighbor distances between heavy chain
#'           sequences in the \code{dist_nearest} column if \code{cross=NULL}. If \code{cross} was 
#'           specified, distances will be added as the \code{cross_dist_nearest} column. 
#'           
#'           Note that distances between light chain sequences are not calculated, even if light chains 
#'           were used for VJL grouping via \code{groupUsingOnlyIGH=FALSE}. Light chain sequences, if any,
#'           will have \code{NA} in the \code{dist_nearest} field.
#'           
#'           Note that the output \code{vCallColumn} and \code{jCallColumn} columns will be converted to 
#'           \code{character} if they were \code{factor} in the input \code{db}.
#'
#' @details
#' 
#' To invoke single-cell mode, both \code{cellIdColumn} and \code{locusColumn} must be supplied. 
#' Otherwise, the function will run under non-single-cell mode.
#' 
#' Under single-cell mode, only heavy chain sequences will be used for calculating nearest neighbor
#' distances. Under non-single-cell mode, all input sequences will be used for calculating nearest
#' neighbor distances, regardless of the values in the \code{locusColumn} field (if present).
#' 
#' For single-cell mode, the input format is the same as that for \link[alakazam]{groupGenes}. 
#' Namely, each row represents a sequence/chain. Sequences/chains from the same cell are linked
#' by a cell ID in the \code{cellIdColumn} field. Under this mode, there is a choice of whether 
#' grouping should be done using only heavy chain (\code{IGH}) sequences only, or using both 
#' heavy chain (\code{IGH}) and light chain (\code{IGK}, \code{IGL}) sequences. This is governed 
#' by \code{groupUsingOnlyIGH}.
#' 
#' If used, values in the \code{locusColumn} column must be one of \code{"IGH"}, \code{"IGK"}, and \code{"IGL"}.
#' 
#' Note that for \code{distToNearest}, a cell with multiple heavy chains is not allowed.
#' 
#' The distance to nearest (heavy chain) neighbor can be used to estimate a threshold for assigning 
#' Ig sequences to clonal groups. A histogram of the resulting vector is often bimodal, with the 
#' ideal threshold being a value that separates the two modes.
#' 
#' The following distance measures are accepted by the \code{model} parameter.
#' 
#' \itemize{
#'   \item \code{"ham"}:          Single nucleotide Hamming distance matrix from \link[alakazam]{getDNAMatrix} 
#'                                with gaps assigned zero distance.
#'   \item \code{"aa"}:           Single amino acid Hamming distance matrix from \link[alakazam]{getAAMatrix}.
#'   \item \code{"hh_s1f"}:       Human single nucleotide distance matrix derived from \link{HH_S1F} with 
#'                                \link{calcTargetingDistance}.
#'   \item \code{"hh_s5f"}:       Human 5-mer nucleotide context distance matix derived from \link{HH_S5F} with 
#'                                \link{calcTargetingDistance}.
#'   \item \code{"mk_rs1nf"}:     Mouse single nucleotide distance matrix derived from \link{MK_RS1NF} with 
#'                                \link{calcTargetingDistance}.
#'   \item \code{"mk_rs5nf"}:     Mouse 5-mer nucleotide context distance matrix derived from \link{MK_RS1NF} with 
#'                                \link{calcTargetingDistance}.
#'   \item \code{"hs1f_compat"}:  Backwards compatible human single nucleotide distance matrix used in 
#'                                SHazaM v0.1.4 and Change-O v0.3.3.
#'   \item \code{"m1n_compat"}:   Backwards compatibley mouse single nucleotide distance matrix used in 
#'                                SHazaM v0.1.4 and Change-O v0.3.3.
#' }
#' 
#' Note on \code{NA}s: if, for a given combination of V gene, J gene, and sequence length,
#' there is only 1  heavy chain sequence (as defined by \code{sequenceColumn}), \code{NA} is 
#' returned instead of a distance (since it has no heavy chain neighbor). If for a given combination 
#' there are multiple heavy chain sequences but only 1 unique one, (in which case every heavy cahin 
#' sequence in this group is the de facto nearest neighbor to each other, thus giving rise to distances 
#' of 0), \code{NA}s are returned instead of zero-distances.
#' 
#' Note on \code{subsample}: Subsampling is performed independently in each VJL group for heavy chain
#' sequences. If \code{subsample} is larger than number of heavy chain sequences in the group, it is 
#' ignored. In other words, subsampling is performed only on groups in which the number of heavy chain 
#' sequences is equal to or greater than \code{subsample}. \code{dist_nearest} has values calculated 
#' using all heavy chain sequences in the group for groups with fewer than \code{subsample} heavy chain
#' sequences, and values calculated using a subset of heavy chain sequences for the larger groups. 
#' To select a value of \code{subsample}, it can be useful to explore the group sizes in \code{db} 
#' (and the number of heavy chain sequences in those groups).
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
#'           from a \link{TargetingModel} object. See \link{HH_S5F}, \link{HH_S1F}, 
#'           \link{MK_RS1NF}, \link[alakazam]{getDNAMatrix}, and \link[alakazam]{getAAMatrix}
#'           for individual model details.
#' 
#' @examples
#' # Subset example data to one sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, sample_id == "-1h")
#' 
#' # Use genotyped V assignments, Hamming distance, and normalize by junction length
#' # First partition based on V and J assignments, then by junction length
#' # Take into consideration ambiguous V and J annotations
#' dist <- distToNearest(db, sequenceColumn="junction", 
#'                       vCallColumn="v_call_genotyped", jCallColumn="j_call",
#'                       model="ham", first=FALSE, VJthenLen=TRUE, normalize="len")
#'                            
#' # Plot histogram of non-NA distances
#' p1 <- ggplot(data=subset(dist, !is.na(dist_nearest))) + 
#'       theme_bw() + 
#'       ggtitle("Distance to nearest: Hamming") + 
#'       xlab("distance") +
#'       geom_histogram(aes(x=dist_nearest), binwidth=0.025, 
#'                      fill="steelblue", color="white")
#' plot(p1)
#' 
#' @export
distToNearest <- function(db, sequenceColumn="junction", vCallColumn="v_call", jCallColumn="j_call", 
                          model=c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", "m1n_compat", "hs1f_compat"), 
                          normalize=c("len", "none"), symmetry=c("avg", "min"),
                          first=TRUE, VJthenLen=TRUE, nproc=1, fields=NULL, cross=NULL, 
                          mst=FALSE, subsample=NULL, progress=FALSE,
                          cellIdColumn=NULL, locusColumn=NULL, 
                          groupUsingOnlyIGH=TRUE, keepVJLgroup=TRUE) {
    
    # Hack for visibility of foreach index variables
    i <- NULL
    
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    
    
    # single-cell mode?
    if ( !is.null(cellIdColumn) & !is.null(locusColumn) ) {
        singleCell <- TRUE
        
        if (!all(db[[locusColumn]] %in% c("IGH", "IGK", "IGL"))) {
            stop("The locus column must be one of {IGH, IGK, IGL}.")
        }
        
    } else {
        singleCell <- FALSE
    } 
    
    columns <- c(sequenceColumn, vCallColumn, jCallColumn, fields, cross)
    columns <- columns[!is.null(columns)]
    
    check <- checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn)) 
    
    # Disallow multiple heavy chains per cell
    if (singleCell) {
        bool <- sapply(unique(db[[cellIdColumn]]), function(x) { 
            return( sum( db[[locusColumn]][db[[cellIdColumn]]==x] == "IGH" )>1 ) 
        } )
        if (any(bool)) {
            stop("Detected multiple heavy chains in cell(s). Each cell must contain only 1 heavy chain.")
        }
    }
    
    
    # Check for invalid characters
    # heavy
    valid_seq <- sapply(db[[sequenceColumn]], allValidChars, getCharsInModel(model)) 
    not_valid_seq <- which(!valid_seq)
    if (length(not_valid_seq) > 0) {
        warning("Invalid sequence characters in the ", sequenceColumn, 
                " column. ", length(not_valid_seq), " sequence(s) removed")
        db <- db[valid_seq, ]
    }
    
    # junction length columns (prep for groupGenes)
    junc_len <- "JUNC_LEN"
    db[[junc_len]] <- stri_length(db[[sequenceColumn]])
    
    # create V+J grouping, or V+J+L grouping
    if (VJthenLen) {
        # 2-stage partitioning using first V+J and then L
        # V+J only first
        # creates $vj_group
        db <- groupGenes(db, v_call=vCallColumn, j_call=jCallColumn, junc_len=NULL,
                         cell_id=cellIdColumn, locus=locusColumn, only_igh=groupUsingOnlyIGH,
                         first=first)
        # L (later)  
        group_cols <- c("vj_group", junc_len)
        
    } else {
        # 1-stage partitioning using V+J+L simultaneously
        # creates $vj_group
        # note that despite the name (VJ), this is based on V+J+L
        db <- groupGenes(db, v_call=vCallColumn, j_call=jCallColumn, junc_len=junc_len,
                         cell_id=cellIdColumn, locus=locusColumn, only_igh=groupUsingOnlyIGH,
                         first=first)
        group_cols <- c("vj_group")
    }
    
    # groups to use
    if (!is.null(fields)) {
        group_cols <- append(group_cols,fields)
    }
    # unique groups
    # not necessary but good practice to force as df and assign colnames
    # (in case group_cols has length 1; which can happen in groupBaseline)
    uniqueGroups <- data.frame(unique(db[, group_cols]), stringsAsFactors=FALSE)
    colnames(uniqueGroups) <- group_cols
    rownames(uniqueGroups) <- NULL
    # indices
    # crucial to have simplify=FALSE 
    # (otherwise won't return a list if uniqueClones has length 1)
    uniqueGroupsIdx <- sapply(1:nrow(uniqueGroups), function(i){
        curGroup <- data.frame(uniqueGroups[i, ], stringsAsFactors=FALSE)
        colnames(curGroup) <- group_cols
        # match for each field
        curIdx <- sapply(group_cols, function(coln){
            db[[coln]]==curGroup[, coln]
        }, simplify=FALSE)
        curIdx <- do.call(rbind, curIdx)
        # intersect to get match across fields 
        curIdx <- which(colSums(curIdx)==length(group_cols))
        # sanity check
        # no NA
        stopifnot( all(!is.na(curIdx)) )
        # index within range of db
        stopifnot( max(curIdx) <= nrow(db) )
        return(curIdx)
    }, simplify=FALSE)

    # Create new column for distance to nearest neighbor
    db$TMP_DIST_NEAREST <- rep(NA, nrow(db))
    db$ROW_ID <- 1:nrow(db)
    
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
    
    # Export groups to the clusters
    if (nproc > 1) { 
        export_functions <- list("db",
                                 "uniqueGroupsIdx", 
                                 "cross",
                                 "mst",
                                 "subsample",
                                 "sequenceColumn", 
                                 "model",
                                 "normalize",
                                 "symmetry",
                                 "nearestDist", 
                                 "HH_S1F_Distance",
                                 "MK_RS1NF_Distance",
                                 "HH_S5F_Distance",
                                 "MK_RS5NF_Distance",
                                 "HS1F_Compat",
                                 "M1N_Compat",
                                 "calcTargetingDistance",
                                 "findUniqSeq",
                                 "pairwise5MerDist",
                                 "nonsquare5MerDist",
                                 "singleCell",
                                 "locusColumn")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    }
    
    
    
    n_groups <- length(uniqueGroupsIdx)
    if (progress) { 
        pb <- progressBar(n_groups) 
    }
    list_db <- foreach(i=1:n_groups, .errorhandling='stop') %dopar% {
        # wrt db
        idx <- uniqueGroupsIdx[[i]]
        
        if (singleCell) {
            # only use IGH
            # wrt idx
            idxBool <- db[[locusColumn]][idx] == "IGH"
        } else {
            idxBool <- rep(TRUE, length(idx))
        }
        
        db_group <- db[idx, ]
        
        crossGroups <- NULL
        if (!is.null(cross)) {
            crossGroups <- db_group %>% dplyr::group_indices(!!!rlang::syms(cross))
        }
        
        arrSeqs <-  db[[sequenceColumn]][idx]
            
        db_group$TMP_DIST_NEAREST[idxBool] <- nearestDist(arrSeqs[idxBool], 
                                                          model=model,
                                                          normalize=normalize,
                                                          symmetry=symmetry,
                                                          crossGroups=crossGroups[idxBool],
                                                          mst=mst,
                                                          subsample=subsample)
        # Update progress
        if (progress) { pb$tick() }
        
        return(db_group)
    }
    
    # Convert list from foreach into a db data.frame
    db <- do.call(rbind, list_db)
    db <- db[order(db$ROW_ID), ]
    
    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    if (!is.null(cross)) {
        db$cross_dist_nearest <- db$TMP_DIST_NEAREST
    } else {
        db$dist_nearest <- db$TMP_DIST_NEAREST
    }
    
    # prepare db for return
    if ((!VJthenLen) && keepVJLgroup) {
        db$vjl_group <- db[["vj_group"]]
    }
    db <- db[, !(names(db) %in% c(junc_len, "vj_group", "ROW_ID", "V1", "J1","TMP_DIST_NEAREST"))]
    
    return(db)
}


#### Distance Threshold Detection ####

#' Find distance threshold
#'
#' \code{findThreshold} automtically determines an optimal threshold for clonal assignment of
#' Ig sequences using a vector of nearest neighbor distances. It provides two alternative methods 
#' using either a Gamma/Gaussian Mixture Model fit (\code{method="gmm"}) or kernel density 
#' fit (\code{method="density"}).
#'
#' @param    distances  numeric vector containing nearest neighbor distances. 
#' @param    method     string defining the method to use for determining the optimal threshold.
#'                      One of \code{"gmm"} or \code{"density"}. See Details for methodological
#'                      descriptions.
#' @param    edge       upper range as a fraction of the data density to rule initialization of 
#'                      Gaussian fit parameters. Default value is 90% of the entries (0.9).
#'                      Applies only when \code{method="density"}. .
#' @param    cross      supplementary nearest neighbor distance vector output from \link{distToNearest} 
#'                      for initialization of the Gaussian fit parameters. 
#'                      Applies only when \code{method="gmm"}. 
#' @param    subsample  maximum number of distances to subsample to before threshold detection.
#' @param    model      allows the user to choose among four possible combinations of fitting curves: 
#'                      \code{"norm-norm"}, \code{"norm-gamma"}, \code{"gamma-norm"}, 
#'                      and \code{"gamma-gamma"}. Applies only when \code{method="gmm"}.
#' @param    cutoff     method to use for threshold selection: the optimal threshold \code{"opt"}, 
#'                      the intersection point of the two fitted curves \code{"intersect"}, or 
#'                      a value defined by user for one of the sensitivity or specificity \code{"user"}.
#'                      Applies only when \code{method="gmm"}.
#' @param    sen        sensitivity required. Applies only when \code{method="gmm"} and \code{cutoff="user"}.
#' @param    spc        specificity required. Applies only when \code{method="gmm"} and \code{cutoff="user"}.
#'                      
#' @param    progress   if \code{TRUE} print a progress bar. 
#' @return   
#' \itemize{
#'   \item \code{"gmm"} method:      Returns a \link{GmmThreshold} object including the  
#'                                   \code{threshold} and the function fit parameters, i.e.
#'                                   mixing weight, mean, and standard deviation of a Normal distribution, or 
#'                                   mixing weight, shape and scale of a Gamma distribution.
#'   \item \code{"density"} method:  Returns a \link{DensityThreshold} object including the optimum 
#'                                   \code{threshold} and the density fit parameters.
#' }
#'
#' @details 
#' \itemize{ 
#'   \item \code{"gmm"}:     Performs a maximum-likelihood fitting procedure, for learning 
#'                           the parameters of two mixture univariate, either Gamma or Gaussian, distributions 
#'                           which fit the bimodal distribution entries. Retrieving the fit parameters, 
#'                           it then calculates the optimum threshold \code{method="optimal"}, where the 
#'                           average of the sensitivity plus specificity reaches its maximum. In addition, 
#'                           the \code{findThreshold} function is also able 
#'                           to calculate the intersection point (\code{method="intersect"}) of the two fitted curves 
#'                           and allows the user to invoke its value as the cut-off point, instead of optimal point.
#'   \item \code{"density"}: Fits a binned approximation to the ordinary kernel density estimate
#'                           to the nearest neighbor distances after determining the optimal
#'                           bandwidth for the density estimate via least-squares cross-validation of 
#'                           the 4th derivative of the kernel density estimator. The optimal threshold
#'                           is set as the minimum value in the valley in the density estimate
#'                           between the two modes of the distribution.
#' }
#' 
#' @seealso  See \link{distToNearest} for generating the nearest neighbor distance vectors.
#'           See \link{plotGmmThreshold} and \link{plotDensityThreshold} for plotting output.
#'      
#' @note 
#' Visually inspecting the resulting distribution fits is strongly recommended when using 
#' either fitting method. Empirical observations imply that the bimodality 
#' of the distance-to-nearest distribution is detectable for a minimum of 1,000 distances.
#' Larger numbers of distances will improve the fitting procedure, although this can come 
#' at the expense of higher computational demands.
#' 
#' @examples
#' \donttest{
#' # Subset example data to one sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, sample_id == "-1h")
#' 
#' # Use nucleotide Hamming distance and normalize by junction length
#' db <- distToNearest(db, sequenceColumn="junction", vCallColumn="v_call",
#'                     jCallColumn="j_call", model="ham", normalize="len", nproc=1)
#'                             
#' # Find threshold using the "gmm" method with optimal threshold
#' output <- findThreshold(db$dist_nearest, method="gmm", model="gamma-gamma", cutoff="opt")
#' plot(output, binwidth=0.02, title=paste0(output@model, "   loglk=", output@loglk))
#' print(output)
#'
#' # Find threshold using the "gmm" method with user defined specificity
#' output <- findThreshold(db$dist_nearest, method="gmm", model="gamma-gamma", 
#'                         cutoff="user", spc=0.99)
#' plot(output, binwidth=0.02, title=paste0(output@model, "   loglk=", output@loglk))
#' print(output)
#'
#' # Find threshold using the "density" method and plot the results
#' output <- findThreshold(db$dist_nearest, method="density")
#' plot(output)
#' print(output)
#' }
#' @export
findThreshold <- function (distances, method=c("density", "gmm"), 
                           edge=0.9, cross=NULL, subsample=NULL,
                           model=c("gamma-gamma", "gamma-norm", "norm-gamma", "norm-norm"),
                           cutoff=c("optimal", "intersect", "user"), sen=NULL, spc=NULL, 
                           progress=FALSE){
    # Check arguments
    method <- match.arg(method)
    model <- match.arg(model)
    cutoff <- match.arg(cutoff)
    
    # Subsample input distances
    if(!is.null(subsample)) {
        subsample <- min(length(distances), subsample)
        distances <- sample(distances, subsample, replace=FALSE)
    }
    
    if (method == "gmm") {
        if (cutoff == "user"){
            if (is.null(sen) & is.null(spc)) {
                cat("Error: one of 'sen' or 'spc' values should be specified.")
                output <- NA
            } else if (!is.null(sen) & !is.null(spc)) {
                cat("Error: only one of 'sen' or 'spc' values can be specified.")
                output <- NA
            } else {
                output <- gmmFit(ent=distances, edge=edge, cross=cross, model=model, cutoff=cutoff, 
                                 sen=sen, spc=spc, progress=progress)
            }
        } else {
            output <- gmmFit(ent=distances, edge=edge, cross=cross, model=model, cutoff=cutoff, 
                             sen=sen, spc=spc, progress=progress)
        }
    } else if (method == "density") {
        output <- smoothValley(distances)
    } else {
        cat("Error: assigned method has not been found.\n")
        output <- NA
    }
    
    return(output)
}


# Find distance threshold with \code{"density"} Method
#
# Infer value of the minimum between the two modes in a bimodal distribution.
#
# @param    distances  numeric vector of distances.
# 
# @return   Returns distance threshold that separates two modes of the input distribution.
#
# @details
# The distance to nearest neighbor can be used to estimate a threshold for assigning Ig
# sequences to clonal groups. A histogram of the resulting vector is often bimodal, 
# with the ideal threshold being a value that separates the two modes. This function takes 
# as input a vector of such distances and infers the ideal threshold.
# 
# @seealso  
# \itemize{
# \item     See \link{distToNearest} for details on generating the input distance vector.
# \item         See \link{gmmFit} for a different threshold inference methodology.
# \item           See \link{findThreshold} to switch between available methods.
#}
# 
# 
# @examples
# # Subset example data to one sample as a demo
# data(ExampleDb, package="alakazam")
# db <- subset(ExampleDb, sample_id == "-1h")
# 
# # Use genotyped V assignments, HS1F model, and normalize by junction length
# dist_hs1f <- distToNearest(db, sequenceColumn="junction", vCallColumn="v_call_genotyped",
#                            jCallColumn="j_call",                
#                            model="hs1f", first=FALSE, normalize="len")
#                  
# # using findThreshold switch
# threshold <- findThreshold(dist_hs1f$dist_nearest, method="density")
# # or
# threshold <- smoothValley(dist_hs1f$dist_nearest)
#                            
# # Plot histogram of non-NA distances
# p1 <- ggplot(data=subset(dist_hs1f, !is.na(dist_nearest))) + theme_bw() + 
#     ggtitle("Distance to nearest: hs1f") + xlab("distance") +
#     geom_histogram(aes(x=dist_nearest), binwidth=0.025, 
#                    fill="steelblue", color="white") + 
#     geom_vline(xintercept=threshold, linetype="dashed")
# plot(p1)
#
# @export
 smoothValley <- function(distances) {
    # Remove NA, NaN, and infinite distances
    distances <- distances[!is.na(distances) & !is.nan(distances) & !is.infinite(distances)]

    # Gaussian distribution bandwidth scale parameter
    # gaussian_scaling <- (1/(4 * pi))^(1/10)
    
    # Ideal bandwidth
    bandwidth <- kedd::h.ucv(unique(distances), 4)$h
    #bandwidth <- kedd::h.ucv(distances, 4)$h
    #bandwidth <- ks::hucv(unique(distances), deriv.order=4)
    
    # Density estimate
    dens <- KernSmooth::bkde(distances, bandwidth=bandwidth, canonical=TRUE)
    #dens <- KernSmooth::bkde(distances, bandwidth=bandwidth)
    xdens <- dens$x
    ydens <- dens$y
    #dens <- ks::kde(distances, h=bandwidth*gaussian_scaling, binned=TRUE)
    #xdens <- dens$eval.points
    #ydens <- dens$estimate
    
    # Find threshold
    tryCatch(threshold <- xdens[which(diff(sign(diff(ydens))) == 2)[1] + 1], 
             error = function(e) {
                 warning('No minimum was found between two modes.')
                 return(NULL) })
    
    results <- new("DensityThreshold",
                   x=distances,
                   bandwidth=bandwidth,
                   xdens=xdens,
                   ydens=ydens,
                   threshold=threshold)
    
    return(results)
}



# Find distance threshold with Gaussian Mixture Method
#
# Fits a bimodal distribution with two Gaussian functions and calculates maximum of the average of the 
# Sensitivity plus Specificity corresponding to the Gaussian distributions.
# 
# @param    ent         numeric vector of distances returned from \link{distToNearest} function.
# @param    edge        upper range (a fraction of the data density) to rule initialization of 
#                       Gaussian fit parameters. Default value is equal to \eqn{90}\% of the entries.
# @param    cross       a supplementary info (numeric vector) invoked from \link{distToNearest} 
#                       function, to support initialization of the Gaussian fit parameters. 
# @param    progress    if \code{TRUE} print progress.
#
# @return   returns an object including optimum "\code{threshold}" cut and the Gaussian fit parameters, 
#           such as mixing proportion ("\code{omega1}" and "\code{omega2}"), mean ("\code{mu1}" and "\code{mu2}"), 
#           and standard deviation ("\code{sigma1}" and "\code{sigma2}"). Returns "\code{NULL}" if no fit has found.         
#
# @seealso  
# \itemize{
# \item     See \link{distToNearest} for details on generating the input distance vector.
# \item     See \link{smoothValley} for a different threshold inference methodology.
# \item     See \link{findThreshold} to switch between available methods.
#}
#
#
# @details This function follows a Gaussian Mixture Model (GMM) procedure, 
#          including the Expectation Maximization (EM) algorithm, for learning the parameters  
#          of two univariate Gaussians which fit the bimodal distribution entries. 
#          Retrieving the fit parameters, it then calculates, analytically, the optimum threshold, 
#          where the average of the Sensitivity plus Specificity reaches its maximum. This threshold 
#          can be then invoked for assigning Ig sequences to clonal groups.
#
# @examples
# # Subset example data to one sample as a demo
# data(ExampleDb, package="alakazam")
# db <- subset(ExampleDb, sample_id == "-1h")
#
# # Use nucleotide Hamming distance and normalize by junction length
# db <- distToNearest(db, sequenceColumn="junction", vCallColumn="v_call_genotyped",
#                     jCallColumn="j_call", model="ham", first=FALSE, normalize="len", nproc=1)
#                             
# # To find the Threshold cut use either findThreshold-switch
# output <- findThreshold(db$dist_nearest, method="gmm", edge=0.9)
# # or 
# output <- gmmFit(db$dist_nearest, edge=0.9) 
gmmFit <- function(ent, edge=0.9, cross=NULL, model, cutoff, sen, spc, progress=FALSE) {
    
    #************* Filter Unknown Data *************#
    ent <- ent[!is.na(ent) & !is.nan(ent) & !is.infinite(ent)]
    if (is.null(cross)) {
        m <- FALSE
    } else {
        m <- mean(cross, na.rm = TRUE) 
    }
    
    #************* Defult edge *************#
    cut <- edge*length(ent)
    
    #************* Define Scan Step For Initializing *************#
    if (ent[which.max(ent)] <= 5) {
        scan_step <- 0.1
    } else {
        scan_step <- 1
    }
    
    #************* Print some info *************#
    if (progress) {
        valley_loc <- 0
        while (1) {
            valley_loc <- valley_loc + scan_step
            if ( length(ent[ent<=valley_loc]) > cut ) break
        }
        n_iter <- ceiling(valley_loc/scan_step)-1
        cat("      STEP> ", "Parameter initialization\n", sep="")
        cat("    VALUES> ", length(ent), "\n", sep="")
        cat("ITERATIONS> ", n_iter, "\n", sep="")
        pb <- progressBar(n_iter)
    }
    
    #*************  set rand seed *************#
    set.seed(NULL)
    
    #*************  define Number of Gaussians *************#
    num_G <- 2
    
    vec.omega1 <- 0; vec.omega2 <- 0
    vec.mu1 <- 0;    vec.mu2 <- 0
    vec.sigma1 <- 0; vec.sigma2 <- 0
    vec.lkhood <- 0
    valley.itr <- 0
    valley_loc <- 0
    nEve <- length(ent)
    
    while (1) {        
        #*************  guess the valley loc *************#
        valley_loc <- valley_loc + scan_step
        if ( length(ent[ent<=valley_loc]) > cut ) break
        
        #*************  Choosing Random Omega *************#
        omega <- runif(1)
        omega <- c(omega, 1.-omega)
        
        #*************  Choosing Random Mean *************#
        mu_int <- mean(ent[ent<=valley_loc])
        mu_int <- c(mu_int, mean(ent[ent>valley_loc]))
        
        #*************  Choosing Random Sigma *************#
        sigma_int <- sd(ent[ent<valley_loc])
        sigma_int <- c(sigma_int, sd(ent[ent>valley_loc]))
        
        #*************  EM Algorithm *************#
        temp_lk <- 0
        itr <- 0
        while (1){
            mu <- 0
            sigma <- 0
            for (j in 1:num_G){
                mu[j] <- mu_int[j]
                sigma[j] <- sigma_int[j]
            }
            
            #*************  E-step Expectation *************#
            resp <- array(0, dim=c(nEve,num_G))
            for(i in 1:nEve){
                for (j in 1:num_G)
                    resp[i,j] <- omega[j]*dnorm(ent[i], mu[j], sigma[j])
                resp[i,] <- resp[i,]/sum(resp[i,])
            }
            
            #*************  M-step Maximization *************#
            for (j in 1:num_G){
                m_c <- sum(resp[,j])
                
                omega[j] <- m_c / nEve
                
                mu[j] <- sum(resp[,j]*ent) 
                mu[j] <- mu[j] / m_c
                
                sigma[j] <- sum(resp[,j]*(ent-mu[j])*(ent-mu[j]))
                sigma[j] <- sigma[j] / m_c
                sigma[j] <- sqrt(sigma[j])
            }
            
            #*************  Log-likelihood calculation *************#
            log_lk <- 0.
            for (i in 1:nEve){
                s <- 0
                for (j in 1:num_G)
                    s <- s + omega[j]*dnorm(ent[i], mu[j], sigma[j])
                log_lk <- log_lk + log(s, base = exp(1))
            }
            log_lk_err <- abs(log_lk - temp_lk)
            itr = itr + 1
            #print(paste0("scaned: ", valley_loc, " itr # ", itr, " -> ", log_lk_err))
            if (is.na(log_lk_err) | is.nan(log_lk_err) | is.infinite(log_lk_err)) break
            if (log_lk_err < 1.e-7) break
            temp_lk <- log_lk;
        }
        
        #************************************************************# 
        #*************  JUST FOR VISUALIZATION PURPOSES *************#
        # print(paste0("scaned: ", valley_loc, " --------> Log-Likelihood: ", log_lk))
        # if (ent[which.min(ent)] >= 0 & ent[which.max(ent)] <= 5) {
        #   h_min <- 0.0
        #   h_max <- 1
        #   dh = 0.02
        # } else {
        #   h_min <- 0.0
        #   h_max <- ent[which.max(ent)]
        #   dh = 1
        # }
        # h <- hist(ent, plot = FALSE, breaks=seq(h_min, h_max, by=dh))
        # plot(h, freq=FALSE, col="steelblue", border="white", xlim=c(h_min, h_max))
        # curve(omega[1]*dnorm(x, mu[1], sigma[1]), add=TRUE, col="darkblue", lwd=2, xlim = c(h_min, h_max))
        # curve(omega[2]*dnorm(x, mu[2], sigma[2]), add=TRUE, col="darkred", lwd=2, xlim = c(h_min, h_max))
        #************************************************************#
        #************************************************************#
        if (!is.na(log_lk_err) & !is.nan(log_lk_err) & !is.infinite(log_lk_err)){
            if (!as.logical(m)){
                valley.itr <- valley.itr + 1
                vec.omega1[valley.itr] <- omega[1]
                vec.omega2[valley.itr] <- omega[2]
                vec.mu1[valley.itr] <- mu[1]
                vec.mu2[valley.itr] <- mu[2]
                vec.sigma1[valley.itr] <- sigma[1]
                vec.sigma2[valley.itr] <- sigma[2]
                vec.lkhood[valley.itr] <- log_lk
            } else if ((mu[1]< m & m < mu[2]) | (mu[2]< m & m < mu[1]) | (mu[1]< m & mu[2]< m) ){
                valley.itr <- valley.itr + 1
                vec.omega1[valley.itr] <- omega[1]
                vec.omega2[valley.itr] <- omega[2]
                vec.mu1[valley.itr] <- mu[1]
                vec.mu2[valley.itr] <- mu[2]
                vec.sigma1[valley.itr] <- sigma[1]
                vec.sigma2[valley.itr] <- sigma[2]
                vec.lkhood[valley.itr] <- log_lk
            }
        }
        
        # Update progress
        if (progress) { pb$tick() }
    }
    
    if (valley.itr != 0) {
        # MaxLoc <- which.max(vec.lkhood)
        MaxLoc <- which.max(abs(vec.lkhood))
        
        omega[1] <- vec.omega1[MaxLoc]; omega[2] <- vec.omega2[MaxLoc]
        mu[1] <- vec.mu1[MaxLoc];       mu[2] <- vec.mu2[MaxLoc]
        sigma[1] <- vec.sigma1[MaxLoc]; sigma[2] <- vec.sigma2[MaxLoc]
        
        # Invoke Gaussians parameters
        omega.gmm <- c(omega[1], omega[2])
        mu.gmm    <- c(mu[1], mu[2]) 
        sigma.gmm <- c(sigma[1], sigma[2]) 
        
        fit_results <- rocSpace(ent=ent, omega.gmm=omega.gmm , mu.gmm=mu.gmm, sigma.gmm=sigma.gmm, 
                               model=model, cutoff=cutoff, sen=sen, spc=spc, progress=progress)
        results <- new("GmmThreshold",
                       x=ent,
                       model=model,
                       cutoff=cutoff,
                       a1=fit_results@a1, 
                       b1=fit_results@b1, 
                       c1=fit_results@c1, 
                       a2=fit_results@a2,  
                       b2=fit_results@b2,  
                       c2=fit_results@c2, 
                       loglk=fit_results@loglk,
                       threshold=fit_results@threshold,
                       sensitivity=fit_results@sensitivity, 
                       specificity=fit_results@specificity,
                       pvalue=fit_results@pvalue)
    } else {    
        print("Error: No fit found")
        results <- NULL
    }
    
    return(results)
}


rocSpace <- function(ent, omega.gmm, mu.gmm, sigma.gmm, model, cutoff, sen, spc, progress=FALSE) {
    func <- model
    bits <- strsplit(func,'-')[[1]]

    # Define mixture Function properties
    if (bits[1] == "norm"){
        func1.0 <- round(omega.gmm[1], digits = 3)                     # -> prob: omega
        func1.1 <- mu.gmm[1]                                           # -> mean: mu
        func1.2 <- sigma.gmm[1]                                        # -> sd: sigma
    } else if (bits[1] == "gamma"){
        func1.0 <- round(omega.gmm[1], digits = 3)                     # -> prob: omega
        func1.1 <- (mu.gmm[1]/sigma.gmm[1])*(mu.gmm[1]/sigma.gmm[1])   # -> shape: k
        func1.2 <- sigma.gmm[1]*sigma.gmm[1]/mu.gmm[1]                 # -> scale: theta
    }
    
    if (bits[2] == "norm"){ 
        func2.1 = mu.gmm[2]                                            # -> mean: mu
        func2.2 = sigma.gmm[2]                                         # -> sd: sigma
    } else if (bits[2] == "gamma"){
        func2.1 <- (mu.gmm[2]/sigma.gmm[2])*(mu.gmm[2]/sigma.gmm[2])   # -> shape: k
        func2.2 <- sigma.gmm[2]*sigma.gmm[2]/mu.gmm[2]                 # -> scale: theta
    }
    
    # Save mixture Function properties
    gmmfunc1.1 <- func1.1
    gmmfunc1.2 <- func1.2
    gmmfunc2.1 <- func2.1
    gmmfunc2.2 <- func2.2
    
    set.seed(NULL)
    # options(warn=-1)
    LOG_LIK<-0
    
    if (progress) {
        cat("      STEP> ", "Fitting ", func, "\n", sep="")
        pb <- progressBar(15)
    }
    for (i in 1:15) {
        #itr<-1
        key<-FALSE
        while (!key){
            # print(paste0(i,":",itr))
            # Fit mixture Functions
            MixModel <- try(suppressWarnings(fitdistr(na.exclude(ent), mixFunction, 
                                     first_curve = bits[1], second_curve = bits[2], 
                                     start=list(omega = func1.0, 
                                                func1.1 = func1.1, func1.2 = func1.2,
                                                func2.1 = func2.1, func2.2 = func2.2), 
                                     lower = c(0.001, 0.001, 0.001, 0.001, 0.001), upper = c(0.999, +Inf, +Inf, +Inf, +Inf))), 
                            silent = TRUE)
            if (inherits(MixModel, "try-error")) {
                func1.0 <- runif(1)
                func1.1 <- abs(gmmfunc1.1 + sample(c(-1,1), 1)*runif(1))
                func1.2 <- abs(gmmfunc1.2 + sample(c(-1,1), 1)*runif(1))
                func2.1 <- abs(gmmfunc2.1 + sample(c(-1,1), 1)*runif(1))
                func2.2 <- abs(gmmfunc2.2 + sample(c(-1,1), 1)*runif(1))
                #itr<-itr+1
                next
            } else if ( (bits[1] == "norm"  & bits[2] == "gamma" & MixModel$estimate[[2]] > MixModel$estimate[[4]] * MixModel$estimate[[5]]) |
                        (bits[1] == "gamma" & bits[2] == "norm"  & MixModel$estimate[[2]] * MixModel$estimate[[3]] > MixModel$estimate[[4]]) |
                        MixModel$estimate[[1]] == 0.001 |
                        MixModel$estimate[[1]] == 0.999) {
                func1.0 <- runif(1)
                func1.1 <- abs(gmmfunc1.1 + sample(c(-1,1), 1)*runif(1))
                func1.2 <- abs(gmmfunc1.2 + sample(c(-1,1), 1)*runif(1))
                func2.1 <- abs(gmmfunc2.1 + sample(c(-1,1), 1)*runif(1))
                func2.2 <- abs(gmmfunc2.2 + sample(c(-1,1), 1)*runif(1))
                # print("here")
                #itr<-itr+1
                next
            } else {
                key<-TRUE
            }
        }
        # print(paste0(func, " fit done. Loglik= ", round(MixModel$loglik, digits = 2)))
        # Invoke fit parameters
        # log_lik <- round(MixModel$loglik, digits = 2)
        log_lik <- round(abs(MixModel$loglik), digits = 2)
        if (log_lik > LOG_LIK){
            LOG_LIK <- log_lik
            
            FUNC1.0 <- MixModel$estimate[[1]]
            FUNC1.1 <- MixModel$estimate[[2]] 
            FUNC1.2 <- MixModel$estimate[[3]]
            
            FUNC2.0 <- 1. - MixModel$estimate[[1]] 
            FUNC2.1 <- MixModel$estimate[[4]] 
            FUNC2.2 <- MixModel$estimate[[5]]
        }
        
        # New fit parameters for next loop
        func1.0 <- runif(1)
        func1.1 <- abs(gmmfunc1.1 + sample(c(-1,1), 1)*runif(1))
        func1.2 <- abs(gmmfunc1.2 + sample(c(-1,1), 1)*runif(1))
        func2.1 <- abs(gmmfunc2.1 + sample(c(-1,1), 1)*runif(1))
        func2.2 <- abs(gmmfunc2.2 + sample(c(-1,1), 1)*runif(1))
        
        # if (i==1 & itr == 1) break
        if (progress) { pb$tick() }
    }
    # options(warn=0)

    # Invoke best fit parameters
    log_lik  <- LOG_LIK
    
    func1.0 <- FUNC1.0
    func1.1 <- FUNC1.1 
    func1.2 <- FUNC1.2
    
    func2.0 <- FUNC2.0 
    func2.1 <- FUNC2.1
    func2.2 <- FUNC2.2
    
    # order fit parameters
    if (bits[1]=="norm" & bits[2]=="norm" & func1.1>func2.1) {
        FUNC0 <- func1.0
        FUNC1 <- func1.1 
        FUNC2 <- func1.2
        
        func1.0 <- func2.0 
        func1.1 <- func2.1
        func1.2 <- func2.2
        
        func2.0 <- FUNC0 
        func2.1 <- FUNC1
        func2.2 <- FUNC2
    } else if (bits[1]=="gamma" & bits[2]=="gamma" & func1.1*func1.2>func2.1*func2.2) {
        FUNC0 <- func1.0
        FUNC1 <- func1.1 
        FUNC2 <- func1.2
        
        func1.0 <- func2.0 
        func1.1 <- func2.1
        func1.2 <- func2.2
        
        func2.0 <- FUNC0 
        func2.1 <- FUNC1
        func2.2 <- FUNC2
    }
    
    # domain [t1,t2] under distribution
    t1<-min(ent)
    t2<-max(ent)
    
    # domain [minInt,maxInt] to search for opt and root
    if (bits[1] == "norm") {
        minInt<-func1.1 
    } else if (bits[1] == "gamma") {
        minInt<-func1.1*func1.2
    }
    
    if (bits[2] == "norm") {
        maxInt<-func2.1 
    } else if (bits[2] == "gamma") {
        maxInt<-func2.1*func2.2
    }    
    
    if (cutoff == "optimal"){
        # Calculate optimum
        opt <- optimize(avgSenSpc, interval = c(minInt, maxInt), tol=1e-8, maximum = TRUE, 
                        t1=t1, t2=t2, 
                        first_curve = bits[1], second_curve = bits[2], 
                        func1.0=func1.0, func1.1=func1.1, func1.2=func1.2, 
                        func2.0=func2.0, func2.1=func2.1, func2.2=func2.2)
        threshold <- opt$maximum
    } else if (cutoff == "intersect") {
        # Calculate intersection
        intxn <- uniroot(intersectPoint, interval = c(minInt, maxInt), tol=1e-8, extendInt="yes",
                         first_curve = bits[1], second_curve = bits[2], 
                         func1.0=func1.0, func1.1=func1.1, func1.2=func1.2, 
                         func2.0=func2.0, func2.1=func2.1, func2.2=func2.2)
        threshold <- intxn$root
    } else if (cutoff == "user") {
        user <- uniroot(userDefineSenSpc, interval = c(t1, t2), tol=1e-8, extendInt="no",
                      t1=t1, t2=t2, 
                      first_curve = bits[1], second_curve = bits[2], 
                      sen = sen, spc = spc,
                      func1.0=func1.0, func1.1=func1.1, func1.2=func1.2, 
                      func2.0=func2.0, func2.1=func2.1, func2.2=func2.2)
        threshold <- user$root
    }
    
    # Calculate Sensitivity and Specificity
    if (bits[1]=="norm") {
        TP = normArea(t1=t1, t2=threshold, omega=func1.0, mu=func1.1, sigma=func1.2)
    } else if (bits[1]=="gamma") {
        TP = gammaArea(t1=t1, t2=threshold, omega=func1.0, k=func1.1, theta=func1.2)
    }
    if (bits[1]=="norm") {
        FN = normArea(t1=threshold, t2=t2, omega=func1.0, mu=func1.1, sigma=func1.2)
    } else if (bits[1]=="gamma") {
        FN = gammaArea(t1=threshold, t2=t2, omega=func1.0, k=func1.1, theta=func1.2)
    }
    
    if (bits[2]=="norm") {
        TN = normArea(t1=threshold, t2=t2, omega=func2.0, mu=func2.1, sigma=func2.2)
    } else if (bits[2]=="gamma") {
        TN = gammaArea(t1=threshold, t2=t2, omega=func2.0, k=func2.1, theta=func2.2)
    }
    if (bits[2]=="norm") {
        FP = normArea(t1=t1, t2=threshold, omega=func2.0, mu=func2.1, sigma=func2.2)
    } else if (bits[2]=="gamma") {
        FP = gammaArea(t1=t1, t2=threshold, omega=func2.0, k=func2.1, theta=func2.2)
    }
    
    sensitivity <- TP/(TP+FN)
    specificity <- TN/(TN+FP)  
    
    # Hartigans dip statistic (HDS) test
    invisible(capture.output(pvalue <- dip.test(ent)$p.value[[1]], type="message"))
    
    fit_results <- new("GmmThreshold",
                       x=numeric(), model=character(), cutoff=character(),
                       a1=func1.0, b1=func1.1, c1=func1.2, 
                       a2=func2.0, b2=func2.1, c2=func2.2,
                       loglk=log_lik, threshold=threshold,
                       sensitivity=sensitivity, specificity=specificity,
                       pvalue=pvalue)
    
    return(fit_results)
}

# Calculates the area (integral) bounded 
# in domain[t1,t2] under Gamma distribution
gammaArea <- function (t1, t2, omega, k, theta){
    trm1 <- pgamma(t1/theta, shape=k, lower.tail=FALSE) * gamma(k)
    trm2 <- pgamma(t2/theta, shape=k, lower.tail=FALSE) * gamma(k)
    
    area <- omega*(trm1 - trm2)/gamma(k)
    
    return(area)
}

# Calculates the area (integral) bounded 
# in domain[t1,t2] under Normal distribution
normArea <- function (t1, t2, omega, mu, sigma){
    erf1 <- (t1-mu)/(sqrt(2)*sigma)
    erf1 <- 2*pnorm(erf1*sqrt(2)) - 1
    
    erf2 <- (t2-mu)/(sqrt(2)*sigma)
    erf2 <- 2*pnorm(erf2*sqrt(2)) - 1
    
    area <- sigma * omega * (-erf1 + erf2) / (2*sigma)
    
    return(area)
}

# find the optimum threshold using 
# optimize function fit
avgSenSpc <- function(t, t1=0, t2=0, first_curve=NULL, second_curve=NULL, 
                      func1.0 = 0, func1.1 = 0, func1.2 = 0,
                      func2.0 = 0, func2.1 = 0, func2.2 = 0) {
    
    if (first_curve == "norm") {
        TP <- normArea(t1=t1, t2=t, omega=func1.0, mu=func1.1, sigma=func1.2)
        FN <- normArea(t1=t, t2=t2, omega=func1.0, mu=func1.1, sigma=func1.2)
    } else if (first_curve == "gamma") {
        TP <- gammaArea(t1=t1, t2=t, omega=func1.0, k=func1.1, theta=func1.2)
        FN <- gammaArea(t1=t, t2=t2, omega=func1.0, k=func1.1, theta=func1.2)
    }
    
    if (second_curve == "norm") {
        FP <- normArea(t1=t1, t2=t, omega=func2.0, mu=func2.1, sigma=func2.2)
        TN <- normArea(t1=t, t2=t2, omega=func2.0, mu=func2.1, sigma=func2.2)
    } else if (second_curve == "gamma") {
        FP <- gammaArea(t1=t1, t2=t, omega=func2.0, k=func2.1, theta=func2.2)
        TN <- gammaArea(t1=t, t2=t2, omega=func2.0, k=func2.1, theta=func2.2)
    }
    
    SEN <- TP/(TP + FN)        
    SPC <- TN/(TN + FP)
    
    return((SEN + SPC)/2)
}

# Intersection Function
intersectPoint <- function(t, first_curve=NULL, second_curve=NULL, 
                         func1.0 = 0, func1.1 = 0, func1.2 = 0,
                         func2.0 = 0, func2.1 = 0, func2.2 = 0) {
    
    if (first_curve == "norm") {
        fit1 <- func1.0*dnorm(t, mean = func1.1, sd = func1.2)
    } else if (first_curve == "gamma") {
        fit1 <- func1.0*dgamma(t, shape = func1.1, scale = func1.2)
    }
    
    if (second_curve == "norm") {
        fit2 <- func2.0*dnorm(t, mean = func2.1, sd = func2.2)
    } else if (second_curve == "gamma") {
        fit2 <- func2.0*dgamma(t, shape = func2.1, scale = func2.2)
    }
    
    return(fit1 - fit2)
}

# useDefineSenSpc
userDefineSenSpc <- function(t, t1=0, t2=0, first_curve=NULL, second_curve=NULL, 
                             sen=NULL, spc=NULL,
                             func1.0=0, func1.1=0, func1.2=0,
                             func2.0=0, func2.1=0, func2.2=0) {
    if (!is.null(sen)) {
        if (first_curve == "norm") {
            TP <- normArea(t1=t1, t2=t, omega=func1.0, mu=func1.1, sigma=func1.2)
            FN <- normArea(t1=t, t2=t2, omega=func1.0, mu=func1.1, sigma=func1.2)
        } else if (first_curve == "gamma") {
            TP <- gammaArea(t1=t1, t2=t, omega=func1.0, k=func1.1, theta=func1.2)
            FN <- gammaArea(t1=t, t2=t2, omega=func1.0, k=func1.1, theta=func1.2)
        }
        threshold <- (TP/(TP+FN)) - sen
    } else if (!is.null(spc)) {
        if (second_curve == "norm") {
            FP <- normArea(t1=t1, t2=t, omega=func2.0, mu=func2.1, sigma=func2.2)
            TN <- normArea(t1=t, t2=t2, omega=func2.0, mu=func2.1, sigma=func2.2)
        } else if (second_curve == "gamma") {
            FP <- gammaArea(t1=t1, t2=t, omega=func2.0, k=func2.1, theta=func2.2)
            TN <- gammaArea(t1=t, t2=t2, omega=func2.0, k=func2.1, theta=func2.2)
        }
        threshold <- (TN/(TN+FP)) - spc
    }
    
    return(threshold)
}

# Mixture Functions
mixFunction <- function(t, first_curve=NULL, second_curve=NULL,
                         omega = 0, 
                         func1.1 = 0, func1.2 = 0,
                         func2.1 = 0, func2.2 = 0) {
    
    if (first_curve == "norm"){
        r <- omega*dnorm(t, mean=func1.1, sd=func1.2)
    } else if (first_curve == "gamma") {
        r <- omega*dgamma(t, shape=func1.1, scale=func1.2)
    } 
    
    if (second_curve == "norm"){
        r <- r + (1-omega)*dnorm(t, mean=func2.1, sd=func2.2)
    } else if (second_curve == "gamma") {
        r <- r + (1-omega)*dgamma(t, shape=func2.1, scale=func2.2)
    } 
}

#' Plot findThreshold results for the gmm method
#' 
#' \code{plotGmmThreshold} plots the results from \code{"gmm"} method of 
#' \link{findThreshold}, including the Gaussian distributions, input nearest neighbor 
#' distance histogram, and threshold selected.
#'
#' @param    data      \link{GmmThreshold} object output by the \code{"gmm"} method 
#'                     of \link{findThreshold}.
#' @param    cross     numeric vector of distances from \link{distToNearest} to draw as a
#'                     histogram below the \code{data} histogram for comparison purposes.
#' @param    xmin      minimum limit for plotting the x-axis. If \code{NULL} the limit will 
#'                     be set automatically.
#' @param    xmax      maximum limit for plotting the x-axis. If \code{NULL} the limit will 
#'                     be set automatically.
#' @param    breaks    number of breaks to show on the x-axis. If \code{NULL} the breaks will 
#'                     be set automatically.
#' @param    binwidth  binwidth for the histogram. If \code{NULL} the binwidth 
#'                     will be set automatically.
#' @param    title     string defining the plot title.
#' @param    size      numeric value for lines in the plot.
#' @param    silent    if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                     object; if \code{FALSE} draw the plot.
#' @param    ...       additional arguments to pass to ggplot2::theme.
#' 
#' @return   A ggplot object defining the plot.
#'
#' @seealso  See \link{GmmThreshold} for the the input object definition and 
#'           \link{findThreshold} for generating the input object. See 
#'           \link{distToNearest} calculating nearest neighbor distances.
#'
#' @examples
#' \donttest{
#' # Subset example data to one sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, sample_id == "-1h")
#'
#' # Use nucleotide Hamming distance and normalize by junction length
#' db <- distToNearest(db, sequenceColumn="junction", vCallColumn="v_call_genotyped",
#'                     jCallColumn="j_call", model="ham", normalize="len", nproc=1)
#' 
#' # To find the threshold cut, call findThreshold function for "gmm" method.
#' output <- findThreshold(db$dist_nearest, method="gmm", model="norm-norm", cutoff="opt")
#' print(output)
#' 
#' # Plot results
#' plotGmmThreshold(output, binwidth=0.02)
#' }
#' @export
plotGmmThreshold <- function(data, cross=NULL, xmin=NULL, xmax=NULL, breaks=NULL, 
                             binwidth=NULL, title=NULL, size=1, silent=FALSE, ...) {
    # Define histogram data.frame and threshold
    xdf <- data.frame(x=data@x)
    
    # Generate curves
    gx <- seq(min(xdf$x), max(xdf$x), by=0.002)
    bits <- strsplit(data@model,'-')[[1]]
    if (bits[1] == "norm") {
        fit1 <- data.frame(x=gx, y=data@a1*dnorm(gx, mean=data@b1, sd=data@c1))
    } else if (bits[1] == "gamma") {
        fit1 <- data.frame(x=gx, y=data@a1*dgamma(gx, shape=data@b1, scale=data@c1))
    }
    
    if (bits[2] == "norm") {
        fit2 <- data.frame(x=gx, y=data@a2*dnorm(gx, mean = data@b2, sd=data@c2))
    } else if (bits[2] == "gamma") {
        fit2 <- data.frame(x=gx, y=data@a2*dgamma(gx, shape = data@b2, scale=data@c2))
    }
    
    # ggplot workaround
    if (is.null(xmin)) { xmin <- NA }
    if (is.null(xmax)) { xmax <- NA }
    
    # Plot distToNearest distribution plus Gaussian fits
    p <- ggplot(xdf, aes_string(x="x")) +
        baseTheme() + 
        xlab("Distance") + 
        ylab("Density") +
        geom_histogram(aes_string(y="..density.."), binwidth=binwidth, 
                       fill="gray40", color="white") +
        geom_line(data=fit1, aes_string(x="x", y="y"), color="darkslateblue", size=size) +
        geom_line(data=fit2, aes_string(x="x", y="y"), color="darkslateblue", size=size) +
        geom_vline(xintercept=data@threshold, color="firebrick", 
                   linetype="longdash", size=size)
    
    # Add cross histogram
    if (!is.null(cross)) {
        cdf <- data.frame(x=cross[is.finite(cross)])
        p <- p + geom_histogram(data=cdf, aes_q(x=~x, y=~-(..density..)), binwidth=binwidth, 
                                fill="gray40", color="white", position="identity") +
            scale_y_continuous(labels=abs)
    }
    
    # Add x limits
    if (is.null(breaks) & (!is.na(xmin) | !is.na(xmax))) {
        p <- p + xlim(xmin, xmax)
    }
    # Set breaks
    if (!is.null(breaks)) {
        p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=breaks),
                                    limits=c(xmin, xmax))
    }
    # Add Title
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    
    # Add additional theme elements
    p <- p + do.call(theme, list(...))
    
    # Plot
    if (!silent) {
        plot(p)
    } else {
        return(p)
    }
}


#' Plot findThreshold results for the density method
#' 
#' \code{plotDensityThreshold} plots the results from \code{"density"} method of 
#' \link{findThreshold}, including the smoothed density estimate, input nearest neighbor 
#' distance histogram, and threshold selected.
#'                           
#' @param    data      \link{DensityThreshold} object output by the \code{"density"} method 
#'                     of \link{findThreshold}.
#' @param    cross     numeric vector of distances from \link{distToNearest} to draw as a
#'                     histogram below the \code{data} histogram for comparison purposes.
#' @param    xmin      minimum limit for plotting the x-axis. If \code{NULL} the limit will 
#'                     be set automatically.
#' @param    xmax      maximum limit for plotting the x-axis. If \code{NULL} the limit will 
#'                     be set automatically.
#' @param    breaks    number of breaks to show on the x-axis. If \code{NULL} the breaks will 
#'                     be set automatically.
#' @param    binwidth  binwidth for the histogram. If \code{NULL} the binwidth 
#'                     will be set automatically to the bandwidth parameter determined by
#'                     \link{findThreshold}.
#' @param    title     string defining the plot title.
#' @param    size      numeric value for the plot line sizes.
#' @param    silent    if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                     object; if \code{FALSE} draw the plot.
#' @param    ...       additional arguments to pass to ggplot2::theme.
#' 
#' @return   A ggplot object defining the plot.
#'
#' @seealso  See \link{DensityThreshold} for the the input object definition and 
#'           \link{findThreshold} for generating the input object. See 
#'           \link{distToNearest} calculating nearest neighbor distances.
#'           
#' @examples
#' \donttest{
#' # Subset example data to one sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, sample_id == "-1h")
#'
#' # Use nucleotide Hamming distance and normalize by junction length
#' db <- distToNearest(db, sequenceColumn="junction", vCallColumn="v_call_genotyped",
#'                     jCallColumn="j_call", model="ham", normalize="len", nproc=1)
#' 
#' # To find the threshold cut, call findThreshold function for "gmm" method.
#' output <- findThreshold(db$dist_nearest, method="density")
#' print(output)
#' 
#' # Plot
#' plotDensityThreshold(output)
#' }
#' @export
plotDensityThreshold <- function(data, cross=NULL, xmin=NULL, xmax=NULL, breaks=NULL, 
                                 binwidth=NULL, title=NULL, size=1, silent=FALSE, ...) {
    # Define plot data.frames
    xdf <- data.frame(x=data@x)
    ddf <- data.frame(x=data@xdens, y=data@ydens)
    ddf <- ddf[ddf$x > 0, ]
    
    # Set binwidth
    if (is.null(binwidth)) { binwidth <- data@bandwidth }
    
    # ggplot workaround
    if (is.null(xmin)) { xmin <- NA }
    if (is.null(xmax)) { xmax <- NA }
    
    # Plot distToNearest distribution plus Gaussian fits
    p <- ggplot(xdf, aes_string(x="x")) +
        baseTheme() +
        xlab("Distance") + 
        ylab("Density") +
        geom_histogram(aes_string(y="..density.."), binwidth=binwidth, 
                       fill="gray40", color="white") +
        geom_line(data=ddf, aes_string(x="x", y="y"), 
                  color="darkslateblue", size=size) +
        geom_vline(xintercept=data@threshold, 
                   color="firebrick", linetype="longdash", size=size)
    
    # Add cross histogram
    if (!is.null(cross)) {
        cdf <- data.frame(x=cross[is.finite(cross)])
        p <- p + geom_histogram(data=cdf, aes_q(x=~x, y=~-(..density..)), binwidth=binwidth, 
                                fill="gray40", color="white", position="identity") +
             scale_y_continuous(labels=abs)
    }
    
    # Add x limits
    if (is.null(breaks) & (!is.na(xmin) | !is.na(xmax))) {
        p <- p + coord_cartesian(xlim = c(xmin, xmax))
    }
    # Set breaks
    if (!is.null(breaks)) {
        p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=breaks),
                                    limits=c(xmin, xmax))
    }
    # Add title
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    
    # Add additional theme elements
    p <- p + do.call(theme, list(...))
    
    # Plot
    if (!silent) {
        plot(p)
    } else {
        return(p)
    }
}

