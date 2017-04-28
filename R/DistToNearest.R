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
#' @slot   omega1       first Gaussain mixing proportion.
#' @slot   omega2       second Gaussain mixing proportion.
#' @slot   mu1          first Gaussian mean.
#' @slot   mu2          second Gaussain mean.
#' @slot   sigma1       first Gaussain standard deviation.
#' @slot   sigma2       second Gaussain standard deviation.
#' @slot   threshold    optimum threshold cut.
#'
#' @seealso      \link{findThreshold}
#'
#' @name         GmmThreshold-class
#' @rdname       GmmThreshold-class
#' @aliases      GmmThreshold
#' @exportClass  GmmThreshold
setClass("GmmThreshold",
         slots=c(x="numeric",
                 omega1="numeric", 
                 omega2="numeric",
                 mu1="numeric",      
                 mu2="numeric",
                 sigma1="numeric", 
                 sigma2="numeric",
                 threshold = "numeric"))


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
    numbOfMutation <- sum(fivemersWithMu)
    
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
pairwise5MerDist <- function(sequences, targetingDistance, 
                             symmetry=c("avg", "min")) {
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
    
    return(dist_mat)
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
nearestDist<- function(sequences, model=c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", 
                                          "hs1f_compat", "m1n_compat"),
                       normalize=c("none", "len", "mut"),
                       symmetry=c("avg", "min"),
                       crossGroups=NULL, mst=FALSE) {
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
        
        # Get distance matrix
        if (model == "ham") {
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=getDNAMatrix(gap=0))
        } else if (model == "aa") {
            seq_uniq <- setNames(alakazam::translateDNA(seq_uniq), seq_uniq)
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=getAAMatrix())
        } else if (model == "hh_s1f") {
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=HH_S1F_Distance)
        } else if (model == "mk_rs1nf") {
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=MK_RS1NF_Distance)
        } else if (model == "hh_s5f") {
            dist_mat <- pairwise5MerDist(seq_uniq, HH_S5F_Distance, symmetry=symmetry)
        } else if (model == "mk_rs5nf") {
            dist_mat <- pairwise5MerDist(seq_uniq, MK_RS5NF_Distance, symmetry=symmetry)
        } else if (model == "hs1f_compat") {
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=HS1F_Compat)
        } else if (model == "m1n_compat") {
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=M1N_Compat)
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
            this_group <- crossGroups[i]
            other_groups <-  which(crossGroups != this_group)
            other_seq <- unique(sequences[other_groups])
            other_idx <- match(other_seq, seq_uniq)
            this_idx <- match(sequences[i], seq_uniq)
            r <- dist_mat[this_idx, other_idx]
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
#' Get non-zero distance of every sequence (as defined by \code{sequenceColumn}) to its 
#' nearest sequence sharing same V gene, J gene, and sequence length.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing nucleotide sequences to compare. 
#'                           Also used to determine sequence length for grouping.
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
#' @param    strip_d         if \code{TRUE} V genes that are paralogous (IGKV3-11 vs IGKVD3-11) 
#'                           are made identical (IGKV3-11). if \code{FALSE} they are kept as 
#' 						     unique V gene variants.
#' @param    nproc           number of cores to distribute the function over.
#' @param    fields          additional fields to use for grouping.
#' @param    cross           columns for grouping to calculate distances across groups 
#'                           (self vs others).
#' @param    mst             if \code{TRUE}, return comma-separated branch lengths from minimum 
#'                           spanning tree.
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
#' there is only 1 sequence (as defined by \code{sequenceColumn}), \code{NA} is returned 
#' instead of a distance (since it has no neighbor). If for a given combination there are 
#' multiple sequences but only 1 unique sequence, (in which case every sequence in this 
#' group is the de facto nearest neighbor to each other, thus giving rise to distances 
#' of 0), \code{NA}s are returned instead of zero-distances.
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
#' db <- subset(ExampleDb, SAMPLE == "-1h")
#' 
#' # Use genotyped V assignments, Hamming distance, and normalize by junction length
#' dist <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", model="ham", 
#'                       first=FALSE, normalize="len")
#'                            
#' # Plot histogram of non-NA distances
#' p1 <- ggplot(data=subset(dist, !is.na(DIST_NEAREST))) + 
#'       theme_bw() + 
#'       ggtitle("Distance to nearest: Hamming") + 
#'       xlab("distance") +
#'       geom_histogram(aes(x=DIST_NEAREST), binwidth=0.025, 
#'                      fill="steelblue", color="white")
#' plot(p1)
#'
#' # Use human 5-mer model
#' dist <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", model="hh_s5f")
#' 
#' @export
distToNearest <- function(db, sequenceColumn="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL", 
                          model=c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", "m1n_compat", "hs1f_compat"), 
                          normalize=c("len", "none"), symmetry=c("avg", "min"), 
                          first=TRUE, strip_d=TRUE, nproc=1, fields=NULL, cross=NULL, mst=FALSE) {
							  
    # Hack for visibility of data.table and foreach index variables
    idx <- yidx <- .I <- NULL
    
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    
    # Check for valid columns
    columns <- c(sequenceColumn, vCallColumn, jCallColumn, fields, cross)
    columns <- columns[!is.null(columns)]
    
    check <- checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn))
    
    # Check for invalid characters
    valid_seq <- sapply(db[[sequenceColumn]], allValidChars, getCharsInModel(model))
    not_valid_seq <- which(!valid_seq)
    if (length(not_valid_seq)>0) {
        warning("Invalid sequence characters in the ", sequenceColumn, " column.",
                length(not_valid_seq)," sequence(s) removed")
        db <- db[valid_seq,]
    }    
    
    # Get targeting distance
    targeting_distance <- if (model == "hs5f") { 
        calcTargetingDistance(HH_S5F) 
        } else { NULL }

    # Parse V and J columns to get gene
    # cat("V+J Column parsing\n")
    if (first) {
        db$V <- getGene(db[[vCallColumn]], strip_d=strip_d)
        db$J <- getGene(db[[jCallColumn]], strip_d=strip_d)
    } else {
        db$V1 <- getGene(db[[vCallColumn]], first=FALSE, strip_d=strip_d)
        db$J1 <- getGene(db[[jCallColumn]], first=FALSE, strip_d=strip_d)
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
    
    # Export groups to the clusters
    if (nproc > 1) { 
        export_functions <- list("db",
                                 "groups", 
                                 "cross",
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
                                 "targeting_distance")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    }
    
   

    list_db <- foreach(idx=iterators::icount(lenGroups), .errorhandling='stop') %dopar% {
        db_group <- db[groups[[idx]], ]
        crossGroups <- NULL
        if (!is.null(cross)) {
            crossGroups <- db_group %>% dplyr::group_indices_(.dots=cross)
        }
        arrSeqs <-  as.vector(unlist(db[groups[[idx]], sequenceColumn]))
        db_group$TMP_DIST_NEAREST <- nearestDist(arrSeqs,
                                                 model=model,
                                                 normalize=normalize,
                                                 symmetry=symmetry,
                                                 crossGroups=crossGroups,
                                                 mst=mst)
        return(db_group)
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


#### Distance Threshold Detection ####

#' Find distance threshold
#'
#' \code{findThreshold} automtically determines and optimal threshold for clonal assignment of
#' Ig sequences using a vector of nearest neighbor distances. It provides two alternative methods 
#' using either a Guassian Mixture Model fit (\code{method="gmm"}) or kernel density 
#' fit (\code{method="density"}).
#'
#' @param    data       numeric vector containing nearest neighbor distances. 
#' @param    method     string defining the method to use for determining the optimal threshold.
#'                      One of \code{"gmm"} or \code{"density"}. See Details for methodological
#'                      descriptions.
#' @param    cutEdge    upper range (a fraction of the data density) to rule initialization of 
#'                      Gaussian fit parameters. Default value is equal to 90% of the entries.
#'                      Applies only to the \code{"gmm"} method.
#' @param    cross      supplementary nearest neighbor distance vector output from \link{distToNearest} 
#'                      for initialization of the Gaussian fit parameters. 
#'                      Applies only to the \code{"gmm"} method.
#' @param    subsample  number of distances to subsample for speeding up bandwidth inference.
#'                      Applies only to the \code{"density"} method. If \code{NULL} no subsampling
#'                      is performed. As bandwith inferrence is computationally expensive, subsampling
#'                      is recommended for large data sets.
#' 
#' @return   
#' \itemize{
#'   \item \code{"gmm"} method:      Returns a \link{GmmThreshold} object including the optimum 
#'                                   \code{threshold} and the Gaussian fit parameters.
#'   \item \code{"density"} method:  Returns a \link{DensityThreshold} object including the optimum 
#'                                   \code{threshold} and the density fit parameters.
#' }
#'
#' @details 
#' \itemize{ 
#'   \item \code{"gmm"}:     Performs a Gaussian Mixture Model (GMM) procedure, 
#'                           including the Expectation Maximization (EM) algorithm, for learning 
#'                           the parameters  of two univariate Gaussians which fit the bimodal 
#'                           distribution entries. Retrieving the fit parameters, it then calculates
#'                           the optimum threshold, where the average of the sensitivity plus 
#'                           specificity reaches its maximum.
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
#' @examples
#' \donttest{
#' # Subset example data to one sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, SAMPLE == "-1h")
#' 
#' # Use nucleotide Hamming distance and normalize by junction length
#' db <- distToNearest(db, model="ham", normalize="len", nproc=1)
#'                             
#' # Find threshold using the "gmm" method
#' output <- findThreshold(db$DIST_NEAREST, method="gmm")
#' print(output)
#' # Plot "gmm" method results
#' plot(output, binwidth=0.02)
#'
#' # Find threshold using the "density" method 
#' output <- findThreshold(db$DIST_NEAREST, method="density")
#' print(output)
#' # Plot "density" method results
#' plot(output)
#' }
#' @export
findThreshold <- function (data, method=c("gmm", "density"), cutEdge=0.9, cross=NULL, 
                           subsample=NULL){
  # Check arguments
  method <- match.arg(method)
  
  if (method == "gmm") {
    output <- gmmFit(data, cutEdge, cross)
  } else if (method == "density") {
    output <- smoothValley(data, subsample)
  } else {
    print("Error: assigned method has not been found")
    output <- NA
  }
  
  return(output)
}


# Find distance threshold with \code{"density"} Method
#
# Infer value of the minimum between the two modes in a bimodal distribution.
#
# @param    distances  numeric vector of distances.
# @param    subsample  number of distances to subsample for speeding up bandwidth inference.
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
# db <- subset(ExampleDb, SAMPLE == "-1h")
# 
# # Use genotyped V assignments, HS1F model, and normalize by junction length
# dist_hs1f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
#                            model="hs1f", first=FALSE, normalize="len")
#                  
# # using findThreshold switch
# threshold <- findThreshold(dist_hs1f$DIST_NEAREST, method="density")
# # or
# threshold <- smoothValley(dist_hs1f$DIST_NEAREST)
#                            
# # Plot histogram of non-NA distances
# p1 <- ggplot(data=subset(dist_hs1f, !is.na(DIST_NEAREST))) + theme_bw() + 
#     ggtitle("Distance to nearest: hs1f") + xlab("distance") +
#     geom_histogram(aes(x=DIST_NEAREST), binwidth=0.025, 
#                    fill="steelblue", color="white") + 
#     geom_vline(xintercept=threshold, linetype="dashed")
# plot(p1)
#
# @export
smoothValley <- function(distances, subsample=NULL) {
    # Remove NA, NaN, and infinite distances
    distances <- distances[!is.na(distances) & !is.nan(distances) & !is.infinite(distances)]
    # Subsample input distances
    if(!is.null(subsample)) {
        subsample <- min(length(distances), subsample)
        distances <- sample(distances, subsample, replace=FALSE)
    }
    # Ideal bandwidth
    h_ucv <- h.ucv(distances, 4)$h
    # Density estimate
    dens <- bkde(distances, bandwidth=h_ucv)
    # Find threshold
    tryCatch(threshold <- dens$x[which(diff(sign(diff(dens$y))) == 2)[1] + 1], 
             error = function(e) {
                 warning('No minimum was found between two modes.')
                 return(NULL) })
    
    results <- new("DensityThreshold",
                   x=distances,
                   bandwidth=h_ucv,
                   xdens=dens$x,
                   ydens=dens$y,
                   threshold=threshold)
    
    return(results)
}



# Find distance threshold with Gaussian Mixture Method
#
# Fits a bimodal distribution with two Gaussian functions and calculates maximum of the average of the 
# Sensitivity plus Specificity corresponding to the Gaussian distributions.
# 
# @param    ent         numeric vector of distances returned from \link{distToNearest} function.
# @param    cutEdge     upper range (a fraction of the data density) to rule initialization of 
#                       Gaussian fit parameters. Default value is equal to \eqn{90}\% of the entries.
# @param    cross       a supplementary info (numeric vector) invoked from \link{distToNearest} 
#                       function, to support initialization of the Gaussian fit parameters. 
#
# @return   returns an objrct including optimum "\code{threshold}" cut and the Gaussian fit parameters, 
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
# db <- subset(ExampleDb, SAMPLE == "-1h")
#
# # Use nucleotide Hamming distance and normalize by junction length
# db <- distToNearest(db, model="ham", first=FALSE, normalize="len", nproc=1)
#                             
# # To find the Threshold cut use either findThreshold-switch
# output <- findThreshold(db$DIST_NEAREST, method="gmm", cutEdge=0.9)
# # or 
# output <- gmmFit(db$DIST_NEAREST, cutEdge=0.9) 
# 
# # Retrieve outputs:
# omega <- c(output@omega1, output@omega2)
# mu <-    c(output@mu1,    output@mu2) 
# sigma <- c(output@sigma1, output@sigma2) 
# threshold <- output@threshold
# 
# # To check the quality of the fit performance and corresponding 
# # threshold location use "plotGmmThreshold" function in shazam. 
# 
# 
# @export
gmmFit <- function(ent, cutEdge=0.9, cross=NULL) {
    
    #************* Filter Unknown Data *************#
    ent <- ent[!is.na(ent) & !is.nan(ent) & !is.infinite(ent)]
    if (is.null(cross)) {
        m <- FALSE
    } else {
        m <- mean(cross, na.rm = TRUE) 
    }
    
    #************* Defult cutEdge *************#
    cut <- cutEdge*length(ent)
    
    #************* Define Scan Step For Initializing *************#
    if (ent[which.max(ent)] <= 5) {
        scan_step <- 0.1
    } else {
        scan_step <- 1
    }
    
    #************* Print some info *************#
    valley_loc <- 0
    while (1) {
        valley_loc <- valley_loc + scan_step
        if ( length(ent[ent<=valley_loc]) > cut ) break
    }
    cat(length(ent), "non-NA entries\n")
    cat("GMM done in", ceiling(valley_loc/scan_step), "iterations\n")
    
    #*************  set rand seed *************#
    set.seed(3000)
    
    #*************  define Number of Gaussians *************#
    num_G <- 2
    
    vec.omega1 <- 0; vec.omega2 <- 0
    vec.mu1 <- 0;    vec.mu2 <- 0
    vec.sigma1 <- 0; vec.sigma2 <- 0
    vec.lkhood <- 0
    valley.itr <- 0
    valley_loc <- 0
    nEve <- length(ent)
    cat("[")
    while (1) {
        cat("=")
        
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
    }
    cat("]\n")
    if (valley.itr != 0) {
        MaxLoc <- which.max(vec.lkhood)
        
        omega[1] <- vec.omega1[MaxLoc]; omega[2] <- vec.omega2[MaxLoc]
        mu[1] <- vec.mu1[MaxLoc];       mu[2] <- vec.mu2[MaxLoc]
        sigma[1] <- vec.sigma1[MaxLoc]; sigma[2] <- vec.sigma2[MaxLoc]  
        threshold <- gmmThreshold (omega, mu, sigma)
        
        #************* Print OutPuts *************#
        # print(paste0("Omega_1= ",    round(omega[1],digits = 6),    ",  Omega_2= ", round(omega[2],digits = 6)))
        # print(paste0("Mu_1= ",       round(mu[1],digits = 6),       ",  Mu_2= ",    round(mu[2],digits = 6)))
        # print(paste0("Sigma_1= ",    round(sigma[1],digits = 6),    ",  Sigma_2= ", round(sigma[2],digits = 6)))
        # print(paste0("Threshold= ", round(threshold,digits = 3)))
        
        #************* Return OutPuts *************#
        results<-new("GmmThreshold",
                     x=ent,
                     omega1=omega[1], 
                     omega2=omega[2],
                     mu1=mu[1],
                     mu2=mu[2],
                     sigma1=sigma[1], 
                     sigma2=sigma[2],
                     threshold=threshold)        
    } else {    
        print("Error: No Gaussian fit found")
        results<-NULL
    }
    
    return(results)
}

# Calculates the area (integral) bounded in domain[a,b] under an univariate Gaussian 
# with parameters of mixing proportion (omega), mean (mu), and standard deviation (sigma).   
nGaussianArea <- function (a, b, omega, mu, sigma){
  erf1 <- (a-mu)/(sqrt(2)*sigma)
  erf1 <- 2*pnorm(erf1*sqrt(2)) - 1
  
  erf2 <- (b-mu)/(sqrt(2)*sigma)
  erf2 <- 2*pnorm(erf2*sqrt(2)) - 1
  
  area <- sigma * omega * (-erf1 + erf2) / (2*sigma)
  return(area)
}


# Calculates, analytically, the optimum threshold, where the average of 
# the Sensitivity and Specificity reaches its maximum. 
gmmThreshold <- function (omega, mu, sigma) {
  
  f <- 1.e-6
  TopRange <- mu[2] + sqrt(2.*sigma[2]*sigma[2]*log(1/f, base = exp(1)))
  
  erf1 <- (TopRange - mu[1])/(sqrt(2)*sigma[1])
  erf1 <- 2*pnorm(erf1*sqrt(2)) - 1
  erf2 <- (mu[1])/(sqrt(2)*sigma[1])
  erf2 <- 2*pnorm(erf2*sqrt(2)) - 1
  A <- 1./(erf1 + erf2)
  A <- A/sigma[1]
  
  
  erf1 <- (TopRange - mu[2])/(sqrt(2)*sigma[2])
  erf1 <- 2*pnorm(erf1*sqrt(2)) - 1
  erf2 <- (mu[2])/(sqrt(2)*sigma[2])
  erf2 <- 2*pnorm(erf2*sqrt(2)) - 1
  B <- 1./(erf1 + erf2)
  B <- B/sigma[2]
  
  
  a <- sigma[2]*sigma[2] - sigma[1]*sigma[1]
  b <- 2.*(mu[2]*sigma[1]*sigma[1] - mu[1]*sigma[2]*sigma[2])
  c <- mu[1]*mu[1]*sigma[2]*sigma[2] 
  c <- c - mu[2]*mu[2]*sigma[1]*sigma[1]
  c <- c - 2.*sigma[1]*sigma[1]*sigma[2]*sigma[2]*log(A/B, base = exp(1))
  
  x1 <- -b + sqrt(b*b - 4.*a*c)
  x1 <- x1/(2.*a)
  
  x2 <- -b - sqrt(b*b - 4.*a*c)
  x2 <- x2/(2.*a)
  
  TP <- nGaussianArea(0, x1, omega[1], mu[1], sigma[1])
  FN <- nGaussianArea(x1, TopRange, omega[1], mu[1], sigma[1])
  FP <- nGaussianArea(0, x1, omega[2], mu[2], sigma[2])
  TN <- nGaussianArea(x1, TopRange, omega[2], mu[2], sigma[2])
  SEN <- TP/(TP+FN)
  SPC <- TN/(TN+FP)
  ave1 <- (SEN+SPC)/2.
  
  TP <- nGaussianArea(0, x2, omega[1], mu[1], sigma[1])
  FN <- nGaussianArea(x2, TopRange, omega[1], mu[1], sigma[1])
  FP <- nGaussianArea(0, x2, omega[2], mu[2], sigma[2])
  TN <- nGaussianArea(x2, TopRange, omega[2], mu[2], sigma[2])
  SEN <- TP/(TP+FN)
  SPC <- TN/(TN+FP)
  ave2 <- (SEN+SPC)/2.
  
  if (ave1 > ave2) {
    threshold <- x1
  } else {
    threshold <- x2
  }
  
  return(threshold)
}


#' Plot findThreshold results for the gmm method
#' 
#' \code{plotGmmThreshold} plots the results from \code{"gmm"} method of 
#' \link{findThreshold}, including the Guassian distributions, input nearest neighbor 
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
#' @param    breaks    number of breaks o show on the x-axis. If \code{NULL} the breaks will 
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
#' db <- subset(ExampleDb, SAMPLE == "-1h")
#'
#' # Use nucleotide Hamming distance and normalize by junction length
#' db <- distToNearest(db, model="ham", normalize="len", nproc=1)
#' 
#' # To find the threshold cut, call findThreshold function for "gmm" method.
#' output <- findThreshold(db$DIST_NEAREST, method="gmm")
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
    
    # Invoke Gaussians parameters
    omega <- c(data@omega1, data@omega2)
    mu <-    c(data@mu1, data@mu2) 
    sigma <- c(data@sigma1, data@sigma2) 
    
    # Generate Gaussian curves
    gx <- seq(min(xdf$x), max(xdf$x), by=0.002)
    g1 <- data.frame(x=gx, y=omega[1]*dnorm(gx, mu[1], sigma[1]))
    g2 <- data.frame(x=gx, y=omega[2]*dnorm(gx, mu[2], sigma[2]))
    
    # ggplot workaround
    if (is.null(xmin)) { xmin <- NA }
    if (is.null(xmax)) { xmax <- NA }
    
    # Plot distToNearest distribution plus Gaussian fits
    p <- ggplot(xdf, aes_string(x="x")) +
        getBaseTheme() + 
        xlab("Distance") + 
        ylab("Density") +
        geom_histogram(aes_string(y="..density.."), binwidth=binwidth, 
                       fill="gray40", color="white") +
        geom_line(data=g1, aes_string(x="x", y="y"), color="darkslateblue", size=size) +
        geom_line(data=g2, aes_string(x="x", y="y"), color="darkslateblue", size=size) +
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
    if (!is.na(xmin) | !is.na(xmax) & is.null(breaks)) {
        p <- p + xlim(xmin, xmax)
    } else if (!is.null(breaks)) {
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
#' @param    breaks    number of breaks o show on the x-axis. If \code{NULL} the breaks will 
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
#' db <- subset(ExampleDb, SAMPLE == "-1h")
#'
#' # Use nucleotide Hamming distance and normalize by junction length
#' db <- distToNearest(db, model="ham", normalize="len", nproc=1)
#' 
#' # To find the threshold cut, call findThreshold function for "gmm" method.
#' output <- findThreshold(db$DIST_NEAREST, method="density")
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
    
    # Set binwidth
    if (is.null(binwidth)) { binwidth <- data@bandwidth }
    
    # ggplot workaround
    if (is.null(xmin)) { xmin <- NA }
    if (is.null(xmax)) { xmax <- NA }
    
    # Plot distToNearest distribution plus Gaussian fits
    p <- ggplot(xdf, aes_string(x="x")) +
        getBaseTheme() +
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
    if (!is.na(xmin) | !is.na(xmax) & is.null(breaks)) {
        p <- p + xlim(xmin, xmax)
    } else if (!is.null(breaks)) {
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