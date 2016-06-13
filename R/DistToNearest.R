# Generates distance to nearest neighbor

#' @include shazam.R
NULL


# Returns a 5-mer sliding window of given sequence
#
# @param   sequence   sequence string
# @return  An array of 5-mer sliding windows
#
# @examples
# window5Mers("ACGTNACGTNACGTN")
window5Mers <- function(sequence){
    n <- stri_length(sequence)
    w <- substr(rep(sequence, n - 4), 1:(n - 4), 5:n)
    
    return(w)
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1            first nucleotide sequence, broken into 5mers.
# @param    seq2            second nucleotide sequence, broken into 5mers.
# @param    targetingModel  targeting model.
# @param    normalize       method of normalization. Default is "none".
#                           "length" = normalize distance by length of junction.
# @param    symmetry        if model is hs5f, distance between seq1 and seq2 is either 
#                           the average (avg) of seq1->seq2 and seq2->seq1 or the 
#                           minimum (min).
# @return   distance between two sequences.
#
# @examples
# seq1 <- c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", 
#          "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
# seq2 <- c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", 
#          "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
#
# shazam:::dist5Mers(seq1, seq2, HS5FModel)
dist5Mers <- function(seq1, seq2, targetingModel, 
                      normalize=c("none", "length", "mutations"),
                      symmetry=c("avg", "min")) {
    # Evaluate choices
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    
    # Get distance from targeting model
    targeting_dist <- calcTargetingDistance(targetingModel)
    
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
    }, error = function(e) {
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


# Given an array of nucleotide sequences, find the pairwise distances
# 
# @param   sequences   character vector of nucleotide sequences.
# @param   model       targeting model.
# @param   normalize   method of normalization. Default is "none".
#                      "length" = normalize distance by length of junction.
# @param   symmetry    if model is hs5f, distance between seq1 and seq2 is either the
#                      average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
# @return  A matrix of pairwise distances between junction sequences.
# 
# @details
# needs method details
# 
# @seealso needs links
# 
# @examples
# # working example
pairwise5MerDist <- function(sequences, targetingModel, 
                             normalize=c("none", "length", "mutations"),
                             symmetry=c("avg", "min")) {
    # Initial checks
    normalize <- match.arg(normalize)
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
                                                  targetingModel,
                                                  normalize=normalize,
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


# Given an array of sequences, find the distance to the closest sequence
#
# @param    sequences      character vector of sequences.
# @param    model          DNA (ham) or amino acid (aa) hamming distance model or
#                          mouse (m1n) or human (hs1f) single nucleotide distance model;
#                          or human 5-mer (hs5f) nucleotide model.
# @param    normalize      method of normalization. Default is "none".
#                          "length" = normalize distance by length of junction.
#                          "mutations" = normalize distance by number of mutations in 
#                          junction.
# @param   targetingModel  targeting model for model="hs5f".
# @param   symmetry        if model is hs5f, distance between seq1 and seq2 is either the
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
# shazam:::nearestDist(sequences, model="ham", normalize="length")
# shazam:::nearestDist(sequences, model="aa", normalize="length")
nearestDist<- function(sequences, model=c("ham", "aa", "m1n", "hs1f", "hs5f"),
                       normalize=c("none", "length", "mutations"),
                       targetingModel=NULL, symmetry=c("avg", "min"),
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
        } else if (model == "m1n") {
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=M1NDistance)
        } else if (model == "hs1f") {
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=HS1FDistance)
        } else if (model == "aa") {
            # Translate sequences
            seq_uniq <- setNames(alakazam::translateDNA(seq_uniq), seq_uniq)
            dist_mat <- pairwiseDist(seq_uniq, dist_mat=getAAMatrix())
        } else if (model == "hs5f") {
            if (is.null(targetingModel)) {
                stop("Must provide targetingModel for model='hs5f'")
            }
            dist_mat <- pairwise5MerDist(seq_uniq, targetingModel, 
                                             normalize, symmetry)
        }

        ## DEBUG
        # cat("\n-> seq_uniq:\n")
        # print(seq_uniq)
        # cat("\n-> dist_mat (raw):\n")
        # print(dist_mat)

        # Normalize distances
        if (normalize == "length") { 
            dist_mat <- dist_mat / seq_length
        } else if (normalize == "mutations") {
            #dist <- dist/sum(strsplit(seq1,"")[[1]] != strsplit(seq2,"")[[1]])
            stop("Sorry! nomalize=mutations is not available.")
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
#'                           \code{"none"} then no normalization if performed.
#' @param    symmetry        if model is hs5f, distance between seq1 and seq2 is either the
#'                           average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
#' @param    first           if \code{TRUE} only the first call of the gene assignments 
#'                           is used. if \code{FALSE} the union of ambiguous gene 
#'                           assignments is used to group all sequences with any 
#'                           overlapping gene calls.
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
    targeting_model <- if (model == "hs5f") { HS5FModel } else { NULL }

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
                                 "HS1FDistance",
                                 "calcTargetingDistance",
                                 "findUniqSeq",
                                 "pairwise5MerDist",
                                 "targeting_model")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    }
    

    list_db <- foreach(idx=iterators::icount(lenGroups), .errorhandling='pass') %dopar% {
        db_group <- db[groups[[idx]], ]
        crossGroups <- NULL
        if (!is.null(cross)) {
            crossGroups <- db_group %>% dplyr::group_indices_(.dots=cross)
        }
        arrSeqs <-  as.vector(unlist(db[groups[[idx]], sequenceColumn]))
        db_group$TMP_DIST_NEAREST <- nearestDist(arrSeqs,
                                                 model=model,
                                                 normalize=normalize,
                                                 targetingModel=targeting_model,
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
