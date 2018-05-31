# Mutation profiling

#' @include Shazam.R
NULL

#### Clonal consensus building functions ####

#' Constructs effective clonal sequences for all clones
#'
#' \code{collapseClones} creates effective input and germline sequences for each clonal 
#' group and appends columns containing the consensus sequences to the input 
#' \code{data.frame}.
#'
#' @param   db                  \code{data.frame} containing sequence data. Required.
#' @param   cloneColumn         \code{character} name of the column containing clonal 
#'                              identifiers. Required.
#' @param   sequenceColumn      \code{character} name of the column containing input 
#'                              sequences. Required. The length of each input sequence should 
#'                              match that of its corresponding germline sequence.
#' @param   germlineColumn      \code{character} name of the column containing germline 
#'                              sequences. Required. The length of each germline sequence 
#'                              should match that of its corresponding input sequence.
#' @param   muFreqColumn        \code{character} name of the column containing mutation
#'                              frequency. Optional. Applicable to the \code{"mostMutated"}
#'                              and \code{"leastMutated"} methods. If not supplied, mutation
#'                              frequency is computed by calling \code{observedMutations}.
#'                              Default is \code{NULL}. See Cautions for note on usage.
#' @param   regionDefinition    \link{RegionDefinition} object defining the regions
#'                              and boundaries of the Ig sequences. Optional. Default is 
#'                              \code{NULL}.
#' @param   method              method for calculating input consensus sequence. Required. 
#'                              One of \code{"thresholdedFreq"}, \code{"mostCommon"}, 
#'                              \code{"catchAll"}, \code{"mostMutated"}, or 
#'                              \code{"leastMutated"}. See "Methods" for details.
#' @param   minimumFrequency    frequency threshold for calculating input consensus sequence.
#'                              Applicable to and required for the \code{"thresholdedFreq"} 
#'                              method. A canonical choice is 0.6. Default is \code{NULL}. 
#' @param   includeAmbiguous    whether to use ambiguous characters to represent positions 
#'                              at which there are multiple characters with frequencies that 
#'                              are at least \code{minimumFrequency} or that are maximal 
#'                              (i.e. ties). Applicable to and required for the 
#'                              \code{"thresholdedFreq"} and \code{"mostCommon"} methods. 
#'                              Default is \code{FALSE}. See "Choosing ambiguous characters" 
#'                              for rules on choosing ambiguous characters.
#' @param   breakTiesStochastic In case of ties, whether to randomly pick a sequence from 
#'                              sequences that fulfill the criteria as consensus. Applicable 
#'                              to and required for all methods except for \code{"catchAll"}. 
#'                              Default is \code{FALSE}. See "Methods" for details. 
#' @param   breakTiesByColumns  A list of the form 
#'                              \code{list(c(col_1, col_2, ...), c(fun_1, fun_2, ...))}, 
#'                              where \code{col_i} is a \code{character} name of a column 
#'                              in \code{db}, and \code{fun_i} is a function to be applied 
#'                              on that column. Currently, only \code{max} and \code{min} 
#'                              are supported. Note that the two \code{c()}'s in \code{list()} 
#'                              are essential (i.e. if there is only 1 column, the list should 
#'                              be of the form \code{list(c(col_1), c(func_1))}. Applicable 
#'                              to and optional for the \code{"mostMutated"} and 
#'                              \code{"leastMutated"} methods. If supplied, \code{fun_i}'s 
#'                              are applied on \code{col_i}'s to help break ties. Default 
#'                              is \code{NULL}. See "Methods" for details. 
#' @param   expandedDb          \code{logical} indicating whether or not to return the 
#'                              expanded \code{db}, containing all the sequences (as opposed
#'                              to returning just one sequence per clone).
#' @param   nproc               Number of cores to distribute the operation over. If the 
#'                              \code{cluster} has already been set earlier, then pass the 
#'                              \code{cluster}. This will ensure that it is not reset.
#'                              
#' 
#' @return   A modified \code{db} with the following additional columns: 
#'           \itemize{
#'             \item \code{CLONAL_SEQUENCE}:  effective sequence for the clone.
#'             \item \code{CLONAL_GERMLINE}:  germline sequence for the clone.
#'             \item \code{CLONAL_SEQUENCE_MUFREQ}:  mutation frequency of 
#'                   \code{CLONAL_SEQUENCE}; only added for the \code{"mostMutated"}
#'                   and \code{"leastMutated"} methods.
#'           }
#'                      
#'           \code{CLONAL_SEQUENCE} is generated with the method of choice indicated 
#'           by \code{method}, and \code{CLONAL_GERMLINE} is generated with the 
#'           \code{"mostCommon"} method, along with, where applicable, user-defined 
#'           parameters such as \code{minimumFrequency}, \code{includeAmbiguous}, 
#'           \code{breakTiesStochastic}, and \code{breakTiesByColumns}.
#'           
#'
#' @section Consensus lengths: For each clone, \code{CLONAL_SEQUENCE} and 
#'          \code{CLONAL_GERMLINE} have the same length. 
#'          
#'          \itemize{
#'                \item For the \code{"thresholdedFreq"}, \code{"mostCommon"}, and 
#'                \code{"catchAll"} methods:
#'          
#'                The length of the consensus sequences is determined by the longest possible
#'                consensus sequence (baesd on \code{inputSeq} and \code{germlineSeq}) and 
#'                \code{regionDefinition@seqLength} (if supplied), whichever is shorter.
#'
#'                Given a set of sequences of potentially varying lengths, the longest possible 
#'                length of their consensus sequence is taken to be the longest length along 
#'                which there is information contained at every nucleotide position across 
#'                majority of the sequences. Majority is defined to be greater than 
#'                \code{floor(n/2)}, where \code{n} is the number of sequences. If the longest 
#'                possible consensus length is 0, there will be a warning and an empty string 
#'                (\code{""}) will be returned. 
#'          
#'                If a length limit is defined by supplying a \code{regionDefinition} via 
#'                \code{regionDefinition@seqLength}, the consensus length will be further 
#'                restricted to the shorter of the longest possible length and 
#'                \code{regionDefinition@seqLength}.
#'          
#'                \item For the \code{"mostMutated"} and \code{"leastMutated"} methods:
#'                
#'                The length of the consensus sequences depends on that of the most/least 
#'                mutated input sequence, and, if supplied, the length limit defined by 
#'                \code{regionDefinition@seqLength}, whichever is shorter. If the germline 
#'                consensus computed using the \code{"mostCommon"} method is longer than 
#'                the most/least mutated input sequence, the germline consensus is trimmed 
#'                to be of the same length as the input consensus.
#'               
#'           }
#'
#' @section Methods: The descriptions below use "sequences" as a generalization of input 
#'          sequences and germline sequences. 
#'          
#'          \itemize{
#'          
#'              \item \code{method="thresholdedFreq"}
#'              
#'                    A threshold must be supplied to the argument \code{minimumFrequency}. At 
#'                    each position along the length of the consensus sequence, the frequency 
#'                    of each nucleotide/character across sequences is tabulated. The 
#'                    nucleotide/character whose frequency is at least (i.e. \code{>=}) 
#'                    \code{minimumFrequency} becomes the consensus; if there is none, the
#'                    consensus nucleotide will be \code{"N"}.
#'                    
#'                    When there are ties (frequencies of multiple nucleotides/characters 
#'                    are at least \code{minimumFrequency}), this method can be deterministic 
#'                    or stochastic, depending on additional parameters.
#'                    
#'                    \itemize{
#'                         \item With \code{includeAmbiguous=TRUE}, ties are resolved 
#'                               deterministically by representing ties using ambiguous 
#'                               characters. See "Choosing ambiguous characters" for how 
#'                               ambiguous characters are chosen.
#'                         \item With \code{breakTiesStochastic=TRUE}, ties are resolved 
#'                               stochastically by randomly picking a character amongst the 
#'                               ties.
#'                         \item When both \code{TRUE}, \code{includeAmbiguous} takes 
#'                               precedence over \code{breakTiesStochastic}.
#'                         \item When both \code{FALSE}, the first character from the ties is 
#'                               taken to be the consensus following the order of \code{"A"}, 
#'                               \code{"T"}, \code{"G"}, \code{"C"}, \code{"N"}, \code{"."}, 
#'                               and \code{"-"}.
#'                    }
#'                    
#'                    Below are some examples looking at a single position based on 5 
#'                    sequences with \code{minimumFrequency=0.6}, 
#'                    \code{includeAmbiguous=FALSE}, and \code{breakTiesStochastic=FALSE}:
#'                    
#'                    \itemize{
#'                         \item If the sequences have \code{"A"}, \code{"A"}, \code{"A"}, 
#'                               \code{"T"}, \code{"C"}, the consensus will be \code{"A"}, 
#'                               because \code{"A"} has frequency 0.6, which is at least 
#'                               \code{minimumFrequency}.
#'                         \item If the sequences have \code{"A"}, \code{"A"}, \code{"T"}, 
#'                               \code{"T"}, \code{"C"}, the consensus will be \code{"N"}, 
#'                               because none of \code{"A"}, \code{"T"}, or \code{"C"} has 
#'                               frequency that is at least \code{minimumFrequency}.
#'                    }
#'          
#'               \item \code{method="mostCommon"}
#'               
#'                     The most frequent nucleotide/character across sequences at each 
#'                     position along the length of the consensus sequence makes up the consensus.
#'                    
#'                     When there are ties (multiple nucleotides/characters with equally 
#'                     maximal frequencies), this method can be deterministic or stochastic, 
#'                     depending on additional parameters. The same rules for breaking ties 
#'                     for \code{method="thresholdedFreq"} apply.
#'                    
#'                     Below are some examples looking at a single position based on 5 
#'                     sequences with \code{includeAmbiguous=FALSE}, and 
#'                     \code{breakTiesStochastic=FALSE}:
#'                     
#'                     \itemize{
#'                          \item If the sequences have \code{"A"}, \code{"A"}, \code{"T"}, 
#'                                \code{"A"}, \code{"C"}, the consensus will be \code{"A"}.
#'                          \item If the sequences have \code{"T"}, \code{"T"}, \code{"C"}, 
#'                                \code{"C"}, \code{"G"}, the consensus will be \code{"T"}, 
#'                                because \code{"T"} is before \code{"C"} in the order of 
#'                                \code{"A"}, \code{"T"}, \code{"G"}, \code{"C"}, \code{"N"}, 
#'                                \code{"."}, and \code{"-"}. 
#'                     }       
#'                     
#'                     
#'               \item \code{method="catchAll"}
#'               
#'                     This method returns a consensus sequence capturing most of the 
#'                     information contained in the sequences. Ambiguous characters are 
#'                     used where applicable. See "Choosing ambiguous characters" for how 
#'                     ambiguous characters are chosen. This method is deterministic and 
#'                     does not involve breaking ties.
#'                     
#'                     Below are some examples for \code{method="catchAll"} looking at a 
#'                     single position based on 5 sequences:
#'                     
#'                     \itemize{
#'                          \item If the sequences have \code{"N"}, \code{"N"}, \code{"N"}, 
#'                                \code{"N"}, \code{"N"}, the consensus will be \code{"N"}.
#'                          \item If the sequences have \code{"N"}, \code{"A"}, \code{"A"}, 
#'                                \code{"A"}, \code{"A"}, the consensus will be \code{"A"}.
#'                          \item If the sequences have \code{"N"}, \code{"A"}, \code{"G"}, 
#'                                \code{"A"}, \code{"A"}, the consensus will be \code{"R"}.
#'                          \item If the sequences have \code{"-"}, \code{"-"}, \code{"."}, 
#'                                \code{"."}, \code{"."}, the consensus will be \code{"-"}.
#'                          \item If the sequences have \code{"-"}, \code{"-"}, \code{"-"}, 
#'                                \code{"-"}, \code{"-"}, the consensus will be \code{"-"}.
#'                          \item If the sequences have \code{"."}, \code{"."}, \code{"."}, 
#'                                \code{"."}, \code{"."}, the consensus will be \code{"."}.
#'                    }
#'                    
#'              \item \code{method="mostMutated"} and \code{method="leastMutated"}
#'              
#'                    These methods return the most/least mutated sequence as the consensus 
#'                    sequence. 
#'                    
#'                    When there are ties (multple sequences have the maximal/minimal mutation
#'                    frequency), this method can be deterministic or stochastic, depending on 
#'                    additional parameters.
#'                    
#'                    \itemize{
#'                         \item With \code{breakTiesStochastic=TRUE}, ties are resolved 
#'                               stochastically by randomly picking a sequence out of 
#'                               sequences with the maximal/minimal mutation frequency.
#'                         \item When \code{breakTiesByColumns} is supplied, ties are resolved
#'                               deterministically. Column by column, a function is applied on 
#'                               the column and sequences with column value matching the 
#'                               functional value are retained, until ties are resolved or 
#'                               columns run out. In the latter case, the first remaining 
#'                               sequence is taken as the consensus.
#'                         \item When \code{breakTiesStochastic=TRUE} and 
#'                               \code{breakTiesByColumns} is also supplied, 
#'                               \code{breakTiesStochastic} takes precedence over 
#'                               \code{breakTiesByColumns}.
#'                         \item When \code{breakTiesStochastic=FALSE} and 
#'                               \code{breakTiesByColumns} is not supplied (i.e. \code{NULL}), 
#'                               the sequence that appears first amongst the ties is taken 
#'                               as the consensus.
#'                    }
#'          
#'          }
#'          
#' 
#' @section Choosing ambiguous characters: 
#'          
#'          Ambiguous characters may be present in the returned consensuses when using the
#'          \code{"catchAll"} method and when using the \code{"thresholdedFreq"} or 
#'          \code{"mostCommon"} methods with \code{includeAmbiguous=TRUE}. 
#'          
#'          The rules on choosing ambiguous characters are as follows:
#'          
#'          \itemize{
#'               \item If a position contains only \code{"N"} across sequences, the consensus 
#'                     at that position is \code{"N"}.
#'               \item If a position contains one or more of \code{"A"}, \code{"T"}, 
#'                     \code{"G"}, or \code{"C"}, the consensus will be an IUPAC character 
#'                     representing all of the characters present, regardless of whether 
#'                     \code{"N"}, \code{"-"}, or \code{"."} is present.
#'               \item If a position contains only \code{"-"} and \code{"."} across sequences, 
#'                     the consensus at thatp osition is taken to be \code{"-"}. 
#'               \item If a position contains only one of \code{"-"} or \code{"."} across 
#'                     sequences, the consensus at that position is taken to be the character 
#'                     present. 
#'          }
#' 
#' @section Cautions: 
#' 
#'          \itemize{
#'               \item   Note that this function does not perform multiple sequence alignment. 
#'                       As a prerequisite, it is assumed that the sequences in 
#'                       \code{sequenceColumn} and \code{germlineColumn} have been aligned 
#'                       somehow. In the case of immunoglobulin repertoire analysis, this 
#'                       usually means that the sequences are IMGT-gapped.
#'               \item   When using the \code{"mostMutated"} and \code{"leastMutated"} methods, 
#'                       if you supply both \code{muFreqColumn} and \code{regionDefinition},
#'                       it is your responsibility to ensure that the mutation frequency in
#'                       \code{muFreqColumn} was calculated with sequence lengths restricted 
#'                       to the \strong{same} \code{regionDefinition} you are supplying. 
#'                       Otherwise, the "most/least mutated" sequence you obtain might not 
#'                       be the most/least mutated given the \code{regionDefinition} supplied, 
#'                       because your mutation frequency was based on a 
#'                       \code{regionDefinition} different from the one supplied.
#'               \item   If you intend to run \code{collapseClones} before 
#'                       building a 5-mer targeting model, you \strong{must} choose 
#'                       parameters such that your collapsed clonal consensuses do 
#'                       \strong{not} include ambiguous characters. This is because the 
#'                       targeting model functions do NOT support ambiguous characters 
#'                       in their inputs.
#'               }
#' 
#' @seealso
#' See \link{IMGT_SCHEMES} for a set of predefined \link{RegionDefinition} objects.
#' 
#' @examples
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d" &
#'                         CLONE %in% c("3100", "3141", "3184"))
#' 
#' # thresholdedFreq method, resolving ties deterministically without using ambiguous characters
#' clones <- collapseClones(db, method="thresholdedFreq", minimumFrequency=0.6,
#'                          includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#'
#' # mostCommon method, resolving ties deterministically using ambiguous characters
#' clones <- collapseClones(db, method="mostCommon", 
#'                          includeAmbiguous=TRUE, breakTiesStochastic=FALSE)
#' 
#' # Make a copy of db that has a mutation frequency column
#' db2 <- observedMutations(db, frequency=TRUE, combine=TRUE)
#' # mostMutated method, resolving ties stochastically
#' clones <- collapseClones(db2, method="mostMutated", muFreqColumn="MU_FREQ", 
#'                          breakTiesStochastic=TRUE, breakTiesByColumns=NULL)
#' # mostMutated method, resolving ties deterministically using additional columns
#' clones <- collapseClones(db2, method="mostMutated", muFreqColumn="MU_FREQ", 
#'                          breakTiesStochastic=FALSE, 
#'                          breakTiesByColumns=list(c("DUPCOUNT"), c(max)))
#' 
#' # catchAll method
#' clones <- collapseClones(db, method="catchAll")
#' 
#' # Build clonal consensus for V-region only
#' clones <- collapseClones(db, method="mostCommon", regionDefinition=IMGT_V)
#' 
#' # Return the same number of rows as the input
#' clones <- collapseClones(db, method="mostCommon", expandedDb=TRUE)
#' 
#' @export
collapseClones <- function(db, 
                           cloneColumn="CLONE", 
                           sequenceColumn="SEQUENCE_IMGT",
                           germlineColumn="GERMLINE_IMGT_D_MASK",
                           muFreqColumn=NULL,
                           regionDefinition=NULL,
                           method=c("mostCommon", "thresholdedFreq", "catchAll", 
                                    "mostMutated", "leastMutated"),
                           minimumFrequency=NULL,
                           includeAmbiguous=FALSE,
                           breakTiesStochastic=FALSE,
                           breakTiesByColumns=NULL,
                           expandedDb=FALSE,
                           nproc=1) {
    # Hack for visibility of foreach index variables
    idx <- NULL
    
    ## DEBUG
    # cloneColumn="CLONE"; sequenceColumn="SEQUENCE_IMGT"; germlineColumn="GERMLINE_IMGT_D_MASK"
    # expandedDb=FALSE; regionDefinition=NULL; method="mostCommon"; nproc=1
    
    #### parameter checks
    
    method <- match.arg(method)
    
    # check minimumFrequency for thresholdedFreq method
    if (method=="thresholdedFreq") {
        if (!is.numeric(minimumFrequency)) {
            stop("minimumFrequency must be a numeric value.")
        } else {
            if ( minimumFrequency<0 | minimumFrequency>1 ) {
                stop("minimumFrequency must be between 0 and 1.")
            }
        }
    }
    
    # check includeAmbiguous & breakTiesStochastic for methods other than catchAll
    if (method %in% c("thresholdedFreq", "mostCommon", "mostMutated", "leastMutated")) {
        if (class(includeAmbiguous) != "logical") {
            stop ("includeAmbiguous must be TRUE or FALSE.")
        }
        if (class(breakTiesStochastic) != "logical") {
            stop ("breakTiesStochastic must be TRUE or FALSE.")
        }
    }
    
    # check breakTiesByColumns and muFreqColumn for methods most/leastMutated
    if (method %in% c("mostMutated", "leastMutated")) {
        
        if (!is.null(breakTiesByColumns)) {
            if (class(breakTiesByColumns) != "list") {
                stop ("breakTiesByColumns must be a list.")
            }
            if (length(breakTiesByColumns) != 2) {
                stop ("breakTiesByColumns must be a nested list of length 2.")
            }
            if ( length(breakTiesByColumns[[1]]) != length(breakTiesByColumns[[2]]) ) {
                stop ("Nested vectors in breakTiesByColumns must have the same lengths.")
            }
            if ( !all(is.character(breakTiesByColumns[[1]])) ) {
                stop ("The first vector in breakTiesByColumns must contain column names.")
            }
            if ( !all( unlist( lapply(breakTiesByColumns[[2]], is.function) ) ) ) {
                stop ("The second vector in breakTiesByColumns must contain functions.")
            }
            if ( !all(breakTiesByColumns[[1]] %in% colnames(db)) ) {
                stop ("All column named included in breakTiesByColumns must be present in db.")
            }
        }
        
        if ( (!is.null(muFreqColumn)) && (!muFreqColumn %in% colnames(db)) ) {
            stop ("If specified, muFreqColumn must be a column present in db.")
        }
    }
    
    # check mutual exclusivitiy
    if (method %in% c("thresholdedFreq", "mostCommon")){
        if (includeAmbiguous & breakTiesStochastic) {
            message("includeAmbiguous and breakTiesStochastic are mutually exclusive. When both TRUE, includeAmbiguous will take precedence.")
        }
        #if ( (!includeAmbiguous) & (!breakTiesStochastic) ) {
        #    message("When both includeAmbiguous and breakTiesStochastic are FALSE, ties are broken in the order of 'A', 'T', 'G', 'C', 'N', '.', and '-'.")
        #}
        if (!is.null(breakTiesByColumns)) {
            message("breakTiesByColumns is ignored when method is thresholdedFreq or mostCommon.")
        }
    }
    
    if (method %in% c("mostMutated", "leastMutated")){
        if (breakTiesStochastic & !is.null(breakTiesByColumns)) {
            message("breakTiesStochastic and breakTiesByColumns are mutually exclusive. When both set, breakTiesStochastic will take precedence.")
        }
        #if ( (!breakTiesStochastic) & is.null(breakTiesByColumns) ) {
        #    message("When breakTiesStochastic is FALSE and breakTiesByColumns is NULL, ties are broken by taking the sequence that appears earlier in the data.frame.")
        #}
        if (includeAmbiguous) {
            message("includeAmbiguous is ignored when method is mostMutated or leastMutated.")
        }
    }
    
    # Check for valid columns
    check <- checkColumns(db, c(cloneColumn, sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    ### Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster <- nproc 
        nproc <- 0
    }
    
    if (class(expandedDb) != "logical") {
        stop ("expandedDb must be TRUE or FALSE.")
    }
    
    # Convert clone identifiers to strings
    db[[cloneColumn]] <- as.character(db[[cloneColumn]])
    
    # get row indices in db for each unique clone
    uniqueClones <- unique(db[[cloneColumn]])
    # crucial to have simplify=FALSE (otherwise won't return a list if uniqueClones has length 1)
    uniqueClonesIdx <- sapply(uniqueClones, function(i){which(db[[cloneColumn]]==i)}, simplify=FALSE)
    
    # if method is most/leastMutated and muFreqColumn not specified,
    # first calculate mutation frequency ($MU_FREQ)
    # IMPORTANT: do this OUTSIDE foreach loop for calcClonalConsensus
    # otherwise will get an error saying muFreqColumn not found in db
    # (something to do with parallelization/foreach)
    if ( (method %in% c("mostMutated", "leastMutated")) & is.null(muFreqColumn) ) {
        message("Calculating observed mutation frequency...")
        db <- observedMutations(db=db, sequenceColumn=sequenceColumn,
                                germlineColumn=germlineColumn, 
                                regionDefinition=regionDefinition,
                                frequency=TRUE, combine=TRUE, 
                                mutationDefinition=NULL, nproc=nproc)
        muFreqColumn <- "MU_FREQ"
    }
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, cpuCount())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else {
        if (nproc != 0) { 
            #cluster <- makeCluster(nproc, type="SOCK") 
            cluster <- parallel::makeCluster(nproc, type= "PSOCK")
        }
        parallel::clusterExport(cluster, 
                                list('db', 'cloneColumn', 'sequenceColumn', 'germlineColumn', 'muFreqColumn',
                                     'regionDefinition', 'method', 'minimumFrequency','includeAmbiguous',
                                     'breakTiesStochastic', 'breakTiesByColumns', 
                                     'calcClonalConsensus', 'calcClonalConsensusHelper', 'breakTiesHelper',
                                     'chars2Ambiguous', 'nucs2IUPAC', 'IUPAC_DNA_2',  'NUCLEOTIDES_AMBIGUOUS',
                                     'uniqueClonesIdx', 'c2s', 's2c'), 
                                envir=environment() )
        registerDoParallel(cluster)
    }
    
    # Printing status to console
    #cat("Collapsing clonal sequences...\n")
    
    # avoid .combine="cbind"!
    # if there is only 1 unique clone, .combind="cbind" will result in a vector (as opposed to
    # a matrix) being returned, which will subsequently result a failure in
    # cons_db$CLONAL_SEQUENCE <- cons_mat[, 1]
    cons_mat <- foreach(idx=1:length(uniqueClonesIdx),
                        .verbose=FALSE, .errorhandling='stop') %dopar% {
                            
                            cloneIdx <- uniqueClonesIdx[[idx]]
                            cloneDb <- db[cloneIdx, ]
                            
                            # collapse clone
                            calcClonalConsensus(db=cloneDb,
                                                sequenceColumn=sequenceColumn,
                                                germlineColumn=germlineColumn,
                                                muFreqColumn=muFreqColumn,
                                                regionDefinition=regionDefinition,
                                                method=method,
                                                minimumFrequency=minimumFrequency,
                                                includeAmbiguous=includeAmbiguous,
                                                breakTiesStochastic=breakTiesStochastic,
                                                breakTiesByColumns=breakTiesByColumns)
                        }
    # using cbind below will give a matrix with columns being clones
    # use rbind to have rows be clones
    # cols: inputCons, germlineCons, inputMuFreq
    cons_mat <- do.call(rbind, cons_mat)
    
    # Stop cluster
    if(nproc > 1) { parallel::stopCluster(cluster) }
    
    # Build return data.frame
    if (expandedDb) { 
        # Fill all rows with the consensus sequence
        clone_index <- match(db[[cloneColumn]], uniqueClones)
        cons_db <- db
        cons_db$CLONAL_SEQUENCE <- unlist(cons_mat[, 1])[clone_index]
        cons_db$CLONAL_GERMLINE <- unlist(cons_mat[, 2])[clone_index]
        
        # assign mutation frequency corresponding to consensus into CLONAL_SEQUENCE_MUFREQ
        if (method %in% c("mostMutated", "leastMutated")) {
            cons_db$CLONAL_SEQUENCE_MUFREQ <- unlist(cons_mat[, 3])[clone_index]
        }
    } else {
        # Return only the first row of each clone
        clone_index <- match(uniqueClones, db[[cloneColumn]])
        cons_db <- db[clone_index, ]
        cons_db$CLONAL_SEQUENCE <- unlist(cons_mat[, 1])
        cons_db$CLONAL_GERMLINE <- unlist(cons_mat[, 2])
        
        # assign mutation frequency corresponding to consensus into CLONAL_SEQUENCE_MUFREQ
        if (method %in% c("mostMutated", "leastMutated")) {
            cons_db$CLONAL_SEQUENCE_MUFREQ <- unlist(cons_mat[, 3])
        }
    }
    
    return(cons_db)
}


# Break ties given additional columns in db and functions to compute on them
#
# @param     idx     vector of indices.
# @param     cols    character vector of colnames. Currently, only columns containing
#                    numeric values are supported/expected.
# @param     funs    list of functions. Currently, only \code{max} and \code{min} are 
#                    supported/expected.
# @param     db      \code{data.frame} containing columns named after \code{cols} with 
#                    corresponding rows for \code{idx}.
#
# @return    a single value from \code{idx}.
# 
# @details   Column by column, \code{breakTiesHelper} calls the corresponding function
#            from \code{funs} on a column in \code{db} and finds the index/indices in 
#            \code{idx} that match(es) the returned value from the function. This stops
#            when only a single matching index is obtained, or columns run out. In the 
#            latter case, the first remaining index is returned.
#
# testing
# expect index 18
# test.idx = c(2,4,18,37,102,76)
# test.db = data.frame(cbind(DUPCOUNT= c(3,5,5,4,5,1),
#                            CONSCOUNT=c(6,6,6,2,3,4), 
#                            ERR=c(0.9, 0.14, 0.12, 0.07, 0.3, 0.5)))
# test.cols = c("DUPCOUNT", "CONSCOUNT", "ERR")
# test.funs = c(max, max, min)
# stopifnot( breakTiesHelper(test.idx, test.cols, test.funs, test.db)==18 )
# # make index 4 and 18 tie for ERR
# # index 4 is expected because it appears before 18
# test.db[3,"ERR"] = 0.14
# stopifnot( breakTiesHelper(test.idx, test.cols, test.funs, test.db)==4 )
#
breakTiesHelper <- function(idx, cols, funs, db) {
    # debug
    # idx=test.idx; cols=test.cols; funs=test.funs; db=test.db
    
    counter <- 1
    while (length(idx)>1 & counter<=length(cols)) {
        cur.col <- cols[counter]
        cur.fun <- funs[[counter]]
        cur.db <- db[[cur.col]]
        
        target <- cur.fun(cur.db)
        tol <- 1e-5 # tolerance
        target.idx <- which( abs(cur.db-target)<=tol ) # wrt idx & db
        
        idx <- idx[target.idx]
        db <- db[target.idx, ]
        counter <- counter+1
    }
    
    if (length(idx)==1) {
        return(idx)
    } else if (length(idx)>1) {
        #print("Failed to resolve ties.") # for testing/debugging
        return(idx[1])
    } else {
        stop("breakTieHelper failed unexpectedly.")
    }
}

# Helper function for calcClonalConsensus
#
# @param   seqs                 a character vector of sequences.
# @param   muFreqColumn         \code{character} name of the column in db containing mutation
#                               frequency. Applicable to and required for the \code{"mostMutated"} 
#                               and \code{"leastmutated"} methods. Default is \code{NULL}. 
# @param   lenLimit             limit on consensus length. 
# @param   mtd                  method to calculate consensus sequence. One of 
#                               \code{"thresholdedFreq"}, \code{"mostCommon"}, \code{"catchAll"},
#                               \code{"mostMutated"}, or \code{"leastMutated"}. See "Methods" under
#                               \link{collapseClones} for details.
# @param   minFreq              frequency threshold for calculating input consensus sequence.
#                               Applicable to and required for the \code{"thresholdedFreq"} method. 
#                               A canonical choice is 0.6. Default is \code{NULL}.
# @param   includeAmbiguous     whether to use ambiguous characters to represent positions at
#                               which there are multiple characters with frequencies that are at least
#                               \code{minimumFrequency} or that are maximal (i.e. ties). Applicable to 
#                               and required for the \code{"thresholdedFreq"} and \code{"mostCommon"} 
#                               methods. Default is \code{FALSE}. See "Choosing ambiguous characters" 
#                               under \link{collapseClones} for rules on choosing ambiguous characters. 
# @param   breakTiesStochastic  In case of ties, whether to randomly pick a sequence from sequences that
#                               fulfill the criteria as consensus. Applicable to and required for all methods
#                               except for \code{"catchAll"}. Default is \code{FALSE}. See "Methods" 
#                               under \link{collapseClones} for details. 
# @param   breakTiesByColumns   A list of the form \code{list(c(col_1, col_2, ...), c(fun_1, fun_2, ...))}, 
#                               where \code{col_i} is a \code{character} name of a column in \code{db},
#                               and \code{fun_i} is a function to be applied on that column. Currently, 
#                               only \code{max} and \code{min} are supported. Note that the two \code{c()}'s
#                               in \code{list()} are essential (i.e. if there is only 1 column, the list 
#                               should be of the form \code{list(c(col_1), c(func_1))}. Applicable to and 
#                               optional for the \code{"mostMutated"} and \code{"leastMutated"} methods. 
#                               If supplied, \code{fun_i}'s are applied on \code{col_i}'s to help break 
#                               ties. Default is \code{NULL}. See "Methods" under \link{collapseClones} 
#                               for details.                                                  
# @param   db                   \code{data.frame} containing sequence data for a single clone. 
#                               Applicable to and required for the \code{"mostMutated"} and 
#                               \code{"leastmutated"} methods. Default is \code{NULL}.
#
# @return  A list containing \code{cons}, which is a character string that is the consensus sequence 
#          for \code{seqs}; and \code{muFreq}, which is the maximal/minimal mutation frequency of
#          the consensus sequence for the \code{"mostMutated"} and \code{"leastMutated"} methods, or
#          \code{NULL} for all other methods.
#
# @details See \link{collapseClones} for detailed documentation on methods and additional parameters.
# 
# @seealso
# \link{calcClonalConsensus} and \link{collapseClones}.
# 
# @examples
# Subset example data
# data(ExampleDb, package="alakazam")
# db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
# 
# Data corresponding to a single clone
# clone <- db[db$CLONE=="3192", ]
# Number of sequences in this clone
# nrow(clone)
# first compute mutation frequency for most/leastMutated methods
# clone = observedMutations(db=clone, frequency=TRUE, combine=TRUE)
# manually create tie
# clone = rbind(clone, clone[which.max(clone$MU_FREQ), ])
# 
# Get consensus input sequence
# thresholdedFreq method, resolve ties deterministically without using ambiguous characters
# consInput1 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn=NULL, lenLimit=NULL,
#                                         mtd="thresholdedFreq", minFreq=0.3,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=FALSE,
#                                         breakTiesByColumns=NULL, db=NULL)$cons
# thresholdedFreq method, resolve ties deterministically using ambiguous characters
# consInput2 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn=NULL, lenLimit=NULL,
#                                         mtd="thresholdedFreq", minFreq=0.3,
#                                         includeAmbiguous=TRUE, 
#                                         breakTiesStochastic=FALSE,
#                                         breakTiesByColumns=NULL, db=NULL)$cons
# thresholdedFreq method, resolve ties stochastically
# consInput3 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn=NULL, lenLimit=NULL,
#                                         mtd="thresholdedFreq", minFreq=0.3,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=TRUE,
#                                         breakTiesByColumns=NULL, db=NULL)$cons
# mostCommon method, resolve ties deterministically without using ambiguous characters
# consInput4 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn=NULL, lenLimit=NULL,
#                                         mtd="mostCommon", minFreq=NULL,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=FALSE,
#                                         breakTiesByColumns=NULL, db=NULL)$cons
# mostCommon method, resolve ties deterministically using ambiguous characters
# consInput5 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn=NULL, lenLimit=NULL,
#                                         mtd="mostCommon", minFreq=NULL,
#                                         includeAmbiguous=TRUE, 
#                                         breakTiesStochastic=FALSE,
#                                         breakTiesByColumns=NULL, db=NULL)$cons
# mostCommon method, resolve ties stochastically
# consInput6 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn=NULL, lenLimit=NULL,
#                                         mtd="mostCommon", minFreq=NULL,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=TRUE,
#                                         breakTiesByColumns=NULL, db=NULL)$cons
# catchAll method
# consInput7 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn=NULL, lenLimit=NULL,
#                                         mtd="catchAll", minFreq=NULL,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=FALSE,
#                                         breakTiesByColumns=NULL, db=NULL)$cons
# mostMutated method, resolve ties stochastically
# consInput8 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn="MU_FREQ", lenLimit=NULL,
#                                         mtd="mostMutated", minFreq=NULL,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=TRUE,
#                                         breakTiesByColumns=NULL, db=clone)$cons
# mostMutated method, resolve ties deterministically using additional columns
# consInput9 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn="MU_FREQ", lenLimit=NULL,
#                                         mtd="mostMutated", minFreq=NULL,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=FALSE,
#                                         breakTiesByColumns=list(c("JUNCTION_LENGTH","DUPCOUNT"), c(max, max)), 
#                                         db=clone)$cons
# consInput10 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                         muFreqColumn="MU_FREQ", lenLimit=NULL,
#                                         mtd="mostMutated", minFreq=NULL,
#                                         includeAmbiguous=FALSE, 
#                                         breakTiesStochastic=FALSE,
#                                         breakTiesByColumns=list(c("DUPCOUNT"), c(max)), 
#                                         db=clone)$cons
# mostMutated method, resolve ties deterministically withou using additional columns
# consInput11 <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT,
#                                          muFreqColumn="MU_FREQ", lenLimit=NULL,
#                                          mtd="mostMutated", minFreq=NULL,
#                                          includeAmbiguous=FALSE, 
#                                          breakTiesStochastic=FALSE,
#                                          breakTiesByColumns=NULL, db=clone)$cons
# 
calcClonalConsensusHelper = function(seqs, muFreqColumn=NULL, lenLimit=NULL,
                                     mtd=c("mostCommon", "thresholdedFreq", "catchAll", "mostMutated", "leastMutated"),
                                     minFreq=NULL,
                                     includeAmbiguous=FALSE,
                                     breakTiesStochastic=FALSE,
                                     breakTiesByColumns=NULL, db=NULL) {
    mtd <- match.arg(mtd)
    
    # check muFreqColumn and get muFreq for most/leastMutated
    if (mtd %in% c("mostMutated", "leastMutated")) {
        if ( is.null(muFreqColumn) ) {
            stop ("muFreqColumn must be specified when method is most/leastMutated.")
        }
        if ( is.null(db) ) {
            stop ("db containing muFreqColumn must be supplied when method is most/leastMutated.")
        }
        if (!muFreqColumn %in% colnames(db)) {
            print(c("Helper", muFreqColumn))
            print(c("Helper", colnames(db)))
            stop ("muFreqColumn must be a column present in db.")
        }
        
        # get muFreq
        muFreq <- db[[muFreqColumn]]
    }
    
    
    numSeqs <- length(seqs)
    
    ##### if only one sequence in clone, return it
    if (numSeqs==1) {
        # restrict length if there is a lenLimit
        if (!is.null(lenLimit)) {
            consensus <- substr(seqs, 1, min(lenLimit, stri_length(seqs)))
        } else {
            # otherwise, return as is
            consensus <- seqs
        }
        
        # return with mutation frequency (if applicable)
        if (mtd %in% c("mostMutated", "leastMutated")) {
            return(list(cons=consensus, muFreq=db[[muFreqColumn]]))
        } else {
            return(list(cons=consensus, muFreq=NULL))
        }
    }
    
    ##### if all sequences are the same, return now
    if (length(unique(seqs))==1) {
        # restrict length if there is a lenLimit
        if (!is.null(lenLimit)) {
            consensus <- substr(seqs[1], 1, min(lenLimit, stri_length(seqs)))
        } else {
            # otherwise, return as is
            consensus <- seqs[1]
        }
        
        # return with mutation frequency (if applicable)
        if (mtd %in% c("mostMutated", "leastMutated")) {
            return(list(cons=consensus, muFreq=db[[muFreqColumn]][1]))
        } else {
            return(list(cons=consensus, muFreq=NULL))
        }
    }
    
    ##### length of longest sequence in seqs
    lenSeqs <- stri_length(seqs)
    lenMax <- max(lenSeqs, na.rm=T)
    
    ##### methods = thresholdedFreq, mostCommon, catchAll
    if (mtd %in% c("thresholdedFreq", "mostCommon", "catchAll")) {
        ##### convert seqs to a matrix
        # if there's no more nucleotide when a seq ends, fill position with NA
        seqsMtx <- matrix(NA, nrow=numSeqs, ncol=lenMax)
        for (i in 1:numSeqs) {
            seqsMtx[i, 1:lenSeqs[i]] <- s2c(seqs[i])
        }
        
        ##### tabulation matrix
        # col: nucleotide position
        # row: A,T,G,C,N,.,-,na (to distinguish from NA)
        tabMtxRownames <- c("A","T","G","C","N",".","-","na")
        tabMtx <- matrix(0, ncol=lenMax, nrow=8, 
                        dimnames=list(tabMtxRownames, NULL))
        ## across seqs, at each nuc position, how many A, T, G, C, N, ., -? 
        # this does not capture NA
        for (j in 1:ncol(seqsMtx)) {
            tab <- table(seqsMtx[, j])
            tabMtx[match(names(tab), tabMtxRownames), j] <- tab
        }
        ## across seqs, at each nuc position, how many NAs?
        numNAs <- colSums(is.na(seqsMtx))
        tabMtx["na", ] <- numNAs
        # sanity check: counts at each nuc pos (colSum) should sum up to number of seqs
        stopifnot( sum( colSums(tabMtx)==numSeqs )  == ncol(tabMtx)  )
        
        ##### only keep positions at which majority of seqs contain information
        ### if there are odd number of n seqs, keep position if it has > floor(n/2) non-NAs
        # e.g. 5 input seqs, >2 non-NA; 2=floor(5/2)
        ### if there are even number of n seqs,  keep position if it has > n/2 non-NAs
        # e.g. 6 input seqs, >3 non-NA; 3=6/2=floor(6/2)
        numNonNAs <- numSeqs - numNAs
        nonNA.keep <- numNonNAs > floor(numSeqs/2)
        # length of longest possible consensus seq
        lenConsensus <- sum(nonNA.keep)
        if (lenConsensus==0) {
            warning("Consensus cannot be produced. Empty string returned.")
            return("")
        }
        ##### if there is a lenLimit, restrict consensus length to 
        # the shorter of longest possible length and lenLimit
        if (!is.null(lenLimit)) {
            lenConsensus <- min(lenConsensus, lenLimit)
        }
        # drop=FALSE so that it works even with lenConsensus of 1
        tabMtx <- tabMtx[, 1:lenConsensus, drop=FALSE]
        
        ### convert absolute count to fraction
        tabMtx <- tabMtx/numSeqs
        # remove "na" row
        # drop=FALSE so that it works even with lenConsensus of 1
        tabMtx <- tabMtx[-which(rownames(tabMtx)=="na"), , drop=FALSE]
        
        if (mtd=="thresholdedFreq") {
            #print(mtd) # for testing
            # use as.matrix so that apply won't break with ncol(tabMtx)=1
            consensus <- apply(as.matrix(tabMtx), 2, function(x){
                idx <- which(x >= minFreq)
                # if no character >= the threshold, assign an N
                if (length(idx)==0) {
                    return("N")
                    # if there is no tie
                } else if (length(idx)==1){
                    return(names(x)[idx])
                    # if there are ties (multiple characters >= the threhold)
                } else if (length(idx)>1) {
                    # ambiguous character allowed
                    if (includeAmbiguous) {
                        return(chars2Ambiguous(tabMtxRownames[idx]))
                        # ambiguous characters not allowed
                    } else {
                        # stochastic
                        if (breakTiesStochastic) {
                            return(names(x)[sample(x=idx, size=1)])
                            # first one is returned
                            # the order is built-in from tabMtxRownames
                        } else {
                            return(names(x)[idx[1]])
                        }
                    }
                }
            })
        } else if (mtd=="mostCommon") { 
            #print(mtd) # for testing
            # use as.matrix so that apply won't break with ncol(tabMtx)=1
            consensus <- apply(as.matrix(tabMtx), 2, function(x){
                max.freq <- max(x)
                tol <- 1e-5 # tolerance
                max.idx <- which( abs(x-max.freq)<=tol )
                
                # if there is no tie
                if (length(max.idx)==1){
                    return(names(x)[max.idx])
                    # if there are ties (multiple characters with maximal frequency)
                } else if (length(max.idx)>1) {
                    # ambiguous character allowed
                    if (includeAmbiguous) {
                        return(chars2Ambiguous(tabMtxRownames[max.idx]))
                        # ambiguous characters not allowed
                    } else {
                        # stochastic
                        if (breakTiesStochastic) {
                            return(names(x)[sample(x=max.idx, size=1)])
                            # first one is returned
                            # the order is built-in from tabMtxRownames
                        } else {
                            return(names(x)[max.idx[1]])
                        }
                    }
                }
            })
        } else if (mtd=="catchAll") {
            #print(mtd) # for testing
            # use as.matrix so that apply won't break with ncol(tabMtx)=1
            consensus <- apply(as.matrix(tabMtx), 2, function(x){
                # all characters that appear at a position across seqs
                nonZeroNucs <- rownames(tabMtx)[x>0]
                # convert characters to (ambiguous) characters
                return(chars2Ambiguous(nonZeroNucs))
            })
        }
        
        # check there is no ambiguous characters if includeAmbiguous if F 
        if ( (mtd=="thresholdedFreq" | mtd=="mostCommon") & !includeAmbiguous ) {
            ambiguous <- NUCLEOTIDES_AMBIGUOUS[!NUCLEOTIDES_AMBIGUOUS %in% 
                                                  c("A","C","G","T","N","-",".")]
            stopifnot( !any(consensus %in% ambiguous) )
        }
        
        # convert from character vector to string
        consensus <- c2s(consensus)
        # sanity check
        stopifnot( stri_length(consensus)==lenConsensus )
    }
    
    ##### methods = mostMutated, leastMutated
    if (mtd %in% c("mostMutated", "leastMutated")) {
        # if there's a lenLimit
        # if a seq is longer than lenLimit, trim it; otherwise, leave it as is
        if (!is.null(lenLimit)) {
            idxLong <- which(lenSeqs > lenLimit)
            seqs[idxLong] <- substr(seqs[idxLong], 1, lenLimit)
        }
        
        ##### get index of seqs that fulfill the criterion
        # muFreq should have been calculated being on sequences with restricted lengths as defined by
        # regionDefinition (which gives rise to lenLimit)
        if (mtd=="mostMutated") {
            #print(mtd) # for testing
            targetMuFreq <- max(muFreq)
        } else if (mtd=="leastMutated") {
            #print(mtd) # for testing
            targetMuFreq <- min(muFreq)
        }
        tol <- 1e-5 # tolerance
        idx <- which( abs(muFreq-targetMuFreq)<=tol )
        
        ##### if there are no ties
        if (length(idx)==1) {
            consensus <- seqs[idx]
            ##### if there are ties
        } else if (length(idx)>1) {
            
            ### stochastic: randomly pick one from idx
            if (breakTiesStochastic) {
                consensus <- seqs[sample(x=idx, size=1)]
                
                ### deterministic: pick one from idx based on breakTiesByColumns    
            } else if (!is.null(breakTiesByColumns)) {
                idx <- breakTiesHelper(idx=idx, cols=breakTiesByColumns[[1]], 
                                      funs=breakTiesByColumns[[2]], db=db[idx, ])
                consensus <- seqs[idx]
                
                ### deterministic: pick first one from idx    
            } else {
                consensus <- seqs[idx[1]]
            }
        }
        
    }
    
    # check length
    if (!is.null(lenLimit)) {
        stopifnot(stri_length(consensus) <= lenLimit)
    }
    
    if (mtd %in% c("mostMutated", "leastMutated")) {
        return(list(cons=consensus, muFreq=targetMuFreq))
    } else {
        return(list(cons=consensus, muFreq=NULL))
    }
}

# Calculate clonal consensus for a single clone
# 
# Given an aligned set of input/observed sequences and an aligned set of germline sequences, 
# generate an input/observed consensus and a germline consensus. 
#
# @param   db                  \code{data.frame} containing sequence data for a single clone. 
#                              Required.
# @param   sequenceColumn      \code{character} name of the column containing input 
#                              sequences. Required. The length of each input sequence should 
#                              match that of its corresponding germline sequence.
# @param   germlineColumn      \code{character} name of the column containing germline 
#                              sequences. Required. The length of each germline sequence should 
#                              match that of its corresponding input sequence.
# @param   muFreqColumn        \code{character} name of the column containing mutation
#                              frequency. Applicable to and required for the \code{"mostMutated"}
#                              and \code{"leastMutated"} methods. Default is \code{NULL}. See 
#                              "Details" for a note of caution.
# @param   regionDefinition    \link{RegionDefinition} object defining the regions and boundaries 
#                              of the Ig sequences. Optional. Default is \code{NULL}.
# @param   method              method for calculating input consensus sequence. Required. One of 
#                              \code{"thresholdedFreq"}, \code{"mostCommon"}, \code{"catchAll"},
#                              \code{"mostMutated"}, or \code{"leastMutated"}. See "Methods" under
#                              \link{collapseClones} for details.
# @param   minimumFrequency    frequency threshold for calculating input consensus sequence.
#                              Applicable to and required for the \code{"thresholdedFreq"} method. 
#                              A canonical choice is 0.6. Default is \code{NULL}. 
# @param   includeAmbiguous    whether to use ambiguous characters to represent positions at
#                              which there are multiple characters with frequencies that are at least
#                              \code{minimumFrequency} or that are maximal (i.e. ties). Applicable to 
#                              and required for the \code{"thresholdedFreq"} and \code{"mostCommon"} 
#                              methods. Default is \code{FALSE}. See "Choosing ambiguous characters" 
#                              under \link{collapseClones} for rules on choosing ambiguous characters. 
# @param   breakTiesStochastic In case of ties, whether to randomly pick a sequence from sequences that
#                              fulfill the criteria as consensus. Applicable to and required for all methods
#                              except for \code{"catchAll"}. Default is \code{FALSE}. See "Methods" 
#                              under \link{collapseClones} for details. 
# @param   breakTiesByColumns  A list of the form \code{list(c(col_1, col_2, ...), c(fun_1, fun_2, ...))}, 
#                              where \code{col_i} is a \code{character} name of a column in \code{db},
#                              and \code{fun_i} is a function to be applied on that column. Currently, 
#                              only \code{max} and \code{min} are supported. Applicable to and optional for
#                              the \code{"mostMutated"} and \code{"leastMutated"} methods. If supplied, 
#                              \code{fun_i}'s are applied on \code{col_i}'s to help break ties. Default is 
#                              \code{NULL}. See "Methods" under \link{collapseClones} for details.                                                             
#                              
# @return  A named list of length 3. "inputCons" and "germlineCons" are the consensus sequences. 
#          The input and germline consensus sequences have the same length. "inputMuFreq" is the 
#          maximal/minimal mutation frequency for the input consensus for the \code{"mostMutated"} 
#          and \code{"leastMutated"} methods, and \code{NULL} for all other methods.
# 
# @details See \link{collapseClones} for detailed documention on methods and additional parameters.
# 
#          Caution: when using the \code{"mostMutated"} and \code{"leastMutated"} methods, if you 
#          supply a \code{regionDefinition}, it is your responsibility to ensure that the mutation 
#          frequency in\code{muFreqColumn} was calculated with sequence lengths restricted to the 
#          \strong{same} \code{regionDefinition} you are supplying. Otherwise, the 
#          "most/least mutated" sequence you obtain might not be the most/least mutated given the 
#          \code{regionDefinition} supplied, because your mutation frequency was based on a 
#          \code{regionDefinition} different from the one supplied.
#                                       
# @seealso
# See \link{collapseClones} for constructing consensus for all clones.
# 
# @examples
# # Subset example data
# data(ExampleDb, package="alakazam")
# db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
# 
# # Data corresponding to a single clone
# clone <- db[db$CLONE=="3192", ]
# # Number of sequences in this clone
# nrow(clone)
# # compute mutation frequency for most/leastMutated methods
# clone <- observedMutations(db=clone, frequency=TRUE, combine=TRUE)
# # manually create a tie
# clone <- rbind(clone, clone[which.max(clone$MU_FREQ), ])
# 
# # Get consensus input and germline sequences
# # thresholdedFreq method, resolve ties deterministically without using ambiguous characters
# cons1 <- calcClonalConsensus(db=clone,
#                              muFreqColumn=NULL, regionDefinition=NULL,
#                              method="thresholdedFreq", 
#                              minimumFrequency=0.3, includeAmbiguous=FALSE,
#                              breakTiesStochastic=FALSE, breakTiesByColumns=NULL)
# # thresholdedFreq method, resolve ties deterministically using ambiguous characters
# cons2 <- calcClonalConsensus(db=clone,
#                              muFreqColumn=NULL, regionDefinition=NULL,
#                              method="thresholdedFreq", 
#                              minimumFrequency=0.3, includeAmbiguous=TRUE,
#                              breakTiesStochastic=FALSE, breakTiesByColumns=NULL)  
# # thresholdedFreq method, resolve ties stochastically
# cons3 <- calcClonalConsensus(db=clone,
#                              muFreqColumn=NULL, regionDefinition=NULL,
#                              method="thresholdedFreq", 
#                              minimumFrequency=0.3, includeAmbiguous=FALSE,
#                              breakTiesStochastic=TRUE, breakTiesByColumns=NULL)  
# # mostCommon method, resolve ties deterministically without using ambiguous characters
# cons4 <- calcClonalConsensus(db=clone,
#                              muFreqColumn=NULL, regionDefinition=NULL,
#                              method="mostCommon", 
#                              minimumFrequency=NULL, includeAmbiguous=FALSE,
#                              breakTiesStochastic=FALSE, breakTiesByColumns=NULL)
# # mostCommon method, resolve ties deterministically  using ambiguous characters
# cons5 <- calcClonalConsensus(db=clone,
#                              muFreqColumn=NULL, regionDefinition=NULL,
#                              method="mostCommon", 
#                              minimumFrequency=NULL, includeAmbiguous=TRUE,
#                              breakTiesStochastic=FALSE, breakTiesByColumns=NULL)
# # mostCommon method, resolve ties stochastically
# cons6 <- calcClonalConsensus(db=clone,
#                              muFreqColumn=NULL, regionDefinition=NULL,
#                              method="mostCommon", 
#                              minimumFrequency=NULL, includeAmbiguous=FALSE,
#                              breakTiesStochastic=TRUE, breakTiesByColumns=NULL)
# # catchAll method
# cons7 <- calcClonalConsensus(db=clone,
#                              muFreqColumn=NULL, regionDefinition=NULL,
#                              method="catchAll", 
#                              minimumFrequency=NULL, includeAmbiguous=FALSE,
#                              breakTiesStochastic=FALSE, breakTiesByColumns=NULL)
# # mostMutated method, resolve ties stochastically
# cons8 <- calcClonalConsensus(db=clone,
#                              muFreqColumn="MU_FREQ", regionDefinition=NULL,
#                              method="mostMutated", 
#                              minimumFrequency=NULL, includeAmbiguous=FALSE,
#                              breakTiesStochastic=TRUE, breakTiesByColumns=NULL)
# # mostMutated method, resolve ties deterministically using additional columns
# cons9 <- calcClonalConsensus(db=clone,
#                              muFreqColumn="MU_FREQ", regionDefinition=NULL,
#                              method="mostMutated", 
#                              minimumFrequency=NULL, includeAmbiguous=FALSE,
#                              breakTiesStochastic=FALSE, 
#                              breakTiesByColumns=list(c("JUNCTION_LENGTH", "DUPCOUNT"), c(max, max)))
# cons10 <- calcClonalConsensus(db=clone,
#                               muFreqColumn="MU_FREQ", regionDefinition=NULL,
#                               method="mostMutated", 
#                               minimumFrequency=NULL, includeAmbiguous=FALSE,
#                               breakTiesStochastic=FALSE, 
#                               breakTiesByColumns=list(c("DUPCOUNT"), c(max)))
# # mostMutated method, resolve ties deterministically without using additional columns
# cons11 <- calcClonalConsensus(db=clone,
#                               muFreqColumn="MU_FREQ", regionDefinition=NULL,
#                               method="mostMutated", 
#                               minimumFrequency=NULL, includeAmbiguous=FALSE,
#                               breakTiesStochastic=FALSE, breakTiesByColumns=NULL)
# @export
calcClonalConsensus <- function(db, 
                                sequenceColumn="SEQUENCE_IMGT", 
                                germlineColumn="GERMLINE_IMGT_D_MASK", 
                                muFreqColumn=NULL,
                                regionDefinition=NULL, 
                                method=c("mostCommon", "thresholdedFreq", "catchAll", "mostMutated", "leastMutated"), 
                                minimumFrequency=NULL, includeAmbiguous=FALSE,
                                breakTiesStochastic=FALSE, breakTiesByColumns=NULL) {
    method <- match.arg(method)
    
    inputSeq <- db[[sequenceColumn]]
    germlineSeq <- db[[germlineColumn]]
    
    # length of seqs in inputSeq and those in germlineSeq should match
    if ( sum(stri_length(inputSeq)==stri_length(germlineSeq)) != length(inputSeq) ) {
        stop("Sequences in inputSeq and germlineSeq have different lengths.")
    }
    
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # length limit from regionDefinition
    if (!is.null(regionDefinition)) {
        lenRegion <- regionDefinition@seqLength
    } else {
        lenRegion <- NULL
    }
    
    ##### get consensus germline sequence (most common)
    # NULL for minFreq and muFreqColumn b/c mostCommon definitely doesn't need them
    germCons <- calcClonalConsensusHelper(seqs=germlineSeq, minFreq=NULL, lenLimit=lenRegion,
                                         mtd="mostCommon", 
                                         includeAmbiguous=includeAmbiguous,
                                         breakTiesStochastic=breakTiesStochastic,
                                         breakTiesByColumns=NULL,
                                         muFreqColumn=NULL, db=NULL)$cons
    
    ##### get consensus observed sequence
    inputConsMuFreq <- calcClonalConsensusHelper(seqs=inputSeq, minFreq=minimumFrequency, lenLimit=lenRegion,
                                                mtd=method, 
                                                includeAmbiguous=includeAmbiguous,
                                                breakTiesStochastic=breakTiesStochastic,
                                                breakTiesByColumns=breakTiesByColumns,
                                                muFreqColumn=muFreqColumn, db=db)
    inputCons <- inputConsMuFreq$cons
    inputMuFreq <- inputConsMuFreq$muFreq
    
    if (method %in% c("mostMutated", "leastMutated")) {
        # possible to have inputCons and germCons of varying lengths
        # germCons (mostCommon) length is "longest possible length" for mostCommon
        # inputCons length is min of length of most/least mutated and lenLimit
        # if different, trim the two to same length
        lenInput <- stri_length(inputCons)
        lenGerm <- stri_length(germCons)
        if (lenInput != lenGerm) {
            minLen <- min(lenInput, lenGerm)
            inputCons <- substr(inputCons, 1, minLen)
            germCons <- substr(germCons, 1, minLen)
        }
    }
    
    # sanity check: length of germCons and inputCons should be the same
    # all methods other than most/leastMutated should expect same lengths of inputCons & germCons
    stopifnot( stri_length(germCons)==stri_length(inputCons) )
    
    return(list("inputCons"=inputCons, "germlineCons"=germCons, "inputMuFreq"=inputMuFreq))
}


#### Mutation counting functions ####

#' Calculate observed numbers of mutations
#'
#' \code{observedMutations} calculates the observed number of mutations for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      \code{character} name of the column containing input 
#'                               sequences. IUPAC ambiguous characters for DNA are 
#'                               supported.
#' @param    germlineColumn      \code{character} name of the column containing 
#'                               the germline or reference sequence. IUPAC ambiguous 
#'                               characters for DNA are supported.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences. If NULL, mutations 
#'                               are counted for entire sequence.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    ambiguousMode       whether to consider ambiguous characters as 
#'                               \code{"either or"} or \code{"and"} when determining and 
#'                               counting the type(s) of mutations. Applicable only if
#'                               \code{sequenceColumn} and/or \code{germlineColumn} 
#'                               contain(s) ambiguous characters. One of 
#'                               \code{c("eitherOr", "and")}. Default is \code{"eitherOr"}.                               
#' @param    frequency           \code{logical} indicating whether or not to calculate
#'                               mutation frequencies. Default is \code{FALSE}.
#' @param    combine             \code{logical} indicating whether for each sequence should
#'                               the mutation counts for the different regions (CDR, FWR) and 
#'                               mutation types be combined and return one value of 
#'                               count/frequency per sequence instead of 
#'                               multiple values. Default is \code{FALSE}.                          
#' @param    nproc               number of cores to distribute the operation over. If the 
#'                               cluster has already been set the call function with 
#'                               \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                               \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with observed mutation counts for each 
#'           sequence listed. The columns names are dynamically created based on the
#'           regions in the \code{regionDefinition}. For example, when using the
#'           \link{IMGT_V} definition, which defines positions for CDR and
#'           FWR, the following columns are added:
#'           \itemize{
#'             \item  \code{MU_COUNT_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{MU_COUNT_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{MU_COUNT_FWR_R}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{MU_COUNT_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           If \code{frequency=TRUE}, R and S mutation frequencies are
#'           calculated over the number of non-N positions in the speficied regions.
#'           \itemize{
#'             \item  \code{MU_FREQ_CDR_R}:  frequency of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{MU_FREQ_CDR_S}:  frequency of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{MU_FREQ_FWR_R}:  frequency of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{MU_FREQ_FWR_S}:  frequency of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           } 
#'           If \code{frequency=TRUE} and \code{combine=TRUE}, the mutations and non-N positions
#'           are aggregated and a single \code{MU_FREQ} value is returned
#'           \itemize{
#'             \item  \code{MU_FREQ}:  frequency of replacement and silent mutations in the 
#'                                      specified region
#'           }     
#'                                  
#' @details
#' Mutation counts are determined by comparing the input sequences (in the column specified 
#' by \code{sequenceColumn}) to the germline sequence (in the column specified by 
#' \code{germlineColumn}). See \link{calcObservedMutations} for more technical details, 
#' \strong{including criteria for which sequence differences are included in the mutation 
#' counts and which are not}.
#' 
#' The mutations are binned as either replacement (R) or silent (S) across the different 
#' regions of the sequences as defined by \code{regionDefinition}. Typically, this would 
#' be the framework (FWR) and complementarity determining (CDR) regions of IMGT-gapped 
#' nucleotide sequences. Mutation counts are appended to the input \code{db} as 
#' additional columns.
#' 
#' 
#' @seealso  
#' \link{calcObservedMutations} is called by this function to get the number of mutations 
#' in each sequence grouped by the \link{RegionDefinition}. 
#' See \link{IMGT_SCHEMES} for a set of predefined \link{RegionDefinition} objects.
#' See \link{expectedMutations} for calculating expected mutation frequencies.
#'           
#' 
#' @examples
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, ISOTYPE == "IgG" & SAMPLE == "+7d")
#'
#' # Calculate mutation frequency over the entire sequence
#' db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
#'                             germlineColumn="GERMLINE_IMGT_D_MASK",
#'                             frequency=TRUE,
#'                             nproc=1)
#'
#' # Count of V-region mutations split by FWR and CDR
#' # With mutations only considered replacement if charge changes
#' db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
#'                             germlineColumn="GERMLINE_IMGT_D_MASK",
#'                             regionDefinition=IMGT_V,
#'                             mutationDefinition=CHARGE_MUTATIONS,
#'                             nproc=1)
#'                      
#' @export
observedMutations <- function(db, 
                              sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              regionDefinition=NULL,
                              mutationDefinition=NULL,
                              ambiguousMode=c("eitherOr", "and"),
                              frequency=FALSE,
                              combine=FALSE,
                              nproc=1) {
    # Hack for visibility of foreach index variable
    idx <- NULL
    
    ambiguousMode <- match.arg(ambiguousMode)
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # Check if mutation count/freq columns already exist
    # and throw overwritting warning
    if (!is.null(regionDefinition)) {
        labels <- regionDefinition@labels
    } else {
        labels <- makeNullRegionDefinition()@labels
    }
    if (frequency == TRUE) {
        if (combine) {
            labels <- "MU_FREQ"
        } else {
            labels <- paste("MU_FREQ_", labels, sep="")
        }
    } else {
        if (combine) {
            labels <- "MU_COUNT"
        } else {
            labels <- paste("MU_COUNT_", labels, sep="")
        }
    }
    
    label_exists <- labels[labels %in% colnames(db)]
    if (length(label_exists)>0) {
        warning(paste0("Columns ", 
                       paste(label_exists, collapse=", "),
                       " exist and will be overwritten")
        )
        db[,label_exists] <- NULL
    }
    
    
    # Check mutation definition
    if (!is.null(mutationDefinition) & !is(mutationDefinition, "MutationDefinition")) {
        stop(deparse(substitute(mutationDefinition)), " is not a valid MutationDefinition object")
    }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster <- nproc 
        nproc <- 0
    }
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, cpuCount())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if (nproc > 1) {        
        cluster <- parallel::makeCluster(nproc, type = "PSOCK")
        parallel::clusterExport(cluster, list('db', 'sequenceColumn', 'germlineColumn', 
                                              'regionDefinition', 'frequency', 'combine',
                                              'ambiguousMode', 
                                              'calcObservedMutations','s2c','c2s','NUCLEOTIDES',
                                              'NUCLEOTIDES_AMBIGUOUS', 'IUPAC2nucs',
                                              'EXPANDED_AMBIGUOUS_CODONS',
                                              'makeNullRegionDefinition', 'mutationDefinition',
                                              'getCodonPos','getContextInCodon','mutationType',
                                              'translateCodonToAminoAcid','AMINO_ACIDS',
                                              'binMutationsByRegion', 'countNonNByRegion'), 
                                envir=environment())
        registerDoParallel(cluster)
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Printing status to console
    #cat("Calculating observed number of mutations...\n")
    
    # Identify all the mutations in the sequences
    numbOfSeqs <- nrow(db)
    observedMutations_list <-
        foreach(idx=iterators::icount(numbOfSeqs)) %dopar% {
            oM <- calcObservedMutations(db[[sequenceColumn]][idx], 
                                        db[[germlineColumn]][idx],
                                        frequency=frequency & !combine,
                                        regionDefinition=regionDefinition,
                                        mutationDefinition=mutationDefinition,
                                        returnRaw = combine,
                                        ambiguousMode=ambiguousMode)
            if (combine) {
                num_mutations <- 0
                if (!all(is.na(oM$pos))) {
                    num_mutations <- sum(oM$pos$R, oM$pos$S)
                }
                if (!frequency) {
                    num_mutations
                } else {
                    num_nonN <- sum(oM$nonN)
                    mu_freq <- num_mutations/num_nonN
                    mu_freq
                }
            } else {
                oM   
            }
        }
    
    # Convert list of mutations to data.frame
    if (combine) {
        labels_length <- 1
    } else if (!is.null(regionDefinition)) {
        labels_length <- length(regionDefinition@labels)
    } else{
        #labels_length=1
        labels_length <- length(makeNullRegionDefinition()@labels)
    }
    observed_mutations <- do.call( rbind, lapply(observedMutations_list, function(x) { 
        length(x) <- labels_length 
        return(x)
    }))
    
    
    sep <- "_"
    if (ncol(observed_mutations) > 1) sep <- "_"
    observed_mutations[is.na(observed_mutations)] <- 0
    if (frequency == TRUE) {
        colnames(observed_mutations) <- gsub("_$","",paste("MU_FREQ", colnames(observed_mutations), sep=sep))
    } else {
        colnames(observed_mutations) <- gsub("_$","",paste("MU_COUNT", colnames(observed_mutations), sep=sep))
    }
    
    # Properly shutting down the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, observed_mutations)
    return(db_new)    
}


#' Count the number of observed mutations in a sequence.
#'
#' \code{calcObservedMutations} determines all the mutations in a given input seqeunce 
#' compared to its germline sequence.
#'
#' @param    inputSeq            input sequence. IUPAC ambiguous characters for DNA are 
#'                               supported.
#' @param    germlineSeq         germline sequence. IUPAC ambiguous characters for DNA 
#'                               are supported.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences. Note, only the part of
#'                               sequences defined in \code{regionDefinition} are analyzed.
#'                               If NULL, mutations are counted for entire sequence.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    ambiguousMode       whether to consider ambiguous characters as 
#'                               \code{"either or"} or \code{"and"} when determining and 
#'                               counting the type(s) of mutations. Applicable only if
#'                               \code{inputSeq} and/or \code{germlineSeq} 
#'                               contain(s) ambiguous characters. One of 
#'                               \code{c("eitherOr", "and")}. Default is \code{"eitherOr"}.                               
#' @param    returnRaw           return the positions of point mutations and their 
#'                               corresponding mutation types, as opposed to counts of 
#'                               mutations across positions. Also returns the number of 
#'                               bases used as the denominator when calculating frequency. 
#'                               Default is \code{FALSE}.                               
#' @param    frequency           \code{logical} indicating whether or not to calculate
#'                               mutation frequencies. The denominator used is the number 
#'                               of bases that are not one of "N", "-", or "." in either 
#'                               the input or the germline sequences. If set, this 
#'                               overwrites \code{returnRaw}. Default is \code{FALSE}.
#'                               
#' @return   For \code{returnRaw=FALSE}, an \code{array} with the numbers of replacement (R) 
#'           and silent (S) mutations. 
#'           
#'           For \code{returnRaw=TRUE}, a list containing 
#'           \itemize{
#'                \item \code{$pos}: A data frame whose columns (\code{position}, \code{R}, 
#'                      \code{S}, and \code{region}) indicate, respecitively, the nucleotide 
#'                      position, the number of R mutations at that position, the number of S 
#'                      mutations at that position, and the region in which that nucleotide
#'                      is in.
#'                \item \code{$nonN}: A vector indicating the number of bases in regions 
#'                      defined by \code{regionDefinition} (excluding non-triplet overhang, 
#'                      if any) that are not one of "N", "-", or "." in either the 
#'                      \code{inputSeq} or \code{germlineSeq}.
#'           }
#'           
#'           For \code{frequency=TRUE}, regardless of \code{returnRaw}, an \code{array} 
#'           with the frequencies of replacement (R) and silent (S) mutations.
#'           
#' @details
#' \strong{Each mutation is considered independently in the germline context}. For illustration,
#' consider the case where the germline is \code{TGG} and the observed is \code{TAC}.
#' When determining the mutation type at position 2, which sees a change from \code{G} to 
#' \code{A}, we compare the codon \code{TGG} (germline) to \code{TAG} (mutation at position
#' 2 independent of other mutations in the germline context). Similarly, when determining 
#' the mutation type at position 3, which sees a change from \code{G} to \code{C}, we 
#' compare the codon \code{TGG} (germline) to \code{TGC} (mutation at position 3 independent 
#' of other mutations in the germline context).
#' 
#' If specified, only the part of \code{inputSeq} defined in \code{regionDefinition} is 
#' analyzed. For example, when using the default \link{IMGT_V} definition, then mutations 
#' in positions beyond 312 will be ignored. Additionally, non-triplet overhang at the 
#' sequence end is ignored.
#' 
#' Only replacement (R) and silent (S) mutations are included in the results. \strong{Excluded}
#' are: 
#' \itemize{
#'      \item Stop mutations
#'      
#'            E.g.: the case where \code{TAGTGG} is observed for the germline \code{TGGTGG}.
#'            
#'      \item Mutations occurring in codons where one or both of the observed and the 
#'            germline involve(s) one or more of "N", "-", or ".".
#'            
#'            E.g.: the case where \code{TTG} is observed for the germline being any one of 
#'            \code{TNG}, \code{.TG}, or \code{-TG}. Similarly, the case where any one of 
#'            \code{TTN}, \code{TT.}, or \code{TT-} is observed for the germline \code{TTG}.
#'            
#' }
#' In other words, a result that is \code{NA} or zero indicates absence of R and S mutations, 
#' not necessarily all types of mutations, such as the excluded ones mentioned above.
#' 
#' \code{NA} is also returned if \code{inputSeq} or \code{germlineSeq} is shorter than 3
#' nucleotides.
#' 
#' @section Ambiguous characters:
#' When there are ambiguous characters present, the user could choose how mutations involving
#' ambiguous characters are counted through \code{ambiguousMode}. The two available modes 
#' are \code{"eitherOr"} and \code{"and"}. 
#' \itemize{
#'      \item With \code{"eitherOr"}, ambiguous characters are each expanded but only 
#'            1 mutation is recorded. When determining the type of mutation, the 
#'            priority for different types of mutations, in decreasing order, is as follows:
#'            no mutation, replacement mutation, silent mutation, and stop mutation. 
#'            
#'            When counting the number of non-N, non-dash, and non-dot positions, each
#'            position is counted only once, regardless of the presence of ambiguous
#'            characters.
#'            
#'            As an example, consider the case where \code{germlineSeq} is \code{"TST"} and
#'            \code{inputSeq} is \code{"THT"}. Expanding \code{"H"} at position 2 in 
#'            \code{inputSeq} into \code{"A"}, \code{"C"}, and \code{"T"}, as well as 
#'            expanding \code{"S"} at position 2 in \code{germlineSeq} into \code{"C"} and 
#'            \code{"G"}, one gets:
#'            
#'            \itemize{
#'                 \item \code{"TCT"} (germline) to \code{"TAT"} (observed): replacement
#'                 \item \code{"TCT"} (germline) to \code{"TCT"} (observed): no mutation
#'                 \item \code{"TCT"} (germline) to \code{"TTT"} (observed): replacement 
#'                 \item \code{"TGT"} (germline) to \code{"TAT"} (observed): replacement 
#'                 \item \code{"TGT"} (germline) to \code{"TCT"} (observed): replacement
#'                 \item \code{"TGT"} (germline) to \code{"TTT"} (observed): replacement
#'            }
#'            
#'            Because "no mutation" takes priority over replacement mutation, the final 
#'            mutation count returned for this example is \code{NA} (recall that only R and 
#'            S mutations are returned). The number of non-N, non-dash, and non-dot 
#'            positions is 3.
#'            
#'      \item With \code{"and"}, ambiguous characters are each expanded and mutation(s)
#'            from all expansions are recorded.
#'            
#'            When counting the number of non-N, non-dash, and non-dot positions, if a 
#'            position contains ambiguous character(s) in \code{inputSeq} and/or 
#'            \code{germlineSeq}, the count at that position is taken to be the total 
#'            number of combinations of germline and observed codons after expansion.
#'            
#'            Using the same example from above, the final result returned for this example
#'            is that there are 5 R mutations at position 2. The number of non-N, non-dash,
#'            and non-dot positions is 8, since there are 6 combinations stemming from 
#'            position 2 after expanding the germline codon (\code{"TST"}) and the observed 
#'            codon (\code{"THT"}).
#' }
#' 
#' @seealso  See \link{observedMutations} for counting the number of observed mutations 
#' in a \code{data.frame}.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' in_seq <- ExampleDb[["SEQUENCE_IMGT"]][100]
#' germ_seq <-  ExampleDb[["GERMLINE_IMGT_D_MASK"]][100]
#' 
#' # Identify all mutations in the sequence
#' ex1_raw <- calcObservedMutations(in_seq, germ_seq, returnRaw=TRUE)
#' # Count all mutations in the sequence
#' ex1_count <- calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE)
#' ex1_freq <- calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE, frequency=TRUE)
#' # Compare this with ex1_count
#' table(ex1_raw$pos$region, ex1_raw$pos$R)[, "1"]
#' table(ex1_raw$pos$region, ex1_raw$pos$S)[, "1"]
#' # Compare this with ex1_freq
#' table(ex1_raw$pos$region, ex1_raw$pos$R)[, "1"]/ex1_raw$nonN
#' table(ex1_raw$pos$region, ex1_raw$pos$S)[, "1"]/ex1_raw$nonN
#' 
#' # Identify only mutations the V segment minus CDR3
#' ex2_raw <- calcObservedMutations(in_seq, germ_seq, 
#'                                 regionDefinition=IMGT_V, returnRaw=TRUE)
#' # Count only mutations the V segment minus CDR3
#' ex2_count <- calcObservedMutations(in_seq, germ_seq, 
#'                                   regionDefinition=IMGT_V, returnRaw=FALSE)
#' ex2_freq <- calcObservedMutations(in_seq, germ_seq, 
#'                                  regionDefinition=IMGT_V, returnRaw=FALSE,
#'                                  frequency=TRUE)
#' # Compare this with ex2_count
#' table(ex2_raw$pos$region, ex2_raw$pos$R)[, "1"]
#' table(ex2_raw$pos$region, ex2_raw$pos$S)[, "1"]                              
#' # Compare this with ex2_freq
#' table(ex2_raw$pos$region, ex2_raw$pos$R)[, "1"]/ex2_raw$nonN     
#' table(ex2_raw$pos$region, ex2_raw$pos$S)[, "1"]/ex2_raw$nonN                                       
#' 
#' # Identify mutations by change in hydropathy class
#' ex3_raw <- calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
#'                                 mutationDefinition=HYDROPATHY_MUTATIONS, 
#'                                 returnRaw=TRUE)
#' # Count mutations by change in hydropathy class
#' ex3_count <- calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
#'                                   mutationDefinition=HYDROPATHY_MUTATIONS, 
#'                                   returnRaw=FALSE)
#' ex3_freq <- calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
#'                                  mutationDefinition=HYDROPATHY_MUTATIONS, 
#'                                  returnRaw=FALSE, frequency=TRUE)
#' # Compre this with ex3_count
#' table(ex3_raw$pos$region, ex3_raw$pos$R)[, "1"]
#' table(ex3_raw$pos$region, ex3_raw$pos$S)[, "1"]
#' # Compare this with ex3_freq
#' table(ex3_raw$pos$region, ex3_raw$pos$R)[, "1"]/ex3_raw$nonN                                        
#' table(ex3_raw$pos$region, ex3_raw$pos$S)[, "1"]/ex3_raw$nonN                                        
#'                                 
#' @export
calcObservedMutations <- function(inputSeq, germlineSeq,
                                  regionDefinition=NULL, mutationDefinition=NULL,
                                  ambiguousMode=c("eitherOr", "and"),
                                  returnRaw=FALSE, frequency=FALSE) {
    
    ambiguousMode <- match.arg(ambiguousMode)
    
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # Check mutation definition
    if (!is.null(mutationDefinition) & !is(mutationDefinition, "MutationDefinition")) {
        stop(deparse(substitute(mutationDefinition)), " is not a valid MutationDefinition object")
    }
    
    # IMPORTANT: convert to uppercase 
    # NUCLEOTIDES, NUCLEOTIDES_AMBIGUOUS are in uppercases only
    inputSeq <- toupper(inputSeq)
    germlineSeq <- toupper(germlineSeq)
    
    # Assign mutation definition
    aminoAcidClasses <- if (is.null(mutationDefinition)) { NULL } else { mutationDefinition@classes }
    
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to ZZZ any other . is not an IMGT gap & will be treated like N
    germlineSeq <- gsub("\\.\\.\\.", "ZZZ", germlineSeq)
    #If there is a single gap left convert it to an N
    germlineSeq <- gsub("\\.", "N", germlineSeq)
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    germlineSeq <- gsub("ZZZ", "...", germlineSeq)
    
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to ZZZ any other . is not an IMGT gap & will be treated like N
    inputSeq <- gsub("\\.\\.\\.", "ZZZ", inputSeq)
    #If there is a single gap left convert it to an N
    inputSeq <- gsub("\\.", "N", inputSeq)
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    inputSeq <- gsub("ZZZ", "...", inputSeq)    
    
    # Trim the input and germline sequence to the shortest
    len_inputSeq <- stri_length(inputSeq)
    len_germlineSeq <- stri_length(germlineSeq)
    
    # If a regionDefinition is passed,
    # then only analyze till the end of the defined length
    if(!is.null(regionDefinition)) {
        rdLength  <- regionDefinition@seqLength
    } else {
        rdLength <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        # Create full sequence RegionDefinition object
        regionDefinition <- makeNullRegionDefinition(rdLength)
    }
    len_shortest <- min(c(len_inputSeq, len_germlineSeq, rdLength), na.rm=TRUE)
    
    c_inputSeq <- s2c(inputSeq)[1:len_shortest]
    c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
    
    # If the sequence and germline (which now should be the same length) is shorter
    # than the rdLength, pad it with Ns
    if(len_shortest<rdLength){
        fillWithNs <- array("N",rdLength-len_shortest)
        c_inputSeq <- c( c_inputSeq, fillWithNs)
        c_germlineSeq <- c( c_germlineSeq, fillWithNs)
    }
    
    # length of c_inputSeq and c_germlineSeq should be multiples of 3; if not, trim
    # at this point, c_inputSeq and c_germlineSeq have the same length
    # this is NECESSARY because otherwise the example below could happen:
    # inputSeq 400..402 (codon 134) is "G  " (no info at 401 and 402);
    # c_inputSeq_codons 400..402 will end up being "G" NA NA,
    # which will be turned by the code strsplit(gsub... into
    # "GNA" "NA" (2 codons!)
    stopifnot(length(c_inputSeq)==length(c_germlineSeq))
    seqLen <- length(c_inputSeq)
    # return NA if seqLen shorter than one complete codon
    # consistent with policy that non-triplet overhang is ignored
    if (seqLen<3) {
        tooShort <- TRUE
    } else {
        tooShort <- FALSE
        # if there's non-triplet overhang, trim/ignore
        if ( (seqLen%%3)!=0 ) {
            c_inputSeq <- c_inputSeq[ 1:(seqLen-(seqLen%%3)) ]
            c_germlineSeq <- c_germlineSeq[ 1:(seqLen-(seqLen%%3)) ]
        }
        stopifnot( (length(c_inputSeq)%%3)==0 )
        stopifnot( (length(c_germlineSeq)%%3)==0 )
    }
    
    mutations_array_raw <- NA
    mutations_array <- setNames(object=rep(NA, length(regionDefinition@labels)), 
                                nm=regionDefinition@labels)
    
    if (!tooShort) {
        # locate mutations
        # germline is one of ATGCN and IUPAC ambiguous characters
        # input is one of ATGCN and IUPAC ambiguous characters
        # character mismatch between germline & input (captures both cases like A vs. T, and W vs. T)
        mutations = ( (c_germlineSeq != c_inputSeq) & 
                          (c_germlineSeq %in% NUCLEOTIDES_AMBIGUOUS[1:14]) & 
                          (c_inputSeq %in% NUCLEOTIDES_AMBIGUOUS[1:14]) ) 
        #print(sum(mutations))
        if (sum(mutations) > 0) {
            # The nucleotide positions of the mutations
            mutations_pos <- which(mutations==TRUE)
            # For every mutations_pos, extract the entire codon from germline
            mutations_pos_codons <- array(sapply(mutations_pos, getCodonPos))
            c_germlineSeq_codons <- c_germlineSeq[mutations_pos_codons]
            # For every mutations_pos, extract the codon from input (without other mutations 
            # at the same codon, if any).
            c_inputSeq_codons <- array(sapply(mutations_pos, function(x) {
                seqP <- c_germlineSeq[getCodonPos(x)]
                seqP[getContextInCodon(x)] <- c_inputSeq[x]
                return(seqP) }))
            # split the string of codons into vector of codons
            c_germlineSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_germlineSeq_codons)), " ")[[1]]
            c_inputSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_inputSeq_codons)), " ")[[1]]
            
            # Determine whether the mutations are R or S
            # a table where rows are R/S/Stop/na, cols are codon positions
            # Count ambiguous characters as "eithe-or" or "and" based on user setting 
            
            # Makes use of the fact that c_germlineSeq_codons and c_inputSeqCodons have
            # the same length
            mutations_array_raw <- sapply(1:length(c_germlineSeq_codons),
                                         function(i){
                                             mutationType(codonFrom=c_germlineSeq_codons[i], 
                                                          codonTo=c_inputSeq_codons[i],
                                                          ambiguousMode=ambiguousMode,
                                                          aminoAcidClasses)
                                         })
            
            # check dimension before assigning nucleotide positions to colnames
            stopifnot(ncol(mutations_array_raw)==length(mutations_pos))
            colnames(mutations_array_raw) <- mutations_pos
            
            # keep only columns in which there are R or S mutations; and keep only R and S rows
            # use drop=FALSE so that matrix won't be collapsed into a vector if there is only 1 TRUE in keep.idx
            keep.idx <- apply(mutations_array_raw, 2, function(x) { any(x[c("R", "S")]>0) } )
            keep.pos <- colnames(mutations_array_raw)[keep.idx]
            mutations_array_raw <- mutations_array_raw[c("R", "S"), keep.idx, drop=FALSE]
            colnames(mutations_array_raw) <- keep.pos
            
            # if none of columns have R or S > 1, dim will be 2x0
            if ( ncol(mutations_array_raw)==0 ) {
                # NA if mutations_array_raw contains all NAs and they have all been removed
                mutations_array_raw <- NA
                mutations_array <- setNames(object=rep(NA, length(regionDefinition@labels)), 
                                            nm=regionDefinition@labels)
            } else {
                # count each mutation type by region
                mutations_array <- binMutationsByRegion(mutations_array_raw, regionDefinition)
            }
        }
    }
    
    # frequency=TRUE overrides returnRaw=FALSE/TRUE
    if (frequency) {
        # avoid is.na(mutations_array_raw) to avoid warning in case mutations_array_raw is a vector
        if (length(mutations_array_raw) == sum(is.na(mutations_array_raw))) {
            return(mutations_array)
        } else {
            # Freq = numb of mutations / numb of non N bases (in both seq and gl)
            denoms <- countNonNByRegion(regDef=regionDefinition, ambiMode=ambiguousMode, 
                                        inputChars=c_inputSeq, germChars=c_germlineSeq,
                                        inputCodons=c_inputSeq_codons, 
                                        germCodons=c_germlineSeq_codons, 
                                        mutPos=mutations_pos)
            mutations_array <- mutations_array/rep(denoms, each=2)
            return(mutations_array)
        }
    }
    
    # return positions of point mutations and their mutation types ("raw")
    if (returnRaw){
        
        if (length(mutations_array_raw) == sum(is.na(mutations_array_raw))) {
            # if mutations_array_raw is NA, or 
            # if mutations_array_raw is empty due to all mutations being "Stop" and hence removed
            # avoid is.na(mutations_array_raw) to avoid warning in case mutations_array_raw is a vector
            
            if (!tooShort) {
                # when input and germline are >=3 nucleotides but there's no mutation
                # c_inputSeq_codons, c_germlineSeq_codons, and mutations_pos won't exist
                # this won't be a problem if ambiguousMode="eitherOr", but would for "and"
                # set inputCodons, germCodons, and mutPos to NULL to work around that
                nonN.denoms <- countNonNByRegion(regDef=regionDefinition, ambiMode=ambiguousMode, 
                                                inputChars=c_inputSeq, germChars=c_germlineSeq,
                                                inputCodons=NULL, 
                                                germCodons=NULL, 
                                                mutPos=NULL)
            } else {
                nonN.denoms <- setNames(object=rep(NA, length(regionDefinition@regions)), 
                                       nm=regionDefinition@regions)
            }
            
            return(list(pos=mutations_array_raw, nonN=nonN.denoms))
        } else {
            
            nonN.denoms <- countNonNByRegion(regDef=regionDefinition, ambiMode=ambiguousMode, 
                                            inputChars=c_inputSeq, germChars=c_germlineSeq,
                                            inputCodons=c_inputSeq_codons, 
                                            germCodons=c_germlineSeq_codons, 
                                            mutPos=mutations_pos)
            
            # df indicating position, mutation type (R or S), and region of each mutation
            rawDf <- data.frame(as.numeric(colnames(mutations_array_raw)))
            rawDf <- cbind(rawDf,
                          mutations_array_raw["R", ],
                          mutations_array_raw["S", ],
                          as.character(regionDefinition@boundaries[as.numeric(colnames(mutations_array_raw))]),
                          stringsAsFactors=F)
            colnames(rawDf) <- c("position", "R", "S", "region")
            return(list(pos=rawDf, nonN=nonN.denoms))
        }
    } else {
        # return counts of each mutation type  
        return(mutations_array)
    }
}


# Aggregate mutations by region
#
# \code{binMutationsByRegion} takes an array of observed mutations (e.g. from 
# \code{calcObservedMutations}) and bins them by the different regions defined in the 
# \code{regionDefinition}.
#
# @param   mutationsArray     \code{array} containing the number of R and S mutations 
#                             at the nucleotide positions where there are mutations.                             
# @param   regionDefinition   \link{RegionDefinition} object defining the regions
#                             and boundaries of the Ig sequences.
# 
# @return An \code{array} of R/S mutations binned across all the unique regions defined
# by \code{regionDefinition}.
# 
# @details
# Note, only the part of sequences defined in \code{regionDefinition} are analyzed.
# For example, if the default \link{IMGT_V} definition is used, then mutations
# in positions beyond 312 will be ignored.
# 
# @seealso  
# See \link{observedMutations} for identifying and counting the 
# number of observed mutations.
# This function is also used in \link{calcObservedMutations}.
# 
# @examples
# # Generate a random mutation array
# numbOfMutPos <- sample(3:10, 1)
# posOfMutations <- sort(sample(330, numbOfMutPos))
# mutations_array <- matrix(0, nrow=2, ncol=numbOfMutPos, dimnames=list(c("R", "S"), posOfMutations))
# mutations_array["R", ] = sample(x=0:10, size=numbOfMutPos, replace=TRUE)
# mutations_array["S", ] = sample(x=0:10, size=numbOfMutPos, replace=TRUE)

# # Random mutations
# binMutationsByRegion(mutations_array, regionDefinition=NULL)
# binMutationsByRegion(mutations_array, regionDefinition=IMGT_V)
binMutationsByRegion <- function(mutationsArray, 
                                 regionDefinition=NULL) {
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # Create full sequence RegionDefinition object 
    # The seqLength will be the largest index of a mutation
    if (is.null(regionDefinition)) {
        regionDefinition <- makeNullRegionDefinition(max(as.numeric(colnames(mutationsArray))))
    }
    
    # get 2 vectors, 1 for R, 1 for S, along length of 1:regionDefinition@seqLength
    # each vector records the number of R/S at each position
    mutatedPositions <- as.numeric(colnames(mutationsArray)) 
    
    mutations_R <- array(NA,  dim=regionDefinition@seqLength)
    mutations_S <- array(NA,  dim=regionDefinition@seqLength)
    mutations_R[mutatedPositions] <- mutationsArray["R", ]
    mutations_S[mutatedPositions] <- mutationsArray["S", ]
    mutations_R <- mutations_R[1:regionDefinition@seqLength]
    mutations_S <- mutations_S[1:regionDefinition@seqLength]
    
    # count number of R/S in each region
    mutations_region_counts <- rep(0, length(regionDefinition@labels))
    names(mutations_region_counts) <- regionDefinition@labels
    for (reg in regionDefinition@regions) {
        mutations_region_counts[paste0(reg, "_R")] <- sum(mutations_R[regionDefinition@boundaries==reg], na.rm=T)
        mutations_region_counts[paste0(reg, "_S")] <- sum(mutations_S[regionDefinition@boundaries==reg], na.rm=T)
    }
    
    return(mutations_region_counts)
}

# Count the number of non-N, non-dash, and non-dot positions
# 
# @param   regDef       regionDefinition
# @param   ambiMode     ambiguousMode
# @param   inputChars   c_inputSeq
# @param   germChars    c_germlineSeq
# @param   inputCodons  c_inputSeq_codons
# @param   germCodons   c_germlineSeq_codons
# @param   mutPos       mutations_pos
# 
# @return  The number of non-N, non-dash, and non-dot characters. Calculation method
#          differs depending on ambiMode being "eitherOr" or "and". By design, when 
#          there is no ambiguous character in the input or germline, the result should be
#          the same regardless of ambiMode.
#
# @details This is a helper function for calcObservedMutations() and is not intended to
#          be called directly. All input arguments are, by design, expected to be 
#          generated as intermediate products during a call to calcObservedMutations().
#          
countNonNByRegion <- function(regDef, ambiMode, inputChars, germChars,
                             inputCodons, germCodons, mutPos) {
    
    regionNames <- unique(sapply(regDef@labels, 
                                function(x) { substr(x, 1, stri_length(x)-2) }))
    
    if (ambiMode=="eitherOr") {
        
        # Subset boundaries to only non-N & non-dash & non-dot bases (in both seq and gl)
        # "which" in next line is ESSENTIAL; otherwise @boundaries won't be truncated
        # e.g. (1:6)[c(T,T,T)] returns 1:6, not 1:3
        boundaries <- regDef@boundaries[
            which(inputChars %in% NUCLEOTIDES_AMBIGUOUS[1:14] &  
                  germChars %in% NUCLEOTIDES_AMBIGUOUS[1:14])]
        
        # number of non-N & non-dash & non-dot bases (in both seq and gl)
        nonN <- sapply(regionNames, function(x) { sum(boundaries==x) })
        
    } else if (ambiMode=="and") {
        
        ### positions where there's no mutation:
        # simply count the positions where both input and germline are 
        # non-N, non-dash, and non-dot
        boundaries.1 <- regDef@boundaries[
            which(inputChars %in% NUCLEOTIDES_AMBIGUOUS[1:14] &  
                  germChars %in% NUCLEOTIDES_AMBIGUOUS[1:14] &
                  (germChars == inputChars))]
        nonN.1 <- sapply(regionNames, function(x) { sum(boundaries.1==x) })
        
        ### positions where there's mutation:
        if ( (!is.null(inputCodons)) & (!is.null(germCodons)) & (!is.null(mutPos)) ) {
            # expand codon with ambiguous character(s) into codons with unambiguous characters
            # calculate the number of possible combinations between input and germline codons
            
            # this makes use of the important fact that each mutation is considered 
            # independently in the germline context
            inputNumExpanded <- sapply(inputCodons, 
                                      function(codon){
                                          length(EXPANDED_AMBIGUOUS_CODONS[[codon]])
                                      })
            germlineNumExpanded <- sapply(germCodons, 
                                         function(codon){
                                             length(EXPANDED_AMBIGUOUS_CODONS[[codon]])
                                         })
            totalNumExpanded <- inputNumExpanded * germlineNumExpanded
            
            # use mutations_pos to capture positions at which R/S is absent (Stop or na instead)
            # such positions would have been omitted from mutations_array_raw or mutations_array
            boundaries.2 <- regDef@boundaries[mutPos]
            # makes use of the fact that inputCodons, germCodons, and 
            # mutPos align exactly
            nonN.2 <- sapply(regionNames, function(x){ sum(totalNumExpanded[boundaries.2==x]) })
        } else {
            nonN.2 <- setNames(object=rep(0, length(regionNames)), nm=regionNames)
        }
        nonN <- nonN.1 + nonN.2
    }
    return(nonN)
}


#### Sliding window approach ####
#' Sliding window approach towards filtering a single sequence
#'
#' \code{slideWindowSeq} determines whether an input sequence contains equal to or more than 
#' a given number of mutations in a given length of consecutive nucleotides (a "window") 
#' when compared to a germline sequence.
#' 
#' @param    inputSeq            input sequence.
#' @param    germlineSeq         germline sequence.
#' @param    mutThresh           threshold on the number of mutations in \code{windowSize} 
#'                               consecutive nucleotides. Must be between 1 and \code{windowSize} 
#'                               inclusive. 
#' @param    windowSize          length of consecutive nucleotides. Must be at least 2.
#'                               
#' @return  \code{TRUE} if there are equal to or more than \code{mutThresh} number of mutations
#'          in any window of \code{windowSize} consecutive nucleotides (i.e. the sequence should
#'          be filtered); \code{FALSE} if otherwise.
#' 
#' @seealso  \link{calcObservedMutations} is called by \code{slideWindowSeq} to identify observed 
#'           mutations. See \link{slideWindowDb} for applying the sliding window approach on a 
#'           \code{data.frame}. See \link{slideWindowTune} for parameter tuning for \code{mutThresh}
#'           and \code{windowSize}.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' in_seq <- ExampleDb[100, "SEQUENCE_IMGT"]
#' germ_seq <-  ExampleDb[100, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Determine if in_seq has 6 or more mutations in 10 consecutive nucleotides
#' slideWindowSeq(inputSeq=in_seq, germlineSeq=germ_seq, mutThresh=6, windowSize=10)
#'                                 
#' @export
slideWindowSeq <- function(inputSeq, germlineSeq, mutThresh, windowSize){
    # identify all R and S mutations in input sequence
    inputMut <- calcObservedMutations(inputSeq=inputSeq, germlineSeq=germlineSeq, returnRaw=T)
    
    # call helper
    return(slideWindowSeqHelper(mutPos=inputMut$pos, mutThresh=mutThresh, windowSize=windowSize))
}


# NOTE: DO NOT MERGE slideWindowSeqHelper with slideWindowSeq (very different input formats)
#       slideWindowTune needs to call slideWindowSeqHelper directly for efficiency
# Helper for sliding window approach towards filtering sequences
#
# @param    mutPos              a \code{data.frame} containing positions and types of point 
#                               mutations as returned in \code{$pos} by 
#                               \code{calcObserverdMutations()} with \code{returnRaw=TRUE}. 
#                               Can be \code{NA}, in which case the returned value will be 
#                               \code{FALSE}.
# @param    mutThresh           threshold on the number of mutations in \code{windowSize} 
#                               consecutive nucleotides. Must be between 1 and \code{windowSize} 
#                               inclusive.
# @param    windowSize          length of consecutive nucleotides. Must be at least 2.
#
# @return   \code{TRUE} if there are equal to or more than \code{mutThresh} number of mutations
#           in any window of \code{windowSize} consecutive nucleotides; \code{FALSE} if otherwise.
#
slideWindowSeqHelper <- function(mutPos, mutThresh, windowSize){
    # check preconditions
    stopifnot(mutThresh >= 1 & mutThresh <= windowSize & windowSize>=2)
    
    if (length(mutPos) == 1 && is.na(mutPos)) {
        # use && instead of & to short-circuit in case length(mutPos)!=1 (otherwise warning)
        return(FALSE)
    } else {
        # general idea:
        # only need to check windows containing mutations (as opposed to every possible window)
        for (i in 1:nrow(mutPos)){
            # get window limits
            lower <- mutPos$position[i]
            upper <- lower + windowSize - 1
            # how many mutations fall within current window
            windowCount <- sum(mutPos[mutPos$position>=lower & mutPos$position<=upper, c("R","S")])
            # return as soon as a window has >= mutThresh mutations
            if (windowCount >= mutThresh) { return(TRUE) }
        }
        
        return(FALSE)
    }
}


#' Sliding window approach towards filtering sequences in a \code{data.frame}
#'
#' \code{slideWindowDb} determines whether each input sequence in a \code{data.frame} 
#' contains equal to or more than a given number of mutations in a given length of 
#' consecutive nucleotides (a "window") when compared to their respective germline 
#' sequence.
#' 
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    mutThresh           threshold on the number of mutations in \code{windowSize} 
#'                               consecutive nucleotides. Must be between 1 and \code{windowSize} 
#'                               inclusive. 
#' @param    windowSize          length of consecutive nucleotides. Must be at least 2.
#'                               
#' @return   a logical vector. The length of the vector matches the number of input sequences in 
#'           \code{db}. Each entry in the vector indicates whether the corresponding input sequence
#'           should be filtered based on the given parameters.
#' 
#' @seealso  See \link{slideWindowSeq} for applying the sliding window approach on a single sequence. 
#'           See \link{slideWindowTune} for parameter tuning for \code{mutThresh} and \code{windowSize}.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' 
#' # Apply the sliding window approach on a subset of ExampleDb
#' slideWindowDb(db = ExampleDb[1:10, ], mutThresh=6, windowSize=10)
#' 
#' @export
slideWindowDb <- function(db, sequenceColumn="SEQUENCE_IMGT", 
                          germlineColumn="GERMLINE_IMGT_D_MASK",
                          mutThresh, windowSize){
    db_filter <- sapply(1:nrow(db), function(i) { slideWindowSeq(inputSeq = db[i, sequenceColumn],
                                                                 germlineSeq = db[i, germlineColumn],
                                                                 mutThresh = mutThresh,
                                                                 windowSize = windowSize)})
    return(db_filter)
}


#' Parameter tuning for sliding window approach
#'
#' Apply \link{slideWindowDb} over a search grid made of combinations of \code{mutThresh} and 
#' \code{windowSize} to help with picking a pair of values for these parameters. Parameter 
#' tuning can be performed by choosing a combination that gives a reasonable number of 
#' filtered/remaining sequences. 
#' 
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    dbMutList           if supplied, this should be a list consisting of \code{data.frame}s 
#'                               returned as \code{$pos} in the nested list produced by 
#'                               \link{calcObservedMutations} with \code{returnRaw=TRUE}; otherwise, 
#'                               \link{calcObservedMutations} is called on columns \code{sequenceColumn}
#'                               and \code{germlineColumn} of \code{db}. Default is \code{NULL}. 
#' @param    mutThreshRange      range of threshold on the number of mutations in \code{windowSize} 
#'                               consecutive nucleotides to try. Must be between 1 and 
#'                               maximum \code{windowSizeRange} inclusive. 
#' @param    windowSizeRange     range of length of consecutive nucleotides to try. The lower end
#'                               must be at least 2.
#' @param    verbose             whether to print out messages indicating current progress. Default
#'                               is \code{TRUE}.              
#'                               
#' @return   a list of logical matrices. Each matrix corresponds to a \code{windowSize} in 
#'           \code{windowSizeRange}. Each column in a matrix corresponds to a \code{mutThresh} in
#'           \code{mutThreshRange}.
#' 
#' @details  If, in a given combination of \code{mutThresh} and \code{windowSize}, \code{mutThresh} 
#'           is greater than \code{windowSize}, \code{NA}s will be returned for that particular
#'           combination. A message indicating that the combination has been "skipped" will be 
#'           printed if \code{verbose=TRUE}.
#'           
#'           If \link{calcObservedMutations} was previously run on \code{db} and saved, supplying
#'           \code{$pos} from the saved result as \code{dbMutList} could save time by skipping a
#'           second call of \link{calcObservedMutations}. This could be helpful especially when 
#'           \code{db} is large.
#' 
#' @seealso  \link{slideWindowDb} is called on \code{db} for tuning. See \link{slideWindowTunePlot} 
#'           for visualization. See \link{calcObservedMutations} for generating \code{dbMutList}.
#' 
#' @examples
#' # Load and subset example data
#' data(ExampleDb, package="alakazam")
#' db <- ExampleDb[1:5, ]
#' 
#' # Try out thresholds of 2-4 mutations in window sizes of 7-9 nucleotides. 
#' # In this case, all combinations are legal.
#' slideWindowTune(db, mutThreshRange=2:4, windowSizeRange=7:9)
#' 
#' # Illegal combinations are skipped, returning NAs.
#' slideWindowTune(db, mutThreshRange=2:4, windowSizeRange=2:4, 
#'                 verbose=FALSE)
#'                                                             
#' # Run calcObservedMutations separately
#' exDbMutList <- sapply(1:5, function(i) {
#'     calcObservedMutations(inputSeq=db[i, "SEQUENCE_IMGT"],
#'                           germlineSeq=db[i, "GERMLINE_IMGT_D_MASK"],
#'                           returnRaw=TRUE)$pos })
#' slideWindowTune(db, dbMutList=exDbMutList, 
#'                 mutThreshRange=2:4, windowSizeRange=2:4)
#'                                                            
#' @export
slideWindowTune <- function(db, sequenceColumn="SEQUENCE_IMGT", 
                            germlineColumn="GERMLINE_IMGT_D_MASK",
                            dbMutList=NULL,
                            mutThreshRange, windowSizeRange, verbose=TRUE){
    # check preconditions
    stopifnot(!is.null(db))
    stopifnot(min(mutThreshRange) >= 1 & 
                  max(mutThreshRange) <= max(windowSizeRange) &
                  min(windowSizeRange) >= 2)
    
    
    # get positions of R/S mutations for sequences in db
    # do this here and then call slideWindowSeqHelper (so it's done only once)
    # instead of calling slideWindowDb which does this every time it is called
    if (is.null(dbMutList)) {
        inputMutList <- sapply(1:nrow(db), 
                              function(i){
                                  calcObservedMutations(inputSeq=db[i, sequenceColumn],
                                                        germlineSeq=db[i, germlineColumn],
                                                        returnRaw=T)$pos})    
    } else {
        if (verbose) {cat("dbMutList supplied; skipped calling calcObservedMutations()\n")}
        inputMutList <- dbMutList
    }
    
    # apply slideWindow on combinations of windowSize and mutThresh
    for (size in windowSizeRange) {
        if (verbose) {cat(paste0("now computing for windowSize = ", size, "\n"))}
        
        for (thresh in mutThreshRange) {
            if (thresh <= size){
                if (verbose) {cat(paste0(">>> mutThresh = ", thresh, "\n"))}
                # apply slideWindow using current pair of parameters
                cur.logical <- unlist(lapply(inputMutList,
                                            slideWindowSeqHelper,
                                            mutThresh = thresh, windowSize = size))
            } else {
                if (verbose) {cat(paste0(">>> mutThresh = ", thresh, " > windowSize = ", 
                                         size, " (skipped)\n"))}
                # NA if skipped
                cur.logical <- rep(NA, nrow(db))
            }
            # store results for each thresh as a column in a logical matrix
            if (thresh == mutThreshRange[1]) {
                cur.mtx <- matrix(data=cur.logical, nrow=length(cur.logical))
            } else {
                cur.mtx <- cbind(cur.mtx, cur.logical)
            }
        }
        colnames(cur.mtx) <- as.character(mutThreshRange)
        
        # store results for each size (and threshes under that size) as a logical matrix in a list
        if (size == windowSizeRange[1]) {
            cur.list <- list(cur.mtx)
        } else {
            cur.list <- c(cur.list, list(cur.mtx))
        }
    }
    names(cur.list) <- as.character(windowSizeRange)
    
    return(cur.list)
}


#' Visualize parameter tuning for sliding window approach
#'
#' Visualize results from \link{slideWindowTune}
#' 
#' @param    tuneList            a list of logical matrices returned by \link{slideWindowTune}.
#' @param    plotFiltered        whether to plot the number of filtered sequences (as opposed to
#'                               the number of remaining sequences). Default is \code{TRUE}.
#' @param    percentage          whether to plot on the y-axis the percentage of filtered sequences
#'                               (as opposed to the absolute number). Default is \code{FALSE}.                             
#' @param    jitter.x            whether to jitter x-axis values. Default is \code{FALSE}.                               
#' @param    jitter.x.amt        amount of jittering to be applied on x-axis values if 
#'                               \code{jitter.x=TRUE}. Default is 0.1.
#' @param    jitter.y            whether to jitter y-axis values. Default is \code{FALSE}.
#' @param    jitter.y.amt        amount of jittering to be applied on y-axis values if 
#'                               \code{jitter.y=TRUE}. Default is 0.1.                               
#' @param    pchs                point types to pass on to \link{plot}.
#' @param    ltys                line types to pass on to \link{plot}.
#' @param    cols                colors to pass on to \link{plot}.                             
#' @param    plotLegend          whether to plot legend. Default is \code{TRUE}.
#' @param    legendPos           position of legend to pass on to \link{legend}. Can be either a
#'                               numeric vector specifying x-y coordinates, or one of 
#'                               \code{"topright"}, \code{"center"}, etc. Default is \code{"topright"}.
#' @param    legendHoriz         whether to make legend horizontal. Default is \code{FALSE}.
#' @param    legendCex           numeric values by which legend should be magnified relative to 1.
#' @param    title               plot main title. Default is NULL (no title)
#' 
#' @details  For each \code{windowSize}, the numbers of sequences filtered or remaining after applying
#'           the sliding window approach are plotted on the y-axis against thresholds on the number of
#'           mutations in a window on the x-axis.
#'           
#'           When plotting, a user-defined \code{amount} of jittering can be applied on values plotted
#'           on either axis or both axes via adjusting \code{jitter.x}, \code{jitter.y}, 
#'           \code{jitter.x.amt} and \code{jitter.y.amt}. This may be help with visually distinguishing
#'           lines for different window sizes in case they are very close or identical to each other. 
#'           If plotting percentages (\code{percentage=TRUE}) and using jittering on the y-axis values 
#'           (\code{jitter.y=TRUE}), it is strongly recommended that \code{jitter.y.amt} be set very
#'           small (e.g. 0.01). 
#'           
#'           \code{NA} for a combination of \code{mutThresh} and \code{windowSize} where 
#'           \code{mutThresh} is greater than \code{windowSize} will not be plotted. 
#' 
#' @seealso  See \link{slideWindowTune} for how to get \code{tuneList}. See \link{jitter} for 
#'           use of \code{amount} of jittering.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' 
#' # Try out thresholds of 2-4 mutations in window sizes of 3-5 nucleotides 
#' # on a subset of ExampleDb
#' tuneList <- slideWindowTune(db = ExampleDb[1:10, ], 
#'                            mutThreshRange = 2:4, windowSizeRange = 3:5,
#'                            verbose = FALSE)
#'
#' # Visualize
#' # Plot numbers of sequences filtered without jittering y-axis values
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
#'                     plotFiltered=TRUE, jitter.y=FALSE)
#'                     
#' # Notice that some of the lines overlap
#' # Jittering could help
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
#'                     plotFiltered=TRUE, jitter.y=TRUE)
#'                     
#' # Plot numbers of sequences remaining instead of filtered
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
#'                     plotFiltered=FALSE, jitter.y=TRUE, 
#'                     legendPos="bottomright")
#'                     
#' # Plot percentages of sequences filtered with a tiny amount of jittering
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
#'                     plotFiltered=TRUE, percentage=TRUE, 
#'                     jitter.y=TRUE, jitter.y.amt=0.01)
#'                                                             
#' @export
slideWindowTunePlot <- function(tuneList, plotFiltered = TRUE, percentage = FALSE,
                               jitter.x = FALSE, jitter.x.amt = 0.1,
                               jitter.y = FALSE, jitter.y.amt = 0.1,
                               pchs = 1, ltys = 2, cols = 1,
                               plotLegend = TRUE, legendPos = "topright", 
                               legendHoriz = FALSE, legendCex = 1, title=NULL){
    
    # invert (!) tuneList if plotting retained sequences
    ylab.part.2 <- "filtered"
    if (!plotFiltered) {
        tuneList <- lapply(tuneList, function(x){!x})
        ylab.part.2 <- "remaining"}
    
    # if number of pchs/ltys/cols provided does not match number of lines expected
    # expand into vector with repeating values (otherwise legend would break)
    if (length(pchs)!=length(tuneList)) {pchs <- rep(pchs, length.out=length(tuneList))}
    if (length(ltys)!=length(tuneList)) {ltys <- rep(ltys, length.out=length(tuneList))}
    if (length(cols)!=length(tuneList)) {cols <- rep(cols, length.out=length(tuneList))}
    
    # tabulate tuneList (and if applicable convert to percentage)
    plotList <- lapply(tuneList, colSums)
    if (percentage) {plotList <- lapply(plotList, function(x){x/nrow(tuneList[[1]])})}
    
    # get x-axis values (i.e. mutThreshRange; colnames of matrix in tuneList with most columns)
    #threshes = as.numeric(colnames(tuneList[[which.max(lapply(lapply(tuneList, colnames), length))]]))
    threshes <- as.numeric(colnames(tuneList[[1]]))
    
    # plot for first window size
    x1 <- threshes
    if (jitter.x) {x1 <- jitter(x1, amount=jitter.x.amt)}
    y1 <- plotList[[1]]
    if (jitter.y) {y1 <- jitter(y1, amount=jitter.y.amt)}
    
    if (percentage) {
        ylab.part.1 <- "Percentage of sequences"
        # ylim
        ylim.padding <- abs(diff(range(plotList, na.rm=T)))*0.01
        ylims <- c(max(0, min(range(plotList, na.rm=T)) - ylim.padding), 
                  min(1, max(range(plotList, na.rm=T)) + ylim.padding) )
        
    } else {
        ylab.part.1 <- "Number of sequences"
        # ylim: non-negative lower limit; upper limit slight above max tabulated sum
        ylims <- c( max(0, min(range(plotList, na.rm=T)) - max(1, jitter.y.amt) ), 
                   max(range(plotList, na.rm=T)) + max(1, jitter.y.amt) )
    }
    
    plot(x1, # mutThreshRange on x-axis
         y1, # tabulated sums in plotList on y-axis
         ylim = ylims,
         # xlim: +/- jitter.x.amt*2 to accommodate for amount of jittering on x-axis
         xlim = c(min(threshes)-jitter.x.amt*2, max(threshes+jitter.x.amt*2)),
         xaxt="n", xlab="Threshold on number of mutations",
         ylab=paste(ylab.part.1, ylab.part.2),
         cex.lab=1.5, cex.axis=1.5, type="b", lwd=1.5,
         pch=pchs[1], lty=ltys[1], col=cols[1])
    axis(side=1, at=threshes, cex.axis=1.5)
    
    # add title
    if (!is.null(title)) {
        title(main=title)
    }
    
    # plot for the rest of the window sizes
    for (i in 1:length(plotList)){
        if (i>=2) {
            
            xi <- threshes
            if (jitter.x) {xi <- jitter(xi, amount=jitter.x.amt)}
            yi <- plotList[[i]]
            if (jitter.y) {yi <- jitter(yi, amount=jitter.y.amt)}
            
            points(xi, yi, type='b', lwd=1.5,
                   pch=pchs[i], lty=ltys[i], col=cols[i])
        }
    }
    
    # add legend
    if (plotLegend) {
        # if legendPos specified as xy coordinates
        if (is.numeric(legendPos) & length(legendPos)==2) {
            legend(x=legendPos[1], y=legendPos[2], 
                   legend = c("Window Size", names(tuneList)),
                   horiz = legendHoriz, cex = legendCex,
                   pch=c(NA, pchs), lty=c(NA, ltys), col=c(NA, cols))
        } else {
            # if legendPos specified as "center", "topright", etc.  
            legend(legendPos, 
                   legend = c("Window Size", names(tuneList)),
                   horiz = legendHoriz, cex = legendCex,
                   pch=c(NA, pchs), lty=c(NA, ltys), col=c(NA, cols))
        }
    }
    
}


#### Expected frequencies calculating functions ####

#' Calculate expected mutation frequencies
#'
#' \code{expectedMutations} calculates the expected mutation frequencies for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      \code{character} name of the column containing input 
#'                               sequences.
#' @param    germlineColumn      \code{character} name of the column containing 
#'                               the germline or reference sequence.
#' @param    targetingModel      \link{TargetingModel} object. Default is \link{HH_S5F}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    nproc               \code{numeric} number of cores to distribute the operation
#'                               over. If the cluster has already been set the call function with 
#'                               \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                               \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with expected mutation frequencies 
#'           for each region defined in \code{regionDefinition}.
#'          
#'           The columns names are dynamically created based on the regions in  
#'           \code{regionDefinition}. For example, when using the \link{IMGT_V}
#'           definition, which defines positions for CDR and FWR, the following columns are
#'           added:  
#'           \itemize{
#'             \item  \code{MU_EXPECTED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{MU_EXPECTED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{MU_EXPECTED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{MU_EXPECTED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           
#' @details
#' Only the part of the sequences defined in \code{regionDefinition} are analyzed. 
#' For example, when using the \link{IMGT_V} definition, mutations in
#' positions beyond 312 will be ignored.
#' 
#' @seealso  
#' \link{calcExpectedMutations} is called by this function to calculate the expected 
#' mutation frequencies. See \link{observedMutations} for getting observed 
#' mutation counts. See \link{IMGT_SCHEMES} for a set of predefined 
#' \link{RegionDefinition} objects.
#' 
#' @examples
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
#'
#' # Calculate expected mutations over V region
#' db_exp <- expectedMutations(db,
#'                             sequenceColumn="SEQUENCE_IMGT",
#'                             germlineColumn="GERMLINE_IMGT_D_MASK",
#'                             regionDefinition=IMGT_V,
#'                             nproc=1)
#' 
#' # Calculate hydropathy expected mutations over V region
#' db_exp <- expectedMutations(db,
#'                            sequenceColumn="SEQUENCE_IMGT",
#'                            germlineColumn="GERMLINE_IMGT_D_MASK",
#'                            regionDefinition=IMGT_V,
#'                            mutationDefinition=HYDROPATHY_MUTATIONS,
#'                            nproc=1)
#'
#' @export
expectedMutations <- function(db, 
                              sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              targetingModel=HH_S5F,
                              regionDefinition=NULL,
                              mutationDefinition=NULL,
                              nproc=1) {
    # Hack for visibility of foreach index variable
    idx <- NULL
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # Check mutation definition
    if (!is.null(mutationDefinition) & !is(mutationDefinition, "MutationDefinition")) {
        stop(deparse(substitute(mutationDefinition)), " is not a valid MutationDefinition object")
    }
    
    # Check if mutation count/freq columns already exist
    # and throw overwritting warning
    if (!is.null(regionDefinition)) {
        labels <- regionDefinition@labels
    } else {
        labels <- makeNullRegionDefinition()@labels
    }
    
    labels <- paste("MU_EXPECTED_", labels, sep="")
    
    label_exists <- labels[labels %in% colnames(db)]
    if (length(label_exists)>0) {
        warning(paste0("Columns ", 
                       paste(label_exists, collapse=", "),
                       " exist and will be overwritten")
        )
        db[label_exists] <- NULL
    }    
    
    # Check targeting model
    if (!is(targetingModel, "TargetingModel")) {
        stop(deparse(substitute(targetingModel)), " is not a valid TargetingModel object")
    }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, cpuCount(), na.rm=T)
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if (nproc > 1) {        
        cluster <- parallel::makeCluster(nproc, type = "PSOCK")
        parallel::clusterExport(cluster, list('db', 'sequenceColumn', 'germlineColumn', 
                                              'regionDefinition','targetingModel',
                                              'calcExpectedMutations','calculateTargeting',
                                              's2c','c2s','NUCLEOTIDES','HH_S5F',
                                              'calculateMutationalPaths','CODON_TABLE'),
                                envir=environment() )
        registerDoParallel(cluster)
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    
    # Printing status to console
    cat("Calculating the expected frequencies of mutations...\n")
    
    # Calculate targeting for each sequence (based on the germline)
    # Should be a 5 x N matrix where N in the number of nucleotides defined by
    # the regionDefinition
    numbOfSeqs <- nrow(db)
    
    targeting_list <-
        foreach (idx=iterators::icount(numbOfSeqs)) %dopar% {
            calcExpectedMutations(germlineSeq=db[[germlineColumn]][idx],
                                  inputSeq=db[[sequenceColumn]][idx],
                                  targetingModel=targetingModel,
                                  regionDefinition=regionDefinition,
                                  mutationDefinition=mutationDefinition)
        }
    
    # Convert list of expected mutation freq to data.frame
    if (is.null(regionDefinition)) {
        labels_length <- length(makeNullRegionDefinition()@labels)
    } else {
        labels_length <- length(regionDefinition@labels)
    }
    expectedMutationFrequencies <- do.call(rbind, lapply(targeting_list, function(x) { 
        length(x) <- labels_length 
        return(x) })) 
    
    expectedMutationFrequencies[is.na(expectedMutationFrequencies)] <- 0
    colnames(expectedMutationFrequencies) <- paste0("MU_EXPECTED_", colnames(expectedMutationFrequencies))
    
    # Properly shutting down the cluster
    if(nproc>1){ parallel::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, expectedMutationFrequencies)
    return(db_new)    
    
}


#' Calculate expected mutation frequencies of a sequence
#'
#' \code{calcExpectedMutations} calculates the expected mutation
#' frequencies of a given sequence. This is primarily a helper function for
#' \link{expectedMutations}. 
#'
#' @param    germlineSeq         germline (reference) sequence.
#' @param    inputSeq            input (observed) sequence. If this is not \code{NULL}, 
#'                               then \code{germlineSeq} will be processed to be the same
#'                               same length as \code{inputSeq} and positions in 
#'                               \code{germlineSeq} corresponding to positions with Ns in 
#'                               \code{inputSeq} will also be assigned an N. 
#' @param    targetingModel      \link{TargetingModel} object. Default is \link{HH_S5F}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#'                               
#' @return   A \code{numeric} vector of the expected frequencies of mutations in the 
#'           regions in the \code{regionDefinition}. For example, when using the default 
#'           \link{IMGT_V} definition, which defines positions for CDR and 
#'           FWR, the following columns are calculated:
#'           \itemize{
#'              \item  \code{MU_EXPECTED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                             CDR2 of the V-segment.
#'              \item  \code{MU_EXPECTED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                             of the V-segment.
#'              \item  \code{MU_EXPECTED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                             FWR2 and FWR3 of the V-segment.
#'              \item  \code{MU_EXPECTED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                             FWR3 of the V-segment.
#'            }
#'           
#' @details
#' \code{calcExpectedMutations} calculates the expected mutation frequencies of a 
#' given sequence and its germline. 
#' 
#' Note, only the part of the sequences defined in \code{regionDefinition} are analyzed. 
#' For example, when using the default \link{IMGT_V} definition, mutations in
#' positions beyond 312 will be ignored.
#' 
#' @seealso  \link{expectedMutations} calls this function.
#' To create a custom \code{targetingModel} see \link{createTargetingModel}.
#' See \link{calcObservedMutations} for getting observed mutation counts.
#' 
#' @examples
#' # Load example data
#' data(ExampleDb, package="alakazam")
#' 
#' # Use first entry in the exampled data for input and germline sequence
#' in_seq <- ExampleDb[["SEQUENCE_IMGT"]][1]
#' germ_seq <-  ExampleDb[["GERMLINE_IMGT_D_MASK"]][1]
#' 
#' # Identify all mutations in the sequence
#' calcExpectedMutations(in_seq, germ_seq)
#' 
#' # Identify only mutations the V segment minus CDR3
#' calcExpectedMutations(in_seq, germ_seq, regionDefinition=IMGT_V)
#' 
#' # Define mutations based on hydropathy
#' calcExpectedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
#'                       mutationDefinition=HYDROPATHY_MUTATIONS)
#' 
#' @export
calcExpectedMutations <- function(germlineSeq,
                                  inputSeq=NULL,
                                  targetingModel=HH_S5F,
                                  regionDefinition=NULL,
                                  mutationDefinition=NULL) {
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # Check mutation definition
    if (!is.null(mutationDefinition) & !is(mutationDefinition, "MutationDefinition")) {
        stop(deparse(substitute(mutationDefinition)), " is not a valid MutationDefinition object")
    }
    
    # Check targeting model
    if (!is(targetingModel, "TargetingModel")) {
        stop(deparse(substitute(targetingModel)), " is not a valid TargetingModel object")
    }
    
    # Mask ambiguous nucleotide characters
    germlineSeq <- gsub("[MRWSYKVHDB]", "N", germlineSeq)
    
    # Assign codon table
    codonTable <- if (is.null(mutationDefinition)) { CODON_TABLE } else { mutationDefinition@codonTable }
    
    # Get targeting
    targeting <- calculateTargeting(germlineSeq=germlineSeq, 
                                    inputSeq=inputSeq,
                                    targetingModel=targetingModel,
                                    regionDefinition=regionDefinition)
    
    # Determine the mutations paths (i.e. determine R and S mutation frequencies)
    mutationalPaths <- calculateMutationalPaths(germlineSeq=c2s(colnames(targeting)), 
                                                regionDefinition=regionDefinition,
                                                codonTable=codonTable)
    
    typesOfMutations <- c("R", "S")
    mutationalPaths[!(mutationalPaths %in% typesOfMutations)] <- NA
    
    if (is.null(regionDefinition)) {
        rdLength <- max(stri_length(inputSeq), stri_length(germlineSeq), na.rm=TRUE)
        regionDefinition <- makeNullRegionDefinition(rdLength)
    }
    listExpectedMutationFrequencies <- list()
    for(region in regionDefinition@regions){
        for(typeOfMutation in typesOfMutations){
            region_mutation <- paste(region, typeOfMutation, sep="_")    
            
            targeting_region <- targeting[1:4, regionDefinition@boundaries %in% region]
            mutationalPaths_region <- mutationalPaths[, regionDefinition@boundaries[1:ncol(mutationalPaths)] %in% region]
            targeting_typeOfMutation_region <- sum(targeting_region[mutationalPaths_region == typeOfMutation], na.rm=TRUE)
            
            listExpectedMutationFrequencies[[region_mutation]] <- targeting_typeOfMutation_region
            
        }
    }
    expectedMutationFrequencies <- unlist(listExpectedMutationFrequencies)
    expectedMutationFrequencies[!is.finite(expectedMutationFrequencies)] <- NA
    expectedMutationFrequencies <- expectedMutationFrequencies / sum(expectedMutationFrequencies, na.rm=TRUE)
    return(expectedMutationFrequencies)    
}


calculateTargeting <- function(germlineSeq,
                               inputSeq=NULL,
                               targetingModel=HH_S5F,
                               regionDefinition=NULL) {
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # Check targeting model
    if (!is(targetingModel, "TargetingModel")) {
        stop(deparse(substitute(targetingModel)), " is not a valid TargetingModel object")
    }
    
    # If an inputSequence is passed then process the germlineSequence
    # to be the same legth, mask germlineSequence with Ns where inputSequence is also N
    # If not needed then  you may skip this step by passing in inputSequence=NULL 
    # (which is default). 
    if(!is.null(inputSeq)){    
        # Trim the input and germline sequence to the shortest
        len_inputSeq <- stri_length(inputSeq)
        len_germlineSeq <- stri_length(germlineSeq)
        # If a regionDefinition is passed,
        # then only analyze till the end of the defined length
        if(!is.null(regionDefinition)){
            length_regionDefinition  <- regionDefinition@seqLength
        } else{
            length_regionDefinition <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        }
        len_shortest <- min( c(len_inputSeq,len_germlineSeq,length_regionDefinition),  na.rm=TRUE)
        
        c_inputSeq <- s2c(inputSeq)[1:len_shortest]
        c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
        
        # If the sequence and germline (which now should be the same length) is shorter
        # than the length_regionDefinition, pad it with Ns
        if(len_shortest < length_regionDefinition){
            fillWithNs <- array("N", length_regionDefinition - len_shortest)
            c_inputSeq <- c(c_inputSeq, fillWithNs)
            c_germlineSeq <- c( c_germlineSeq, fillWithNs)
        }
        
        # Mask germline with Ns where input sequence has Ns
        c_germlineSeq[c_inputSeq == "N" |  !c_inputSeq %in% c(NUCLEOTIDES[1:5], ".") ] <- "N"    
        s_germlineSeq <- c2s(c_germlineSeq)
    } else {
        s_germlineSeq <- germlineSeq
    }
    
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to ZZZ any other . is not an IMGT gap & will be treated like N
    gaplessSeq <- gsub("\\.\\.\\.", "ZZZ", s_germlineSeq)
    #If there is a single gap left convert it to an N
    gaplessSeq <- gsub("\\.", "N", gaplessSeq)
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    gaplessSeq <- gsub("ZZZ", "...", gaplessSeq)

    # Vector of seq    
    c_germlineSeq <- s2c(gaplessSeq)
    # Matrix to hold targeting values for each position in c_germlineSeq
    germlineSeqTargeting <- matrix(NA, 
                                   ncol=stri_length(gaplessSeq), 
                                   nrow=length(NUCLEOTIDES[1:5]),
                                   dimnames=list(NUCLEOTIDES[1:5], c_germlineSeq))
    
    # Now remove the IMGT gaps so that the correct 5mers can be made to calculate
    # targeting. e.g.
    # GAGAAA......TAG yields: "GAGAA" "AGAAA" "GAAAT" "AAATA" "AATAG"
    # (because the IMGT gaps are NOT real gaps in sequence!!!)
    gaplessSeq <- gsub("\\.\\.\\.", "", gaplessSeq)
    gaplessSeqLen <- stri_length(gaplessSeq)
    
    #Slide through 5-mers and look up targeting
    gaplessSeq <- paste("NN", gaplessSeq, "NN", sep="")
    gaplessSeqLen <- stri_length(gaplessSeq)
    pos <- 3:(gaplessSeqLen - 2)
    subSeq <- substr(rep(gaplessSeq, gaplessSeqLen - 4), (pos - 2), (pos + 2))
    germlineSeqTargeting_gapless <- targetingModel@targeting[, subSeq]
    #     germlineSeqTargeting_gapless <- sapply(subSeq, function(x) { 
    #         targetingModel@targeting[, x] })

    germlineSeqTargeting[, c_germlineSeq != "."] <- germlineSeqTargeting_gapless  
    
    # Set self-mutating targeting values to be NA
    mutatingToSelf <- colnames(germlineSeqTargeting)
    mutatingToSelf[!(mutatingToSelf %in% NUCLEOTIDES[1:5])] <- "N"
    #     # TODO: What's with this <<- business?
    #     # TODO: I think this is assigning NA to all self-mutations, which are already NA
    #     sapply(1:ncol(germlineSeqTargeting), function(pos) { germlineSeqTargeting[mutatingToSelf[pos], pos] <<- NA })
    
    germlineSeqTargeting[!is.finite(germlineSeqTargeting)] <- NA
    return(germlineSeqTargeting)
}

calculateMutationalPaths <- function(germlineSeq,
                                     inputSeq=NULL,
                                     regionDefinition=NULL,
                                     codonTable=NULL) {    
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    # Set codon table if required
    if (is.null(codonTable)) { codonTable <- CODON_TABLE }
    
    # If an inputSequence is passed then process the germlineSequence
    # to be the same length, mask germlineSequence with Ns where inputSequence is also N
    # If this function is being called after running calculateTargeting you may skip
    # this step by passing in inputSequence=NULL (which is default). This way you save
    # some processing time.
    if(!is.null(inputSeq)){    
        # Trim the input and germline sequence to the shortest
        len_inputSeq <- stri_length(inputSeq)
        len_germlineSeq <- stri_length(germlineSeq)
        # If a regionDefinition is passed,
        # then only analyze till the end of the defined length
        if(!is.null(regionDefinition)){
            length_regionDefinition  <- regionDefinition@seqLength
        } else{
            length_regionDefinition <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        }
        len_shortest <- min( c(len_inputSeq,len_germlineSeq,length_regionDefinition),  na.rm=TRUE)
        
        c_inputSeq <- s2c(inputSeq)[1:len_shortest]
        c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
        
        # If the sequence and germline (which now should be the same length) is shorter
        # than the length_regionDefinition, pad it with Ns
        if(len_shortest<length_regionDefinition){
            fillWithNs <- array("N",length_regionDefinition-len_shortest)
            c_inputSeq <- c(c_inputSeq, fillWithNs)
            c_germlineSeq <- c(c_germlineSeq, fillWithNs)
        }
        
        # Mask germline with Ns where input sequence has Ns
        c_germlineSeq[c_inputSeq == "N" |  !c_inputSeq %in% c(NUCLEOTIDES[1:5], ".") ] = "N"    
        s_germlineSeq <- c2s(c_germlineSeq)
    } else {
        s_germlineSeq <- germlineSeq
        c_germlineSeq <- s2c(s_germlineSeq)
    }
    
    
    s_germlineSeq_len <- stri_length(s_germlineSeq)    
    vecCodons <- sapply({1:(s_germlineSeq_len/3)}*3 - 2, function(x) { substr(s_germlineSeq, x, x + 2) })
    vecCodons[!vecCodons %in% colnames(codonTable)] <- "NNN"
    matMutationTypes <- matrix(codonTable[, vecCodons], nrow=4, byrow=F,
                              dimnames=list(NUCLEOTIDES[1:4], c_germlineSeq[ 1:(3 * length(vecCodons))]))
    
    return(matMutationTypes)
}

#### Additional helper functions ####

# Convert one or more nucleotide characters to IUPAC code 
# for incomplete nucleic acid specification
# 
# @param   nucs     a character vector of nucleotides. One or more of 
#                   \code{c("A", "C", "G", "T")}.
# 
# @return  a single character from the IUPAC ambiguous code.
#
nucs2IUPAC <- function(nucs) {
    # input nucleotides must be one of the characters allowed
    legal <- c("A", "C", "G", "T")
    if (sum(! nucs %in% legal)>0) {
        stop("Input nucleotides must be one of A, C, G, or T.")
    }
    
    # sort by alphabetical order (important)
    nucs <- sort(unique(nucs))
    # concatenate
    nucs <- c2s(nucs)
    
    # convert
    return(IUPAC_DNA_2[nucs])
}

# Convert one or more characters including dash and dots to ambiguous characters
# 
# @param   chars     a character vector of nucleotides. One or more of 
#                    \code{c("A", "C", "G", "T", "N", "-", ".")}.
# 
# @return  a single IUPAC character or "-" or "."
#
chars2Ambiguous <- function(chars) {
    # chars must all be unique
    stopifnot(length(unique(chars)) == length(chars))
    
    # input characters must be one of the characters allowed
    legal <- c("A", "C", "G", "T", "N", "-", ".")
    if (sum(! chars %in% legal) > 0) {
        stop("Input characters must be one of A, C, G, T, N, - (dash), or . (dot)")
    }
    
    # if any of A, T, G, C, N appears
    if (any(chars %in% c("A", "C", "G", "T", "N"))) {
        
        # ignore - and .
        idx.dash.dot <- which(chars == "-" | chars == ".")
        if (length(idx.dash.dot)>0) {
            chars <- chars[-idx.dash.dot]
        }
        
        # if only N appears
        if (sum(chars=="N") == length(chars)) {
            return("N")
        } else {
            # otherwise, if there are any of A, T, G, C
            # remove N
            # e.g. AGN would be treated as AG (R)
            # e.g. ATGN would be treated as AGT (D)
            # e.g. ATGCN would be treated as ACGT (N)
            idx.N <- which(chars == "N")
            if (length(idx.N) > 0) {
                chars <- chars[-idx.N]
            } 
            return(nucs2IUPAC(chars))
        }
        
    } else {
        # otherwise, if only one or both of - and . appear(s)    
        # if both - and . appear, return -
        if (sum(chars %in% c("-", ".")) == 2) {
            return("-")
        } else {
            # if only - or . appears, return that
            return(chars)
        }
    }
    
}

# Convert IUPAC incomplete nucleic acid to one or more characters
#
# @param   code       a single IUPAC character.
# @param   excludeN   if \code{TRUE}, do not translate when \code{code} 
#                     is \code{N}. Default is \code{TRUE}.
# @return  a character vector of nucleotides. One or more of 
#          \code{c("A", "C", "G", "T")}.
#
IUPAC2nucs <- function(code, excludeN=TRUE) {
    # input character must be one of IUPAC codes
    if (! code %in% names(IUPAC_DNA) ) {
        stop("Input character must be one of IUPAC DNA codes.")
    }
    
    # convert
    if (code == "N" & excludeN) {
        return(code)
    } else {
        return(IUPAC_DNA[[code]])
    }
}

# Given a nuclotide position, returns the codon number
# e.g. nuc 86  = codon 29
getCodonNumb <- function(nucPos){
    return( ceiling(nucPos/3) )
}

# Given a codon, returns all the nuc positions that make the codon
getCodonNucs <- function(codonNumb){
    getCodonPos(codonNumb*3)
}

# Given a nucleotide postions return the position in the codon
getContextInCodon <- function(nucPos){
    return((nucPos - 1)%%3 + 1 )
}

# Given a nuclotide position, returns the pos of the 3 nucs that made the codon
# e.g. nuc 86 is part of nucs 85,86,87
getCodonPos <- function(nucPos) {
    codonNum <-  (ceiling(nucPos / 3)) * 3
    return ((codonNum - 2):codonNum)
}

# Translate codon to amino acid
translateCodonToAminoAcid <- function(Codon) {
    return (AMINO_ACIDS[Codon])
}

# Given two codons, tells you if the mutation is R or S (based on your definition)
#
# @param   codonFrom         starting codon. IUPAC ambiguous characters are allowed.
# @param   codonTo           ending codon.  IUPAC ambiguous characters are allowed.
# @param   ambiguousMode     whether to consider ambiguous characters as "either or"
#                            or "and" when determining (and counting) the type(s) of 
#                            mutations. Applicable only if \code{codonFrom} and/or 
#                            \code{codonTo} contains ambiguous characters. One of 
#                            \code{c("eitherOr", "and")}. Default is \code{"eitherOr"}.
# @param   aminoAcidClasses  vector of amino acid trait classes.
#                            if NULL then R or S is determined by amino acid identity
# @return  A vector with entries named by mutation type, including "R" (replacement), 
#          "S" (silent), "Stop" (stop) or "na" (input codons are identical or include NA).
#          Each entry indicates the count of its corresponding type of mutation.
#
# @details When there are ambiguous characters in \code{codonFrom} and/or \code{codonTo}:
#          \itemize{
#               \item  If \code{ambiguousMode="eitherOr"}, ambiguous characters will each 
#                      be expanded but only 1 mutation will be recorded. The priority for 
#                      different types of mutations, in decreasing order, is as follows:
#                      no mutation ("na"), replacement mutation ("R"), silent mutation ("S"),
#                      and stop mutation ("Stop"). The returned vector will have exactly one
#                      entry with a count of 1 and 0 in all other entries.
#               \item  If \code{ambiguousMode="and"}, ambiguous characters will each be 
#                      expanded and mutation(s) from each expansion will be recorded. 
#                      Each entry in the returned vector could potentially be greater than 1.
#          }
#
# @examples
# # Without classes
# mutationType("TTT", "TTC")
# mutationType("TTT", "TTA")
# mutationType("TTT", "TGA")
# mutationType("TGG", "TWG")
#
# # With classes
# classes <- HYDROPATHY_MUTATIONS@classes
# mutationType("TTT", "TTC", aminoAcidClasses=classes)
# mutationType("TTT", "TTA", aminoAcidClasses=classes)
# mutationType("TTT", "TCT", aminoAcidClasses=classes)
# mutationType("TTT", "TGA", aminoAcidClasses=classes)
# 
mutationType <- function(codonFrom, codonTo, 
                         ambiguousMode=c("eitherOr", "and"),
                         aminoAcidClasses=NULL) {
    # codonFrom="TTT"; codonTo="TTA"
    # codonFrom="TTT"; codonTo="TGA"
    
    ambiguousMode <- match.arg(ambiguousMode)
    
    # placeholder for tabulation
    tab <- setNames(object=rep(0, 4), nm=c("R", "S", "Stop", "na"))
    
    if (grepl(pattern="[-.]", x=codonFrom) | grepl(pattern="[-.]", x=codonTo)) {
        # "na"
        tab[4] <- 1
    } else {
        codonFrom.all <- EXPANDED_AMBIGUOUS_CODONS[[codonFrom]]
        codonTo.all <- EXPANDED_AMBIGUOUS_CODONS[[codonTo]]
        
        for (cur.codonFrom in codonFrom.all) {
            for (cur.codonTo in codonTo.all) {
                
                # if codons are the same, there is no mutation; count as NA
                if (cur.codonFrom == cur.codonTo) {
                    # "na"
                    tab[4] <- tab[4] + 1
                } else {
                    # Translate codons
                    cur.aaFrom <- translateCodonToAminoAcid(cur.codonFrom)
                    cur.aaTo <- translateCodonToAminoAcid(cur.codonTo)
                    
                    # If any codon is NA then return NA
                    if (any(is.na(c(codonFrom, codonTo, cur.aaFrom, cur.aaTo)))) { 
                        # "na"
                        tab[4] <- tab[4] + 1
                    } else if (any(c(cur.aaFrom, cur.aaTo) == "*")) {
                        # If any amino acid is Stop then return "Stop"
                        tab[3] <- tab[3] + 1
                    } else if (is.null(aminoAcidClasses)) {
                        # Check for exact identity if no amino acid classes are specified
                        mutation <- if (cur.aaFrom == cur.aaTo) { "S" } else { "R" }
                        tab[mutation] <- tab[mutation]+1
                    } else {
                        # Check for amino acid class identity if classes are specified
                        mutation <- if (aminoAcidClasses[cur.aaFrom] == 
                                        aminoAcidClasses[cur.aaTo]) { "S" } else { "R" }
                        tab[mutation] <- tab[mutation]+1
                    }
                }
            }
        }
        
        # if there's ambiguous char in observed or germline
        if ((length(codonFrom.all) > 1) | (length(codonTo.all) > 1)) {
            if (ambiguousMode=="eitherOr") {
                if (tab[4] > 0) { # "na"
                    tab <- setNames(object=c(0, 0, 0, 1), nm=c("R", "S", "Stop", "na"))
                } else if (tab[2] > 0) { # "S"
                    tab <- setNames(object=c(0, 1, 0, 0), nm=c("R", "S", "Stop", "na"))
                } else if (tab[1] > 0) { # "R"
                    tab <- setNames(object=c(1, 0, 0, 0), nm=c("R", "S", "Stop", "na"))
                } else {
                    tab <- setNames(object=c(0, 0, 1, 0), nm=c("R", "S", "Stop", "na"))
                }
                stopifnot(sum(tab) == 1)
            } else {
                stopifnot(sum(tab) >= 1)
            }
        } else {
            # no need to do anything if there isn't ambiguous char in observed or germline
            # there should be only 1 mutation 
            stopifnot(sum(tab) == 1)
        }
    }
    
    return(tab)
}

