# Targeting models

#' @include Shazam.R
#' @include Core.R
NULL

#### Data ####

#' Uniform 5-mer null targeting model.
#'
#' A null 5-mer model of somatic hypermutation targeting where all substitution, mutability
#' and targeting rates are uniformly distributed.
#'
#' @format   A \link{TargetingModel} object.
#' 
#' @seealso  See \link{HH_S5F} and \link{HKL_S5F} for the human 5-mer targeting models; and 
#'           \link{MK_RS5NF} for the mouse 5-mer targeting model.
"U5N"


#' Human heavy chain, silent, 1-mer, functional substitution model.
#'
#' 1-mer substitution model of somatic hypermutation based on analysis of silent mutations
#' in functional heavy chain Ig sequences from Homo sapiens.
#'
#' @format   A 4x4 matrix of nucleotide substitution rates. The rates are normalized,
#'           therefore each row sums up to 1.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#' }
#' 
#' @seealso  See \link{HKL_S1F} for the human light chain 1-mer substitution model and 
#'           \link{MK_RS1NF} for the mouse light chain 1-mer substitution model.
#'  
#' @note     \code{HH_S1F} replaces \code{HS1FDistance} in versions of SHazaM prior to 0.1.5.
"HH_S1F"


#' Human kappa and lambda chain, silent, 1-mer, functional substitution model.
#'
#' 1-mer substitution model of somatic hypermutation based on analysis of silent mutations
#' in functional kappa and lambda light chain Ig sequences from Homo sapiens.
#' 
#' @format   A 4x4 matrix of nucleotide substitution rates. The rates are normalized,
#'           therefore each row sums up to 1. 
#' 
#' @references
#' \enumerate{
#'   \item  Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,
#'   Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
#'   Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
#'   Immunology, 197(9), 3566-3574.
#'  }
#' 
#' @seealso  See \link{HH_S1F} for the human heavy chain 1-mer substitution model and 
#'           \link{MK_RS1NF} for the mouse light chain 1-mer substitution model.
#'           
#' @note     Reported in Table III in Cui et al, 2016.
"HKL_S1F"


#' Mouse kappa chain, replacement and silent, 1-mer, non-functional substitution model.
#'
#' 1-mer substitution model of somatic hypermutation based on analysis of replacement and
#' silent mutations in non-functional kappa light chain Ig sequences from NP-immunized Mus
#' musculus.
#'           
#' @format   A 4x4 matrix of nucleotide substitution rates. The rates are normalized,
#'           therefore each row sums up to 1. 
#'           
#' @references
#' \enumerate{
#'   \item  Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,
#'   Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
#'   Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
#'   Immunology, 197(9), 3566-3574.
#'  }
#' 
#' @seealso  See \link{HH_S1F} for the human heavy chain 1-mer substitution model and 
#'           \link{HKL_S1F} for the human light chain 1-mer substitution model.
#'
#' @note     \code{MK_RS1NF} replaces \code{M1NDistance} from versions of SHazaM prior to 0.1.5.
"MK_RS1NF"


#' Human heavy chain, silent, 5-mer, functional targeting model.
#'
#' 5-mer model of somatic hypermutation targeting based on analysis of silent mutations
#' in functional heavy chain Ig sequences from Homo sapiens.
#'
#' @format   A \link{TargetingModel} object.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'  
#' @seealso  See \link{HH_S1F} for the 1-mer substitution matrix from the same 
#'           publication; \link{HKL_S5F} for the human light chain 5-mer targeting model; 
#'           \link{MK_RS5NF} for the mouse 5-mer targeting model; and \link{U5N} for the 
#'           uniform 5-mer null targeting model.
"HH_S5F"


#' Human kappa and lambda light chain, silent, 5-mer, functional targeting model.
#'
#' 5-mer model of somatic hypermutation targeting based on analysis of silent mutations
#' in functional kappa and lambda light chain Ig sequences from Homo sapiens.
#'
#' @format A \link{TargetingModel} object.
#' 
#' @references
#' \enumerate{
#'   \item  Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,
#'   Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
#'   Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
#'   Immunology, 197(9), 3566-3574.
#'  }
#'  
#' @seealso  See \link{HH_S5F} for the human heavy chain 5-mer targeting model; 
#'           \link{MK_RS5NF} for the mouse kappa light chain 5-mer targeting model; 
#'           and \link{U5N} for the uniform 5-mer null targeting model.
"HKL_S5F"


#' Mouse kappa light chain, replacement and silent, 5-mer, non-functional targeting model.
#'
#' 5-mer model of somatic hypermutation targeting based on analysis of replacement and
#' silent mutations in non-functional kappa light chain Ig sequences from NP-immunized 
#' Mus musculus.
#'
#' @format \link{TargetingModel} object.
#' 
#' @references
#' \enumerate{
#'   \item  Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,
#'   Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
#'   Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
#'   Immunology, 197(9), 3566-3574.
#'  }
#'  
#' @seealso  See \link{MK_RS1NF} for the 1-mer substitution matrix from the same
#'           publication; \link{HH_S5F} for the human heavy chain silent 5-mer 
#'           functional targeting model; \link{HKL_S5F} for the human light chain 
#'           silent 5-mer functional targeting model; and \link{U5N} for the 
#'           uniform 5-mer null targeting model.
"MK_RS5NF"


#### Classes ####

#' S4 class defining a mutability model
#' 
#' \code{MutabilityModel} defines a data structure for the 5-mer motif-based SHM targeting
#' mutability model. 
#' 
#' @slot    .Data          numeric vector containing 5-mer mutability estimates
#' @slot    source         character vector annotating whether the mutability was
#'                         inferred or directly measured.
#' @slot    numMutS        a number indicating the number of silent mutations used for 
#'                         estimating mutability
#' @slot    numMutR        a number indicating the number of replacement mutations used 
#'                         for estimating mutability  
#'                         
#' @name         MutabilityModel-class
#' @rdname       MutabilityModel-class
#' @aliases      MutabilityModel
#' @exportClass  MutabilityModel
MutabilityModel <- setClass("MutabilityModel", 
                            slots=c(source="character",
                                    numMutS="numeric",
                                    numMutR="numeric"),
                            contains="numeric")


#' S4 class defining a targeting matrix
#' 
#' \code{TargetingMatrix} defines a data structure for just the targeting matrix 
#' (as opposed to the entire \code{TargetingModel})
#' 
#' @slot    .Data          matrix.
#' @slot    numMutS        number indicating the number of silent mutations used for 
#'                         estimating mutability.
#' @slot    numMutR        number indicating the number of replacement mutations used 
#'                         for estimating mutability. 
#'                         
#' @name         TargetingMatrix-class
#' @rdname       TargetingMatrix-class
#' @aliases      TargetingMatrix
#' @exportClass  TargetingMatrix
TargetingMatrix <- setClass("TargetingMatrix",
                            slots=c(numMutS="numeric",
                                    numMutR="numeric"),
                            contains="matrix")


#' S4 class defining a targeting model
#' 
#' \code{TargetingModel} defines a common data structure for mutability, substitution and
#' targeting of immunoglobulin (Ig) sequencing data in a 5-mer microsequence context.
#' 
#' @slot     name          Name of the model.
#' @slot     description   Description of the model and its source data.
#' @slot     species       Genus and species of the source sequencing data.
#' @slot     date          Date the model was built.
#' @slot     citation      Publication source.
#' @slot     substitution  Normalized rates of the center nucleotide of a given 5-mer 
#'                         mutating to a different nucleotide. The substitution model 
#'                         is stored as a 5x3125 matrix of rates. Rows define
#'                         the mutated nucleotide at the center of each 5-mer, one of 
#'                         \code{c("A", "C", "G", "T", "N")}, and columns define the 
#'                         complete 5-mer of the unmutated nucleotide sequence.
#' @slot     mutability    Normalized rates of a given 5-mer being mutated. The 
#'                         mutability model is stored as a numeric vector of length 3125 
#'                         with mutability rates for each 5-mer. Note that "normalized" 
#'                         means that the mutability rates for the 1024 5-mers that 
#'                         contain no "N" at any position sums up to 1 (as opposed to 
#'                         the entire vector summing up to 1).
#' @slot     targeting     Rate matrix of a given mutation ocurring, defined as 
#'                         \eqn{mutability * substitution}. The targeting model 
#'                         is stored as a 5x3125 matrix. Rows define
#'                         the mutated nucleotide at the center of each 5-mer, one of 
#'                         \code{c("A", "C", "G", "T", "N")}, and columns define the complete 5-mer 
#'                         of the unmutated nucleotide sequence.
#' @slot    numMutS        number indicating the number of silent mutations used for 
#'                         estimating mutability.
#' @slot    numMutR        number indicating the number of replacement mutations used 
#'                         for estimating mutability.
#' 
#' @seealso  See \link{createTargetingModel} building models from sequencing data.
#'           
#' @name         TargetingModel-class
#' @rdname       TargetingModel-class
#' @aliases      TargetingModel
#' @exportClass  TargetingModel
setClass("TargetingModel", 
         slots=c(name="character",
                 description="character",
                 species="character",
                 date="character",
                 citation="character",
                 mutability="numeric",
                 substitution="matrix",
                 targeting="matrix",
                 numMutS="numeric",
                 numMutR="numeric"),
         prototype=list(name="name",
                        description="description",
                        species="species",
                        date="2000-01-01",
                        citation="citation",
                        mutability=numeric(3125),
                        substitution=matrix(0, 5, 3125),
                        targeting=matrix(0, 5, 3125),
                        numMutS=as.numeric(NA),
                        numMutR=as.numeric(NA)))


#### Methods ####

#' @param    x    \code{MutabilityModel} object.
#' 
#' @rdname   MutabilityModel-class
#' @aliases  MutabilityModel-method
#' @export
setMethod("print", c(x="MutabilityModel"),
          function(x) { vec <- x@.Data; names(vec) <- names(x); print(vec) })

#' @param    x    \code{MutabilityModel} object.
#' 
#' @rdname   MutabilityModel-class
#' @aliases  MutabilityModel-method
#' @export
setMethod("as.data.frame", c(x="MutabilityModel"),
          function(x) { data.frame(motif=names(x), mutability=x, source=x@source[names(x)]) })

#' @param    x    \code{TargetingModel} object.
#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotMutability}.
#' 
#' @rdname   TargetingModel-class
#' @aliases  TargetingModel-method
#' @export
setMethod("plot", c(x="TargetingModel", y="missing"),
          function(x, y, ...) { plotMutability(x, ...) })


#### Model building functions #####

#' Builds a substitution model
#'
#' \code{createSubstitutionMatrix} builds a 5-mer nucleotide substitution model by counting 
#' the number of substitution mutations occuring in the center position for all 5-mer 
#' motifs.
#'
#' @param    db                data.frame containing sequence data.
#' @param    model             type of model to create. The default model, "s", 
#'                             builds a model by counting only silent mutations. \code{model="s"}
#'                             should be used for data that includes functional sequences.
#'                             Setting \code{model="rs"} creates a model by counting both 
#'                             replacement and silent mutations and may be used on fully 
#'                             non-functional sequence data sets.
#' @param    sequenceColumn    name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn    name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn       name of the column containing the V-segment allele call.
#' @param    multipleMutation  string specifying how to handle multiple mutations occuring 
#'                             within the same 5-mer. If \code{"independent"} then multiple 
#'                             mutations within the same 5-mer are counted indepedently. 
#'                             If \code{"ignore"} then 5-mers with multiple mutations are 
#'                             excluded from the total mutation tally.
#' @param    returnModel       string specifying what type of model to return; one of
#'                             \code{c("5mer", "1mer", "1mer_raw")}. If \code{"5mer"} 
#'                             (the default) then a 5-mer nucleotide context model is 
#'                             returned. If \code{"1mer"} or \code{"1mer_raw"} then a single 
#'                             nucleotide substitution matrix (no context) is returned;
#'                             where \code{"1mer_raw"} is the unnormalized version of the 
#'                             \code{"1mer"} model. Note, neither 1-mer model may be used
#'                             as input to \link{createMutabilityMatrix}.
#' @param    minNumMutations   minimum number of mutations required to compute the 5-mer 
#'                             substitution rates. If the number of mutations for a 5-mer
#'                             is below this threshold, its substitution rates will be 
#'                             estimated from neighboring 5-mers. Default is 50. 
#'                             Not required if \code{numMutationsOnly=TRUE}. 
#' @param    numMutationsOnly  when \code{TRUE}, return counting information on the number
#'                             of mutations for each 5-mer, instead of building a substitution
#'                             matrix. This option can be used for parameter tuning for 
#'                             \code{minNumMutations} during preliminary analysis. 
#'                             Default is \code{FALSE}. Only applies when \code{returnModel} 
#'                             is set to \code{"5mer"}. The \code{data.frame} returned when
#'                             this argument is \code{TRUE} can serve as the input for
#'                             \link{minNumMutationsTune}.                                                       
#' 
#' @return   For \code{returnModel = "5mer"}: 
#' 
#'           When \code{numMutationsOnly} is \code{FALSE}, a 4x1024 matrix of column 
#'           normalized substitution rates for each 5-mer motif with row names defining 
#'           the center nucleotide, one of \code{c("A", "C", "G", "T")}, and column names 
#'           defining the 5-mer nucleotide sequence. 
#'           
#'           When \code{numMutationsOnly} is 
#'           \code{TRUE}, a 1024x4 data frame with each row providing information on 
#'           counting the number of mutations for a 5-mer. Columns are named 
#'           \code{fivemer.total}, \code{fivemer.every}, \code{inner3.total}, and
#'           \code{inner3.every}, corresponding to, respectively,
#'           the total number of mutations when counted as a 5-mer, 
#'           whether there is mutation to every other base when counted as a 5-mer,
#'           the total number of mutations when counted as an inner 3-mer, and
#'           whether there is mutation to every other base when counted as an inner 3-mer.
#'           
#'           For \code{returnModel = "1mer"} or \code{"1mer_raw"}:
#'           a 4x4 normalized or un-normalized 1-mer substitution matrix respectively.
#' 
#' @details  \strong{Caution: The targeting model functions do NOT support ambiguous 
#'           characters in their inputs. You MUST make sure that your input and germline
#'           sequences do NOT contain ambiguous characters (especially if they are
#'           clonal consensuses returned from \code{collapseClones}).}
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'
#' @seealso  \link{extendSubstitutionMatrix}, \link{createMutabilityMatrix}, 
#'           \link{createTargetingMatrix}, \link{createTargetingModel},
#'           \link{minNumMutationsTune}.
#' 
#' @examples
#' \donttest{
#' # Subset example data to one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")[1:25,]
#'
#' # Count the number of mutations per 5-mer
#' subCount <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
#'                                      germlineColumn="germline_alignment_d_mask",
#'                                      vCallColumn="v_call",
#'                                      model="s", multipleMutation="independent",
#'                                      returnModel="5mer", numMutationsOnly=TRUE)
#'
#' # Create model using only silent mutations
#' sub <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
#'                                 germlineColumn="germline_alignment_d_mask",
#'                                 vCallColumn="v_call",
#'                                 model="s", multipleMutation="independent",
#'                                 returnModel="5mer", numMutationsOnly=FALSE,
#'                                 minNumMutations=20)
#' }
#' 
#' @export
createSubstitutionMatrix <- function(db, model=c("s", "rs"), 
                                     sequenceColumn="sequence_alignment",
                                     germlineColumn="germline_alignment_d_mask",
                                     vCallColumn="v_call",
                                     multipleMutation=c("independent", "ignore"),
                                     returnModel=c("5mer", "1mer", "1mer_raw"),
                                     minNumMutations=50,
                                     numMutationsOnly=FALSE)  {
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    returnModel <- match.arg(returnModel)
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn, vCallColumn))
    if (check != TRUE) { stop(check) }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # Check validity of input sequences 
    # (MUST NOT CONTAIN AMBIGUOUS CHARACTERS -- not supported)
    bool_obsv <- checkAmbiguousExist(db[[sequenceColumn]])
    bool_germ <- checkAmbiguousExist(db[[germlineColumn]])
    if (any(bool_obsv | bool_germ)) {
        stop("Ambiguous characters are not supported in input sequences.")
    }
    
    # Setup
    nuc_chars <- NUCLEOTIDES[1:4]
    nuc_words <- seqinr::words(4, nuc_chars)
    
    # Define v_families (heavy or light chain) to only those found in the data
    v_families <- getFamily(db[[vCallColumn]])
    
    # Define empty return list of lists
    substitutionMatrix <- matrix(0, ncol=4, nrow=4, dimnames=list(nuc_chars, nuc_chars))
    substitutionList <- list()    
    for(v_fam in unique(v_families)) {
        substitutionList[[v_fam]] <- list()
        for(word in nuc_words){
            substitutionList[[v_fam]][[word]] <- substitutionMatrix
        }
    }
    
    # Remove IMGT gaps in the germline & input sequences
    matInputCollapsed <- removeCodonGaps(db[, c(sequenceColumn, germlineColumn)])
    # TODO: Unnecessary conversion
    db[[sequenceColumn]] <- matInputCollapsed[, 1]
    db[[germlineColumn]] <- matInputCollapsed[, 2]
    
    # Get mutations
    mutations <- listObservedMutations(db, sequenceColumn=sequenceColumn, 
                                       germlineColumn=germlineColumn,
                                       multipleMutation=multipleMutation,
                                       model=model)
    
    if (model == "s") { # Silent model
        for(index in 1:length(mutations)) {
            cSeq <-  s2c(db[[sequenceColumn]][index])
            cGL  <-  s2c(db[[germlineColumn]][index])
            indexMutation <- mutations[[index]]
            v_fam <- v_families[index]
            
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions<=VLENGTH]
            positions <- positions[!is.na(positions)]
            for( position in  positions){
                wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
                codonNucs <- getCodonPos(position)
                codonGL <- cGL[codonNucs]
                codonSeq <- cSeq[codonNucs]
                muCodonPos <- {position-1}%%3+1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if (!any(codonGL=="N") & !any(codonSeq=="N")) {
                    codonPermutate <- matrix(rep(codonGL,3),ncol=3,byrow=T)
                    codonPermutate[,muCodonPos] <- canMutateTo(glAtMutation)[-4]
                    codonPermutate <- apply(codonPermutate,1,paste,collapse="")
                    codonPermutate <- matrix( c( codonPermutate, rep(c2s(codonGL),3) ), ncol=2, byrow=F)
                    # not intended to be used where input sequences have 
                    # ambiguous characters; it assumes that only 1 entry (r/s/stop/na) from
                    # mutationType is non-zero/1
                    muType <- mutationTypeOptimized(codonPermutate)
                    if (!length(grep("N",wrd))) {
                        if (sum(muType=="s") == length(muType) ){
                            substitutionList[[v_fam]][[wrd]][glAtMutation,seqAtMutation] <- (substitutionList[[v_fam]][[wrd]][glAtMutation,seqAtMutation] + 1)
                        }
                    }
                }
            }
        }
    } else if (model == "rs") { # RS model (All mutations)
        for (index in 1:length(mutations)) {
            cSeq <-  s2c(db[[sequenceColumn]][index])
            cGL  <-  s2c(db[[germlineColumn]][index])
            indexMutation <- mutations[[index]]
            v_fam <- v_families[index]
            
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions<=VLENGTH]
            positions <- positions[!is.na(positions)]
            for( position in  positions){
                wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
                codonNucs <- getCodonPos(position)
                codonGL <- cGL[codonNucs]
                codonSeq <- cSeq[codonNucs]
                muCodonPos <- {position-1}%%3+1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if( !any(codonGL=="N") & !any(codonSeq=="N") ){
                    if(!length(grep("N",wrd))){
                        substitutionList[[v_fam]][[wrd]][glAtMutation,seqAtMutation] <- substitutionList[[v_fam]][[wrd]][glAtMutation, seqAtMutation] + 1
                    }
                }
            }
        }
    }
    
    
    # Convert substitutionList to listSubstitution to facilitate the aggregation of mutations
    arrNames <- c(outer(unique(v_families), nuc_words, paste, sep = "_"))
    listSubstitution <- array(0, dim=c(length(arrNames), 4, 4), 
                              dimnames=list(arrNames, nuc_chars, nuc_chars))
    
    for(v_fam in unique(v_families)){
        listSubstitution[paste(v_fam, nuc_words, sep="_"), , ] <- t(sapply(nuc_words, 
                                                                           function(word) { substitutionList[[v_fam]][[word]] }))
    }
    
    # Aggregate mutations from all V families
    M <- list()
    subMat1mer <- matrix(0, 4, 4) # a single substitution matrix for all fivemers
    listSubNames <- sapply(dimnames(listSubstitution)[[1]], 
                           function(x) { strsplit(x, "_", fixed=TRUE)[[1]] })
    
    .sumSub <- function(i, n) {
        x <- listSubstitution[listSubNames[2, ] == n, i, ]
        if(is.null(dim(x))) {
            return (x)
        } else {
            return (colSums(x))
        }
    }
    for (nuc_word in nuc_words) {
        # Sums mutations from all families
        M[[nuc_word]] <- t(sapply(1:4, .sumSub, n=nuc_word))
        rownames(M[[nuc_word]]) <- nuc_chars
        subMat1mer <- subMat1mer + M[[nuc_word]]
    }
    
    # Return 1-mer substitution model; this output cannot be used for createMutabilityMatrix
    if (returnModel == "1mer") {
        subMat1merNorm <- t(apply(subMat1mer, 1, function(x){x/sum(x)}))
        return (subMat1merNorm)
    } else if (returnModel == "1mer_raw") {
        return (subMat1mer)
    } 
    
    ##### for a given 5mer, count number of mutations
    # fivemer=M; FIVEMER="CCATT"
    .simplifivemer <- function(fivemer, FIVEMER, Thresh=50, count=F) {
      # center
      Nuc=substr(FIVEMER,3,3)
      # neighbors
      Nei=paste(substr(FIVEMER,1,2),substr(FIVEMER,4,5),collapse="",sep="")
      
      ### using 5mer
      # aggregate mutations
      FIVE.5 <- fivemer[[Nei]][Nuc,]
      # count total number of mutations for a given 5mer
      fivemer.total <- sum(FIVE.5)
      # are there mutations to every other base?
      fivemer.every <- ( sum(FIVE.5==0)==1 )
      
      ### using inner 3mer
      # aggregate mutations from 5-mers with the same inner 3-mer
      FIVE.3 <- FIVE.5
      for(i in 1:4){
        for(j in 1:4){
          MutatedNeighbor=paste(nuc_chars[i],substring(Nei,2,3),nuc_chars[j],collapse="",sep="")
          FIVE.3=FIVE.3+fivemer[[MutatedNeighbor]][Nuc,]
        }
      }
      # count total number of mutations for inner 3mer
      inner3.total <- sum(FIVE.3)
      # are there mutations to every other base?
      inner3.every <- ( sum(FIVE.3==0)==1 )
      
      ### using 1mer
      FIVE.1 <- FIVE.5
      MutatedNeighbors <- seqinr::words(4, nuc_chars)
      for (MutatedNeighbor in MutatedNeighbors) {
        FIVE.1=FIVE.1+fivemer[[MutatedNeighbor]][Nuc,]
      }
      
      if (!count) {
        # For a 5mer, if the total number of mutations is greater than Thresh, 
        # and if there are mutations to every other base, compute for the 5mer 
        if ( fivemer.total > Thresh & fivemer.every ){
            return(FIVE.5)
        } else if ( inner3.total > Thresh & inner3.every ) { 
        # Otherwise aggregate mutations from 5-mers with the same inner 3-mer
            return(FIVE.3)
        } else {
        # If the total number of mutations is still not enough, 
        # aggregate mutations from all 5-mers (i.e., use 1-mer model)
            return(FIVE.1)
        } 
      }
      
      return(data.frame(fivemer.total, fivemer.every,
                        inner3.total, inner3.every, 
                        stringsAsFactors=F))
    }
    
    # either construct 5mer substition matrix, and normalize (numMutationsOnly = F)
    if (!numMutationsOnly) {
      substitutionModel <- sapply(seqinr::words(5, nuc_chars), 
                                  function(x) { .simplifivemer(M, x, 
                                                               Thresh = minNumMutations,
                                                               count = numMutationsOnly) }, simplify=T)
      # Assign A->A, C->C, G->G, T->T to NA
      center_nuc <- gsub("..([ACGT])..", "\\1", colnames(substitutionModel))
      for (i in 1:length(center_nuc)) {
        substitutionModel[center_nuc[i], i] <- NA
      }
      
      # Normalize by column
      substitutionModel <- apply(substitutionModel, 2, function(x) { x / sum(x, na.rm=TRUE) })
      substitutionModel[!is.finite(substitutionModel)] <- NA    
    } else {
      # or count number of mutations (numMutationsOnly = T), return data frame
      # need to set simplify to F in sapply() and then use bind_rows; otherwise 
      # every entry in df would be a list
      substitutionModel <- sapply(seqinr::words(5, nuc_chars), 
                                  function(x) { .simplifivemer(M, x, 
                                                               Thresh = minNumMutations,
                                                               count = numMutationsOnly) }, simplify=F)
      substitutionModel <- dplyr::bind_rows(substitutionModel)
      rownames(substitutionModel) <- seqinr::words(5, nuc_chars)
    }
    
    return(substitutionModel)
}

#' Parameter tuning for minNumMutations
#' 
#' \code{minNumMutationsTune} helps with picking a threshold value for \code{minNumMutations}
#' in \link{createSubstitutionMatrix} by tabulating the number of 5-mers for which 
#' substitution rates would be computed directly or inferred at various threshold values.
#' 
#' @param subCount                 \code{data.frame} returned by \link{createSubstitutionMatrix}
#'                                 with \code{numMutationsOnly=TRUE}.
#' @param minNumMutationsRange     a number or a vector indicating the value or range of values
#'                                 of \code{minNumMutations} to try.
#' 
#' @return      A 3xn \code{matrix}, where n is the number of trial values of \code{minNumMutations}
#'              supplied in \code{minNumMutationsRange}. Each column corresponds to a value
#'              in \code{minNumMutationsRange}. The rows correspond to the number of 5-mers
#'              for which substitution rates would be computed directly using the 5-mer itself 
#'              (\code{"5mer"}), using its inner 3-mer (\code{"3mer"}), and using the central 
#'              1-mer (\code{"1mer"}), respectively.
#' 
#' @details     At a given threshold value of \code{minNumMutations}, for a given 5-mer,
#'              if the total number of mutations is greater than the threshold and there
#'              are mutations to every other base, substitution rates are computed directly
#'              for the 5-mer using its mutations. Otherwise, mutations from 5-mers with 
#'              the same inner 3-mer as the 5-mer of interest are aggregated. If the number 
#'              of such mutations is greater than the threshold and there are mutations to 
#'              every other base, these mutations are used for inferring the substitution 
#'              rates for the 5-mer of interest; if not, mutations from all 5-mers with the 
#'              same center nucleotide are aggregated and used for inferring the substitution
#'              rates for the 5-mer of interest (i.e. the 1-mer model).
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso     See argument \code{numMutationsOnly} in \link{createSubstitutionMatrix} 
#'              for generating the required input \code{data.frame} \code{subCount}. 
#'              See argument \code{minNumMutations} in \link{createSubstitutionMatrix}
#'              for what it does.  
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")
#'
#' # Count the number of mutations per 5-mer
#' subCount <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
#'                                      germlineColumn="germline_alignment_d_mask",
#'                                      vCallColumn="v_call",
#'                                      model="s", multipleMutation="independent",
#'                                      returnModel="5mer", numMutationsOnly=TRUE)
#' 
#' # Tune minNumMutations
#' minNumMutationsTune(subCount, seq(from=10, to=80, by=10))
#'                                       
#' @export
minNumMutationsTune <- function(subCount, minNumMutationsRange) {
  stopifnot( nrow(subCount)==1024 & ncol(subCount)==4 )

  tuneMtx <- sapply(minNumMutationsRange, 
                     function(thresh) {
                       method.count <- c(# as 5mer
                                        sum( subCount$fivemer.total > thresh & 
                                             subCount$fivemer.every ),
                                        # as inner 3mer
                                        sum( !(subCount$fivemer.total > thresh & 
                                               subCount$fivemer.every) &
                                             (subCount$inner3.total > thresh & 
                                              subCount$inner3.every) ),
                                        # as 1mer
                                        sum( !(subCount$fivemer.total > thresh & 
                                               subCount$fivemer.every) & 
                                             !(subCount$inner3.total > thresh & 
                                               subCount$inner3.every) )
                                       )
                       names(method.count) <- c("5mer", "3mer", "1mer")
                       stopifnot( sum(method.count)==1024 )
                       return(method.count)
                     })
  colnames(tuneMtx) <- minNumMutationsRange
  return(tuneMtx)
}

#' Builds a mutability model
#'
#' \code{createMutabilityMatrix} builds a 5-mer nucleotide mutability model by counting 
#' the number of mutations occuring in the center position for all 5-mer motifs.
#'
#' @param    db                  data.frame containing sequence data.
#' @param    substitutionModel   matrix of 5-mer substitution rates built by 
#'                               \link{createSubstitutionMatrix}. Note, this model will
#'                               only impact mutability scores when \code{model="s"}
#'                               (using only silent mutations).
#' @param    model               type of model to create. The default model, "s", 
#'                               builds a model by counting only silent mutations. \code{model="s"}
#'                               should be used for data that includes functional sequences.
#'                               Setting \code{model="rs"} creates a model by counting both 
#'                               replacement and silent mutations and may be used on fully 
#'                               non-functional sequence data sets.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn         name of the column containing the V-segment allele call.
#' @param    multipleMutation    string specifying how to handle multiple mutations occuring 
#'                               within the same 5-mer. If \code{"independent"} then multiple 
#'                               mutations within the same 5-mer are counted indepedently. 
#'                               If \code{"ignore"} then 5-mers with multiple mutations are 
#'                               excluded from the total mutation tally.
#' @param    minNumSeqMutations  minimum number of mutations in sequences containing each 5-mer
#'                               to compute the mutability rates. If the number is smaller 
#'                               than this threshold, the mutability for the 5-mer will be 
#'                               inferred. Default is 500. Not required if 
#'                               \code{numSeqMutationsOnly=TRUE}.   
#' @param    numSeqMutationsOnly when \code{TRUE}, return only a vector counting the number of 
#'                               observed mutations in sequences containing each 5-mer. This 
#'                               option can be used for parameter tuning for \code{minNumSeqMutations} 
#'                               during preliminary analysis using \link{minNumSeqMutationsTune}. 
#'                               Default is \code{FALSE}.                              
#'
#' @return   When \code{numSeqMutationsOnly} is \code{FALSE}, a \code{MutabilityModel} containing a
#'           named numeric vector of 1024 normalized mutability rates for each 5-mer motif with names 
#'           defining the 5-mer nucleotide sequence.
#'           
#'           When \code{numSeqMutationsOnly} is \code{TRUE}, a named numeric
#'           vector of length 1024 counting the number of observed mutations in sequences containing 
#'           each 5-mer.
#' 
#' @details  \strong{Caution: The targeting model functions do NOT support ambiguous 
#'           characters in their inputs. You MUST make sure that your input and germline
#'           sequences do NOT contain ambiguous characters (especially if they are 
#'           clonal consensuses returned from \code{collapseClones}).}
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso  \link{MutabilityModel}, \link{extendMutabilityMatrix}, \link{createSubstitutionMatrix}, 
#'           \link{createTargetingMatrix}, \link{createTargetingModel},
#'           \link{minNumSeqMutationsTune}
#' 
#' @examples
#' \donttest{
#' # Subset example data to 50 sequences of one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")[1:50,]
#'
#' # Create model using only silent mutations
#' sub_model <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
#'                                       germlineColumn="germline_alignment_d_mask",
#'                                       vCallColumn="v_call",model="s")
#' mut_model <- createMutabilityMatrix(db, sub_model, model="s", 
#'                                     sequenceColumn="sequence_alignment",
#'                                     germlineColumn="germline_alignment_d_mask",
#'                                     vCallColumn="v_call",
#'                                     minNumSeqMutations=200,
#'                                     numSeqMutationsOnly=FALSE)
#'                                     
#' # View top 5 mutability estimates
#' head(sort(mut_model, decreasing=TRUE), 5)
#' 
#' # View the number of S mutations used for estimating mutabilities
#' mut_model@numMutS
#' 
#' # Count the number of mutations in sequences containing each 5-mer
#' mut_count <- createMutabilityMatrix(db, sub_model, model="s", 
#'                                     sequenceColumn="sequence_alignment",
#'                                     germlineColumn="germline_alignment_d_mask",
#'                                     vCallColumn="v_call",
#'                                     numSeqMutationsOnly=TRUE)
#' }
#' 
#' @export
createMutabilityMatrix <- function(db, substitutionModel, model=c("s", "rs"),
                                   sequenceColumn="sequence_alignment", 
                                   germlineColumn="germline_alignment_d_mask",
                                   vCallColumn="v_call",
                                   multipleMutation=c("independent", "ignore"),
                                   minNumSeqMutations=500, 
                                   numSeqMutationsOnly=FALSE) {
    # substitutionModel=sub_model; model="s"; sequenceColumn="sequence_alignment"; germlineColumn="germline_alignment_d_mask"
    # vCallColumn="v_call"; multipleMutation="ignore"; minNumSeqMutations=10
    
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn, vCallColumn))
    if (check != TRUE) { stop(check) }
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    
    # Check validity of input sequences 
    # (MUST NOT CONTAIN AMBIGUOUS CHARACTERS -- not supported)
    bool_obsv <- checkAmbiguousExist(db[[sequenceColumn]])
    bool_germ <- checkAmbiguousExist(db[[germlineColumn]])
    if (any(bool_obsv | bool_germ)) {
        stop("Ambiguous characters are not supported in input sequences.")
    }
    
    # Check that the substitution model is valid
    if (any(dim(substitutionModel) != c(4, 1024))) {
        stop ("Please supply a valid 5-mer substitutionModel.")
    }

    # Set constants for function
    nuc_chars <- NUCLEOTIDES[1:4]
    
    # Remove IMGT gaps in the germline & input sequences
    matInputCollapsed <- removeCodonGaps(db[, c(sequenceColumn, germlineColumn)])
    # TODO: Unnecessary conversion
    db[[sequenceColumn]] <- matInputCollapsed[, 1]
    db[[germlineColumn]] <- matInputCollapsed[, 2]
    
    # Count mutations
    # TODO: this could be listMutations() instead, and skip the conversion from matInputCollapsed back to a data.frame
    mutations <- listObservedMutations(db, sequenceColumn=sequenceColumn, 
                                       germlineColumn=germlineColumn,
                                       multipleMutation=multipleMutation,
                                       model=model)
    
    # Foreground Count: Count the number of observed mutations for each 5-mer
    template <- rep(0, 1024)
    names(template) <- seqinr::words(5, nuc_chars)
    COUNT <- list()
    for(index in 1:length(mutations)){
        COUNT[[index]] <- template
        indexMutation <- mutations[[index]]
        if(!sum(is.na(indexMutation))){
            cSeq <-  s2c(db[[sequenceColumn]][index])
            cGL  <-  s2c(db[[germlineColumn]][index])
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions <= VLENGTH]
            positions <- positions[!is.na(positions)]
            for (position in  positions){
                wrd5 <- substr(db[[germlineColumn]][index], position - 2, position + 2)
                if(!grepl("[^ACGT]", wrd5) & nchar(wrd5) == 5){
                    codonNucs <- getCodonPos(position)
                    codonGL <- cGL[codonNucs]
                    codonSeq <- cSeq[codonNucs]
                    muCodonPos <- {position - 1} %% 3 + 1
                    #seqAtMutation <- codonSeq[muCodonPos]
                    glAtMutation <- codonGL[muCodonPos]
                    if (!any(codonGL %in% c("N", "-", ".")) & !any(codonSeq %in% c("N", "-", "."))) {
                        if (!length(grep("N", wrd5))) {
                            COUNT[[index]][wrd5]<- COUNT[[index]][wrd5] + 1;
                        }
                    }
                }
            }
        }
    }

    # Define sum of rates for nucleotide sets from substitution model
    # Two character sets
    wrd2Index <- combn(1:4, 2)
    wrd2Sums <- t(apply(wrd2Index, 2, function(x) colSums(substitutionModel[x, ], na.rm=TRUE)))
    rownames(wrd2Sums) <- apply(wrd2Index, 2, function(x) paste(nuc_chars[x], collapse=""))
    # Three character sets
    wrd3Index <- combn(1:4, 3)
    wrd3Sums <- t(apply(wrd3Index, 2, function(x) colSums(substitutionModel[x, ], na.rm=TRUE)))
    rownames(wrd3Sums) <- apply(wrd3Index, 2, function(x) paste(nuc_chars[x], collapse=""))
    # Merge single character, two character and three character sets
    substitutionSums <- rbind(substitutionModel, wrd2Sums, wrd3Sums)
        
    # Replace dots with Ns
    sSeqVec <- gsub("\\.", "N", db[[sequenceColumn]])
    sGermVec <- gsub("\\.", "N", db[[germlineColumn]])
    
    # Define template for 5-mer sums by position
    countTemplate <- matrix(0, VLENGTH, 1024, dimnames=list(1:VLENGTH, names(template)))
    
    # Background Count: Count the number of occurrences of each 5-mer
    BG_COUNT <- list()
    for (index in 1:length(mutations)) {
        tmpCounts <- countTemplate
        sGL <- sGermVec[index]
        cSeq <-  s2c(sSeqVec[index])
        cGL  <-  s2c(sGL)[1:VLENGTH]
        positions <- 3:(length(cGL) - 2)
        for (pos in  positions) {
            wrd5 <- substr(sGL, pos - 2, pos + 2)
            if (!grepl("[^ACGT]", wrd5) & nchar(wrd5) == 5) {
                codonNucs <- getCodonPos(pos)
                codonGL <- cGL[codonNucs]
                codonSeq <- cSeq[codonNucs]
                muCodonPos <- (pos - 1) %% 3 + 1
                glAtMutation <- codonGL[muCodonPos]
                if (!any(codonGL %in% c("N", "-")) & !any(codonSeq %in% c("N", "-"))) {
                    # Determine mutation types for NUCLEOTIDES[1:4]
                    muType <- CODON_TABLE[1:4 + 4*(muCodonPos - 1), stri_flatten(codonGL)]

                    # Set characters that meet mutation criteria
                    if (model == "s") {
                        muChars <- nuc_chars[1:4][nuc_chars[1:4] != glAtMutation & muType == "s"]
                    } else { 
                        muChars <- nuc_chars[1:4][nuc_chars[1:4] != glAtMutation]
                    }

                    # Update counts
                    if (length(muChars) > 0) {
                        #cat(stri_flatten(muChars), substitutionSums[stri_flatten(muChars), wrd5], "\n")
                        tmpCounts[pos, wrd5] <- substitutionSums[stri_flatten(muChars), wrd5]
                    }
                }
            }
        }
        BG_COUNT[[index]] <- colSums(tmpCounts)
        BG_COUNT[[index]][BG_COUNT[[index]] == 0] <- NA
    }
    
    Mutability <- list()
    for(i in 1:length(mutations)) {
        mut_mat <- COUNT[[i]] / BG_COUNT[[i]]
        mut_mat <- mut_mat / sum(mut_mat, na.rm=TRUE)
        mut_mat[!is.finite(mut_mat)] <- NA
        wgt_mat <- length(mutations[[i]])
        Mutability[[i]] <- list(mut_mat, wgt_mat)
    }
    
    # total counts of mutations
    # each list item is the total S and R mutation counts in a seq
    mutationsTotalLst <- lapply(mutations, function(m){ 
        return( c(S=sum(m=="s", na.rm=T), R=sum(m=="r", na.rm=T)) ) 
    })
    # total S and R mutation counts across seqs
    mutationsTotalRS <- Reduce("+", mutationsTotalLst)
    
    # Aggregate mutability
    MutabilityMatrix <- sapply(Mutability, function(x) x[[1]])
    MutabilityWeights <- sapply(Mutability, function(x) x[[2]])
    Mutability_Mean <- apply(MutabilityMatrix, 1, weighted.mean, w=MutabilityWeights, na.rm=TRUE)
    Mutability_Mean[!is.finite(Mutability_Mean)] <- NA
    Mutability_Mean[Mutability_Mean == 0] <- NA
    
    # Filter out 5-mers with low number of observed mutations in the sequences
    NumSeqMutations <- sapply(1:1024,function(i)sum(MutabilityWeights[!is.na(MutabilityMatrix[i,])])) 
    names(NumSeqMutations) <- names(Mutability_Mean)
    if (numSeqMutationsOnly) {return(NumSeqMutations)}
    
    Mutability_Mean[NumSeqMutations <= minNumSeqMutations] <- NA
    
    # Infer mutability for missing 5-mers
    .fillHot <-function(FIVEMER,mutability){
        if(FIVEMER%in%names(mutability))if(!is.na(mutability[[FIVEMER]]))if(mutability[[FIVEMER]]>=0.0)return(mutability[[FIVEMER]])
        Nuc=substr(FIVEMER,3,3)
        #Nei=paste(substr(FIVEMER,1,2),substr(FIVEMER,4,5),collapse="",sep="")
        FIVE=0
        COUNT=0
        
        # For A/T, infer mutability using the 3-mer model. 
        if(Nuc%in%c("A","T")){
            for(i in 1:3){
                for(j in 1:3){
                    MutatedNeighbor=paste(canMutateTo(substr(FIVEMER,1,1))[i],substr(FIVEMER,2,4),canMutateTo(substr(FIVEMER,5,5))[j],collapse="",sep="")
                    if(!is.na(mutability[[MutatedNeighbor]])){
                        FIVE=FIVE+mutability[[MutatedNeighbor]]
                        COUNT=COUNT+1
                    }
                }
            }
            return(FIVE/COUNT)
        }
        # For G, infer using 5-mers with the same downstream nucleotides 
        if(Nuc=="G"){
            for(i in 1:3){
                for(j in 1:3){
                    MutatedNeighbor=paste(canMutateTo(substr(FIVEMER,1,1))[i],canMutateTo(substr(FIVEMER,2,2))[j],substr(FIVEMER,3,5),collapse="",sep="")
                    if(!is.na(mutability[[MutatedNeighbor]])){
                        FIVE=FIVE+mutability[[MutatedNeighbor]]
                        COUNT=COUNT+1
                    }
                }
            }
            return(FIVE/COUNT)
        }
        
        # For C, infer using 5-mers with the same upstream nucleotides 
        if(Nuc=="C"){
            for(i in 1:3){
                for(j in 1:3){
                    MutatedNeighbor=paste(substr(FIVEMER,1,3),canMutateTo(substr(FIVEMER,4,4))[i],canMutateTo(substr(FIVEMER,5,5))[j],collapse="",sep="")
                    if(!is.na(mutability[[MutatedNeighbor]])){
                        FIVE=FIVE+mutability[[MutatedNeighbor]]
                        COUNT=COUNT+1
                    }
                }
            }
            return(FIVE/COUNT)
        }
    }
    
    Mutability_Mean_Complete <-sapply(words(5, nuc_chars), .fillHot, mutability = Mutability_Mean)
    
    for(i in names(which(is.na(Mutability_Mean_Complete)))){
        Mutability_Mean_Complete[i]<- .fillHot(i,mutability=Mutability_Mean_Complete)
    }
    for(i in names(which((Mutability_Mean_Complete)<1e-6))){
        Mutability_Mean_Complete[i]<- .fillHot(i,mutability=Mutability_Mean_Complete)
    }
    # If the neighboring 5-mers still don't have enough mutations, use 0 instead. 
    if (length(is.na(Mutability_Mean_Complete)) > 0) {
        warning("Insufficient number of mutations to infer some 5-mers. Filled with 0. ")
        Mutability_Mean_Complete[is.na(Mutability_Mean_Complete)] <- 0 
    }
    
    
    # Normalize
    Mutability_Mean_Complete <- Mutability_Mean_Complete / sum(Mutability_Mean_Complete, na.rm=TRUE)

    # Define whether the 5-mer mutability is measured or inferred
    mut_names <- names(Mutability_Mean_Complete)
    mut_source <- setNames(rep("Measured", length(mut_names)), mut_names)
    mut_source[mut_names %in% names(which(is.na(Mutability_Mean)))] <- "Inferred"
  
    # Return MutabilityModel
    mut_model <- MutabilityModel(Mutability_Mean_Complete,
                                 source=mut_source,
                                 numMutS=mutationsTotalRS[["S"]],
                                 numMutR=mutationsTotalRS[["R"]])
    return(mut_model)
}


#' Parameter tuning for minNumSeqMutations
#' 
#' \code{minNumSeqMutationsTune} helps with picking a threshold value for \code{minNumSeqMutations}
#' in \link{createMutabilityMatrix} by tabulating the number of 5-mers for which 
#' mutability would be computed directly or inferred at various threshold values.
#' 
#' @param mutCount                  a \code{vector} of length 1024 returned by 
#'                                  \link{createMutabilityMatrix} with \code{numSeqMutationsOnly=TRUE}.
#' @param minNumSeqMutationsRange   a number or a vector indicating the value or the range of values 
#'                                  of \code{minNumSeqMutations} to try.
#' 
#' @return      A 2xn \code{matrix}, where n is the number of trial values of \code{minNumSeqMutations}
#'              supplied in \code{minNumSeqMutationsRange}. Each column corresponds to a value
#'              in \code{minNumSeqMutationsRange}. The rows correspond to the number of 5-mers
#'              for which mutability would be computed directly (\code{"measured"}) and inferred
#'              (\code{"inferred"}), respectively.
#'              
#' @details     At a given threshold value of \code{minNumSeqMutations}, for a given 5-mer,
#'              if the total number of mutations is greater than the threshold, mutability 
#'              is computed directly. Otherwise, mutability is inferred.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso     See argument \code{numSeqMutationsOnly} in \link{createMutabilityMatrix} 
#'              for generating the required input \code{vector} \code{mutCount}. 
#'              See argument \code{minNumSeqMutations} in \link{createMutabilityMatrix}
#'              for what it does.  
#' 
#' @examples
#' \donttest{
#' # Subset example data to one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")
#' set.seed(112)
#' db <- dplyr::slice_sample(db, n=75)
#' # Create model using only silent mutations
#' sub <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
#'                                 germlineColumn="germline_alignment_d_mask",
#'                                 vCallColumn="v_call", 
#'                                 model="s", multipleMutation="independent",
#'                                 returnModel="5mer", numMutationsOnly=FALSE,
#'                                 minNumMutations=20)
#'
#' # Count the number of mutations in sequences containing each 5-mer
#' mutCount <- createMutabilityMatrix(db, substitutionModel = sub,
#'                                    sequenceColumn="sequence_alignment",
#'                                    germlineColumn="germline_alignment_d_mask",
#'                                    vCallColumn="v_call",
#'                                    model="s", multipleMutation="independent",
#'                                    numSeqMutationsOnly=TRUE)
#' 
#' # Tune minNumSeqMutations
#' minNumSeqMutationsTune(mutCount, seq(from=100, to=300, by=50))
#' }                                      
#' @export
minNumSeqMutationsTune <- function(mutCount, minNumSeqMutationsRange) {
  stopifnot( length(mutCount) == 1024 )
  
  tuneMtx <- sapply(minNumSeqMutationsRange, 
                     function(thresh) {
                       method.count <- c( sum(mutCount > thresh),
                                         sum(mutCount <= thresh) )
                       names(method.count) <- c("measured", "inferred")
                       stopifnot( sum(method.count)==1024 )
                       return(method.count)
                     })
  colnames(tuneMtx) <- minNumSeqMutationsRange
  return(tuneMtx)
}


#' Extends a substitution model to include Ns.
#' 
#' \code{extendSubstitutionMatrix} extends a 5-mer nucleotide substitution model 
#' with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.
#'
#' @param    substitutionModel  matrix of 5-mers substitution counts built by 
#'                              \link{createSubstitutionMatrix}.
#' 
#' @return   A 5x3125 matrix of normalized substitution rate for each 5-mer motif with 
#'           rows names defining the center nucleotide, one of \code{c("A", "C", "G", "T", "N")}, 
#'           and column names defining the 5-mer nucleotide sequence.
#' 
#' @seealso  \link{createSubstitutionMatrix}, \link{extendMutabilityMatrix}
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")
#'
#' # Create model using only silent mutations
#' sub_model <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
#'                                       germlineColumn="germline_alignment_d_mask",
#'                                       vCallColumn="v_call",model="s")
#' ext_model <- extendSubstitutionMatrix(sub_model)
#' 
#' @export
extendSubstitutionMatrix <- function(substitutionModel) {
    # TODO: fix order so Ns are at the end? (c(input_names, words not in input_names))
    
    # Define old and new column/row names
    input_names <- colnames(substitutionModel)
    nuc_chars <- NUCLEOTIDES[1:5]
    nuc_5mers <- seqinr::words(5, alphabet=nuc_chars)
    
    # Define empty extended matrix with Ns
    extend_mat <- matrix(NA, nrow=length(nuc_chars), ncol=length(nuc_5mers), 
                         dimnames=list(nuc_chars, nuc_5mers))
    
    # Extend matrix with Ns
    for (mer in nuc_5mers) {
        if (mer %in% input_names) {
            extend_mat[, mer] <- c(substitutionModel[, mer], "N"=NA)
        } else {
            mer_char <- s2c(mer)
            n_index <- grep("N", mer_char)
            if (any(n_index == 3)) {
                extend_mat[, mer] <- NA
            } else {
                mer_char[n_index] <- "."
                mer_str <- c2s(mer_char)
                mer_index <- grep(mer_str, input_names)
                extend_mat[, mer] <- c(apply(substitutionModel[, mer_index], 1, mean, na.rm=TRUE), "N"=NA)
            }
        }
    }
    
    # Normalize
    #extend_mat <- apply(extend_mat, 2, function(x) { x/sum(x, na.rm=TRUE) })
    extend_mat[!is.finite(extend_mat)] <- NA
    
    return (extend_mat)
}


#' Extends a mutability model to include Ns.
#' 
#' \code{extendMutabilityMatrix} extends a 5-mer nucleotide mutability model 
#' with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.
#'
#' @param    mutabilityModel  vector of 5-mer mutability rates built by 
#'                            \link{createMutabilityMatrix}.
#' 
#' @return   A \code{MutabilityModel} containing a 3125 vector of normalized 
#'           mutability rates for each 5-mer motif with names defining the 5-mer 
#'           nucleotide sequence. Note that "normalized" means that the mutability 
#'           rates for the 1024 5-mers that contain no "N" at any position sums up 
#'           to 1 (as opposed to the entire vector summing up to 1). 
#'           
#'           If the input \code{mutabilityModel} is of class \code{MutabilityModel}, 
#'           then the output \code{MutabilityModel} will carry over the input 
#'           \code{numMutS} and \code{numMutR} slots.
#' 
#' @seealso  \link{createMutabilityMatrix}, \link{extendSubstitutionMatrix}, 
#'           \link{MutabilityModel}
#' 
#' @examples
#' \donttest{
#' # Subset example data to one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")
#' set.seed(112)
#' db <- dplyr::slice_sample(db, n=75)
#'
#' # Create model using only silent mutations and ignore multiple mutations
#' sub_model <- createSubstitutionMatrix(db, model="s", sequenceColumn="sequence_alignment",
#'                                       germlineColumn="germline_alignment_d_mask",
#'                                       vCallColumn="v_call")
#' mut_model <- createMutabilityMatrix(db, sub_model, model="s", 
#'                                     sequenceColumn="sequence_alignment",
#'                                     germlineColumn="germline_alignment_d_mask",
#'                                     vCallColumn="v_call")
#' ext_model <- extendMutabilityMatrix(mut_model)
#' }
#' 
#' @export
extendMutabilityMatrix <- function(mutabilityModel) {
    # TODO: fix order so Ns are at the end? (c(input_names, words not in input_names))
    
    # Define old and new column/row names
    input_names <- names(mutabilityModel)
    nuc_chars <- NUCLEOTIDES[1:5]
    nuc_5mers <- seqinr::words(5, alphabet=nuc_chars)
    
    # Define empty extended matrix with Ns
    extend_mat <- array(NA, dim=length(nuc_5mers), dimnames=list(nuc_5mers))
    
    # Extend matrix with Ns
    for(mer in nuc_5mers) {
        #cat(mer,"\n")
        if (mer %in% input_names) {
            extend_mat[mer] <- mutabilityModel[mer]
        } else {
            mer_char <- s2c(mer)
            n_index <- grep("N", mer_char)
            if (any(n_index == 3)) {
                extend_mat[mer] <- NA
            } else {
                mer_char[n_index] <- "."
                mer_str <- c2s(mer_char)
                mer_index <- grep(mer_str, input_names)
                extend_mat[mer] <- mean(mutabilityModel[mer_index], na.rm=TRUE)
            }
        }
    }
    
    # Normalize    
    #extend_mat <- extend_mat / sum(extend_mat, na.rm=TRUE)
    extend_mat[!is.finite(extend_mat)] <- NA
    
    # Carry over @numMutS and @numMutR, if any
    if (all(c("numMutS", "numMutR") %in% slotNames(mutabilityModel))) {
        mut_s <- mutabilityModel@numMutS
        mut_r <- mutabilityModel@numMutR
    } else {
        mut_s <- as.numeric(NA) 
        mut_r <- as.numeric(NA)
    }
    
    # Carry over @source
    if ("source" %in% slotNames(mutabilityModel)) {
        mut_names <- names(extend_mat)
        mut_source <- setNames(rep("Extended", length(mut_names)), mut_names)
        mut_source[names(mutabilityModel@source)] <- mutabilityModel@source
    } else {
        mut_source <- as.character(NA)  
    }
    
    # Return extended MutabilityModel
    extend_model <- MutabilityModel(extend_mat,
                                    source=mut_source,
                                    numMutS=mut_s,
                                    numMutR=mut_r)

    return(extend_model)
}
 

#' Calculates a targeting rate matrix
#' 
#' \code{createTargetingMatrix} calculates the targeting model matrix as the
#' combined probability of mutability and substitution.
#'
#' @param    substitutionModel  matrix of 5-mers substitution rates built by 
#'                              \link{createSubstitutionMatrix} or 
#'                              \link{extendSubstitutionMatrix}.
#' @param    mutabilityModel    vector of 5-mers mutability rates built by 
#'                              \link{createMutabilityMatrix} or 
#'                              \link{extendMutabilityMatrix}.
#' 
#' @return   A \code{TargetingMatrix} with the same dimensions as the input \code{substitutionModel} 
#'           containing normalized targeting probabilities for each 5-mer motif with 
#'           row names defining the center nucleotide and column names defining the 
#'           5-mer nucleotide sequence. 
#'           
#'           If the input \code{mutabilityModel} is of class \code{MutabilityModel}, then the output 
#'           \code{TargetingMatrix} will carry over the input \code{numMutS} and \code{numMutR} slots.
#' 
#' @details
#' Targeting rates are calculated by multiplying the normalized mutability rate by the 
#' normalized substitution rates for each individual 5-mer.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data.
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso  \link{createSubstitutionMatrix}, \link{extendSubstitutionMatrix}, 
#'           \link{createMutabilityMatrix}, \link{extendMutabilityMatrix}, 
#'           \link{TargetingMatrix}, \link{createTargetingModel}
#' 
#' @examples
#' \donttest{
#' # Subset example data to 50 sequences, of one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")[1:50,]
#'
#' # Create 4x1024 models using only silent mutations
#' sub_model <- createSubstitutionMatrix(db, model="s", sequenceColumn="sequence_alignment",
#'                                       germlineColumn="germline_alignment_d_mask",
#'                                       vCallColumn="v_call")
#' mut_model <- createMutabilityMatrix(db, sub_model, model="s",
#'                                     sequenceColumn="sequence_alignment",
#'                                     germlineColumn="germline_alignment_d_mask",
#'                                     vCallColumn="v_call")
#' 
#' # Extend substitution and mutability to including Ns (5x3125 model)
#' sub_model <- extendSubstitutionMatrix(sub_model)
#' mut_model <- extendMutabilityMatrix(mut_model)
#' 
#' # Create targeting model from substitution and mutability
#' tar_model <- createTargetingMatrix(sub_model, mut_model)
#' }
#' 
#' @export
createTargetingMatrix <- function(substitutionModel, mutabilityModel) {
    # Calculate targeting
    tar_mat <- sweep(substitutionModel, 2, mutabilityModel, `*`)
    
    # Normalize
    #tar_mat <- tar_mat / sum(tar_mat, na.rm=TRUE)
    tar_mat[!is.finite(tar_mat)] <- NA
    
    # carry over @numMutS and @numMutR, if any
    if (all(c("numMutS", "numMutR") %in% slotNames(mutabilityModel))) {
        tar_mat <- TargetingMatrix(tar_mat,
                                   numMutS=mutabilityModel@numMutS,
                                   numMutR=mutabilityModel@numMutR)
    } else {
        tar_mat <- TargetingMatrix(tar_mat,
                                   # NA is of class logical by default
                                   numMutS=as.numeric(NA), 
                                   numMutR=as.numeric(NA))
    }
    
    return(tar_mat)
}


#' Creates a TargetingModel
#' 
#' \code{createTargetingModel} creates a 5-mer \code{TargetingModel}.
#'
#' @param    db                  data.frame containing sequence data.
#' @param    model               type of model to create. The default model, "s", 
#'                               builds a model by counting only silent mutations. \code{model="s"}
#'                               should be used for data that includes functional sequences.
#'                               Setting \code{model="rs"} creates a model by counting both 
#'                               replacement and silent mutations and may be used on fully 
#'                               non-functional sequence data sets.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn         name of the column containing the V-segment allele calls.
#' @param    multipleMutation    string specifying how to handle multiple mutations occuring 
#'                               within the same 5-mer. If \code{"independent"} then multiple 
#'                               mutations within the same 5-mer are counted indepedently. 
#'                               If \code{"ignore"} then 5-mers with multiple mutations are 
#'                               excluded from the otal mutation tally.
#' @param    minNumMutations     minimum number of mutations required to compute the 5-mer 
#'                               substitution rates. If the number of mutations for a 5-mer
#'                               is below this threshold, its substitution rates will be 
#'                               estimated from neighboring 5-mers. Default is 50.   
#' @param    minNumSeqMutations  minimum number of mutations in sequences containing each 5-mer
#'                               to compute the mutability rates. If the number is smaller 
#'                               than this threshold, the mutability for the 5-mer will be 
#'                               inferred. Default is 500.   
#' @param    modelName           name of the model.
#' @param    modelDescription    description of the model and its source data.
#' @param    modelSpecies        genus and species of the source sequencing data.
#' @param    modelDate           date the model was built. If \code{NULL} the current date
#'                               will be used.
#' @param    modelCitation       publication source.
#' 
#' @return   A \link{TargetingModel} object.
#' 
#' @details  \strong{Caution: The targeting model functions do NOT support ambiguous 
#'           characters in their inputs. You MUST make sure that your input and germline
#'           sequences do NOT contain ambiguous characters (especially if they are
#'           clonal consensuses returned from \code{collapseClones}).}
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data.
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso  See \link{TargetingModel} for the return object. 
#'           See \link{plotMutability} plotting a mutability model.
#'           See \link{createSubstitutionMatrix}, \link{extendSubstitutionMatrix}, 
#'           \link{createMutabilityMatrix}, \link{extendMutabilityMatrix} and 
#'           \link{createTargetingMatrix} for component steps in building a model.
#' 
#' @examples
#' \donttest{
#' # Subset example data to one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")[1:80,]
#' 
#' # Create model using only silent mutations and ignore multiple mutations
#' model <- createTargetingModel(db, model="s", sequenceColumn="sequence_alignment",
#'                               germlineColumn="germline_alignment_d_mask",
#'                               vCallColumn="v_call", multipleMutation="ignore")
#' 
#' # View top 5 mutability estimates
#' head(sort(model@mutability, decreasing=TRUE), 5)
#' 
#' # View number of silent mutations used for estimating mutability
#' model@numMutS
#' }
#' 
#' @export
createTargetingModel <- function(db, model=c("s", "rs"), sequenceColumn="sequence_alignment",
                                 germlineColumn="germline_alignment_d_mask",
                                 vCallColumn="v_call",
                                 multipleMutation=c("independent", "ignore"),
                                 minNumMutations=50, minNumSeqMutations=500,
                                 modelName="", modelDescription="", modelSpecies="", 
                                 modelCitation="", modelDate=NULL) {
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn, vCallColumn))
    if (check != TRUE) { stop(check) }
    
    # Check validity of input sequences 
    # (MUST NOT CONTAIN AMBIGUOUS CHARACTERS -- not supported)
    bool_obsv <- checkAmbiguousExist(db[[sequenceColumn]])
    bool_germ <- checkAmbiguousExist(db[[germlineColumn]])
    if (any(bool_obsv | bool_germ)) {
        stop("Ambiguous characters are not supported in input sequences.")
    }
    
    # Set date
    if (is.null(modelDate)) { modelDate <- format(Sys.time(), "%Y-%m-%d") }

    # Create models
    sub_mat<- createSubstitutionMatrix(db, model=model, 
                                       sequenceColumn=sequenceColumn,
                                       germlineColumn=germlineColumn,
                                       vCallColumn=vCallColumn,
                                       multipleMutation=multipleMutation,
                                       minNumMutations=minNumMutations,
                                       returnModel="5mer")
    mut_mat <- createMutabilityMatrix(db, sub_mat, model=model,
                                      sequenceColumn=sequenceColumn,
                                      germlineColumn=germlineColumn,
                                      vCallColumn=vCallColumn,
                                      multipleMutation=multipleMutation,
                                      minNumSeqMutations=minNumSeqMutations)

    # Extend 5-mers with Ns
    sub_mat <- extendSubstitutionMatrix(sub_mat)
    mut_mat <- extendMutabilityMatrix(mut_mat)
    
    # Make targeting matrix
    tar_mat <- createTargetingMatrix(sub_mat, mut_mat) 
    
    # Define TargetingModel object
    model_obj <- new("TargetingModel",
                     name=modelName,
                     description=modelDescription,
                     species=modelSpecies,
                     date=modelDate,
                     citation=modelCitation,
                     substitution=sub_mat,
                     mutability=mut_mat,
                     targeting=tar_mat@.Data,
                     numMutS=mut_mat@numMutS,
                     numMutR=mut_mat@numMutR)

    return(model_obj)
}


#' Calculate total mutability
#' 
#' \code{calculateMutability} calculates the total (summed) mutability for a set of sequences 
#' based on a 5-mer nucleotide mutability model.
#'
#' @param    sequences           character vector of sequences.
#' @param    model               \link{TargetingModel} object with mutation likelihood information.
#' @param    progress            if \code{TRUE} print a progress bar.
#' 
#' @return   Numeric vector with a total mutability score for each sequence.
#' 
#' @examples
#' \donttest{
#' # Subset example data to one isotype and sample as a demo
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")
#'
#' # Calculate mutability of germline sequences using \link{HH_S5F} model
#' mutability <- calculateMutability(sequences=db[["germline_alignment_d_mask"]], model=HH_S5F)
#' }
#' 
#' @export
calculateMutability <- function(sequences, model=HH_S5F, progress=FALSE) {
    # Initialize variables
    alphb <- seqinr::s2c("ACGTN")
    model_kmer <- names(model@mutability)
    model_rates <- as.vector(model@mutability)
    sequences <- toupper(sequences)
    sequences <- gsub("\\.", "N", sequences)
    
    # Mutability calculation
    mutability <- vector(mode="numeric", length=length(sequences))
    if (progress) { 
        pb <- progressBar(length(sequences)) 
    }
    for (s in 1:length(sequences)) {
        kmer <- seqinr::count(seqinr::s2c(sequences[s]), wordsize=5, alphabet=alphb)
        seq_kmer <- names(kmer)
        seq_counts <- as.vector(kmer)
        
        index <- match(seq_kmer, model_kmer)
        mutability[s] <- sum(seq_counts*model_rates[index], na.rm=TRUE)
        
        if (progress) { pb$tick() }
    }
    
    return(mutability)
}



# Create model and rescale mutabilities
# model <- createTargetingModel(db, model="s", multipleMutation="ignore")
# mut <- rescaleMutability(model)
rescaleMutability <- function(model, mean=1.0) {
    if (is(model, "TargetingModel")) {
        model <- model@mutability
    } else if (!is(model, "vector")) {
        stop("Input must be either a mutability vector or TargetingModel object.")
    }
    
    # TODO:  perhaps this is not so useful.  or perhaps it should be relative to max(model).
    # Renormalize
    rescaled <- model / sum(model, na.rm=T) * sum(!is.na(model)) * mean
    rescaled[!is.finite(rescaled)] <- NA
    
    return(rescaled)
}


# Remove in-frame IMGT gaps
# 
# @param    matInput  Nx2 matrix with input and germline sequences in each column
# @return   A two column matrix with "..." codons removed.
removeCodonGaps <- function(matInput) {
    # Function to return valid codon sets
    # i = position, x = codon list 1, y = codon list 2
    .f1 <- function(i, x, y) {
        xcod <- x[i]
        ycod <- y[i]
        if (xcod != "..." & ycod != "...") {
            c(xcod, ycod)
        } else {
            c("", "")
        }
    }

    # Function to parse sequences
    # z = vector of 2 sequences
    .f2 <- function(z) {
        # Split strings into 3-mers
        cods <- stri_extract_all_regex(c(z[1], z[2]), ".{3}")
        cmat <- sapply(1:length(cods[[1]]), .f1, x=cods[[1]], y=cods[[2]])
        apply(cmat, 1, paste, collapse="")
    }
    
    matCollapsed <- t(apply(matInput, 1, .f2))
    
    return(matCollapsed)
}


#' Make a degenerate 5-mer substitution model based on a 1-mer substitution model
#'
#' \code{makeDegenerate5merSub} populates substitution rates from a 1-mer substitution model
#' into 5-mers with corresponding central 1-mers.
#'
#' @param    sub1mer             a 4x4 matrix containing (normalized) substitution rates.
#'                               Row names should correspond to nucleotides to mutate from.
#'                               Column names should correspond to nucleotides to mutate into.
#'                               Nucleotides should include "A", "T", "G", and "C" 
#'                               (case-insensitive).
#' @param    extended            whether to return the unextended (\code{extended=FALSE}) or 
#'                               extended (\code{extended=TRUE}) 5-mer substitution model. 
#'                               Default is \code{FALSE}.
#'
#' @return   For \code{extended=FALSE}, a 4x1024 matrix. For \code{extended=TRUE}, a 5x3125 
#'           matrix.
#'
#' @details  As a concrete example, consider a 1-mer substitution model in which substitution
#'           rates from "A" to "T", "G", and "C" are, respectively, 0.1, 0.6, and 0.3. In the 
#'           resultant degenerate 5-mer substitution model, all the 5-mers (columns) that have 
#'           an "A" as their central 1-mer would have substitution rates (rows) of 0.1, 0.6, and 
#'           0.3 to "T", "G", and "C" respectively. 
#'           
#'           When \code{extended=TRUE}, \code{extendSubstitutionMatrix} is called to extend
#'           the 4x1024 substitution matrix.
#'  
#' @seealso  See \link{makeAverage1merSub} for making a 1-mer substitution model by taking
#'           the average of a 5-mer substitution model. See \link{extendSubstitutionMatrix}
#'           for extending the substitution matrix.
#' 
#' @examples
#' # Make a degenerate 5-mer model (4x1024) based on HKL_S1F (4x4)
#' # Note: not to be confused with HKL_S5F@substitution, which is non-degenerate
#' degenerate5merSub <- makeDegenerate5merSub(sub1mer = HKL_S1F)
#' 
#' # Look at a few 5-mers
#' degenerate5merSub[, c("AAAAT", "AACAT", "AAGAT", "AATAT")]
#' 
#' @export
makeDegenerate5merSub <- function(sub1mer, extended=FALSE) {
    # make sure that rownames and colnames of sub1mer are uppercase
    rownames(sub1mer) <- toupper(rownames(sub1mer))
    colnames(sub1mer) <- toupper(colnames(sub1mer))
    
    # create 5-mer labels using ATGC
    nuc_chars <- NUCLEOTIDES[1:4]
    nuc_words <- seqinr::words(5, nuc_chars)
    
    # get center positions of 5mers
    nuc_centers <- sapply(nuc_words, function(x){seqinr::s2c(x)[3]})
    
    # initiate 5-mer substitution matrix (4x1024)
    sub5mer <- matrix(NA, nrow=4, ncol=length(nuc_words),
                     dimnames=list(nuc_chars, nuc_words))
    
    # assign values from 1-mer model to 5-mer model
    for (from in rownames(sub1mer)) {
        for (to in colnames(sub1mer)) {
            if (from != to) { # if statement keeps diagonals as NA
                colIndex <- which(nuc_centers == from)
                sub5mer[to, colIndex] <- sub1mer[from, to]
            }
        }
    }
    stopifnot(dim(sub5mer) == c(4, 1024))
    
    # if extended=TRUE, extend
    if (extended) {
        sub5mer <- extendSubstitutionMatrix(sub5mer)
        stopifnot(dim(sub5mer) == c(5, 3125))
    }
    
    return(sub5mer)
}

#' Make a degenerate 5-mer mutability model based on a 1-mer mutability model
#'
#' \code{makeDegenerate5merMut} populates mutability rates from a 1-mer mutability model
#' into 5-mers with corresponding central 1-mers.
#'
#' @param    mut1mer             a named vector of length 4 containing (normalized) 
#'                               mutability rates. Names should correspond to nucleotides, 
#'                               which should include "A", "T", "G", and "C" 
#'                               (case-insensitive).
#' @param    extended            whether to return the unextended (\code{extended=FALSE}) or 
#'                               extended (\code{extended=TRUE}) 5-mer mutability model. 
#'                               Default is \code{FALSE}.
#'
#' @return   For \code{extended=FALSE}, a vector of length 1024. The vector returned is 
#'           normalized. For \code{extended=TRUE}, a vector of length 3125. 
#'
#' @details  As a concrete example, consider a 1-mer mutability model in which mutability
#'           rates of "A", "T", "G", and "C" are, respectively, 0.14, 0.23, 0.31, and 0.32. 
#'           In the resultant degenerate 5-mer mutability model, all the 5-mers that have 
#'           an "A" as their central 1-mer would have mutability rate of 0.14/256, where 256 is
#'           the number of such 5-mers. 
#'           
#'           When \code{extended=TRUE}, \code{extendMutabilityMatrix} is called to extend the
#'           mutability vector of length 1024 into a vector of length 3125.
#'  
#' @seealso  See \link{makeAverage1merMut} for making a 1-mer mutability model by 
#'           taking the average of a 5-mer mutability model. See 
#'           \link{extendMutabilityMatrix} for extending the mutability vector.
#' 
#' @examples
#' # Make a degenerate 5-mer model (length of 1024) based on a 1-mer model
#' example1merMut <- c(A=0.2, T=0.1, C=0.4, G=0.3)
#' degenerate5merMut <- makeDegenerate5merMut(mut1mer = example1merMut)
#' 
#' # Look at a few 5-mers
#' degenerate5merMut[c("AAAAT", "AACAT", "AAGAT", "AATAT")]
#' 
#' # Normalized
#' sum(degenerate5merMut)
#' 
#' @export
makeDegenerate5merMut <- function(mut1mer, extended=FALSE) {
    # make sure that names of mut1mer are uppercase
    names(mut1mer) <- toupper(names(mut1mer))
    
    # create 5-mer labels using ATGCN
    nuc_chars <- NUCLEOTIDES[1:4]
    nuc_words <- seqinr::words(5, nuc_chars)
    
    # get center positions of 5mers
    nuc_centers <- sapply(nuc_words, function(x){seqinr::s2c(x)[3]})
    
    # initiate 5-mer mutability vector (length of 3125)
    mut5mer <- rep(NA, length=length(nuc_words))
    names(mut5mer) <- nuc_words
    
    # assign values from 1-mer model to 5-mer model
    for (center in names(mut1mer)) {
        index <- which(nuc_centers == center)
        mut5mer[index] <- mut1mer[center]
    } 
    stopifnot(length(mut5mer) == 1024)
    
    # normalize
    mut5mer <- mut5mer / sum(mut5mer, na.rm=T)
    
    # if extended=TRUE, extend
    if (extended) {
        mut5mer <- extendMutabilityMatrix(mut5mer)
        stopifnot(length(mut5mer) == 3125)
    }
    
    return(mut5mer)
}


#' Make a 1-mer substitution model by averaging over a 5-mer substitution model
#'
#' \code{makeAverage1merSub} averages substitution rates in a 5-mer substitution model
#' to derive a 1-mer substitution model.
#'
#' @param    sub5mer             a 4x1024 matrix such as that returned by 
#'                               \code{createSubstitutionMatrix} and that returned by
#'                               \code{makeDegenerate5merSub} with \code{extended=FALSE}.
#'                               Column names should correspond to 5-mers containing the 
#'                               central 1-mer to mutate from. Row names should correspond to 
#'                               nucleotides to mutate into. Nucleotides should include 
#'                               "A", "T", "G", and "C" (case-insensitive).
#'
#' @return   A 4x4 matrix with row names representing nucleotides to mutate from and column
#'           names representing nucleotides to mutate into. Rates are normalized by row. 
#'
#' @details  For example, the substitution rate from "A" to "T" in the resultant 1-mer model
#'           is derived by averaging the substitution rates into a "T" of all the 5-mers that 
#'           have an "A" as their central 1-mer. 
#'  
#' @seealso  See \link{makeDegenerate5merSub} for making a degenerate 5-mer substitution 
#'           model based on a 1-mer substitution model. 
#' 
#' @examples
#' # Make a degenerate 5-mer model (4x1024) based on HKL_S1F (4x4)
#' degenerate5merSub <- makeDegenerate5merSub(sub1mer = HKL_S1F)
#' 
#' # Now make a 1-mer model by averaging over the degenerate 5-mer model
#' # Expected to get back HKL_S1F
#' makeAverage1merSub(sub5mer = degenerate5merSub)
#' 
#' @export
makeAverage1merSub <- function(sub5mer) {
    stopifnot(dim(sub5mer) == c(4, 1024))
    
    # make sure that rownames and colnames of sub5mer are uppercase
    rownames(sub5mer) <- toupper(rownames(sub5mer))
    colnames(sub5mer) <- toupper(colnames(sub5mer))
    
    # get 5-mers and center positions of 5-mers
    nuc_words <- colnames(sub5mer)
    nuc_centers <- sapply(nuc_words, function(x){seqinr::s2c(x)[3]})
    
    # create 1-mer labels using ATGC
    nuc_chars <- NUCLEOTIDES[1:4]
    
    # initiate 1-mer substitution matrix (4x4)
    sub1mer <- matrix(NA, nrow=length(nuc_chars), ncol=length(nuc_chars),
                     dimnames=list(nuc_chars, nuc_chars))
    
    # assign values from 5-mer model to 1-mer model
    for (from in rownames(sub1mer)) {
        for (to in colnames(sub1mer)) {
            if (from != to) { # if statement keeps diagonals as NA
                colIndex <- which(nuc_centers == from)
                sub1mer[from, to] <- mean(sub5mer[to, colIndex], na.rm=T)
            }
        }
    }
    stopifnot(dim(sub1mer) == c(4, 4))
    
    # normalize
    # tricky: apply transposes result; use t() to transpose back 
    sub1mer <- t(apply(sub1mer, 1, function(x){x/sum(x, na.rm=T)}))
    
    return(sub1mer)
}


#' Make a 1-mer mutability model by averaging over a 5-mer mutability model
#'
#' \code{makeAverage1merMut} averages mutability rates in a 5-mer mutability model
#' to derive a 1-mer mutability model.
#'
#' @param    mut5mer             a named vector of length 1024 such as that returned by 
#'                               \code{createMutabilityMatrix} and that returned by
#'                               \code{makeDegenerate5merMut} with \code{extended=FALSE}.
#'                               Names should correspond to 5-mers made up of "A", "T", 
#'                               "G", and "C" (case-insensitive). \code{NA} values are 
#'                               allowed.
#'
#' @return   A named vector of length 4 containing normalized mutability rates.
#'
#' @details  For example, the mutability rate of "A" in the resultant 1-mer model
#'           is derived by averaging the mutability rates of all the 5-mers that 
#'           have an "A" as their central 1-mer, followed by normalization.
#'  
#' @seealso  See \link{makeDegenerate5merMut} for making a degenerate 5-mer mutability
#'           model based on a 1-mer mutability model. 
#' 
#' @examples
#' # Make a degenerate 5-mer model (length of 1024) based on a 1-mer model
#' example1merMut <- c(A=0.2, T=0.1, C=0.4, G=0.3)
#' degenerate5merMut <- makeDegenerate5merMut(mut1mer = example1merMut)
#'  
#' # Now make a 1-mer model by averaging over the degenerate 5-mer model
#' # Expected to get back example1merMut
#' makeAverage1merMut(mut5mer = degenerate5merMut)
#' 
#' @export
makeAverage1merMut <- function(mut5mer) {
    stopifnot(length(mut5mer) == 1024)
    
    # make sure that names mut5mer are uppercase
    names(mut5mer) <- toupper(names(mut5mer))
    
    # get 5-mers and center positions of 5-mers
    nuc_words <- names(mut5mer)
    nuc_centers <- sapply(nuc_words, function(x){seqinr::s2c(x)[3]})
    
    # create 1-mer labels using ATGC
    nuc_chars <- NUCLEOTIDES[1:4]
    
    # initiate 1-mer mutability vector (length 4)
    mut1mer <- rep(NA, length=length(nuc_chars))
    names(mut1mer) <- nuc_chars
    
    # assign values from 5-mer model to 1-mer model
    for (center in names(mut1mer)) {
        index <- which(nuc_centers == center)
        mut1mer[center] <- mean(mut5mer[index], na.rm=T)
    }
    stopifnot(length(mut1mer) == 4)
    
    # normalize
    mut1mer <- mut1mer / sum(mut1mer, na.rm=T)
    
    return(mut1mer)
}

#### Distance Functions ####

#' Calculates a 5-mer distance matrix from a TargetingModel object
#' 
#' \code{calcTargetingDistance} converts either the targeting rates in a \code{TargetingModel}
#'  model to a matrix of 5-mer to single-nucleotide mutation distances, or the substitution 
#'  rates in a 1-mer substitution model to a symmetric distance matrix.
#' 
#' @param    model     \link{TargetingModel} object with mutation likelihood information, or
#'                     a 4x4 1-mer substitution matrix normalized by row with rownames and 
#'                     colnames consisting of "A", "T", "G", and "C".
#' @param    places    decimal places to round distances to.
#'                                                
#' @return   For input of \link{TargetingModel}, a matrix of distances for each 5-mer motif with 
#'           rows names defining the center nucleotide and column names defining the 5-mer 
#'           nucleotide sequence. For input of 1-mer substitution matrix, a 4x4 symmetric distance
#'           matrix. 
#'           
#' @details
#' The targeting model is transformed into a distance matrix by:
#' \enumerate{
#'   \item  Converting the likelihood of being mutated \eqn{p=mutability*substitution} to 
#'          distance \eqn{d=-log10(p)}.
#'   \item  Dividing this distance by the mean of the distances.
#'   \item  Converting all infinite, no change (e.g., A->A), and NA distances to 
#'          zero.
#' }
#' 
#' The 1-mer substitution matrix is transformed into a distance matrix by:
#' \enumerate{
#'   \item  Symmetrize the 1-mer substitution matrix.
#'   \item  Converting the rates to distance \eqn{d=-log10(p)}.
#'   \item  Dividing this distance by the mean of the distances.
#'   \item  Converting all infinite, no change (e.g., A -> A), and NA distances to 
#'          zero.
#' }
#' 
#' @seealso  See \link{TargetingModel} for this class of objects and
#'           \link{createTargetingModel} for building one.
#' 
#' @examples
#' # Calculate targeting distance of HH_S5F
#' dist <- calcTargetingDistance(HH_S5F)
#' 
#' # Calculate targeting distance of HH_S1F
#' dist <- calcTargetingDistance(HH_S1F)
#' 
#' @export
calcTargetingDistance <- function(model, places=2) {
  # if model is TargetingModel object, assume 5-mer targeting model
  # extract targeting matrix
  if (is(model, "TargetingModel")) {
    input <- "5mer"
    model <- model@targeting
  } else if (is(model, "matrix") & all(dim(model) == c(4, 4))) {
    # if model is a matrix, assume 1-mer substitution matrix and symmetrize
    input <- "1mer"
    model <- symmetrize(model)
  } else {
    # anything else: error
    stop("Input must be either a 4x4 targeting matrix or TargetingModel object.")
  }
  
  # Take negative log10 of all probabilities
  model_dist <- -log10(model)
  # Center distances on 1 (divide by mean)
  model_dist <- model_dist/mean(model_dist, na.rm=T)
  # Set infinite values to NA
  model_dist[!is.finite(model_dist)] <- NA
  
  # TODO: the indexing of A-->A etc positions can probably be done smarter/faster
  # Assign no-change entries to distance 0
  if (input == "5mer") {
      center_nuc <- gsub("..([ACGTN])..", "\\1", colnames(model_dist))
      for (i in 1:length(center_nuc)) {
          model_dist[center_nuc[i], i] <- 0
      }
      # Assign 0 to N and 5mers with N in center position 
      model_dist[,center_nuc=="N"] <- 0
      model_dist["N",] <- 0
  } else if (input == "1mer") {
      diag(model_dist) <- 0
      model_dist <- rbind(model_dist, matrix(0, 3, 4))
      model_dist <- cbind(model_dist, matrix(0, 7, 3))
      colnames(model_dist)[5:7] <- rownames(model_dist)[5:7] <- c("N", "-", ".")
  }
  
  # Round
  model_dist <- round(model_dist, places)
  
  return(model_dist)
}

# (From G Yaari)
# Symmetrize a non-symmetric, 1-mer substitution matrix
# 
# \code{symmetrize} makes a 1-mer substitution matrix symmetric by minimizing the 
# rss between it and a symmetric matrix.
# 
# @param    sub1mer     a 4x4 matrix normalized by row. Each row sums up to 1 and 
#                       reflects substitution probabilities for each nucleotide. 
#                       Rownames and colnames are "A","C","G", and "T".
#
# @return   a 4x4 symmetrix matrix with \code{NA}s on the diagonal.
# 
# @details  The input 1-mer substitution matrix is approximated by a symmetric matrix 
#           both with respect to time (e.g. C->T = T->C), and with respect to strand 
#           (e.g. C->T = G->A). The symmetric matrix has three free parameters that 
#           are estimated by minimizing the sum of squares between this matrix and 
#           the input matrix. The fitted matrix was normalized to ensure that each 
#           row sums up to 1.
symmetrize <- function(sub1mer) {
  rownames(sub1mer) <- toupper(rownames(sub1mer))
  colnames(sub1mer) <- toupper(colnames(sub1mer))
  
  # check that rows sum up to 1
  stopifnot(all.equal(rowSums(sub1mer), rep(1, 4), 
                      check.names=FALSE, tolerance=0.055))
  
  .minDist <- function(Pars, subMtx) {
    (Pars[1] - subMtx["A", "C"])^2 + (Pars[1] - subMtx["C", "A"])^2 + 
    (Pars[1] - subMtx["G", "T"])^2 + (Pars[1] - subMtx["T", "G"])^2 +
    (Pars[2] - subMtx["A", "G"])^2 + (Pars[2] - subMtx["G", "A"])^2 + 
    (Pars[2] - subMtx["C", "T"])^2 + (Pars[2] - subMtx["T", "C"])^2 +
    (Pars[3] - subMtx["A", "T"])^2 + (Pars[3] - subMtx["T", "A"])^2 + 
    (Pars[3] - subMtx["C", "G"])^2 + (Pars[3] - subMtx["G", "C"])^2
  }
  
  pars <- optim(par=rep(0, 3), fn=.minDist, subMtx=sub1mer)$par
  pars <- pars/sum(pars)
  symmetric_substitution <- sub1mer
  symmetric_substitution["A", 2:4] <- pars
  symmetric_substitution["C", c(1, 4, 3)] <- pars
  symmetric_substitution["G", c(4, 1, 2)] <- pars
  symmetric_substitution["T", c(3, 2, 1)] <- pars
  
  # NAs on diagonal instead of 0 so that calcTargetingDistance works with 1-mer model
  diag(symmetric_substitution) <- NA
  return(symmetric_substitution)
}


#' Write targeting model distances to a file
#' 
#' \code{writeTargetingDistance} writes a 5-mer targeting distance matrix 
#' to a tab-delimited file.
#' 
#' @param    model     \link{TargetingModel} object with 
#'                     mutation likelihood information.
#' @param    file      name of file to write.
#'                                                
#' @return   NULL
#'           
#' @details
#' The targeting distance write as a tab-delimited 5x3125 matrix. Rows define the mutated 
#' nucleotide at the center of each 5-mer, one of \code{c("A", "C", "G", "T", "N")}, 
#' and columns define the complete 5-mer of the unmutated nucleotide sequence. 
#' \code{NA} values in the distance matrix are replaced with distance 0.
#'    
#' @seealso  Takes as input a \link{TargetingModel} object and calculates  
#'           distances using \link{calcTargetingDistance}.
#' 
#' @examples
#' \dontrun{
#' # Write HS5F targeting model to working directory as hs5f.tab
#' writeTargetingDistance(HH_S5F, "hh_s5f.tsv") 
#' }
#' 
#' @export
writeTargetingDistance <- function(model, file) {
    to_write <- as.data.frame(calcTargetingDistance(model))
    to_write[is.na(to_write)] <- 0
    write.table(to_write, file, quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)
}


#### Plotting functions ####

#' Plot mutability probabilities
#' 
#' \code{plotMutability} plots the mutability rates of a \code{TargetingModel}.
#' 
#' @param    model        \link{TargetingModel} object or vector containing normalized 
#'                        mutability rates.
#' @param    nucleotides  vector of center nucleotide characters to plot.
#' @param    mark         vector of 5-mer motifs to highlight in the plot. If \code{NULL}
#'                        only highlight classical hot and cold spot motifs.
#' @param    style        type of plot to draw. One of:
#'                        \itemize{
#'                          \item \code{"hedgehog"}:  circular plot showing higher mutability
#'                                                    scores further from the circle. The 5-mer
#'                                                    is denoted by the values of the inner 
#'                                                    circle. The 5-mer is read from the most interior 
#'                                                    position of the 5-mer (5') to most exterior position 
#'                                                    (3'), with the center nucleotide in the center ring.
#'                                                    Note, the order in which the 5-mers are plotted is
#'                                                    different for nucleotides \code{c("A", "C")} and 
#'                                                    \code{c("G", "T")}.
#'                          \item \code{"bar"}:       bar plot of mutability similar to the 
#'                                                    \code{hedgehog} style with the most 5' positions
#'                                                    of each 5-mer at the base of the plot.
#'                        }
#' @param    size         numeric scaling factor for lines and text in the plot.
#' @param    silent       if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                        objects; if \code{FALSE} draw the plot.
#' @param    ...          additional arguments to pass to ggplot2::theme.
#'                                                
#' @return   A named list of ggplot objects defining the plots, with names defined by the 
#'           center nucleotide for the plot object.
#'    
#' @seealso  Takes as input a \link{TargetingModel} object. 
#'           See \link{createTargetingModel} for model building.
#' 
#' 
#' @examples
#' # Plot one nucleotide in circular style
#' plotMutability(HH_S5F, "C")
#' 
#' # Plot two nucleotides in barchart style
#' plotMutability(HH_S5F, c("G", "T"), style="bar")
#' 
#' @export
plotMutability <- function(model, nucleotides=c("A", "C", "G", "T"), mark=NULL,
                           style=c("hedgehog", "bar"), size=1, silent=FALSE, 
                           ...) {
    # model=HH_S5F
    # nucleotides=c("C")
    # nucleotides=c("A", "C", "G", "T")
    # style="hedgehog"
    # size=1
    # silent=FALSE
    
    # Check input
    nucleotides <- toupper(nucleotides)
    style <- match.arg(style)
    
    if (is(model, "TargetingModel")) {
        model <- model@mutability
    } else if (!is(model, "vector")) {
        stop("Input must be either a mutability vector or TargetingModel object.")
    }
    
    # Set base plot settings
    base_theme <- theme_bw() +
        theme(panel.spacing=grid::unit(0, "lines"),
              panel.background=element_blank()) +
        theme(axis.text=element_text(margin=grid::unit(0, "lines"))) +
        theme(text=element_text(size=10*size),
              title=element_text(size=10*size),
              legend.spacing=grid::unit(0, "lines"),
              legend.background=element_blank())

    # Scaling and layout parameters
    score_offset <- 0
    score_scale <- 15
    text_offset <- -5.6
    
    # Set guide colors
    motif_colors <- setNames(c("#000000", "#4daf4a", "#e41a1c", "#094d85", "#999999"),
                             c("Marked", "WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))
    dna_colors <- setNames(c("#7bce77", "#ff9b39", "#f04949", "#5796ca", "#c4c4c4"),
                           c("A", "C", "G", "T", "N"))
    
    # Build data.frame of mutability scores
    mut_scores <- model[!grepl("N", names(model))]
    mut_scores[!is.finite(mut_scores)] <- 0
    mut_words <- names(mut_scores)
    mut_positions <- as.data.frame(t(sapply(mut_words, seqinr::s2c)))
    colnames(mut_positions) <- paste0("pos", 1:ncol(mut_positions))
    mut_df <- data.frame(word=mut_words, 
                         score=mut_scores, 
                         mut_positions)

    # Add hot and cold-spot motif information
    mut_df$motif <- "Neutral"
    mut_df$motif[grepl("(.[AT]A..)|(..T[AT].)", mut_df$word, perl=TRUE)] <- "WA/TW"
    mut_df$motif[grepl("([AT][GA]C..)|(..G[CT][AT])", mut_df$word, perl=TRUE)] <- "WRC/GYW"
    mut_df$motif[grepl("([CG][CT]C..)|(..G[GA][CG])", mut_df$word, perl=TRUE)] <- "SYC/GRS"
    if (is.null(mark)) {
        mut_df$motif <- factor(mut_df$motif, levels=c("WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))
    } else {
        mut_df$motif[mut_df$word %in% mark] <- "Marked"
        mut_df$motif <- factor(mut_df$motif, levels=c("Marked", "WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))
    }
    
    # Subset to nucleotides of interest
    mut_df <- mut_df[mut_df$pos3 %in% nucleotides, ]

    # Functions to transform and revert score for plotting
    score_max <- max(mut_df$score, na.rm=TRUE)
    .transform_score <- function(x) { x / score_max * score_scale + score_offset }
    .invert_score <- function(x) { (x - score_offset) / score_scale * score_max }
    
    # Rescale scores for plotting
    mut_df$score <- .transform_score(mut_df$score)
    
    # Build plots for each center nucleotide
    plot_list <- list()
    for (center_nuc in nucleotides) {
        # center_nuc <- "C"
        # Subset data to center nucleotide
        sub_df <- mut_df[mut_df$pos3 == center_nuc, ]
        
        # Order 5-mers by positions, with reversed order if center nucleotide is G or T
        if (center_nuc %in% c("A", "C")) {
            sub_df <- dplyr::arrange(sub_df, !!!rlang::syms(c("pos1", "pos2", "pos4", "pos5")))
            sub_df$x <- 1:nrow(sub_df)            
        } else if (center_nuc %in% c("G", "T")) {
            sub_df <- dplyr::arrange(sub_df, !!!rlang::syms(c("pos5", "pos4", "pos2", "pos1")))
            sub_df$x <- 1:nrow(sub_df)
        } else {
            stop("Invalid nucleotide choice")
        }
        
        # Melt 5-mer position data
        sub_melt <- sub_df %>% 
            tidyr::gather("pos", "char", !!!rlang::syms(colnames(mut_positions))) %>% 
            select("x", "pos", "char")
        #sub_melt$pos <- factor(sub_melt$pos, levels=mut_names)
        #sub_melt$pos <- as.numeric(sub_melt$pos)
        sub_melt$pos <- as.numeric(gsub("pos", "", sub_melt$pos))

        # Define nucleotide text and rectangle position data
        sub_text <- list()
        for (i in 1:5) {
            nuc_rle <- rle(sub_melt$char[sub_melt$pos == i])

            # Set rectangle x limits
            rect_max <- cumsum(nuc_rle$lengths)
            rect_min <- rect_max - diff(c(0, rect_max))
            
            # Set text position
            if (length(rect_max) > 1) { 
                text_x <- rect_max - diff(c(0, rect_max)) / 2
            } else { 
                text_x <- rect_max / 2 
            }
        
            tmp_df <- data.frame(text_x=text_x, 
                                 text_y=i,
                                 text_label=factor(nuc_rle$values, levels=names(dna_colors)),
                                 rect_min=rect_min,
                                 rect_max=rect_max)
            
            sub_text[[i]] <- tmp_df
        }

        # Define text and rectangle positions for inner circle
        sub_melt$pos <- sub_melt$pos + text_offset
        sub_text <- lapply(sub_text, function(x) { dplyr::mutate(x, text_y=!!rlang::sym("text_y") + !!rlang::sym("text_offset")) })
        sub_rect <- dplyr::bind_rows(sub_text) %>%
            mutate(rect_width=rect_max - rect_min,
                    ymin=!!rlang::sym("text_y") - 0.5,
                    ymax=!!rlang::sym("text_y") + 0.5)
        
        # Use only colors for motifs present in sub_df
        sub_colors <- motif_colors[names(motif_colors) %in% sub_df$motif]
        
        # Define base plot object
        p1 <- ggplot(sub_df) + 
            base_theme + 
            #ggtitle(paste0("NN", center_nuc, "NN")) +
            xlab("") +
            ylab("") + 
            scale_color_manual(name="Motif", values=c(sub_colors, dna_colors), 
                               breaks=names(sub_colors)) +
            scale_fill_manual(name="", values=c(sub_colors, dna_colors), guide="none") +
            geom_rect(data=sub_rect, 
                      mapping=aes(xmin=rect_min, xmax=rect_max, ymin=ymin, ymax=ymax, 
                                         fill=text_label, color=text_label), 
                      size=0.5*size, alpha=1, show.legend=FALSE) +
            #geom_tile(data=sub_rect, 
            #          mapping=aes_string(x="text_x", y="text_y", width="rect_width", fill="text_label"), 
            #          size=0, alpha=0.7, show.legend=FALSE) +
            #geom_tile(data=sub_melt, mapping=aes_string(x="x", y="pos", fill="char"), size=0, alpha=0.7,
            #          show.legend=FALSE) +
            geom_text(data=sub_text[[3]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                      color="black", hjust=0.5, vjust=0.5, size=3*size, fontface=2)
        
        # Add 5-mer nucleotide labels
        if (center_nuc %in% c("A", "C")) {
            p1 <- p1 + geom_text(data=sub_text[[1]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                                 color="black", hjust=0.5, vjust=0.5, size=2*size) +
                geom_text(data=sub_text[[2]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                          color="black", hjust=0.5, vjust=0.5, size=2*size)
        } else if (center_nuc %in% c("G", "T")) {
            p1 <- p1 + geom_text(data=sub_text[[4]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                                 color="black", hjust=0.5, vjust=0.5, size=2*size) +
                geom_text(data=sub_text[[5]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                          color="black", hjust=0.5, vjust=0.5, size=2*size)
        }

        # Add style options and mutability scores
        if (style == "hedgehog") {
            y_limits <- c(text_offset - 1, score_scale + score_offset)
            #orient_x <- sub_text[[3]]$text_x[1]
            #orient_y <- text_offset - 1
            p1 <- p1 + theme(plot.margin=grid::unit(c(0, 0, 0, 0), "lines"),
                             panel.grid=element_blank(), 
                             panel.border=element_blank(),
                             axis.title=element_blank(),
                             axis.text=element_blank(), 
                             axis.ticks=element_blank(),
                             legend.direction="horizontal",
                             legend.justification=c(0.5, 1),
                             legend.position=c(0.5, 1)) +
                scale_x_continuous(expand=c(0, 0)) +
                scale_y_continuous(limits=y_limits, expand=c(0, 0)) +
                coord_polar(theta="x") +
                geom_segment(data=sub_df, mapping=aes(x=x, xend=x, yend=score, color=motif), 
                             y=score_offset, size=0.75*size, position=position_nudge(x = -0.5)) +
                guides(color=guide_legend(override.aes=list(linetype=1, size=2*size)))
              
        } else if (style == "bar") {
            y_breaks <- seq(score_offset, score_scale + score_offset, 1)
            y_limits <- c(text_offset + 0.5, score_scale + score_offset)
            p1 <- p1 + theme(plot.margin=grid::unit(c(1, 1, 1, 1), "lines"),
                             panel.grid=element_blank(), 
                             panel.border=element_rect(color="black"),
                             axis.text.x=element_blank(), 
                             axis.ticks.x=element_blank(),
                             legend.position="top") +
                ylab("Mutability") +
                scale_x_continuous(expand=c(0, 1)) +
                scale_y_continuous(limits=y_limits, breaks=y_breaks, expand=c(0, 0.5),
                                   labels=function(x) scales::scientific(.invert_score(x))) +
                geom_bar(data=sub_df, mapping=aes(x=x, y=score, fill=motif, color=motif), 
                         stat="identity", position=position_nudge(x = -0.5), size=0, width=0.7) +
              guides(color=guide_legend(override.aes=list(fill=sub_colors, linetype=0)))
        }

        # Add additional theme elements
        p1 <- p1 + do.call(theme, list(...))
        
        # Add plot to list
        plot_list[[center_nuc]] <- p1
    }

    
    # Plot
    if (!silent) { 
        do.call(gridPlot, args=c(plot_list, ncol=length(plot_list))) 
    }
    
    invisible(plot_list)
}


#' Visualize parameter tuning for minNumMutations and minNumSeqMutations
#'
#' Visualize results from \link{minNumMutationsTune} and \link{minNumSeqMutationsTune}
#' 
#' @param    tuneMtx             a \code{matrix} or a \code{list} of matrices produced by either 
#'                               \link{minNumMutationsTune} or \link{minNumSeqMutationsTune}.
#'                               In the case of a list, it is assumed that each matrix corresponds
#'                               to a sample and that all matrices in the list were produced using
#'                               the same set of trial values of \code{minNumMutations} or 
#'                               \code{minNumSeqMutations}.
#' @param    thresh              a number or a vector of indicating the value or the range of values
#'                               of \code{minNumMutations} or \code{minNumSeqMutations} to plot. 
#'                               Should correspond to the columns of \code{tuneMtx}.
#' @param    criterion           one of \code{"5mer"}, \code{"3mer"}, \code{"1mer"}, or \code{"3mer+1mer"} 
#'                               (for \code{tuneMtx} produced by \link{minNumMutationsTune}), or either 
#'                               \code{"measured"} or \code{"inferred"} (for \code{tuneMtx} produced by 
#'                               \link{minNumSeqMutationsTune}).                
#' @param    pchs                point types to pass on to \link{plot}.
#' @param    ltys                line types to pass on to \link{plot}.
#' @param    cols                colors to pass on to \link{plot}.                             
#' @param    plotLegend          whether to plot legend. Default is \code{TRUE}. Only applicable 
#'                               if \code{tuneMtx} is a named list with names of the matrices 
#'                               corresponding to the names of the samples.
#' @param    legendPos           position of legend to pass on to \link{legend}. Can be either a
#'                               numeric vector specifying x-y coordinates, or one of 
#'                               \code{"topright"}, \code{"center"}, etc. Default is \code{"topright"}.
#' @param    legendHoriz         whether to make legend horizontal. Default is \code{FALSE}.
#' @param    legendCex           numeric values by which legend should be magnified relative to 1.
#' 
#' @details  For \code{tuneMtx} produced by \link{minNumMutationsTune}, for each sample, depending on
#'           \code{criterion}, the numbers of 5-mers for which substitution rates are directly computed
#'           (\code{"5mer"}), inferred based on inner 3-mers (\code{"3mer"}), inferred based on 
#'           central 1-mers (\code{"1mer"}), or inferred based on inner 3-mers and central 1-mers
#'           (\code{"3mer+1mer"}) are plotted on the y-axis against values of \code{minNumMutations} 
#'           on the x-axis.
#' 
#'           For \code{tuneMtx} produced by \link{minNumSeqMutationsTune}, for each sample, depending on
#'           \code{criterion}, the numbers of 5-mers for which mutability rates are directly measured
#'           (\code{"measured"}) or inferred (\code{"inferred"}) are plotted on the y-axis against values
#'           of \code{minNumSeqMutations} on the x-axis.
#'           
#'           Note that legends will be plotted only if \code{tuneMtx} is a supplied as a named \code{list}
#'           of matrices, ideally with names of each \code{matrix} corresponding to those of the samples 
#'           based on which the matrices were produced, even if \code{plotLegend=TRUE}.
#' 
#' @seealso  See \link{minNumMutationsTune} and \link{minNumSeqMutationsTune} for generating 
#'           \code{tuneMtx}. 
#' 
#' @examples
#' \donttest{
#' # Subset example data to one isotype and 200 sequences
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHA")
#' set.seed(112)
#' db <- dplyr::slice_sample(db, n=50)
#' 
#' tuneMtx = list()
#' for (i in 1:length(unique(db$sample_id))) {
#'     # Get data corresponding to current sample
#'     curDb = db[db[["sample_id"]] == unique(db[["sample_id"]])[i], ]
#'     
#'     # Count the number of mutations per 5-mer
#'     subCount = createSubstitutionMatrix(db=curDb, model="s", 
#'                                         sequenceColumn="sequence_alignment",
#'                                         germlineColumn="germline_alignment_d_mask",
#'                                         vCallColumn="v_call",
#'                                         multipleMutation="independent",
#'                                         returnModel="5mer", numMutationsOnly=TRUE)
#'     
#'     # Tune over minNumMutations = 5..50
#'     subTune = minNumMutationsTune(subCount, seq(from=5, to=50, by=5))
#'     
#'     tuneMtx = c(tuneMtx, list(subTune))
#' }
#'
#' # Name tuneMtx after sample names 
#' names(tuneMtx) = unique(db[["sample_id"]])
#' 
#' # plot with legend for both samples for a subset of minNumMutations values
#' plotTune(tuneMtx, thresh=c(5, 15, 25, 40), criterion="3mer",
#'          pchs=16:17, ltys=1:2, cols=2:3, 
#'          plotLegend=TRUE, legendPos=c(25, 30))
#' 
#' # plot for only 1 sample for all the minNumMutations values (no legend)
#' plotTune(tuneMtx[[1]], thresh=seq(from=5, to=50, by=5), criterion="3mer")
#' }
#' 
#' @export
plotTune <- function(tuneMtx, thresh, 
                    criterion=c("5mer", "3mer", "1mer", "3mer+1mer", 
                                "measured", "inferred"), 
                    pchs = 1, ltys = 2, cols = 1,
                    plotLegend = TRUE, legendPos = "topright", 
                    legendHoriz = FALSE, legendCex = 1) {
  
  stopifnot(length(criterion)==1)
  stopifnot(is.matrix(tuneMtx) | is.list(tuneMtx))
  
  ### extract plot data into plotMtx
  # if tuneMtx is just a matrix
  if (!is.list(tuneMtx)) {
    if (criterion!="3mer+1mer") {
      plotMtx <- matrix(tuneMtx[criterion, as.character(thresh)], nrow=1)
    } else {
      plotMtx <- matrix(colSums(tuneMtx[c("3mer", "1mer"), as.character(thresh)]), nrow=1)
    }
  } else {
  # if tuneMtx is a named list of matrices (e.g. corresponding to multiple samples)  
    if (criterion!="3mer+1mer") {
      plotMtx <- do.call(base::rbind, 
                        lapply(tuneMtx, 
                               function(mtx){mtx[criterion, as.character(thresh)]}))
    } else {
      plotMtx <- do.call(base::rbind, 
                        lapply(tuneMtx, 
                               function(mtx){colSums(mtx[c("3mer", "1mer"), 
                                                         as.character(thresh)])}))
    }
    rownames(plotMtx) <- names(tuneMtx)
  }
  # sanity check: there should not be any NA
  stopifnot(!any(is.na(plotMtx)))
  
  ### if number of pchs/ltys/cols provided does not match number of samples expected
  # expand into vector with repeating values (otherwise legend would break)
  if (length(pchs)!=nrow(plotMtx)) {pchs <- rep(pchs, length.out=nrow(plotMtx))}
  if (length(ltys)!=nrow(plotMtx)) {ltys <- rep(ltys, length.out=nrow(plotMtx))}
  if (length(cols)!=nrow(plotMtx)) {cols <- rep(cols, length.out=nrow(plotMtx))}
  
  
  ### axis labels
  if (criterion %in% c("5mer", "3mer", "1mer", "3mer+1mer")) {
    xlab.name <- "Minimum # mutations per 5-mer to\ndirectly compute 5-mer substitution rates"
    # cannot use switch because variable names cannot start with number
    ylab.name <- "# 5-mers for which substitution rates are\n"
    if (criterion=="5mer") {
      ylab.name <- paste(ylab.name, "directly computed")
    } else if (criterion=="3mer") {
      ylab.name <- paste(ylab.name, "inferred based on inner 3-mers")
    } else if (criterion=="1mer") {
      ylab.name <- paste(ylab.name, "inferred based on central 1-mers")
    } else if (criterion=="3mer+1mer") {
      ylab.name <- paste(ylab.name, "inferred based on 3- and 1-mers")
    }
  } else if (criterion %in% c("measured", "inferred")) {
    xlab.name <- "Minimum # mutations in sequences containing each 5-mer\nto directly compute mutability"
    ylab.name <- paste("# 5-mers for which mutability is", criterion)
  }
  
  ### plot
  # bottom, left, top, right
  par(mar=c(6, 6, 4, 2) + 0.1)
  for (i in 1:nrow(plotMtx)) {
    if (i==1) {
      plot(x=thresh, y=plotMtx[i, ], 
           ylim=range(plotMtx),
           xaxt="n", xlab="", ylab="",
           cex.axis=1.5, type="b", lwd=1.5,
           pch=pchs[i], lty=ltys[i], col=cols[i])
      axis(side=1, at=thresh, cex.axis=1.5)
      mtext(side=1, text=xlab.name, line=4, cex=1.2)
      mtext(side=2, text=ylab.name, line=3, cex=1.2)
    } else {
      points(x=thresh, y=plotMtx[i, ],
             type="b", lwd=1.5,
             pch=pchs[i], lty=ltys[i], col=cols[i])
    }
  }
  
  ### legend (even if plotLegend=T, only plotted if tuneMtx is a named list)
  if ( !is.null(rownames(plotMtx)) & plotLegend ) {
    # if legendPos specified as xy coordinates
    if (is.numeric(legendPos) & length(legendPos)==2) {
      legend(x=legendPos[1], y=legendPos[2], 
             legend = c("Sample", rownames(plotMtx)),
             horiz = legendHoriz, cex = legendCex,
             pch=c(NA, pchs), lty=c(NA, ltys), col=c(NA, cols))
    } else {
    # if legendPos specified as "center", "topright", etc.  
      legend(legendPos, 
             legend = c("Sample", rownames(plotMtx)),
             horiz = legendHoriz, cex = legendCex,
             pch=c(NA, pchs), lty=c(NA, ltys), col=c(NA, cols))
    }
  }
  
}

#### Original BASELINe functions ####

# Given a nuc, returns the other 3 nucs it can mutate to
canMutateTo <- function(nuc) {
    NUCLEOTIDES[1:4][NUCLEOTIDES[1:4] != nuc]
}


# Compute the mutations types
# matOfCodons: nx2; n=pairs of codons; 1st col=codonTo, 2nd col=codonFrom
# NOTE: this function is not intended to be used where input sequences have 
#       ambiguous characters; it assumes that only 1 entry (r/s/stop/na) from
#       mutationType is non-zero/1
mutationTypeOptimized <- function(matOfCodons) {
    # mutType: 4xn; rows: r/s/stop/na
    mutType <- apply(matOfCodons, 1, function(x) { mutationType(x[2], x[1]) })
    idx <- apply(mutType, 2, function(y){which(y>0)[1]})
    mutType <- rownames(mutType)[idx]
    mutType[which(mutType=="na")] <- NA
    return(mutType)
}


# row 1 = GL
# row 2 = Seq
# in_matrix <- matrix(c(c("A","A","A","C","C","C"), c("A","G","G","C","C","A")), 2 ,6, byrow=T)
# analyzeMutations2NucUri(in_matrix)
analyzeMutations2NucUri <- function(in_matrix) {
    if(ncol(in_matrix) > VLENGTH) {
        paramGL <- in_matrix[2,1:VLENGTH]
        paramSeq <- in_matrix[1,1:VLENGTH]
    } else {
        paramGL <- in_matrix[2,]
        paramSeq <- in_matrix[1,]
    }
    #mutations = apply(rbind(paramGL,paramSeq), 2, function(x){!x[1]==x[2]})
    mutations_val <- paramGL != paramSeq
    if (any(mutations_val)) {
        mutationPos <- {1:length(mutations_val)}[mutations_val]
        #mutationPos = mutationPos[sapply(mutationPos, function(x){!any(paramSeq[getCodonPos(x)]=="N")})]
        length_mutations =length(mutationPos)
        mutationInfo <- rep(NA,length_mutations)
        if (any(mutationPos)) {
            pos<- mutationPos
            pos <- pos[!is.na(pos)]
            pos_array <- array(sapply(pos,getCodonPos))
            codonGL <- paramGL[pos_array]
            codonGL[is.na(codonGL)] <- "N"
            codonSeq <- sapply(pos,function(x){
                seqP <- paramGL[getCodonPos(x)]
                muCodonPos <- {x-1}%%3+1
                seqP[muCodonPos] <- paramSeq[x]
                return(seqP)
            })
            codonSeq[is.na(codonSeq)] <- "N"
            GLcodons <-  apply(matrix(codonGL,length_mutations,3,byrow=TRUE),1,c2s)
            Seqcodons <-   apply(codonSeq,2,c2s)
            mutationInfo <- apply(rbind(GLcodons , Seqcodons),2,function(x){
                # not intended to be used where input sequences have 
                # ambiguous characters; it assumes that only 1 entry (r/s/stop/na) from
                # mutationType is non-zero/1
                mutType <- mutationType(c2s(x[1]),c2s(x[2]))
                mutType <- names(mutType)[which(mutType>0)]
                if (mutType=="na") {mutType=NA}
                return(mutType)
                })
            names(mutationInfo) <- mutationPos
        }
        if (any(!is.na(mutationInfo))) {
            return(mutationInfo[!is.na(mutationInfo)])
        } else {
            return(NA)
        }
        
        
    } else {
        return (NA)
    }
}



# List mutations
listMutations <- function(seqInput, seqGL, multipleMutation, model) {
    #if( is.na(c(seqInput, seqGL)) ) return(array(NA,4))
    if (is.na(seqInput) | is.na(seqGL)) { return(NA) }
    seqI <- s2c(seqInput)
    seqG <- s2c(seqGL)
    matIGL <- matrix(c(seqI, seqG), ncol=length(seqI), nrow=2, byrow=T)
    mutations <- analyzeMutations2NucUri(matIGL)
    mutations <- mutations[!is.na(mutations)]
    #positions <- as.numeric(names(mutations))
    # mutations <- mutations[positions <= VLENGTH]
    
    #remove the nucleotide mutations in the codons with multiple mutations
    if (multipleMutation == "ignore") {
       mutationCodons <- getCodonNumb(as.numeric(names(mutations)))
       tableMutationCodons <- table(mutationCodons)
       codonsWithMultipleMutations <- as.numeric(names(tableMutationCodons[tableMutationCodons>1]))
       mutations <- mutations[!(mutationCodons %in% codonsWithMultipleMutations)]
    }
    if (model == "s") {
       mutations <- mutations[mutations == "s"]
    }
    if (length(mutations) > 0) {
        return(mutations)
    } else {
        return(NA)
    }
}


# List the numbers of observed mutations
#
# This lists the observed number of mutation.
#
# @param   db  a data.frame of the DB file.
# @param   sequenceColumn  The name of the sequence column.
# @param   germlineColumn  The name of the germline column.
# 
# @return  list of mutations in each clone
listObservedMutations <- function(db, sequenceColumn="sequence_alignment", 
                                  germlineColumn="germline_alignment_d_mask",
                                  multipleMutation=c("independent", "ignore"),
                                  model = c("rs", "s"))  {
    
    # Make sure the columns specified exist 
    if (!(sequenceColumn %in% names(db))) {
        stop("The sequence column", sequenceColumn, "was not found.")
    } 
    if (!(germlineColumn %in% names(db))) {
        stop("The germline column", germlineColumn, "was not found.")
    } 
    
    mutations <- mapply(listMutations, db[[sequenceColumn]], db[[germlineColumn]], 
                        multipleMutation, model, USE.NAMES=FALSE, SIMPLIFY = F)
    return(mutations)
}


#### Testing functions ####

# Function to make dummy data for testing targetting functions
# 
# @param   nseq  number of sequences
# @param   nmut  number of mutations per sequence
# @param   nmer  number of 5-mers per sequence (sequence length = 5 * nmer)
#
# @return  a data.frame with columns sequence_id, sequence_alignment, germline_alignment_d_mask, v_call.
# 
# @examples
# db <- makeTargetingTestDb(500)
makeTargetingTestDb <- function(nseq=10, nmut=40, nmers=50) {
    nuc_chars <- c("A", "C", "G", "T")
    
    .mut <- function(x, n) {
       i <- sample(1:nchar(x), n) 
       y <- seqinr::s2c(x)
       y[i] <- sapply(y[i], function(z) sample(nuc_chars[nuc_chars != z], 1))
       return(seqinr::c2s(y))
    }
    
    seq <- apply(replicate(nseq, sample(seqinr::words(5, nuc_chars), nmers)), 2, paste, collapse="")
    germ <- sapply(seq, .mut, n=nmut)
    db <- data.frame(sequence_id=paste0("SEQ", 1:nseq),
                     sequence_alignment=seq,
                     germline_alignment_d_mask=germ,
                     v_call="Homsap IGHV3-66*02 F", stringsAsFactors=FALSE)
    rownames(db) <- NULL
    
    return(db)
}
