# Selection analysis using BASELINe

#' @include RegionDefinitions.R
#' @include Shazam.R
NULL

#### Classes ####

#' S4 class defining a BASELINe (selection) object
#' 
#' \code{Baseline} defines a common data structure the results of selection
#' analysis using the BASELINe method.
#' 
#' @slot    description         \code{character} providing general information regarding the 
#'                              sequences, selection analysis and/or object.
#' @slot    db                  \code{data.frame} containing annotation information about 
#'                              the sequences and selection results.
#' @slot    regionDefinition    \link{RegionDefinition} object defining the regions
#'                              and boundaries of the Ig sequences.
#' @slot    testStatistic       \code{character} indicating the statistical framework 
#'                              used to test for selection. For example, \code{"local"} or 
#'                              \code{"focused"}.                           
#' @slot    regions             \code{character} vector defining the regions the BASELINe 
#'                              analysis was carried out on. For \code{"cdr"} and \code{"fwr"} 
#'                              or \code{"cdr1"}, \code{"cdr2"}, \code{"cdr3"}, etc.
#' @slot    numbOfSeqs          \code{matrix} of dimensions \code{r x c} containing the number of 
#'                              sequences or PDFs in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @slot    binomK              \code{matrix} of dimensions \code{r x c} containing the number of 
#'                              successes in the binomial trials in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @slot    binomN              \code{matrix} of dimensions \code{r x c} containing the total 
#'                              number of trials in the binomial in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @slot    binomP              \code{matrix} of dimensions \code{r x c} containing the probability 
#'                              of success in one binomial trial in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @slot    pdfs                \code{list} of matrices containing PDFs with one item for each 
#'                              defined region (e.g. \code{cdr} and \code{fwr}). Matrices have dimensions
#'                              \code{r x c} dementions, where:\cr
#'                              \code{r} = number of rows = number of sequences or groups. \cr
#'                              \code{c} = number of columns = length of the PDF (default 4001).
#' @slot    stats               \code{data.frame} of BASELINe statistics, 
#'                              including: mean selection strength (mean Sigma), 95\% confidence 
#'                              intervals, and p-values with positive signs for the presence of 
#'                              positive selection and/or p-values with negative signs for the
#'                              presence of negative selection.
#'                          
#' @name         Baseline-class
#' @rdname       Baseline-class
#' @aliases      Baseline
#' @exportClass  Baseline
#' @seealso      See \link{summarizeBaseline} for more information on \code{@stats}.
setClass("Baseline", 
         slots=c(description="character",
                 db="data.frame",
                 regionDefinition="RegionDefinition",
                 testStatistic="character",
                 regions="character",
                 numbOfSeqs="matrix",
                 binomK="matrix",
                 binomN="matrix",
                 binomP="matrix",
                 pdfs="list",
                 stats="data.frame"))


#### Methods #####

#' @param    x    \code{Baseline} object.
#' @param    y    name of the column in the \code{db} slot of \code{baseline} 
#'                containing primary identifiers.
#' @param    ...  arguments to pass to \link{plotBaselineDensity}.
#' 
#' @rdname   Baseline-class
#' @aliases  Baseline-method
#' @export
setMethod("plot", c(x="Baseline", y="character"),
          function(x, y, ...) { plotBaselineDensity(x, y, ...) })


#' @param    object  \code{Baseline} object.
#' @param    nproc   number of cores to distribute the operation over.
#' 
#' @rdname   Baseline-class
#' @aliases  Baseline-method
#' @export
setMethod("summary", c(object="Baseline", nproc=integer()),
          function(object, nproc=1) { summarizeBaseline(object, returnType="df", nproc=nproc) })


#### Accessory functions #####

#' Creates a Baseline object
#' 
#' \code{createBaseline} creates and initialize a \code{Baseline} object. 
#' 
#' @param   description         \code{character} providing general information regarding the 
#'                              sequences, selection analysis and/or object.
#' @param   db                  \code{data.frame} containing annotation information about 
#'                              the sequences and selection results.
#' @param   regionDefinition    \link{RegionDefinition} object defining the regions
#'                              and boundaries of the Ig sequences.
#' @param   testStatistic       \code{character} indicating the statistical framework 
#'                              used to test for selection. For example, \code{"local"} or 
#'                              \code{"focused"} or \code{"imbalanced"}.                           
#' @param   regions             \code{character} vector defining the regions the BASELINe 
#'                              analysis was carried out on. For \code{"cdr"} and \code{"fwr"} 
#'                              or \code{"cdr1"}, \code{"cdr2"}, \code{"cdr3"}, etc. If \code{NULL}
#'                              then regions will be determined automatically from \code{regionDefinition}.
#' @param   numbOfSeqs          \code{matrix} of dimensions \code{r x c} containing the number of 
#'                              sequences or PDFs in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @param   binomK              \code{matrix} of dimensions \code{r x c} containing the number of 
#'                              successes in the binomial trials in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @param   binomN              \code{matrix} of dimensions \code{r x c} containing the total 
#'                              number of trials in the binomial in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @param   binomP              \code{matrix} of dimensions \code{r x c} containing the probability 
#'                              of success in one binomial trial in each region, where:\cr
#'                              \code{r} = number of rows = number of groups or sequences.\cr
#'                              \code{c} = number of columns = number of regions.
#' @param   pdfs                \code{list} of matrices containing PDFs with one item for each 
#'                              defined region (e.g. \code{cdr} and \code{fwr}). Matrices have dimensions
#'                              \code{r x c} dementions, where:\cr
#'                              \code{r} = number of rows = number of sequences or groups. \cr
#'                              \code{c} = number of columns = length of the PDF (default 4001).
#' @param   stats               \code{data.frame} of BASELINe statistics, 
#'                              including: mean selection strength (mean Sigma), 95\% confidence 
#'                              intervals, and p-values with positive signs for the presence of 
#'                              positive selection and/or p-values with negative signs for the
#'                              presence of negative selection.
#'                              
#' @return   A \code{Baseline} object.
#' 
#' @details
#' Create and initialize a \code{Baseline} object. 
#' 
#' The \code{testStatistic} indicates the statistical framework used to test for selection. 
#' For example,
#' \itemize{
#'   \item   \code{local} = CDR_R / (CDR_R + CDR_S).
#'   \item   \code{focused} = CDR_R / (CDR_R + CDR_S + FWR_S).
#'   \item   \code{immbalance} = CDR_R + CDR_s / (CDR_R + CDR_S + FWR_S + FWR_R)
#' }
#' For \code{focused} the \code{regionDefinition} must only contain two regions. If more 
#' than two regions are defined, then the \code{local} test statistic will be used.
#' For further information on the frame of these tests see Uduman et al. (2011).
#' 
#' @seealso  See \link{Baseline} for the return object.
#' 
#' @references
#' \enumerate{
#'   \item  Hershberg U, et al. Improved methods for detecting selection by mutation 
#'            analysis of Ig V region sequences. 
#'            Int Immunol. 2008 20(5):683-94.
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data.
#'            Front Immunol. 2013 4(November):358.
#'  }
#'  
#' @examples
#' # Creates an empty Baseline object
#' createBaseline()
#' 
#' @export
createBaseline <- function(description="",
                           db=data.frame(),
                           regionDefinition=createRegionDefinition(),
                           testStatistic="",
                           regions=NULL,
                           numbOfSeqs=matrix(),
                           binomK=matrix(),
                           binomN=matrix(),
                           binomP=matrix(),
                           pdfs=list(),
                           stats=data.frame()) {
    
    if (is.null(regionDefinition)) {
        regionDefinition <- makeNullRegionDefinition()
    }
    # Get regions if not passing in               
    if (is.null(regions)) {
        regions <- regionDefinition@regions
    }
    # Define empty stats data.frame if not passed in
    if (nrow(stats) == 0) {
        stats <- data.frame(group=character(),
                            region=character(),
                            baseline_sigma=character(),
                            baseline_ci_lower=character(),
                            baseline_ci_upper=character(),
                            baseline_ci_pvalue=character(),
                            stringsAsFactors=FALSE) 
    }
    
    # Define RegionDefinition object
    baseline <- new("Baseline",
                    description=description,
                    db=as.data.frame(db),
                    regionDefinition=regionDefinition,
                    testStatistic=testStatistic,
                    regions=regionDefinition@regions,
                    numbOfSeqs=numbOfSeqs,
                    binomK=binomK,
                    binomN=binomN,
                    binomP=binomP,
                    pdfs=pdfs,
                    stats=as.data.frame(stats))
    
    return(baseline)
}


#' Edit the Baseline object
#' 
#' \code{editBaseline} edits a field in a \code{Baseline} object.
#'
#' @param   baseline  \code{Baseline} object to be edited.
#' @param   field     name of the field in the \code{Baseline} object to be edited.
#' @param   value     value to set the \code{field}.
#' 
#' @return   A \code{Baseline} object with the field of choice updated.
#' 
#' @seealso  See \link{Baseline} for the input and return object.
#'
#' @examples
#' \donttest{
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHG" & sample_id == "+7d")
#' 
#' # Make Baseline object
#' baseline <- calcBaseline(db, 
#'                          sequenceColumn="sequence_alignment",
#'                          germlineColumn="germline_alignment_d_mask", 
#'                          testStatistic="focused",
#'                          regionDefinition=IMGT_V,
#'                          targetingModel=HH_S5F,
#'                          nproc=1)
#'                          
#' # Edit the field "description"
#' baseline <- editBaseline(baseline, field="description", 
#'                          value="+7d IGHG")
#' }
#' 
#' @export
editBaseline <- function(baseline, field, value) {
    if (!match(field, slotNames(baseline))) { 
        stop(field, " is not part of the Baseline object.")
    }
    slot(baseline, field) <- value
    
    return(baseline)
}


#### Calculation functions ####


# Helper function for calcBaseline
#
# @param   observed
# @param   expected
# @param   region
# @param   testStatistic
# @param   regionDefinition
# 
# @return  A modified \link{Baseline} object with the BASELINe probability 
#          density function calculated for the regions defined in the \code{regionDefinition}.
calcBaselineHelper <- function(observed,
                                expected,
                                region,
                                testStatistic="local",
                                regionDefinition=NULL) {
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    
    if (is.null(regionDefinition)) {
        regions <- makeNullRegionDefinition()@regions
    } else {
        regions <- regionDefinition@regions
    }
    
    # Evaluate argument choices
    testStatistic <- match.arg(testStatistic, c("local", "focused", "imbalanced"))
    
    # If there are more than two regions (e.g. CDR and FWR then you cannot perform the focused test)
    if (testStatistic=="focused" & length(regions)!=2) {
        testStatistic="local"    
    }    
    
    # local test statistic
    if (testStatistic == "local") { 
        obsX_Index <- grep( paste0("mu_count_", region,"_r"),  names(observed) )
        # important to have "_" after region
        # otherwise this might happen (leading to bugs in results):
        # region = codon_1
        # expect grep to find only codon_1_S and codon_1_R
        # in fact, however, codon_10_S, codon_10_R, codon_101_S, codon_101_R are matched
        obsN_Index <- grep( paste0("mu_count_", region, "_"),  names(observed) )
        
        expX_Index <- grep( paste0("mu_expected_", region,"_r"),  names(expected) )
        # important to have "_" after region
        expN_Index <- grep( paste0("mu_expected_", region, "_"),  names(expected) )       
    }
    
    # focused test statistic
    if (testStatistic == "focused") { 
        obsX_Index <- grep( paste0("mu_count_", region,"_r"),  names(observed) )
        obsN_Index <- 
            grep( 
                paste0( 
                    "mu_count_", region, "|", 
                    "mu_count_", regions[regions!=region], "_s"
                ),
                names(observed) 
            )
        
        expX_Index <- grep( paste0("mu_expected_", region,"_r"),  names(expected) )
        expN_Index <- 
            grep( 
                paste0( 
                    "mu_expected_", region, "|", 
                    "mu_expected_",  regions[regions!=region], "_s"
                ),
                names(expected) 
            )        
    }     
    
    # imbalanced test statistic
    if (testStatistic == "imbalanced") { 
        obsX_Index <- grep( paste0("mu_count_", region),  names(observed) )
        obsN_Index <- grep( "mu_count_",names(observed))  
        
        expX_Index <- grep( paste0("mu_expected_", region),  names(expected) )
        expN_Index <-   grep( "mu_expected_",names(expected)) 
        
    }
    
    
    obsX <- sum(as.numeric( observed[obsX_Index] ))
    obsN <- sum(as.numeric(observed[obsN_Index]), na.rm=T )
    
    expP <-
        as.numeric( 
            sum(expected[expX_Index]) / 
                sum( expected[expN_Index], na.rm=T )
        )
    
    
    return( c( calcBaselineBinomialPdf( x=obsX, n=obsN, p=expP ), obsX, obsN, expP ) )
}

# Calculate the BASELINe probability function in a binomial framework.
calcBaselineBinomialPdf <- function (x=3, 
                                     n=10, 
                                     p=0.33,
                                     CONST_i=CONST_I,
                                     max_sigma=20,
                                     length_sigma=4001) {
    if(n!=0){
        sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)
        sigma_1<-log({CONST_i/{1-CONST_i}}/{p/{1-p}})
        index<-min(n,60)

        y <- dbeta(CONST_i, x+BAYESIAN_FITTED[index], n+BAYESIAN_FITTED[index]-x)*(1-p)*p*exp(sigma_1)/({1-p}^2+2*p*{1-p}*exp(sigma_1)+{p^2}*exp(2*sigma_1))

        if (!sum(is.na(y))) {
            tmp <- approx(sigma_1, y, sigma_s)$y
            return(tmp / sum(tmp) / (2 * max_sigma / (length_sigma - 1)))
        } else {
            return(NA)
        }
    } else {
        return(NA)
    }
}


#' Group BASELINe PDFs
#' 
#' \code{groupBaseline} convolves groups of BASELINe posterior probability density 
#' functions (PDFs) to get combined PDFs for each group.
#'
#' @param    baseline   \code{Baseline} object containing the \code{db} and the 
#'                      BASELINe posterior probability density functions 
#'                      (PDF) for each of the sequences, as returned by
#'                      \link{calcBaseline}.
#' @param    groupBy    The columns in the \code{db} slot of the \code{Baseline}
#'                      object by which to group the sequence PDFs.
#' @param    nproc      number of cores to distribute the operation over. If 
#'                      \code{nproc} = 0 then the \code{cluster} has already been
#'                      set and will not be reset.
#' 
#' @return   A \link{Baseline} object, containing the modified \code{db} and the BASELINe 
#'           posterior probability density functions (PDF) for each of the groups.
#'           
#' @details
#' While the selection strengths predicted by BASELINe perform well on average, 
#' the estimates for individual sequences can be highly variable, especially when the 
#' number of mutations is small. 
#' 
#' To overcome this, PDFs from sequences grouped by biological or experimental relevance,
#' are convolved to from a single PDF for the selection strength. For example, sequences
#' from each sample may be combined together, allowing you to compare selection  across 
#' samples. This is accomplished through a fast numerical convolution technique.
#'               
#' @seealso  To generate the \link{Baseline} object see \link{calcBaseline}.
#'           To calculate BASELINe statistics, such as the mean selection strength
#'           and the 95\% confidence interval, see \link{summarizeBaseline}.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134. 
#'            (Corrections at http://selection.med.yale.edu/baseline/correction/)
#'  }
#' 
#' @examples  
#' \donttest{
#' # Subset example data from alakazam
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call %in% c("IGHM", "IGHG"))
#' 
#' # Collapse clones
#' db <- collapseClones(db, cloneColumn="clone_id",
#'                      sequenceColumn="sequence_alignment",
#'                      germlineColumn="germline_alignment_d_mask",
#'                      method="thresholdedFreq", minimumFrequency=0.6,
#'                      includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#'
#' # Calculate BASELINe
#' baseline <- calcBaseline(db, 
#'                          sequenceColumn="clonal_sequence",
#'                          germlineColumn="clonal_germline", 
#'                          testStatistic="focused",
#'                          regionDefinition=IMGT_V,
#'                          targetingModel=HH_S5F,
#'                          nproc=1)
#'                          
#' # Group PDFs by sample
#' grouped1 <- groupBaseline(baseline, groupBy="sample_id")
#' sample_colors <- c("-1h"="steelblue", "+7d"="firebrick")
#' plotBaselineDensity(grouped1, idColumn="sample_id", colorValues=sample_colors, 
#'                     sigmaLimits=c(-1, 1))
#'  
#' # Group PDFs by both sample (between variable) and isotype (within variable)
#' grouped2 <- groupBaseline(baseline, groupBy=c("sample_id", "c_call"))
#' isotype_colors <- c("IGHM"="darkorchid", "IGHD"="firebrick", 
#'                     "IGHG"="seagreen", "IGHA"="steelblue")
#' plotBaselineDensity(grouped2, idColumn="sample_id", groupColumn="c_call",
#'                     colorElement="group", colorValues=isotype_colors,
#'                     sigmaLimits=c(-1, 1))
#' 
#' # Collapse previous isotype (within variable) grouped PDFs into sample PDFs
#' grouped3 <- groupBaseline(grouped2, groupBy="sample_id")
#' sample_colors <- c("-1h"="steelblue", "+7d"="firebrick")
#' plotBaselineDensity(grouped3, idColumn="sample_id", colorValues=sample_colors,
#'                     sigmaLimits=c(-1, 1))
#' }
#' @export
groupBaseline <- function(baseline, groupBy, nproc=1) {
    # Hack for visibility of foreach index variables
    i <- NULL
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, cpuCount())
    
    # Get indices of unique combinations of field(s) specified by groupBy
    # unique groups
    # crucial to use data.frame and assign colnames (esp. when groupBy has length 1)
    uniqueGroups <- data.frame(unique(baseline@db[, groupBy]))
    colnames(uniqueGroups) <- groupBy 
    rownames(uniqueGroups) <- NULL
    # indices
    # crucial to have simplify=FALSE 
    # (otherwise won't return a list if uniqueClones has length 1)
    uniqueGroupsIdx <- sapply(1:nrow(uniqueGroups), function(i){
        curGroup <- data.frame(uniqueGroups[i, ])
        colnames(curGroup) <- groupBy
        # match for each field
        curIdx <- sapply(groupBy, function(coln){
            baseline@db[, coln]==curGroup[, coln]
        }, simplify=FALSE)
        curIdx <- do.call(rbind, curIdx)
        # intersect to get match across fields 
        curIdx <- which(colSums(curIdx)==length(groupBy))
    }, simplify=FALSE)
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    baseline@db <- data.frame()
    gc()
    if (nproc > 1){        
        cluster <- parallel::makeCluster(nproc, type = "PSOCK")
        parallel::clusterExport( cluster, list('baseline', 'uniqueGroupsIdx',
                                               'break2chunks', 'PowersOfTwo', 
                                               'convolutionPowersOfTwo', 
                                               'convolutionPowersOfTwoByTwos', 
                                               'weighted_conv', 
                                               'calculate_bayesGHelper', 
                                               'groupPosteriors', 'fastConv'), 
                                 envir=environment() )
        registerDoParallel(cluster, cores=nproc)
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    
    # Print status to console
    cat("Grouping BASELINe probability density functions...\n")
    
    # Number of total groups
    numbOfTotalGroups <- length(uniqueGroupsIdx)
    list_pdfs <- list()
    regions <- baseline@regions
    
    # Initialize numbOfSeqs
    # This holds the number of non NA sequences
    numbOfSeqs <- matrix(NA, 
                         ncol=length(baseline@regions), 
                         nrow=numbOfTotalGroups,
                         dimnames=list(1:numbOfTotalGroups, regions))
    templateBinom <- numbOfSeqs
    
    # For every region (e.g. CDR, FWR etc.)
    for (region in regions) {
        
        # Group (convolute) all the PDFS and get one single PDF
        list_region_pdfs  <-
            foreach(i=1:numbOfTotalGroups) %dopar% {
                idx <- uniqueGroupsIdx[[i]]
                
                # Get a matrix (r=numb of sequences/groups * c=4001(i,e. the length of the PDFS))
                # Care was taken to make sure that @pdfs[[region]] should be maintained
                # as a matrix regardless of the number of input sequences (even for a
                # single-sequence input)
                # Thus matrix_GroupPdfs should be expected to be maintained as a matrix as
                # opposed a numeric vector
                matrix_GroupPdfs <- (baseline@pdfs[[region]])[idx, , drop=FALSE]
                stopifnot(is(matrix_GroupPdfs, "matrix"))
                
                # A list version of 
                list_GroupPdfs <- 
                    lapply( 1:nrow(matrix_GroupPdfs), 
                            function(rowIndex) {
                                rowVals <- matrix_GroupPdfs[rowIndex, ]
                                if( !all(is.na(rowVals)) ) { matrix_GroupPdfs[rowIndex, ] }
                            })
                rm(matrix_GroupPdfs)
                gc()
                # Determine the number of sequences that went into creating each of the PDFs
                # If running groupBaseline for the first time after calcBaseline, then
                # each PDF should have a numbOfSeqs=1. 
                numbOfSeqs_region <- baseline@numbOfSeqs[idx, region]
                numbOfSeqs_region <- numbOfSeqs_region[numbOfSeqs_region > 0]
                if(any(numbOfSeqs_region>0)) { 
                    names(numbOfSeqs_region) <- 1:length(numbOfSeqs_region) 
                }
                
                list_GroupPdfs <- list_GroupPdfs[!unlist(lapply(list_GroupPdfs, function(x) { any(is.na(x)) }))]
                list_GroupPdfs <- Filter(Negate(is.null), list_GroupPdfs)
                numbOfNonNASeqs <- length(list_GroupPdfs)
                
                # If all the PDFs in the group are NAs, return a PDF of NAs
                if (length(list_GroupPdfs) == 0) {
                    return(c(rep(NA, 4001), 0))
                }
                
                # If all the PDFs in the group have a numbOfSeqs=1 then
                # call groupPosteriors, which groups PDFs with equal weight
                if (sum(numbOfSeqs_region) == length(numbOfSeqs_region)) { 
                    return(c(groupPosteriors(list_GroupPdfs), numbOfNonNASeqs ) )
                }
                
                # If all the PDFs in the group different numbOfSeqs then call 
                # combineWeightedPosteriors, which groups PDFs weighted by the number of seqs/PDFs
                # that went into creating those PDFs
                if (sum(numbOfSeqs_region) > length(numbOfSeqs_region)) {
                    
                    # sort by number of items
                    len_numbOfSeqs_region <- length(numbOfSeqs_region)
                    sorted_numbOfSeqs_region <- sort(numbOfSeqs_region)
                    rm(numbOfSeqs_region)
                    gc()
                    sorted_list_GroupPdfs <- list()
                    for(newIndex in 1:len_numbOfSeqs_region){
                        sorted_list_GroupPdfs[[newIndex]] <-  list_GroupPdfs[[ as.numeric(names(sorted_numbOfSeqs_region)[newIndex]) ]]
                    }
                    
                    # Group all the PDFs that are created with the equal numbers of seqs/PDFs (i.e. of equal weight)                  
                    repeat {
                        # Count the numb of PDFs with the same weights
                        table_sorted_numbOfSeqs_region <- table(sorted_numbOfSeqs_region)
                        # Weight of interest (the first in the list)
                        pdfWeight <- names(table_sorted_numbOfSeqs_region[table_sorted_numbOfSeqs_region>1])[1]
                        if(is.na(pdfWeight)) { 
                            break
                        }
                        # The corresponding idexes of these PDFs with the same weight
                        indexesOfWeight <- which(sorted_numbOfSeqs_region==pdfWeight)
                        # Convolute these PDFs together
                        list_sameWeightPdfs <- sorted_list_GroupPdfs[indexesOfWeight]
                        updatedPdf <- groupPosteriors(list_sameWeightPdfs)
                        rm(list_sameWeightPdfs)
                        gc()
                        # The new updated weights for this convoluted PDF
                        updatedWeight <- as.numeric(pdfWeight) * length(indexesOfWeight)
                        
                        # remove these from sorted_numbOfSeqs_region & sorted_list_GroupPdfs
                        sorted_numbOfSeqs_region  <- sorted_numbOfSeqs_region[-indexesOfWeight]
                        sorted_list_GroupPdfs <- sorted_list_GroupPdfs[-indexesOfWeight]
                        rm(indexesOfWeight)
                        gc()
                        
                        # add the convoluted PDF and its new weight
                        newLength <- length(sorted_numbOfSeqs_region)+1
                        sorted_numbOfSeqs_region[newLength] <- updatedWeight
                        sorted_list_GroupPdfs[[newLength]] <- updatedPdf
                        rm(updatedWeight)
                        rm(updatedPdf)
                        gc()
                        
                        # sort by number of items
                        len_sorted_numbOfSeqs_region <- length(sorted_numbOfSeqs_region)
                        sorted_numbOfSeqs_region <- sort(sorted_numbOfSeqs_region)
                        names(sorted_numbOfSeqs_region) <- as.character(1:len_sorted_numbOfSeqs_region)
                        list_GroupPdfs <- sorted_list_GroupPdfs
                        sorted_list_GroupPdfs <- list()
                        for(newIndex in 1:len_numbOfSeqs_region){
                            sorted_list_GroupPdfs[[newIndex]] <-  list_GroupPdfs[[ as.numeric(names(sorted_numbOfSeqs_region)[newIndex]) ]]
                        }
                        
                        table_sorted_numbOfSeqs_region <- table(sorted_numbOfSeqs_region)
                        
                        if(sum(table_sorted_numbOfSeqs_region>1)>0){
                            break
                        }
                    }
                    
                    #return( c( groupPosteriors(sorted_list_GroupPdfs), 10 ) )
                    
                    # Do pairwise grouping of PDFs based on weight
                    # 1. sort by weights
                    # 2. group the lowest two weighted PDFs
                    # 3. resort, and repeat till you get one PDFs
                    if(length(list_GroupPdfs)>1){
                        repeat{
                            
                            updatedPdf <- combineWeightedPosteriors(list_GroupPdfs[[1]], 
                                                                    sorted_numbOfSeqs_region[1], 
                                                                    list_GroupPdfs[[2]], 
                                                                    sorted_numbOfSeqs_region[2])
                            updatedWeight <- sorted_numbOfSeqs_region[1] + sorted_numbOfSeqs_region[2]
                            # remove these from sorted_numbOfSeqs_region & sorted_list_GroupPdfs
                            sorted_numbOfSeqs_region  <- sorted_numbOfSeqs_region[-c(1,2)]
                            sorted_list_GroupPdfs <- sorted_list_GroupPdfs[-c(1,2)]
                            
                            # add the convoluted PDF and its new weight
                            newLength <- length(sorted_numbOfSeqs_region)+1
                            sorted_numbOfSeqs_region[newLength] <- updatedWeight
                            rm(updatedWeight)
                            sorted_list_GroupPdfs[[newLength]] <- updatedPdf
                            rm(updatedPdf)
                            gc()
                            # sort by number of items
                            len_sorted_numbOfSeqs_region <- length(sorted_numbOfSeqs_region)
                            sorted_numbOfSeqs_region <- sort(sorted_numbOfSeqs_region)
                            names(sorted_numbOfSeqs_region) <- as.character(1:len_sorted_numbOfSeqs_region)
                            list_GroupPdfs <- sorted_list_GroupPdfs
                            sorted_list_GroupPdfs <- list()
                            for(newIndex in 1:len_numbOfSeqs_region){
                                sorted_list_GroupPdfs[[newIndex]] <-  list_GroupPdfs[[ as.numeric(names(sorted_numbOfSeqs_region)[newIndex]) ]]
                            }
                            
                            if(length(list_GroupPdfs)==1){
                                break
                            }
                        }
                    }
                    return( c( list_GroupPdfs[[1]], as.numeric(sorted_numbOfSeqs_region) ) )
                }
                
            }
        
        # Convert the list of the region's PDFs into a matrix                
        matrix_region_pdfs <- 
            do.call(rbind, lapply(list_region_pdfs, 
                                  function(x) { 
                                      length(x) <- 4002 
                                      return(x)
                                  }))
        
        # Normalize and save PDF matrix
        # Hardcode normalization to max_sigma=20 and sigma_length=4001
        pdf_norm <- 2*20 / 4000
        pdf_mat <- matrix_region_pdfs[, 1:4001, drop=FALSE]
        list_pdfs[[region]] <- pdf_mat / rowSums(pdf_mat, na.rm=TRUE) / pdf_norm
        # Save regions
        numbOfSeqs[, region] <- matrix_region_pdfs[, 4002]
    }
    
    #colnames(numbOfSeqs) <- paste0("NUMB_SEQUENCES_", colnames(numbOfSeqs))
    
    # Create the db, which will now contain the group information
    stopifnot(is.data.frame(uniqueGroups))
    db <- uniqueGroups
    
    # Create a Baseline object with the above results to return
    baseline <- createBaseline(description="",
                               db=as.data.frame(db),
                               regionDefinition=baseline@regionDefinition,
                               testStatistic=baseline@testStatistic,
                               regions=regions,
                               numbOfSeqs=numbOfSeqs,
                               binomK=templateBinom,
                               binomN=templateBinom,
                               binomP=templateBinom,
                               pdfs=list_pdfs)
    
    # Calculate BASELINe stats and update slot
    baseline <- summarizeBaseline(baseline)
    
    # Stop cluster
    if(nproc > 1) { parallel::stopCluster(cluster) }
    
    return(baseline)
}


#' Calculate BASELINe summary statistics
#'
#' \code{summarizeBaseline} calculates BASELINe statistics such as the mean selection 
#' strength (mean Sigma), the 95\% confidence intervals and p-values for the presence of
#' selection.
#'
#' @param    baseline    \code{Baseline} object returned by \link{calcBaseline} containing 
#'                       annotations and BASELINe posterior probability density functions 
#'                       (PDFs) for each sequence.
#' @param    returnType  One of \code{c("baseline", "df")} defining whether
#'                       to return a \code{Baseline} object ("baseline") with an updated
#'                       \code{stats} slot or a data.frame ("df") of summary statistics.
#' @param    nproc       number of cores to distribute the operation over. If 
#'                       \code{nproc} = 0 then the \code{cluster} has already been
#'                       set and will not be reset.
#' 
#' @return   Either a modified \code{Baseline} object or data.frame containing the 
#'           mean BASELINe selection strength, its 95\% confidence intervals, and 
#'           a p-value for the presence of selection.
#' 
#' @details  The returned p-value can be either positive or negative. Its magnitude 
#'           (without the sign) should be interpreted as per normal. Its sign indicates 
#'           the direction of the selection detected. A positive p-value indicates positive
#'           selection, whereas a negative p-value indicates negative selection.
#'                     
#' @seealso  See \link{calcBaseline} for generating \code{Baseline} objects and
#'           \link{groupBaseline} for convolving groups of BASELINe PDFs.
#'           
#' @references
#' \enumerate{
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
#' }
#' 
#' @examples
#' \donttest{
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHG")
#' 
#' # Collapse clones
#' db <- collapseClones(db, cloneColumn="clone_id",
#'                      sequenceColumn="sequence_alignment",
#'                      germlineColumn="germline_alignment_d_mask",
#'                      method="thresholdedFreq", minimumFrequency=0.6,
#'                      includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#'                      
#' # Calculate BASELINe
#' baseline <- calcBaseline(db, 
#'                          sequenceColumn="clonal_sequence",
#'                          germlineColumn="clonal_germline", 
#'                          testStatistic="focused",
#'                          regionDefinition=IMGT_V,
#'                          targetingModel=HH_S5F,
#'                          nproc = 1)
#' 
#' # Grouping the PDFs by the sample annotation
#' grouped <- groupBaseline(baseline, groupBy="sample_id")
#' 
#' # Get a data.frame of the summary statistics
#' stats <- summarizeBaseline(grouped, returnType="df")
#' }                     
#' @export
summarizeBaseline <- function(baseline, returnType=c("baseline", "df"), nproc=1) {
    # Hack for visibility of foreach index variable
    idx <- NULL
    
    # Check arguments
    returnType <- match.arg(returnType)
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, cpuCount())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if (nproc > 1){        
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        parallel::clusterExport(cluster, list('baseline',
                                              'baselineSigma',
                                              'baselineCI',
                                              'baselinePValue'), 
                                envir=environment())
        registerDoParallel(cluster, cores=nproc)
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Printing status to console
    cat("Calculating BASELINe statistics...\n")
    
    # Calculate stats for each sequence/group
    numbOfTotalSeqs <- nrow(baseline@db)
    regions <- baseline@regions
    db <- baseline@db
    if ("sequence_id" %in% colnames(db)) { 
        db <- subset(db, select="sequence_id") 
    } else if ("SEQUENCE_ID" %in% colnames(db)) {
        db <- subset(db, select="SEQUENCE_ID") 
    }
    list_stats <-
        foreach(idx=iterators::icount(numbOfTotalSeqs)) %dopar% {
            df_baseline_seq <- data.frame()
            db_seq <- data.frame(db[idx, ])
            names(db_seq) <- names(db)
            for (region in regions) {
                
                # care was taken to make sure that @pdfs[[region]] should be maintained
                # as a matrix regardless of the number of input sequences (even for a
                # single-sequence input)
                stopifnot(is(baseline@pdfs[[region]], "matrix"))
                baseline_pdf <- baseline@pdfs[[region]][idx, ]
                
                baseline_ci <- baselineCI(baseline_pdf)
                df_baseline_seq_region <- 
                    data.frame(db_seq,
                               region=factor(region, levels=regions),
                               baseline_sigma=baselineSigma(baseline_pdf),
                               baseline_ci_lower=baseline_ci[1],
                               baseline_ci_upper=baseline_ci[2],
                               baseline_ci_pvalue=baselinePValue(baseline_pdf))
                df_baseline_seq <- dplyr::bind_rows(df_baseline_seq, df_baseline_seq_region)
            }
            df_baseline_seq[,1] <- as.vector(unlist(df_baseline_seq[,1]))
            df_baseline_seq[,2] <- as.vector(unlist(df_baseline_seq[,2]))
            return(df_baseline_seq)
        }
    
    # Stop cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    # Convert list of BASELINe stats into a data.frame
    stats <- as.data.frame(dplyr::bind_rows(list_stats))
    
    if (returnType == "df") {
        return(stats)    
    } else if (returnType == "baseline") {
        # Append stats to baseline object
        return(editBaseline(baseline, field="stats", stats))
    } else {
        return(NULL)
    }
}


#' Two-sided test of BASELINe PDFs
#' 
#' \code{testBaseline} performs a two-sample signifance test of BASELINe 
#' posterior probability density functions (PDFs).
#'
#' @param    baseline   \code{Baseline} object containing the \code{db} and grouped 
#'                      BASELINe PDFs returned by \link{groupBaseline}.
#' @param    groupBy    string defining the column in the \code{db} slot of the 
#'                      \code{Baseline} containing sequence or group identifiers.
#' 
#' @return   A data.frame with test results containing the following columns:
#'           \itemize{
#'             \item  \code{region}:  sequence region, such as \code{cdr} and \code{fwr}.
#'             \item  \code{test}:    string defining the groups be compared. The
#'                                    string is formated as the conclusion associated with the
#'                                    p-value in the form \code{GROUP1 != GROUP2}. Meaning,
#'                                    the p-value for rejection of the null hypothesis that 
#'                                    GROUP1 and GROUP2 have equivalent distributions.
#'             \item  \code{pvalue}:  two-sided p-value for the comparison.
#'             \item  \code{fdr}:     FDR corrected \code{pvalue}.
#'           }
#'           
#' @seealso  To generate the \link{Baseline} input object see \link{groupBaseline}.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134. 
#'            (Corretions at http://selection.med.yale.edu/baseline/correction/)
#'  }
#' 
#' @examples
#' \donttest{
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call %in% c("IGHM", "IGHG", "IGHA"))
#'
#' # Collapse clones
#' db <- collapseClones(db, cloneColumn="clone_id",
#'                      sequenceColumn="sequence_alignment",
#'                      germlineColumn="germline_alignment_d_mask",
#'                      method="thresholdedFreq", minimumFrequency=0.6,
#'                      includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#'                      
#' # Calculate BASELINe
#' baseline <- calcBaseline(db, 
#'                          sequenceColumn="clonal_sequence",
#'                          germlineColumn="clonal_germline", 
#'                          testStatistic="focused",
#'                          regionDefinition=IMGT_V,
#'                          targetingModel=HH_S5F,
#'                          nproc=1)
#' 
#' # Group PDFs by the isotype
#' grouped <- groupBaseline(baseline, groupBy="c_call")
#' 
#' # Visualize isotype PDFs
#' plot(grouped, "c_call")
#' 
#' # Perform test on isotype PDFs
#' testBaseline(grouped, groupBy="c_call")
#' }
#' @export
testBaseline <- function(baseline, groupBy) {
    ## DEBUG
    # baseline=grouped; groupBy="sample_id"
    
    # Get test groups
    groups <- as.character(baseline@db[[groupBy]])
    if (length(groups) < 2) {
        stop('The ', groupBy, ' column does not contain at least two groups.')
    }
    pair_indices <- combn(1:length(groups), 2, simplify=F)
    pair_names <- combn(groups, 2, simplify=F)
    test_names <- sapply(pair_names, paste, collapse=" != ")
    
    # Run tests
    test_list <- list()
    for (n in baseline@regions) {
        d <- baseline@pdfs[[n]]
        p <- sapply(pair_indices, function(x) { baseline2DistPValue(d[x[1], ], d[x[2], ])})
        test_list[[n]] <- data.frame(test=test_names, pvalue=p)
    }
    test_df <- bind_rows(test_list, .id="region")
    test_df$fdr <- p.adjust(test_df$pvalue, method="fdr")
    
    return(test_df)
}
    

# Calculate mean sigma of a BASELINe PDF
#
# @param   base          BASLINe PDF vector.
# @param   max_sigma     maximum sigma score.
# @param   length_sigma  length of the PDF vector.
# @return  Mean sigma.
baselineSigma <- function(base, max_sigma=20, length_sigma=4001) {
    # Return NA on bad input
    if (any(is.na(base))) { 
        return(NA) 
    }
    
    sigma_s <- seq(-max_sigma, max_sigma, length.out=length_sigma)
    norm <- sum(base, na.rm=TRUE)
    sigma_mean <- base %*% sigma_s / norm
    
    return(sigma_mean)
}


# Calculate confidence interval for BASELINe PDF
#
# @param   base          BASLINe PDF vector.
# @param   low           lower CI percentile.
# @param   up            upper CI percentile.
# @param   max_sigma     maximum sigma score.
# @param   length_sigma  length of the PDF vector.
# @return  A vector of \code{c(lower, upper)} confidence bounds.
baselineCI <- function (base, low=0.025, up=0.975, max_sigma=20, length_sigma=4001){
    # Return NA on bad input
    if (any(is.na(base))) { 
        return(c(NA, NA)) 
    }
    
    sigma_s <- seq(-max_sigma, max_sigma, length.out=length_sigma)
    cdf <- cumsum(base)
    cdf <- cdf / cdf[length(cdf)]
    intervalLow <- findInterval(low, cdf)
    fractionLow <- (low - cdf[intervalLow])/(cdf[intervalLow + 1] - cdf[intervalLow])
    intervalUp <- findInterval(up,cdf)
    fractionUp <- (up - cdf[intervalUp]) / (cdf[intervalUp] - cdf[intervalUp - 1])
    sigmaLow <- sigma_s[intervalLow] + fractionLow*(sigma_s[intervalLow + 1] - sigma_s[intervalLow])
    sigmaUp <- sigma_s[intervalUp] + fractionUp*(sigma_s[intervalUp + 1] - sigma_s[intervalUp])
    
    return(c(sigmaLow, sigmaUp))
}


# Calculate a p-value that the given BASELINe PDF differs from zero
#
# @param   base          BASLINe PDF vector.
# @param   max_sigma     maximum sigma score.
# @param   length_sigma  length of the PDF vector.
# @return  A p-value. The returned p-value can be either positive or negative. 
#          Its magnitude (without the sign) should be interpreted as per normal. 
#          Its sign indicate the direction of the selection detected. A positive 
#          p-value indicates positive selection, whereas a negative p-value 
#          indicates negative selection.
baselinePValue <- function (base, length_sigma=4001, max_sigma=20){
    # note: since there isn't a null distribution, this "p-value" isn't a p-value in the 
    # conventional sense
    if (!any(is.na(base))) {
        # normalization factor
        #norm <- (length_sigma - 1) / 2 / max_sigma
        # sums up to 100 for default setting (sigma_s from -20 to 20 with length 4001)
        norm <- sum(base, na.rm=TRUE) 
        
        # compute Pr(selection strength < 0):
        # sum up density for sigma from min_sigma up to and right before 0 (area under curve), plus
        # + binomial correction (density at sigma=0 divided by 2); 
        # normalized
        pvalue <- ( sum(base[1:((length_sigma - 1) / 2)]) + 
                    base[((length_sigma + 1) / 2)] / 2 ) / norm
        
        # from Fig 4 caption of Detecting selection in immunoglobulin sequences by Uduman et al. 2011
        # "Note that P values less than zero are indicative of negative selection."
        
        # 1) if Pr(selection strength < 0) <= 0.5, return Pr(selection strength < 0)
        #    this will be positive, and serves as the "p-value" for positive selection
        # 2) if Pr(selection strength < 0) > 0.5, return -Pr(selection strength>0)
        #    this will be negative, and serves as the "p-value" for negative selection
        #    (negative sign highlights the fact that selection is negative)
        if (pvalue > 0.5) { pvalue <- -(1 - pvalue) } 
    } else {
        pvalue <- NA
    }
    
    return(pvalue)
}


# Compute p-value of two BASELINe PDFs
#
# @param   base1  first selection PDF; must be a numeric vector
# @param   base2  second selection PDF; must be a numeric vector
# @return  Two-sided p-value that base1 and base2 differ.
baseline2DistPValue <-function(base1, base2) {
    # NOTE: make sure to supply 2 vectors (not 1-row data.frames) when
    #       calling this function directly
    
    ## Debug
    # base1=grouped@pdfs[["CDR"]][1, ]; base2=grouped@pdfs[["FWR"]][1, ]
    
    # Get lengths
    len1 <- length(base1)
    len2 <- length(base2)
    
    # Check input
    if (len1 != len2) {
        stop("base1 and base2 must be the same length.")
    }
    # NA if all values in pdfs are NA
    if (sum(is.na(base1))==len1 | sum(is.na(base2))==len2) {
        return(NA)
    }

    # Determine p-value
    if (len1 > 1) {
        # Normalize
        base1 <- base1 / sum(base1, na.rm=TRUE)
        base2 <- base2 / sum(base2, na.rm=TRUE)
        # Calculate p-value
        cum2 <- cumsum(base2) - base2/2
        pvalue <- sum(base1*cum2)
        if (pvalue > 0.5) { pvalue <- 1 - pvalue }
    } else {
        pvalue <- NA
    }
    
    return(pvalue)
}  


#### Plotting functions ####

#' Plots BASELINe probability density functions
#' 
#' \code{plotBaselineDensity} plots the probability density functions resulting from selection 
#' analysis using the BASELINe method.
#'
#' @param    baseline       \code{Baseline} object containing selection probability 
#'                          density functions.
#' @param    idColumn       name of the column in the \code{db} slot of \code{baseline} 
#'                          containing primary identifiers.
#' @param    groupColumn    name of the column in the \code{db} slot of \code{baseline} 
#'                          containing secondary grouping identifiers. If \code{NULL}, 
#'                          organize the plot only on values in \code{idColumn}.
#' @param    colorElement   one of \code{c("id", "group")} specifying whether the 
#'                          \code{idColumn} or \code{groupColumn} will be used for color coding. 
#'                          The other entry, if present, will be coded by line style.
#' @param    colorValues    named vector of colors for entries in \code{colorElement}, with 
#'                          names defining unique values in the \code{colorElement} column and values
#'                          being colors. Also controls the order in which values appear on the
#'                          plot. If \code{NULL} alphabetical ordering and a default color palette 
#'                          will be used.
#' @param    facetBy        one of \code{c("region", "group")} specifying which category to facet the
#'                          plot by, either values in \code{groupColumn} ("group") or regions
#'                          defined in the \code{regions} slot of the \code{baseline} object ("region").
#'                          If this is set to "group", then the region will behave as the \code{groupColumn}
#'                          for purposes of the \code{colorElement} argument.
#' @param    title          string defining the plot title.
#' @param    subsetRegions  character vector defining a subset of regions to plot, correspoding 
#'                          to the regions for which the \code{baseline} data was calculated. If
#'                          \code{NULL} all regions in \code{baseline} are plotted.
#' @param    sigmaLimits    numeric vector containing two values defining the \code{c(lower, upper)}
#'                          bounds of the selection scores to plot.
#' @param    style          type of plot to draw. One of:
#'                          \itemize{
#'                            \item \code{"density"}:  plots a set of curves for each probability 
#'                                                     density function in \code{baseline}, 
#'                                                     with colors determined by values in the
#'                                                     \code{colorElement} column.
#'                                                     Faceting is determined by the 
#'                                                     \code{facetBy} argument.
#'                          }
#' @param    sizeElement    one of \code{c("none", "id", "group")} specifying whether the lines in the
#'                          plot should be all of the same size (\code{none}) or have their sizes depend on 
#'                          the values in \code{id} or \code{code}.
#' @param    size           numeric scaling factor for lines, points and text in the plot.
#' @param    silent         if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                          object; if \code{FALSE} draw the plot.
#' @param    ...            additional arguments to pass to ggplot2::theme.
#' 
#' @return   A ggplot object defining the plot.
#'
#' @seealso  Takes as input a \link{Baseline} object returned from \link{groupBaseline}.
#' 
#' @examples
#' \donttest{
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call %in% c("IGHM", "IGHG"))
#' 
#' # Collapse clones
#' db <- collapseClones(db, cloneColumn="clone_id",
#'                      sequenceColumn="sequence_alignment",
#'                      germlineColumn="germline_alignment_d_mask",
#'                      method="thresholdedFreq", minimumFrequency=0.6,
#'                      includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#'                      
#' # Calculate BASELINe
#' baseline <- calcBaseline(db, 
#'                          sequenceColumn="clonal_sequence",
#'                          germlineColumn="clonal_germline", 
#'                          testStatistic="focused",
#'                          regionDefinition=IMGT_V,
#'                          targetingModel=HH_S5F,
#'                          nproc=1)
#'  
#' # Grouping the PDFs by the sample and isotype annotations
#' grouped <- groupBaseline(baseline, groupBy=c("sample_id", "c_call"))
#' 
#' # Plot density faceted by region with custom isotype colors
#' isotype_colors <- c("IGHM"="darkorchid", "IGHD"="firebrick", 
#'                     "IGHG"="seagreen", "IGHA"="steelblue")
#' plotBaselineDensity(grouped, "sample_id", "c_call", colorValues=isotype_colors, 
#'                     colorElement="group", sigmaLimits=c(-1, 1))
#'
#' # Facet by isotype instead of region
#' sample_colors <- c("-1h"="steelblue", "+7d"="firebrick")
#' plotBaselineDensity(grouped, "sample_id", "c_call", facetBy="group",
#'                     colorValues=sample_colors, sigmaLimits=c(-1, 1))
#' }
#' 
#' @export
plotBaselineDensity <- function(baseline, idColumn, groupColumn=NULL, colorElement=c("id", "group"), 
                                colorValues=NULL, title=NULL, subsetRegions=NULL, sigmaLimits=c(-5, 5), 
                                facetBy=c("region", "group"), style=c("density"), 
                                sizeElement=c("none", "id", "group"), size=1, 
                                silent=FALSE, ...) {
    ## DEBUG
    # baseline=grouped
    # idColumn="sample_id"; groupColumn="c_call"; subsetRegions=NULL; sigmaLimits=c(-5, 5)
    # facetBy="region"; style="density"; size=1; silent=FALSE
    
    # Check input
    colorElement <- match.arg(colorElement)
    style <- match.arg(style)
    facetBy <- match.arg(facetBy)
    sizeElement <- match.arg(sizeElement)
    
    # Set base plot settings
    base_theme <- theme_bw() +
        theme(panel.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(color="black", size=0.5)) +
        theme(strip.background=element_rect(fill="white", color="black", size=0.5))
    
    if (style == "density") {
        # Check for proper grouping
        if (any(duplicated(baseline@db[, c(idColumn, groupColumn)]))) {
            stop("More than one unique annotation set per summary statistic. Rerun groupBaseline to combine data.")
        }
        
        # Subset to regions of interest
        dens_names <- baseline@regions        
        if (!is.null(subsetRegions)) {
            dens_names <- dens_names[dens_names %in% subsetRegions]
        }
        dens_list <- baseline@pdfs[dens_names]
        
        # Get row and column names for PDF matrices
        group_df <- subset(baseline@db, select=c(idColumn, groupColumn))
        group_df$GROUP_COLLAPSE <- apply(subset(group_df, select=c(idColumn, groupColumn)), 
                                         1, paste, collapse=",")
        col_names <- seq(-20, 20, length.out=ncol(dens_list[[1]]))
        
        # Update column and rownames for PDF matrices and subset to Sigma in -5:5
        for (i in 1:length(dens_list)) {
            rownames(dens_list[[i]]) <- group_df$GROUP_COLLAPSE
            colnames(dens_list[[i]]) <- col_names
            dens_list[[i]] <- dens_list[[i]][, col_names >= sigmaLimits[1] & 
                                               col_names <= sigmaLimits[2],
                                             drop=FALSE]
        }
        
        # Melt density matrices
        melt_list <- list()
        for (n in dens_names) {
            tmp_df <- as.data.frame(dens_list[[n]])
            tmp_df$GROUP_COLLAPSE <- rownames(dens_list[[n]])
            gather_cols <- names(tmp_df)[names(tmp_df) != "GROUP_COLLAPSE"]
            melt_list[[n]] <- tidyr::gather(tmp_df, "SIGMA", "DENSITY", tidyselect::all_of(gather_cols), convert=TRUE)
        }
        dens_df <- dplyr::bind_rows(melt_list, .id="region")
        
        # Assign id and group columns to density data.frame
        dens_df[, idColumn] <- group_df[match(dens_df$GROUP_COLLAPSE, group_df$GROUP_COLLAPSE), 
                                        idColumn]
        if (!is.null(groupColumn)) {
            dens_df[, groupColumn] <- group_df[match(dens_df$GROUP_COLLAPSE, group_df$GROUP_COLLAPSE), 
                                               groupColumn]
        }    
        
        # Set secondary grouping and faceting columns
        if (facetBy == "group") { 
            secondaryColumn <- "region"
            facetColumn <- groupColumn
        } else if (facetBy == "region") {
            secondaryColumn <- groupColumn
            facetColumn <- "region" 
        }
        
        # Apply color order
        if (!is.null(colorValues)) {
            if (colorElement == "id") {
                dens_df[, idColumn] <- factor(dens_df[, idColumn], levels=names(colorValues))
            } else {
                dens_df[, groupColumn] <- factor(dens_df[, groupColumn], levels=names(colorValues))
            }
        }
        
        # Apply line width
        dens_df[, "size"] <- factor(1)
        
        if (sizeElement=="id") {
            dens_df[, "size"] <- factor(dens_df[, idColumn])
        } else if (sizeElement == "group" ) {
            dens_df[, "size"] <- factor(dens_df[, groupColumn])
        }
        
        size_values <- 1:length(levels(dens_df[,"size"]))
        size_names <- levels(dens_df[, "size"])        
        size_values <- size*size_values/max(size_values)
        names(size_values) <- size_names
        
        # Plot probability density curve
        p1 <- ggplot(dens_df, aes_string(x="SIGMA", y="DENSITY")) +
            base_theme + 
            xlab(expression(Sigma)) +
            ylab("Density") +
            geom_line(aes(size=size))
        # Add line
        if (colorElement == "id" & is.null(secondaryColumn)) {
            p1 <- p1 + aes_string(color=idColumn)
        } else if (colorElement == "id" & !is.null(secondaryColumn)) {
            p1 <- p1 + aes_string(color=idColumn, linetype=secondaryColumn)
        } else if (colorElement == "group") {
            p1 <- p1 + aes_string(color=secondaryColumn, linetype=idColumn)
        } else {
            stop("Incompatible arguments for groupColumn, colorElement and facetBy")
        }
        # Add colors
        if (!is.null(colorValues)) {
            p1 <- p1 + scale_color_manual(values=colorValues)
        }   
        # Add title
        if (!is.null(title)) {
            p1 <- p1 + ggtitle(title)
        }
        # Add facet
        if (is.null(facetColumn)) {
            stop("Cannot facet by group if groupColumn=NULL")
        } else {
            p1 <- p1 + facet_grid(paste(facetColumn, "~ ."))
        }
    }
    
    # Add additional theme elements
    p1 <- p1 + 
        scale_size_manual(breaks=names(size_values), values=size_values)
    
    if (sizeElement == "none") {
        p1 <- p1 +
            guides( 
                size="none", 
                colour = guide_legend(override.aes = list(size = size_values))
                )        
        if (length(unique(c(groupColumn,idColumn)))>1) {
            p1 <- p1 +
                guides (
                    linetype = guide_legend(override.aes = list(size = size_values))
                )
        }
    } else if (sizeElement == colorElement) {
        p1 <- p1 +
            guides( 
                size="none", 
                colour = guide_legend(override.aes = list(size = size_values)))
    } else {
        p1 <- p1 +
            guides( 
                size="none", 
                linetype = guide_legend(override.aes = list(size = size_values))
            ) 
    }
    
      
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { 
        plot(p1)
    }
    
    invisible(p1)
}


#' Plots BASELINe summary statistics
#' 
#' \code{plotBaselineSummary} plots a summary of the results of selection analysis 
#' using the BASELINe method.
#'
#' @param    baseline       either a data.frame returned from \link{summarizeBaseline}
#'                          or a \code{Baseline} object returned from \link{groupBaseline}
#'                          containing selection probability density functions and summary 
#'                          statistics.
#' @param    idColumn       name of the column in \code{baseline} containing primary identifiers. 
#'                          If the input is a \code{Baseline} object, then this will be a column
#'                          in the \code{stats} slot of \code{baseline}.
#' @param    groupColumn    name of the column in \code{baseline} containing secondary grouping 
#'                          identifiers. If the input is a \code{Baseline} object, then this will 
#'                          be a column in the \code{stats} slot of \code{baseline}.
#' @param    groupColors    named vector of colors for entries in \code{groupColumn}, with 
#'                          names defining unique values in the \code{groupColumn} and values
#'                          being colors. Also controls the order in which groups appear on the
#'                          plot. If \code{NULL} alphabetical ordering and a default color palette 
#'                          will be used. Has no effect if \code{facetBy="group"}.
#' @param    subsetRegions  character vector defining a subset of regions to plot, correspoding 
#'                          to the regions for which the \code{baseline} data was calculated. If
#'                          \code{NULL} all regions in \code{baseline} are plotted.
#' @param    facetBy        one of c("group", "region") specifying which category to facet the
#'                          plot by, either values in \code{groupColumn} ("group") or regions
#'                          defined in \code{baseline} ("region"). The data that is not used
#'                          for faceting will be color coded.
#' @param    title          string defining the plot title.
#' @param    style          type of plot to draw. One of:
#'                          \itemize{
#'                            \item \code{"summary"}:  plots the mean and confidence interval for
#'                                                     the selection scores of each value in 
#'                                                     \code{idColumn}. Faceting and coloring
#'                                                     are determine by values in \code{groupColumn}
#'                                                     and regions defined in \code{baseline}, 
#'                                                     depending upon the \code{facetBy} argument.
#'                          }
#' @param    size           numeric scaling factor for lines, points and text in the plot.
#' @param    silent         if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                          object; if \code{FALSE} draw the plot.
#' @param    ...            additional arguments to pass to ggplot2::theme.
#' 
#' @return   A ggplot object defining the plot.
#'
#' @seealso  Takes as input either a \link{Baseline} object returned by \link{groupBaseline} 
#'           or a data.frame returned from \link{summarizeBaseline}.
#' 
#' @examples
#' \donttest{
#' # Subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call %in% c("IGHM", "IGHG"))
#' 
#' # Collapse clones
#' db <- collapseClones(db, cloneColumn="clone_id",
#'                      sequenceColumn="sequence_alignment",
#'                      germlineColumn="germline_alignment_d_mask",
#'                      method="thresholdedFreq", minimumFrequency=0.6,
#'                      includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#'                      
#' # Calculate BASELINe
#' baseline <- calcBaseline(db, 
#'                          sequenceColumn="clonal_sequence",
#'                          germlineColumn="clonal_germline", 
#'                          testStatistic="focused",
#'                          regionDefinition=IMGT_V,
#'                          targetingModel=HH_S5F,
#'                          nproc=1)
#'  
#' # Grouping the PDFs by sample and isotype annotations
#' grouped <- groupBaseline(baseline, groupBy=c("sample_id", "c_call"))
#' 
#' # Plot mean and confidence interval by region with custom group colors
#' isotype_colors <- c("IGHM"="darkorchid", "IGHD"="firebrick", 
#'                     "IGHG"="seagreen", "IGHA"="steelblue")
#' plotBaselineSummary(grouped, "sample_id", "c_call", 
#'                     groupColors=isotype_colors)
#' 
#' # Facet by group instead of region
#' plotBaselineSummary(grouped, "sample_id", "c_call", facetBy="group")
#' }
#' 
#' @export
plotBaselineSummary <- function(baseline, idColumn, groupColumn=NULL, groupColors=NULL, 
                                subsetRegions=NULL, facetBy=c("region", "group"), 
                                title=NULL, style=c("summary"), size=1, silent=FALSE, ...) {
    # Check arguments
    style <- match.arg(style)
    facetBy <- match.arg(facetBy)
    
    # Check input object
    if (is(baseline, "Baseline")) {
        stats_df <- baseline@stats
    } else if (is(baseline, "data.frame")) {
        stats_df <- baseline
    } else {
        stop("Input must be either a data.frame or Baseline object.")
    }
    
    # Check for required columns
    baseline_cols <- c("region", 
                       "baseline_sigma", 
                       "baseline_ci_lower", 
                       "baseline_ci_pvalue")
    if (!(all(baseline_cols %in% names(stats_df)))) {
        stop("Input must contain columns defined by summarizeBaseline.")
    }
    
    # Check for proper grouping
    if (any(duplicated(stats_df[, c(idColumn, groupColumn, "region")]))) {
        stop("More than one unique annotation set per summary statistic. Rerun groupBaseline to combine data.")
    }
    
    # Subset to regions of interest
    if (!is.null(subsetRegions)) {
        stats_df <- stats_df[stats_df$region %in% subsetRegions, ]
    }
    
    # Set base plot settings
    base_theme <- theme_bw() +
        theme(panel.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(color="black", size=0.5)) +
        theme(strip.background=element_rect(fill="white", color="black", size=0.5)) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(), 
              axis.ticks.x=element_blank()) +
        theme(legend.position="top") +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
    
    if (style == "summary") { 
        # Plot mean and confidence intervals
        stats_df <- stats_df[!is.na(stats_df$baseline_sigma), ]
        if (!is.null(groupColumn) & !is.null(groupColors)) {
            stats_df[,groupColumn] <- factor(stats_df[,groupColumn], levels=names(groupColors))
        }
        p1 <- ggplot(stats_df, aes_string(x=idColumn, y="baseline_sigma", ymax=max("baseline_sigma"))) +
            base_theme + 
            xlab("") +
            ylab(expression(Sigma)) +
            geom_hline(yintercept=0, size=1*size, linetype=2, color="grey") +
            geom_point(size=3*size, position=position_dodge(0.6)) +
            geom_errorbar(aes_string(ymin="baseline_ci_lower", ymax="baseline_ci_upper"), 
                          width=0.2, size=0.5*size, alpha=0.8, position=position_dodge(0.6))
        if (!is.null(title)) {
            p1 <- p1 + ggtitle(title)
        }    
        if (is.null(groupColumn) & facetBy == "region") {
            p1 <- p1 + facet_grid(region ~ .)
        } else if (!is.null(groupColumn) & !is.null(groupColors) & facetBy == "region") {
            #groupColors <- factor(groupColors, levels=groupColors)
            p1 <- p1 + scale_color_manual(name=groupColumn, values=groupColors) +
                aes_string(color=groupColumn) + facet_grid(region ~ .)
        } else if (!is.null(groupColumn) & is.null(groupColors) & facetBy == "region") {
            p1 <- p1 + aes_string(color=groupColumn) + facet_grid(region ~ .)
        } else if (!is.null(groupColumn) & facetBy == "group") {
            p1 <- p1 + scale_color_manual(name="Region", values=REGION_PALETTE) +
                aes_string(color="region") + facet_grid(paste(groupColumn, "~ ."))
        } else {
            stop("Cannot facet by group if groupColumn=NULL")
        }
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { 
        plot(p1)
    } else {
        invisible(p1)
    }
}


#### Original BASELINe functions ####

##Covolution
break2chunks<-function(G=1000){
    base<-2^round(log(sqrt(G),2),0)
    return(c(rep(base,floor(G/base)-1),base+G-(floor(G/base)*base)))
}

PowersOfTwo <- function(G=100){
    exponents <- array()
    i = 0
    while(G > 0){
        i=i+1
        exponents[i] <- floor( log2(G) )
        G <- G-2^exponents[i]
    }
    return(exponents)
}

convolutionPowersOfTwo <- function( cons, length_sigma=4001 ){
    G = ncol(cons)
    if(G>1){
        for(gen in log(G,2):1){
            ll<-seq(from=2,to=2^gen,by=2)
            sapply(ll,function(l){cons[,l/2]<<-weighted_conv(cons[,l],cons[,l-1],length_sigma=length_sigma)})
        }
    }
    return( cons[,1] )
}

convolutionPowersOfTwoByTwos <- function( cons, length_sigma=4001,G=1 ){
    if(length(ncol(cons))) G<-ncol(cons)
    groups <- PowersOfTwo(G)
    matG <- matrix(NA, ncol=length(groups), nrow=length(cons)/G )
    startIndex = 1
    for( i in 1:length(groups) ){
        stopIndex <- 2^groups[i] + startIndex - 1
        if(stopIndex!=startIndex){
            matG[,i] <- convolutionPowersOfTwo( cons[,startIndex:stopIndex], length_sigma=length_sigma )
            startIndex = stopIndex + 1
        }
        else {
            if(G>1) matG[,i] <- cons[,startIndex:stopIndex]
            else matG[,i] <- cons
            #startIndex = stopIndex + 1
        }
    }
    return( list( matG, groups ) )
}

weighted_conv <- function(x, y, w=1, m=100, length_sigma=4001){
    lx<-length(x)
    ly<-length(y)
    if({lx<m}| {{lx*w}<m}| {{ly}<m}| {{ly*w}<m}){
        if(w<1){
            y1<-approx(1:ly,y,seq(1,ly,length.out=m))$y
            x1<-approx(1:lx,x,seq(1,lx,length.out=m/w))$y
            lx<-length(x1)
            ly<-length(y1)
        }
        else {
            y1<-approx(1:ly,y,seq(1,ly,length.out=m*w))$y
            x1<-approx(1:lx,x,seq(1,lx,length.out=m))$y
            lx<-length(x1)
            ly<-length(y1)
        }
    }
    else{
        x1<-x
        y1<-approx(1:ly,y,seq(1,ly,length.out=floor(lx*w)))$y
        ly<-length(y1)
    }
    tmp<-approx(x=1:(lx+ly-1),y=convolve(x1,rev(y1),type="open"),xout=seq(1,lx+ly-1,length.out=length_sigma))$y
    tmp[tmp<=0] = 0
    return(tmp/sum(tmp))
}


combineWeightedPosteriors <- function(PDF1, NumberOfSeq1, PDF2, NumberOfSeq2, length_sigma=4001){
    #return(list(calculate_bayesGHelper(list(cbind(PDF1,PDF2),c(log(NumberOfSeq1,2),log(NumberOfSeq2,2))),
    #                                   length_sigma=length_sigma),NumberOfSeq1+NumberOfSeq2))
    return( calculate_bayesGHelper( list( cbind(PDF1,PDF2),
                                          c(log(NumberOfSeq1,2),
                                            log(NumberOfSeq2,2))
    ),
    length_sigma=length_sigma
    ))
}

calculate_bayesGHelper <- function( listMatG,length_sigma=4001 ){
    matG <- listMatG[[1]]
    groups <- listMatG[[2]]
    i = 1
    resConv <- matG[,i]
    denom <- 2^groups[i]
    if(length(groups)>1){
        while( i<length(groups) ){
            i = i + 1
            resConv <- weighted_conv(resConv, matG[,i], w= {{2^groups[i]}/denom} ,length_sigma=length_sigma)
            #cat({{2^groups[i]}/denom},"\n")
            denom <- denom + 2^groups[i]
        }
    }
    return(resConv)
}

# Given a list of PDFs, returns a convoluted PDF
groupPosteriors <- function( listPosteriors, max_sigma=20, length_sigma=4001 ,Threshold=2 ){
    listPosteriors = listPosteriors[ !is.na(listPosteriors) ]
    Length_Postrior<-length(listPosteriors)
    if(Length_Postrior>1 & Length_Postrior<=Threshold){
        cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
        listMatG <- convolutionPowersOfTwoByTwos(cons,length_sigma=length_sigma)
        y<-calculate_bayesGHelper(listMatG,length_sigma=length_sigma)
        return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }else if(Length_Postrior==1) return(listPosteriors[[1]])
    else  if(Length_Postrior==0) return(NA)
    else {
        cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
        y = fastConv(cons,max_sigma=max_sigma, length_sigma=length_sigma )
        return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }
}

fastConv<-function(cons, max_sigma=20, length_sigma=4001){
    chunks<-break2chunks(G=ncol(cons))
    if(ncol(cons)==3) chunks<-2:1
    index_chunks_end <- cumsum(chunks)
    index_chunks_start <- c(1,index_chunks_end[-length(index_chunks_end)]+1)
    index_chunks <- cbind(index_chunks_start,index_chunks_end)
    
    case <- sum(chunks!=chunks[1])
    if(case==1) End <- max(1,((length(index_chunks)/2)-1))
    else End <- max(1,((length(index_chunks)/2)))
    
    firsts <- sapply(1:End,function(i){
        indexes<-index_chunks[i,1]:index_chunks[i,2]
        convolutionPowersOfTwoByTwos(cons[ ,indexes])[[1]]
    })
    if(case==0){
        result<-calculate_bayesGHelper( convolutionPowersOfTwoByTwos(firsts) )
    }else if(case==1){
        last<-list(calculate_bayesGHelper(
            convolutionPowersOfTwoByTwos( cons[ ,index_chunks[length(index_chunks)/2,1]:index_chunks[length(index_chunks)/2,2]] )
        ),0)
        result_first<-calculate_bayesGHelper(convolutionPowersOfTwoByTwos(firsts))
        result<-calculate_bayesGHelper(
            list(
                cbind(
                    result_first,last[[1]]),
                c(log(index_chunks_end[length(index_chunks)/2-1],2),log(index_chunks[length(index_chunks)/2,2]-index_chunks[length(index_chunks)/2,1]+1,2))
            )
        )
    }
    return(as.vector(result))
}



#' Calculate the BASELINe PDFs (including for regions that include CDR3 and FWR4)
#' 
#' \code{calcBaseline} calculates the BASELINe posterior probability density 
#' functions (PDFs) for sequences in the given Change-O \code{data.frame}.
#'
#' @param   db                  \code{data.frame} containing sequence data and annotations.
#' @param   sequenceColumn      \code{character} name of the column in \code{db} 
#'                              containing input sequences.
#' @param   germlineColumn      \code{character} name of the column in \code{db} 
#'                              containing germline sequences.
#' @param   testStatistic       \code{character} indicating the statistical framework 
#'                              used to test for selection. One of 
#'                              \code{c("local", "focused", "imbalanced")}.
#' @param   regionDefinition    \link{RegionDefinition} object defining the regions
#'                              and boundaries of the Ig sequences.
#' @param   targetingModel      \link{TargetingModel} object. Default is  \link{HH_S5F}.
#' @param   mutationDefinition  \link{MutationDefinition} object defining replacement
#'                              and silent mutation criteria. If \code{NULL} then 
#'                              replacement and silent are determined by exact 
#'                              amino acid identity. Note, if the input data.frame 
#'                              already contains observed and expected mutation frequency 
#'                              columns then mutations will not be recalculated and this
#'                              argument will be ignored.
#' @param   calcStats           \code{logical} indicating whether or not to calculate the 
#'                              summary statistics \code{data.frame} stored in the 
#'                              \code{stats} slot of a \link{Baseline} object.
#' @param   nproc               number of cores to distribute the operation over. If 
#'                              \code{nproc=0} then the \code{cluster} has already been
#'                              set and will not be reset.
#' @param   cloneColumn         \code{character} name of the column in \code{db} 
#'                              containing clonal identifiers. Relevant only for 
#'                              when regionDefinition includes CDR and FWR4 (else
#'                              this value can be \code{NULL})
#' @param   juncLengthColumn          \code{character} name of the column in \code{db} 
#'                              containing the junction length. Relevant only for 
#'                              when regionDefinition includes CDR and FWR4 (else
#'                              this value can be \code{NULL})  
#' @return  A \link{Baseline} object containing the modified \code{db} and BASELINe 
#'          posterior probability density functions (PDF) for each of the sequences.
#'           
#' @details 
#' Calculates the BASELINe posterior probability density function (PDF) for 
#' sequences in the provided \code{db}. 
#'          
#' \strong{Note}: Individual sequences within clonal groups are not, strictly speaking, 
#' independent events and it is generally appropriate to only analyze selection 
#' pressures on an effective sequence for each clonal group. For this reason,
#' it is strongly recommended that the input \code{db} contains one effective 
#' sequence per clone. Effective clonal sequences can be obtained by calling 
#' the \link{collapseClones} function.
#'                   
#' If the \code{db} does not contain the 
#' required columns to calculate the PDFs (namely mu_count & mu_expected)
#' then the function will:
#'   \enumerate{
#'   \item  Calculate the numbers of observed mutations.
#'   \item  Calculate the expected frequencies of mutations and modify the provided 
#'          \code{db}. The modified \code{db} will be included as part of the 
#'          returned \code{Baseline} object.
#' }
#'          
#' The \code{testStatistic} indicates the statistical framework used to test for selection. 
#' E.g.
#' \itemize{
#'   \item   \code{local} = CDR_R / (CDR_R + CDR_S).
#'   \item   \code{focused} = CDR_R / (CDR_R + CDR_S + FWR_S).
#'   \item   \code{imbalanced} = CDR_R + CDR_S / (CDR_R + CDR_S + FWR_S + FRW_R).
#' }
#' For \code{focused} the \code{regionDefinition} must only contain two regions. If more 
#' than two regions are defined the \code{local} test statistic will be used.
#' For further information on the frame of these tests see Uduman et al. (2011).
#'                              
#' @references
#' \enumerate{
#'   \item  Hershberg U, et al. Improved methods for detecting selection by mutation 
#'            analysis of Ig V region sequences. 
#'            Int Immunol. 2008 20(5):683-94.
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data.
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso See \link{Baseline} for the return object.
#'          See \link{groupBaseline} and \link{summarizeBaseline} for further processing.
#'          See \link{plotBaselineSummary} and \link{plotBaselineDensity} for plotting results.
#' 
#' @examples
#' # Load and subset example data
#' data(ExampleDb, package="alakazam")
#' db <- subset(ExampleDb, c_call == "IGHG" & sample_id == "+7d")
#' 
#' # Collapse clones
#' db <- collapseClones(db, cloneColumn="clone_id", 
#'                      sequenceColumn="sequence_alignment",
#'                      germlineColumn="germline_alignment_d_mask",
#'                      method="thresholdedFreq", minimumFrequency=0.6,
#'                      includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#'  
#' # Calculate BASELINe
#' baseline <- calcBaseline(db, 
#'                          sequenceColumn="clonal_sequence",
#'                          germlineColumn="clonal_germline", 
#'                          testStatistic="focused",
#'                          regionDefinition=IMGT_V,
#'                          targetingModel=HH_S5F,
#'                          nproc=1)
#'                          
#' @export
calcBaseline <- function(db, 
                         sequenceColumn = "clonal_sequence", 
                         germlineColumn = "clonal_germline", 
                         testStatistic = c("local","focused", "imbalanced"), 
                         regionDefinition = NULL, 
                         targetingModel = HH_S5F, 
                         mutationDefinition = NULL,
                         calcStats = FALSE, 
                         nproc = 1,
                         # following are relevant only when regionDefinition includes CDR3 and FWR4:
                         cloneColumn = NULL,
                         juncLengthColumn = NULL) {
    
    # Hack for visibility of foreach index variable
    idx <- NULL
    
    # Evaluate argument choices
    testStatistic <- match.arg(testStatistic)
    
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    regionDefinitionName <- ""
    if (!is.null(regionDefinition)) {
        regionDefinitionName <- regionDefinition@name
    }
    
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
    
    # Convert sequence columns to uppercase
    db <- toupperColumns(db, c(sequenceColumn, germlineColumn))
    db$tmp_baseline_row_id <- 1:nrow(db)
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, cpuCount())
    # nproc_arg will be passed to any function that has the nproc argument
    # If the cluster is already being set by the parent function then 
    # this will be set to 'cluster', that way the child function does not close
    # the connections and reset the cluster.
    #nproc_arg <- nproc
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if (nproc > 1) {        
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        parallel::clusterExport(cluster, list('db',
                                              'sequenceColumn', 'germlineColumn', 
                                              'cloneColumn', 'juncLengthColumn', 'makeRegion',
                                              'testStatistic', 'regionDefinition',
                                              'targetingModel', 'mutationDefinition','calcStats',
                                              'break2chunks', 'PowersOfTwo', 
                                              'convolutionPowersOfTwo', 
                                              'convolutionPowersOfTwoByTwos', 
                                              'weighted_conv', 
                                              'calculate_bayesGHelper', 
                                              'groupPosteriors', 'fastConv',
                                              'calcBaselineHelper',
                                              'c2s', 's2c', 'words', 'translate',
                                              'calcBaselineBinomialPdf','CONST_I',
                                              'BAYESIAN_FITTED',
                                              'calcObservedMutations','NUCLEOTIDES',
                                              'NUCLEOTIDES_AMBIGUOUS', 'IUPAC2nucs',
                                              'makeNullRegionDefinition',
                                              'getCodonPos','getContextInCodon',
                                              'mutationType',
                                              'AMINO_ACIDS','binMutationsByRegion',
                                              'collapseMatrixToVector','calcExpectedMutations',
                                              'calculateTargeting','HH_S5F','calculateMutationalPaths',
                                              'CODON_TABLE','IMGT_V_BY_REGIONS'), 
                                envir=environment() )    
        registerDoParallel(cluster, cores=nproc)
        #nproc_arg <- cluster
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # If db does not contain the required columns to calculate the PDFs (namely mu_count 
    # & mu_expected mutations), then the function will:
    #          1. Calculate the numbers of observed mutations
    #          2. Calculate the expected frequencies of mutations    
    # After that BASELINe prob. densities can be calcualted per sequence. 
    if (is.null(regionDefinition)) {
        rd_labels <- makeNullRegionDefinition()@labels
        observedColumns <- paste0("mu_count_", rd_labels)
        expectedColumns <- paste0("mu_expected_", rd_labels)
    } else {
        observedColumns <- paste0("mu_count_", regionDefinition@labels)
        expectedColumns <- paste0("mu_expected_", regionDefinition@labels)
    }
    
    if (!all(c(observedColumns, expectedColumns) %in% colnames(db))) {
        
        # If the germlineColumn & sequenceColumn are not found in the db error and quit
        if (!all(c(sequenceColumn, germlineColumn) %in% colnames(db))) {
            stop(paste0("Both ", sequenceColumn, " & ", germlineColumn, 
                        " columns need to be present in the db"))
        }
        
        message(paste0("calcBaseline will calculate observed and expected mutations for ",
                       sequenceColumn," using ", germlineColumn, " as a reference."))
        
        
        # Calculate the numbers of observed mutations
        db <- observedMutations(db,
                                sequenceColumn=sequenceColumn,
                                germlineColumn=germlineColumn,
                                regionDefinition=regionDefinition,
                                mutationDefinition=mutationDefinition,
                                frequency=FALSE, combine=FALSE,
                                nproc=0,
                                cloneColumn=cloneColumn,
                                juncLengthColumn=juncLengthColumn)
        
        # Calculate the expected frequencies of mutations
        db <- expectedMutations(db,
                                sequenceColumn=sequenceColumn,
                                germlineColumn=germlineColumn,
                                regionDefinition=regionDefinition,
                                targetingModel=targetingModel,
                                mutationDefinition=mutationDefinition,
                                nproc=0,
                                cloneColumn=cloneColumn,
                                juncLengthColumn=juncLengthColumn)
    } else {
        message(paste0("calcBaseline will use existing observed and expected mutations, in the fields: ",
                       paste0(observedColumns, sep="", collapse=", ")," and ", 
                       paste0(expectedColumns, sep="", collapse=", ")))
    }
    
    # Calculate PDFs for each sequence
    
    # Print status to console if not using extended regions definitions
    cat("Calculating BASELINe probability density functions...\n")
    
    # Number of sequences (used in foreach)
    totalNumbOfSequences <- nrow(db)
    # The column indexes of the mu_count_ and mu_expected_
    cols_observed <- grep( paste0("mu_count_"),  colnames(db) ) 
    cols_expected <- grep( paste0("mu_expected_"),  colnames(db) ) 
    

    # Exporting additional environment variables and functions needed to run foreach 
    if ( nproc>1 ) {
        parallel::clusterExport( 
            cluster, list('cols_observed', 'cols_expected','calcBaselineHelper'), 
            envir=environment() 
        )
        registerDoParallel(cluster)
    }
    
    list_pdfs <- list()
    list_numbOfSeqs <- list()
    list_k <- list()
    list_n <- list()
    list_p <- list()
    
    if (is.null(regionDefinition)) {
        regions <- makeNullRegionDefinition()@regions   
    } else {
        regions <- regionDefinition@regions
    }
    
    # For every region (e.g. CDR, FWR etc.)
    for (region in regions) {
        
        # Foreach returns a list of PDFs
        list_region_pdfs <- 
            foreach(idx=iterators::icount(totalNumbOfSequences)) %dopar% {  
                
                rd <- regionDefinition
                if (regionDefinitionName %in% c("IMGT_VDJ_BY_REGIONS","IMGT_VDJ")) {
                    ## Prepare extended region definition
                      rd <- makeRegion(juncLength = db[[juncLengthColumn]][idx],
                                       sequenceImgt = db[[sequenceColumn]][idx],
                                       regionDefinition=regionDefinition)
                }
                
                calcBaselineHelper( 
                    observed = db[cols_observed][idx,],
                    expected = db[cols_expected][idx,],
                    region = region,
                    testStatistic = testStatistic,
                    regionDefinition = rd
                )
            }
        
        # Count the number of non NA PDFs 
        list_numbOfSeqs[[region]] <- rep(1,totalNumbOfSequences)
        #is.na(list_region_pdfs)] <- 0
        
        # Convert the list of the region's PDFs into a matrix                
        
        mat_pdfs_binom <- 
            do.call( rbind, 
                     lapply( 
                         list_region_pdfs, 
                         function(x) { 
                             length(x) <- 4004 
                             return(x)
                         }
                     )
            )
        #cat(class(mat_pdfs_binom), "\n") # for debugging
        #cat(dim(mat_pdfs_binom), "\n") # for debugging
        
        # IMPORTANT: if input has a single sequence, mat_pdfs_binom (1-row) gets coerced 
        # into a numeric vector without matrix(..., nrow=nrow(mat_pdfs_binom))
        list_pdfs[[region]] <- matrix(mat_pdfs_binom[, 1:4001], nrow=nrow(mat_pdfs_binom))
        #cat(class(list_pdfs[[region]]), "\n") # for debugging
        stopifnot(is(list_pdfs[[region]], "matrix"))
        
        list_k[[region]] <- mat_pdfs_binom[, 4002]
        list_n[[region]] <- mat_pdfs_binom[, 4003]
        list_p[[region]] <- mat_pdfs_binom[, 4004]
        list_numbOfSeqs[[region]][is.na(list_k[[region]])] <- 0
    }
    
    # WIP - check Milca's funcittion, hotw it handles hte new data struct for rextender regions
    # Template for values for the regions
    mat_template <- matrix( NA, 
                            ncol=length(regions), 
                            nrow=totalNumbOfSequences,
                            dimnames=list(1:totalNumbOfSequences, regions)
    )
    
    # numbOfSeqs
    # This holds the number of non NA sequences
    numbOfSeqs <- mat_template
    for(region in regions){
        numbOfSeqs[,region] <-   list_numbOfSeqs[[region]]
    }
    
    # binomK
    # This holds the number of exact successin in the binomial trials
    binomK <- mat_template
    for(region in regions){
        binomK[,region] <-   list_k[[region]]
    }
    
    # binomN
    # This holds the total numbers trials in the binomial
    binomN <- mat_template
    for(region in regions){
        binomN[,region] <-   list_n[[region]]
    }
    
    # binomP
    # This holds the prob of successin in the binomial trials
    binomP <- mat_template
    for(region in regions){
        binomP[,region] <-   list_p[[region]]
    }
    
    # If regionDefinition include CDR3 and FWR4 - then the regionDefinition is different for each clone,
    # In this case - regionDefinition will be NULL.
    
    # Create a Baseline object with the above results to return
    baseline <- createBaseline(description="",
                               db=as.data.frame(db) %>% 
                                   arrange(!!rlang::sym("tmp_baseline_row_id")) %>% 
                                   select(-!!rlang::sym("tmp_baseline_row_id")),
                               regionDefinition=regionDefinition,
                               testStatistic=testStatistic,
                               regions=regions,
                               numbOfSeqs=numbOfSeqs,
                               binomK=binomK,
                               binomN=binomN,
                               binomP=binomP,
                               pdfs=list_pdfs )
    
    # Calculate BASELINe stats and update slot
    if (calcStats==TRUE) {
        baseline <- summarizeBaseline(baseline)
    }
    
    # Stop cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    return(baseline)
    
}


