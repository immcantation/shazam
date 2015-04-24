#' @include RegionDefinitions.R
#' @include shm.R
NULL

#### Classes ####

#' S4 class defining a BASELINe selection object
#' 
#' \code{BASELINe} defines a common data structure for defining the BASELINe selection
#' results
#' 
#' @slot    id                              Unique identifier.
#' @slot    description                     Description of the sequence.
#' @slot    inputSequence                   Nucleotide sequence being tested for selection.
#' @slot    germlineSequence                The germline or reference sequence.
#' @slot    regionDefinition                \code{\link{RegionDefinition}} object defining the regions
#'                                          and boundaries of the Ig sequences. Note, only the part of
#'                                          sequences defined in \code{regionDefinition} are analyzed.
#'                                          Any mutations outside the definition will be ignored. E.g.
#'                                          If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                                          used, then mutations in positions greater than 312 will not
#'                                          be counted.
#' @slot    observedMutations               The number of mutations observed in the sequence
#' @slot    expectedMutationFrequencies     The expected frequencies of mutations, based on the
#'                                          \code{germlineSequence}.
#' @slot    numbOfSequences                 An \code{array} of the number of non NA sequences (i.e. 
#'                                          sequences with mutations), for each region (defined in the
#'                                          \code{defineRegions}).
#' @slot    regions                         The levels of the boundaries (e.g. CDR & FWR).
#' @slot    labels                          The labels for the boundary/mutations combinations.
#'                                          e.g. CDR_R CDR_S FWR_R, FWR_S.
#' @slot    probDensity                     A \code{list} (one for each sequence in the provided \code{db}
#'                                          of \code{list(s)} (one for each region) contianing the probability
#'                                          density function(s).
#'                          
#' @name BASELINe
#' @export
setClass("BASELINe", 
         slots=c(id="character",
                 description="character",
                 inputSequence="character",
                 germlineSequence="character",
                 observedMutations="array",
                 expectedMutationFrequencies="array",
                 regionDefinition="RegionDefinition",
                 regions="character",
                 labels="character",
                 probDensity="list")
)

#### BASELINe building functions #####

#' Creates a BASELINe object
#' 
#' \code{createBASELINe} creates a \code{RegionDefinition}.
#'
#' 
#' @param   id                              Unique identifier.
#' @param   description                     Description of the sequence.
#' @param   inputSequence                   Nucleotide sequence being tested for selection.
#' @param   germlineSequence                The germline or reference sequence.
#' @param   regionDefinition                \code{\link{RegionDefinition}} object defining the regions
#'                                          and boundaries of the Ig sequences. Note, only the part of
#'                                          sequences defined in \code{regionDefinition} are analyzed.
#'                                          Any mutations outside the definition will be ignored. E.g.
#'                                          If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                                          used, then mutations in positions greater than 312 will not
#'                                          be counted.
#' @param   observedMutations               The number of mutations observed in the sequence
#' @param   expectedMutationFrequencies     The expected frequencies of mutations, based on the
#'                                          \code{germlineSequence}.
#' @param   probDensity                     A \code{list} (one for each sequence in the provided \code{db}
#'                                          of \code{list(s)} (one for each region) contianing the probability
#'                                          density function(s).
#' 
#' @return   A \code{BASELINe} object.
#' 
#' @seealso  See \code{\link{BASELINe}} for the return object.
#' 
#' @examples
#' library(shm)
#' 
#' 
#' @export
createBASELINe <- function(id="",
                           description="",
                           inputSequence="",
                           germlineSequence="",
                           regionDefinition=IMGT_V_NO_CDR3,
                           observedMutations=array(NA,1),
                           expectedMutationFrequencies=array(NA,1),
                           probDensity=list()){
    
    
    if (class(observedMutations) == "data.frame") {
        observedMutations <- array( as.numeric(observedMutations), 
                                    dimnames=list(colnames(observedMutations)) )
    }
    
    if (class(expectedMutationFrequencies) == "data.frame") {
        expectedMutationFrequencies <- array( as.numeric(expectedMutationFrequencies), 
                                              dimnames=list(colnames(expectedMutationFrequencies)) )
    }
    
    # Define RegionDefinition object
    seqBASELINe <- new("BASELINe",
                       id=id,
                       description=description,
                       inputSequence=inputSequence,
                       observedMutations=observedMutations,
                       expectedMutationFrequencies=expectedMutationFrequencies,
                       regionDefinition=regionDefinition,
                       regions=regionDefinition@regions,
                       labels=regionDefinition@labels,
                       probDensity=probDensity)
    
    return(seqBASELINe)
}


# Edit the BASELINe object
# 
# \code{editBASELINe} edits a \code{BASELINe}.
#
# @param   seqBASELINe     The \code{BASELINe} S4 object to be edited.
# @param   field_name      Name of the field in the \code{BASELINe} S4 object to be edited.
# @param   value           The value to set the \code{field_anme}.
# 
# @return   A \code{BASELINe} object.
# 
# @seealso  See \code{\link{BASELINe}} for the return object.
editBASELINe <- function ( seqBASELINe,
                           field_name,
                           value ) {
    if ( !match(field_name, slotNames(seqBASELINe)) ) { stop("field_name not part of BASELINe object!")}
    slot(seqBASELINe, field_name) = value
    return(seqBASELINe)
}





#### BASELINe selection calculating functions ####

#' Return BASELINe probability density functions for each entry in the DB
#'
#' \code{getBASELINePDF} calculates 
#'
#' @param       db                  data.frame containing sequence data.
#' @param       testStatistic       The statistical framework used to test for selection.
#'                                  \code{local} = CDR_R / (CDR_R + CDR_S)
#'                                  \code{focused} = CDR_R / (CDR_R + CDR_S + FWR_S)
#'                                  See Uduman et al. (2011) for further information.
#' @param       groupBy             Aggregate BASELINe probability density functions from 
#' @param       sequenceColumn      name of the column containing sample/input sequences.
#' @param       germlineColumn      name of the column containing germline sequences.
#' @param       regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                                  and boundaries of the Ig sequences. Note, only the part of
#'                                  sequences defined in \code{regionDefinition} are analyzed.
#'                                  Any mutations outside the definition will be ignored. E.g.
#'                                  If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                                  used, then mutations in positions greater than 312 will not
#'                                  be counted.
#' @param       nproc               number of cores to distribute the operation over. If 
#'                                  \code{nproc} = 0 then the \code{cluster} has already been
#'                                  set and will not be reset.
#' 
#' @return      A \code{list} of length two. The first element contains the modified \code{db}
#'              and the second element is a \code{list} of \code{\link{BASELINe}} objects
#'              (elements match the sequences in \code{db}).
#'           
#' @details     \code{getBASELINe} calculates the BASELINe probability density functions for each
#'              sequence in the \code{db.}.
#' 
#' @seealso     To calculate BASELINe statistics, such as the mean selection strength
#'              and the 95\% confidence interval, see \code{\link{getBASELINeStats}}.
#' 
#' 
#' @references
#' \enumerate{
#'   \item  Gur Yaari; Mohamed Uduman; Steven H. Kleinstein. Quantifying selection 
#'          in high-throughput Immunoglobulin sequencing data sets. Nucleic Acids Res.
#'           2012 May 27. 
#'  \item   Mohamed Uduman; Gur Yaari; Uri Hershberg; Mark J. Shlomchik; Steven H. 
#'          Kleinstein. Detecting selection in immunoglobulin sequences. Nucleic Acids 
#'          Res. 2011 Jul;39(Web Server issue):W499-504.  
#'  \item   Hershberg U, Uduman M, Shlomchik MJ, Kleinstein SH. Improved methods for
#'          detecting selection by mutation analysis of Ig V region sequences. Int Immunol.
#'          2008 May;20(5):683-94. doi: 10.1093/intimm/dxn026. Epub 2008 Apr 7. 
#' }
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#'                      
#' @export
getBASELINePDF <- function( db,
                            sequenceColumn="SEQUENCE_IMGT",
                            germlineColumn="GERMLINE_IMGT_D_MASK",
                            testStatistic=c("local","focused"),
                            regionDefinition=IMGT_V_NO_CDR3,
                            nproc=1 ) {
    
    # Evaluate argument choices
    testStatistic <- match.arg(testStatistic, c("local", "focused"))
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- makeCluster(nproc, type = "SOCK")
        clusterExport( cluster, list('db', 'sequenceColumn', 'germlineColumn', 'regionDefinition'), envir=environment() )
        clusterEvalQ( cluster, library("shm") )
        registerDoSNOW(cluster)
    } else if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    
    # If db does not contain the correct OBSERVED & MUTATIONS columns then first
    # create them. After that BASELINe prob. densities can be calcualted per seqeunce.
    observedColumns <- paste0("OBSERVED_", regionDefinition@labels)
    expectedColumns <- paste0("EXPECTED_", regionDefinition@labels)
    
    if( !all( c(observedColumns,expectedColumns) %in% colnames(db) ) ) {
        
        # checking the validity of input parameters
        if ( !all( c(!is.null(sequenceColumn), !is.null(germlineColumn)) ) ) {
            stop("The OBSERVED and EXPECTED values are not in DB. Please specify the 'sequenceColumn', 'germlineColumn'")
        }
        
        db <- getClonalConsensus(db, 
                                 cloneColumn="CLONE", 
                                 sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 collapseByClone=TRUE,nproc=0)
        
        db <- getObservedMutations(db,
                                   sequenceColumn="CLONAL_CONSENSUS_SEQUENCE",
                                   germlineColumn="GERMLINE_IMGT_D_MASK",
                                   regionDefinition=IMGT_V_NO_CDR3,
                                   nproc=0)
        
        
        db <- getExpectedMutationFrequencies(db,
                                             sequenceColumn="CLONAL_CONSENSUS_SEQUENCE",
                                             germlineColumn="GERMLINE_IMGT_D_MASK",
                                             regionDefinition=IMGT_V_NO_CDR3,
                                             nproc=0)
    }
    
    # Compute the BASELINe prob. density values for all the seqeunces in the DB
    # and for each of the regions, defined in the regionDefinition
    
    # Printing status to console
    cat("Calculating BASELINe probability density functions...\n")
    
    numbOfSeqs <- nrow(db)
    db_BASELINe <-
        foreach( i=icount(numbOfSeqs) ) %dopar% {
            seqBASELINe <- createBASELINe( id = db[i,"SEQUENCE_ID"],
                                           observedMutations = db[i, observedColumns],
                                           expectedMutationFrequencies = db[i, expectedColumns] )
            return(
                calculateBASELINePDF( seqBASELINe, 
                                      testStatistic=testStatistic,
                                      regionDefinition=regionDefinition) 
            )
        }
    
    
    # Properly shutting down the cluster
    if(nproc>1){ stopCluster(cluster) }
    
    return(list("db"=db, "list_BASELINe"=db_BASELINe))    
}


#' Return DB of BASELINe stats
#'
#' \code{getBASELINeStats} calculates 
#'
#' @param       db              \code{data.frame} containing sequence data.
#' @param       list_BASELINe   \code{list} of \code{\link{BASELINe}} objects, length matching
#'                              the number of entires in the \code{db}. (This is returned by 
#'                              \code{\link{getBASELINePDF}})
#' @param       nproc           number of cores to distribute the operation over. If 
#'                              \code{nproc} = 0 then the \code{cluster} has already been
#'                              set and will not be reset.
#' 
#' @return      A modified \code{db} data.frame with the BASELINe selection values and 95%
#'              confidence intervals. The columns names are dynamically created based on the
#'              regions in the \code{regionDefinition}. E.g. For the default
#'              \code{\link{IMGT_V_NO_CDR3}} definition, which defines positions for CDR and
#'              FWR, the following columns are added:  
#'              \itemize{
#'                  \item  \code{BASELINe_CDR}:  BASELINe Sigma value for the CDR
#'                  \item  \code{BASELINe_CDR_95CI_LOWER}: 95% CI for the CDR
#'                  \item  \code{BASELINe_FWR}:  BASELINe Sigma value for the FWR
#'                  \item  \code{BASELINe_FWR_95CI_LOWER}: 95% CI for the FWR
#'              }
#'           
#' @details     \code{getBASELINeStats} calculates 
#' 
#' @seealso     See \code{\link{getBASELINePDF}}.
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#'                      
#' @export
getBASELINeStats <- function ( db,
                               list_BASELINe,
                               nproc=1 ) {
    
    
    # Checking if the nrows of db is the same as the length of list_BASELINe
    if ( length(list_BASELINe)!= nrow(db) ) {
        stop("The number of entires in db and list_BASELINe do not match.")
    }
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- makeCluster(nproc, type = "SOCK")
        clusterExport( cluster, list( 'db', 'list_BASELINe', 
                                      'calculateBASELINeStats',
                                      'calculateBASELINeSigma',
                                      'calculateBASELINeCI',
                                      'calculateBASELINePvalue' ), 
                       envir=environment() )
        clusterEvalQ( cluster, library("shm") )
        registerDoSNOW(cluster)
    } else if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Printing status to console
    cat("Calculating BASELINe statistics...\n")
    
    numbOfSeqs <- nrow(db)
    BASELINe_stats_list <-
        foreach( i=icount(numbOfSeqs) ) %dopar% {
            seqBASELINe <- list_BASELINe[[i]]
            seqBASELINe_stats <- list()
            for ( region in seqBASELINe@regions ) {
                seqBASELINe_stats[[region]] <- calculateBASELINeStats( seqBASELINe@probDensity[[region]] )
            }
            seqBASELINe_stats <- unlist( seqBASELINe_stats )
            names(seqBASELINe_stats) <- gsub("\\.","_",names(seqBASELINe_stats))
            return( seqBASELINe_stats )
        }
    
    # Convert list of BASELINe stats into a data.frame
    BASELINe_stats <- do.call(rbind, BASELINe_stats_list)
    
    # Bind the stats to db
    db_new <- cbind(db, BASELINe_stats)
    return(db_new)   
}

#' Calculate BASELINe
#'
#' \code{calculateBASELINePDF} calculates 
#'
#' @param   seqBASELINe         \code{\link{BASELINe}} object representing the sequence
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences. Note, only the part of
#'                              sequences defined in \code{regionDefinition} are analyzed.
#'                              Any mutations outside the definition will be ignored. E.g.
#'                              If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                              used, then mutations in positions greater than 312 will not
#'                              be counted.
#' 
#' @return  A modified \code{\link{BASELINe}} object with the BASELINe probability 
#'          density function calculated for the regions defined in the \code{regionDefinition}.
#'           
#' @details
#' \code{getBASELINe} calculates 
#' 
#' @seealso  
#' See also
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#'                      
#' @export
calculateBASELINePDF  <- function( seqBASELINe,
                                   testStatistic=c("local", "focused"),
                                   regionDefinition=IMGT_V_NO_CDR3 ) {
    
    # Evaluate argument choices
    testStatistic <- match.arg(testStatistic, c("local", "focused"))
    
    #If there are more than two regions (e.g. CDR and FWR then you cannot perform the focused test)
    if (testStatistic=="focused" & length(seqBASELINe@regions)!=2) {
        testStatistic="local"    
    }
    
    list_BASELINE_ProbDensity <- list()
    for (region in seqBASELINe@regions) {
        
        # local test statistic
        if (testStatistic == "local") { 
            obsX_Index <- grep( paste0("OBSERVED_", region,"_R"),  names(seqBASELINe@observedMutations) )
            obsN_Index <- grep( paste0("OBSERVED_", region),  names(seqBASELINe@observedMutations) )
            
            expX_Index <- grep( paste0("EXPECTED_", region,"_R"),  names(seqBASELINe@expectedMutationFrequencies) )
            expN_Index <- grep( paste0("EXPECTED_", region),  names(seqBASELINe@expectedMutationFrequencies) )       
        }
        
        # focused test statistic
        if (testStatistic == "focused") { 
            obsX_Index <- grep( paste0("OBSERVED_", region,"_R"),  names(seqBASELINe@observedMutations) )
            obsN_Index <- 
                grep( 
                    paste0( "OBSERVED_", region, "|", 
                            "OBSERVED_", seqBASELINe@regions[seqBASELINe@regions!=region], "_S"
                    ),
                    names(seqBASELINe@observedMutations) 
                )
            
            expX_Index <- grep( paste0("EXPECTED_", region,"_R"),  names(seqBASELINe@expectedMutationFrequencies) )
            expN_Index <- 
                grep( 
                    paste0( "EXPECTED_", region, "|", 
                            "EXPECTED_", seqBASELINe@regions[seqBASELINe@regions!=region], "_S"
                    ),
                    names(seqBASELINe@expectedMutationFrequencies) 
                )
            
        }     
        
        obsX <- seqBASELINe@observedMutations[obsX_Index]
        obsN <- sum( seqBASELINe@observedMutations[obsN_Index], na.rm=T )
        
        expP <- 
            seqBASELINe@expectedMutationFrequencies[expX_Index] / 
            sum( seqBASELINe@expectedMutationFrequencies[expN_Index], na.rm=T )
        
        list_BASELINE_ProbDensity[[region]] <- calculateBASELINeBinomialPDF( x=obsX, n=obsN, p=expP)
    }
    
    seqBASELINe <- editBASELINe(seqBASELINe, "probDensity", list_BASELINE_ProbDensity)
    
    return(seqBASELINe)
}




# Calculate the BASELINe probability function in a
# binomial framework.
calculateBASELINeBinomialPDF <- function ( x=3, 
                                           n=10, 
                                           p=0.33,
                                           CONST_i=CONST_I,
                                           max_sigma=20,
                                           length_sigma=4001 ) {
    if(n!=0){
        sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)
        sigma_1<-log({CONST_i/{1-CONST_i}}/{p/{1-p}})
        index<-min(n,60)
        y<- dbeta(CONST_i,x+BAYESIAN_FITTED[index],n+BAYESIAN_FITTED[index]-x)*(1-p)*p*exp(sigma_1)/({1-p}^2+2*p*{1-p}*exp(sigma_1)+{p^2}*exp(2*sigma_1))
        if(!sum(is.na(y))){
            tmp<-approx(sigma_1,y,sigma_s)$y
            tmp/sum(tmp)/{2*max_sigma/{length_sigma-1}}
        }else{
            return(NA)
        }
    }else{
        return(NA)
    }
}


# Given a BASELIne PDF calculate mean sigma
calculateBASELINeSigma <- function ( baseline_pdf,
                                     max_sigma=20,
                                    length_sigma=4001 ) {
    
    if ( length(baseline_pdf)!=length_sigma) { return(NA) }
    
    sigma_s <- seq(-max_sigma, max_sigma, length.out=length_sigma)
    norm = {length_sigma-1}/2/max_sigma
    return( (baseline_pdf%*%sigma_s/norm)  )
}


# Given a BASELIne PDF calculate Confidence Interval
calculateBASELINeCI <- function ( baseline_pdf,
                                  low=0.025,
                                  up=0.975,
                                  max_sigma=20,
                                  length_sigma=4001 ){
    
    if ( length(baseline_pdf)!=length_sigma ) { return( c(NA,NA) ) }
    
    sigma_s <- seq(-max_sigma, max_sigma, length.out=length_sigma)
    cdf <- cumsum(baseline_pdf)
    cdf <- cdf/cdf[length(cdf)]
    intervalLow <- findInterval(low,cdf)
    fractionLow <- (low - cdf[intervalLow])/(cdf[intervalLow+1]-cdf[intervalLow])
    intervalUp <- findInterval(up,cdf)
    fractionUp <- (up - cdf[intervalUp])/(cdf[intervalUp]-cdf[intervalUp-1])
    sigmaLow <- sigma_s[intervalLow]+fractionLow*(sigma_s[intervalLow+1]-sigma_s[intervalLow])
    sigmaUp <- sigma_s[intervalUp]+fractionUp*(sigma_s[intervalUp+1]-sigma_s[intervalUp])
    return( c(sigmaLow,sigmaUp) )
}

# Given a BASELIne PDF calculate P value
calculateBASELINePvalue <- function ( baseline_pdf, 
                                      length_sigma=4001, 
                                      max_sigma=20 ){
    if ( length(baseline_pdf)>1 ) {
        norm = {length_sigma-1}/2/max_sigma
        pVal = {sum(baseline_pdf[1:{{length_sigma-1}/2}]) + baseline_pdf[{{length_sigma+1}/2}]/2}/norm
        if(pVal>0.5){
            pVal = pVal-1
        }
        return(pVal)
    }else{
        return(NA)
    }
}


# Given a BASELIne PDF calculate Mean, Confidence Interval (lower & upper) and P value
calculateBASELINeStats <- function ( baseline_pdf,
                                     low=0.025,
                                     up=0.975,
                                     max_sigma=20, 
                                     length_sigma=4001 ){
    
    # if NA (i.e. length of baseline_pdf is 1)
    if ( length(baseline_pdf)==1  ) { return(rep(NA,4)) }
    
    
    baselineSigma <- calculateBASELINeSigma( baseline_pdf=baseline_pdf, 
                                             max_sigma=max_sigma, 
                                             length_sigma=length_sigma )
    
    
    baselineCI <- calculateBASELINeCI( baseline_pdf=baseline_pdf,
                                       low=low,
                                       up=up,
                                       max_sigma=max_sigma,
                                       length_sigma=length_sigma )
    
    baselinePvalue <- calculateBASELINePvalue( baseline_pdf=baseline_pdf,
                                               max_sigma=max_sigma,
                                               length_sigma=length_sigma )
    
    return( c( "Sigma"=baselineSigma, 
               "CI_Lower"=baselineCI[1], 
               "CI_Upper"=baselineCI[2],
               "Pvalue"=baselinePvalue 
               ) 
            )
}


